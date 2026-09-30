// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Code for EveryBeam beam calculations.
//!
//! EveryBeam (<https://git.astron.nl/RD/EveryBeam>) is a C++ library that
//! provides beam models for many telescopes (e.g. SKA-Low, LOFAR, MWA,
//! OSKAR-simulated arrays). Telescopes are described by a measurement set.
//!
//! EveryBeam needs the absolute time and the J2000 (RA, Dec) of each
//! direction; see [`BeamTime`]. EveryBeam can only calculate beam responses on
//! the CPU; if GPU modelling is used, EveryBeam beam responses are calculated
//! on the CPU and copied to the GPU.
//!
//! # Polarisation conventions
//!
//! As documented in `FluxDensity::to_inst_stokes`, hyperdrive orders its
//! polarisations with X as East-West and Y as North-South, i.e. the sky basis
//! is (East, North).
//!
//! Without normalisation (or with "amplitude" normalisation, which is a scalar),
//! EveryBeam's Jones matrices have rows corresponding to the telescope's feeds
//! and columns corresponding to the (North, East) sky basis (with respect to
//! the J2000 celestial pole). Here, the columns are swapped, and the rows are
//! left alone (so that they continue to match the order of the feeds in the
//! visibilities). Note that in this case, the basis of the rows can depend on
//! EveryBeam's element model.
//!
//! With "full" (or "preapplied") normalisation, EveryBeam left-multiplies the
//! response by the inverse of the response at the beam centre, so both the rows
//! and the columns are in the (North, East) sky basis; the response is the
//! identity at the beam centre. Here, both the rows and columns are swapped, so
//! that they are ordered (East, North). For arrays with East-West and
//! North-South dipoles (e.g. SKA-Low, MWA), the rows then correspond to the
//! (X, Y) feeds, up to a sign. This is the default.

mod ffi;
#[cfg(test)]
mod tests;

use std::{
    ffi::{c_char, CStr, CString},
    os::unix::ffi::OsStrExt,
    path::{Path, PathBuf},
    ptr::NonNull,
    sync::Arc,
};

use log::debug;
use marlu::{AzEl, Jones, RADec};
use ndarray::prelude::*;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::{Beam, BeamError, BeamTime, BeamType};

#[cfg(any(feature = "cuda", feature = "hip"))]
use super::{BeamGpu, DevicePointer, GpuFloat};

/// The size of the buffers used to receive error messages from the shim.
const ERR_LEN: usize = 1024;

/// The number of directions given to EveryBeam per call. Each call is
/// independent and is given to a rayon thread.
const DIRECTION_CHUNK_SIZE: usize = 64;

/// Options for EveryBeam. Any options that aren't set use EveryBeam's defaults.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct EveryBeamOptions {
    /// The element response model, e.g. "default", "hamaker", "lobes",
    /// "oskar_dipole", "oskar_spherical_wave", "skala40_wave", "skalow_feko".
    pub element_response_model: Option<String>,

    /// The beam mode, e.g. "full", "array_factor", "element".
    pub beam_mode: Option<String>,

    /// The beam normalisation mode, e.g. "none", "full", "amplitude",
    /// "preapplied".
    pub beam_normalisation_mode: Option<String>,

    /// The path to element-response coefficients (telescope dependent).
    pub coeff_path: Option<PathBuf>,

    /// The data column to use for LOFAR "preapplied" beam information.
    pub data_column_name: Option<String>,

    /// The measurement set field ID to use for the beam pointing.
    pub field_id: usize,

    /// Use the subband (i.e. reference) frequency for the station
    /// beamformer, rather than each channel's frequency (LOFAR only).
    pub use_subband_frequency: bool,

    /// Interpolate the beam over frequency (MWA only).
    pub frequency_interpolation: bool,
}

/// An owned handle to a telescope loaded by EveryBeam.
struct Telescope(NonNull<ffi::eb_telescope>);

// The shim only calls const methods on the telescope, and each call to
// `eb_point_responses` uses its own (non-shared) point-response object.
unsafe impl Send for Telescope {}
unsafe impl Sync for Telescope {}

impl Drop for Telescope {
    fn drop(&mut self) {
        unsafe { ffi::eb_free(self.0.as_ptr()) }
    }
}

/// A wrapper around a telescope loaded by EveryBeam that implements the
/// [`Beam`] trait.
pub(crate) struct EveryBeam {
    inner: Arc<EveryBeamInner>,
}

struct EveryBeamInner {
    telescope: Telescope,
    num_stations: usize,
    ms: PathBuf,
    options: EveryBeamOptions,
}

fn opt_cstring(s: Option<&str>) -> Result<Option<CString>, BeamError> {
    s.map(|s| CString::new(s).map_err(|_| BeamError::EveryBeam(format!("Invalid string '{s}'"))))
        .transpose()
}

fn cstring_ptr(s: &Option<CString>) -> *const c_char {
    s.as_ref().map(|s| s.as_ptr()).unwrap_or(std::ptr::null())
}

fn err_to_string(err: &[c_char]) -> String {
    unsafe { CStr::from_ptr(err.as_ptr()) }
        .to_string_lossy()
        .into_owned()
}

/// Convert an EveryBeam Jones matrix with columns in the (North, East) sky
/// basis to hyperdrive's convention (columns in the (East, North) sky basis).
/// The rows are unchanged.
fn swap_columns(j: Jones<f64>) -> Jones<f64> {
    Jones::from([j[1], j[0], j[3], j[2]])
}

/// Convert an EveryBeam Jones matrix with rows and columns in the (North, East)
/// sky basis to hyperdrive's convention (rows and columns in the (East, North)
/// sky basis).
fn swap_rows_and_columns(j: Jones<f64>) -> Jones<f64> {
    Jones::from([j[3], j[2], j[1], j[0]])
}

/// Does this EveryBeam normalisation mode normalise the beam by
/// left-multiplying by an inverse Jones matrix (which changes the basis of the
/// rows)? `None` means EveryBeam's default, which is no normalisation.
fn is_matrix_normalisation(mode: Option<&str>) -> bool {
    match mode {
        None => false,
        Some(m) => {
            let m = m.to_lowercase().replace('_', "");
            matches!(m.as_str(), "full" | "preapplied" | "preappliedorfull")
        }
    }
}

impl EveryBeam {
    /// Load the telescope described by the measurement set `ms`. If
    /// `num_tiles` is supplied, it is checked against the number of stations
    /// in the measurement set.
    pub(crate) fn new(
        ms: &Path,
        num_tiles: Option<usize>,
        options: EveryBeamOptions,
    ) -> Result<EveryBeam, BeamError> {
        debug!("Loading EveryBeam telescope from {}", ms.display());
        let ms_c = CString::new(ms.as_os_str().as_bytes())
            .map_err(|_| BeamError::EveryBeam(format!("Invalid path {}", ms.display())))?;
        let element_response_model = opt_cstring(options.element_response_model.as_deref())?;
        let beam_mode = opt_cstring(options.beam_mode.as_deref())?;
        let beam_normalisation_mode = opt_cstring(options.beam_normalisation_mode.as_deref())?;
        let coeff_path = opt_cstring(
            options
                .coeff_path
                .as_ref()
                .map(|p| p.to_str().ok_or_else(|| BeamError::EveryBeam(format!("Invalid path {}", p.display()))))
                .transpose()?,
        )?;
        let data_column_name = opt_cstring(options.data_column_name.as_deref())?;
        let eb_options = ffi::eb_options {
            element_response_model: cstring_ptr(&element_response_model),
            beam_mode: cstring_ptr(&beam_mode),
            beam_normalisation_mode: cstring_ptr(&beam_normalisation_mode),
            coeff_path: cstring_ptr(&coeff_path),
            data_column_name: cstring_ptr(&data_column_name),
            use_channel_frequency: (!options.use_subband_frequency).into(),
            frequency_interpolation: options.frequency_interpolation.into(),
        };

        let mut err = [0 as c_char; ERR_LEN];
        let telescope =
            unsafe { ffi::eb_load(ms_c.as_ptr(), &eb_options, err.as_mut_ptr(), ERR_LEN) };
        let telescope = NonNull::new(telescope).ok_or_else(|| {
            BeamError::EveryBeam(format!(
                "Couldn't load telescope from {}: {}",
                ms.display(),
                err_to_string(&err)
            ))
        })?;
        let telescope = Telescope(telescope);
        let num_stations = unsafe { ffi::eb_num_stations(telescope.0.as_ptr()) };
        let beam = EveryBeam {
            inner: Arc::new(EveryBeamInner {
                telescope,
                num_stations,
                ms: ms.to_path_buf(),
                options,
            }),
        };
        debug!("EveryBeam telescope has {num_stations} stations");

        if let Some(num_tiles) = num_tiles {
            if num_tiles != num_stations {
                return Err(BeamError::StationCountMismatch {
                    ms: num_stations,
                    tiles: num_tiles,
                });
            }
        }

        Ok(beam)
    }
}

impl EveryBeamInner {
    /// Calculate beam responses for (RA, Dec) directions, at all specified
    /// frequencies and stations. The results are written into `results`, which
    /// must have dimensions (stations, freqs, directions). The responses are in
    /// hyperdrive's polarisation convention.
    pub(crate) fn calc_jones_radec(
        &self,
        radecs: &[RADec],
        freqs_hz: &[f64],
        stations: &[usize],
        time: BeamTime,
        mut results: ArrayViewMut3<Jones<f64>>,
    ) -> Result<(), BeamError> {
        assert_eq!(results.dim(), (stations.len(), freqs_hz.len(), radecs.len()));
        if let Some(&s) = stations.iter().find(|&&s| s >= self.num_stations) {
            return Err(BeamError::BadTileIndex {
                got: s,
                max: self.num_stations - 1,
            });
        }
        if radecs.is_empty() || freqs_hz.is_empty() || stations.is_empty() {
            return Ok(());
        }

        let (ras, decs): (Vec<f64>, Vec<f64>) = radecs.iter().map(|r| (r.ra, r.dec)).unzip();
        let time_mjd_s = time.epoch.to_mjd_utc_seconds();
        let num_directions = radecs.len();
        let station_stride = freqs_hz.len() * num_directions;
        let freq_stride = num_directions;

        // Each chunk of directions writes to disjoint parts of `results`.
        struct SendPtr(*mut f64);
        unsafe impl Send for SendPtr {}
        unsafe impl Sync for SendPtr {}
        let out = SendPtr(
            results
                .as_slice_mut()
                .expect("is contiguous")
                .as_mut_ptr()
                .cast::<f64>(),
        );

        ras.par_chunks(DIRECTION_CHUNK_SIZE)
            .zip(decs.par_chunks(DIRECTION_CHUNK_SIZE))
            .enumerate()
            .try_for_each(|(i_chunk, (ras, decs))| {
                let out = &out;
                let mut err = [0 as c_char; ERR_LEN];
                let status = unsafe {
                    ffi::eb_point_responses(
                        self.telescope.0.as_ptr(),
                        time_mjd_s,
                        self.options.field_id,
                        ras.as_ptr(),
                        decs.as_ptr(),
                        ras.len(),
                        freqs_hz.as_ptr(),
                        freqs_hz.len(),
                        stations.as_ptr(),
                        stations.len(),
                        out.0.add(8 * i_chunk * DIRECTION_CHUNK_SIZE),
                        station_stride,
                        freq_stride,
                        err.as_mut_ptr(),
                        ERR_LEN,
                    )
                };
                if status == 0 {
                    Ok(())
                } else {
                    Err(BeamError::EveryBeam(err_to_string(&err)))
                }
            })?;

        if is_matrix_normalisation(self.options.beam_normalisation_mode.as_deref()) {
            results.mapv_inplace(swap_rows_and_columns);
        } else {
            results.mapv_inplace(swap_columns);
        }
        Ok(())
    }

    /// Like [`EveryBeam::calc_jones_radec`], but with [`AzEl`] directions.
    /// `time.lst_rad` and `latitude_rad` are used to convert the directions to
    /// (RA, Dec).
    fn calc_jones_azel(
        &self,
        azels: &[AzEl],
        freqs_hz: &[f64],
        stations: &[usize],
        latitude_rad: f64,
        time: Option<BeamTime>,
        results: ArrayViewMut3<Jones<f64>>,
    ) -> Result<(), BeamError> {
        let time = time.ok_or(BeamError::NeedsTime)?;
        let radecs: Vec<RADec> = azels
            .iter()
            .map(|azel| azel.to_hadec(latitude_rad).to_radec(time.lst_rad))
            .collect();
        self.calc_jones_radec(&radecs, freqs_hz, stations, time, results)
    }
}

impl Beam for EveryBeam {
    fn get_beam_type(&self) -> BeamType {
        BeamType::EveryBeam
    }

    fn get_num_tiles(&self) -> usize {
        self.inner.num_stations
    }

    fn get_dipole_delays(&self) -> Option<ArcArray<u32, Dim<[usize; 2]>>> {
        None
    }

    fn get_ideal_dipole_delays(&self) -> Option<[u32; 16]> {
        None
    }

    fn get_dipole_gains(&self) -> Option<ArcArray<f64, Dim<[usize; 2]>>> {
        None
    }

    fn get_beam_file(&self) -> Option<&Path> {
        Some(&self.inner.ms)
    }

    fn calc_jones(
        &self,
        azel: AzEl,
        freq_hz: f64,
        tile_index: Option<usize>,
        latitude_rad: f64,
        time: Option<BeamTime>,
    ) -> Result<Jones<f64>, BeamError> {
        let mut result = [Jones::default()];
        self.calc_jones_array_inner(
            &[azel],
            freq_hz,
            tile_index,
            latitude_rad,
            time,
            &mut result,
        )?;
        Ok(result[0])
    }

    fn calc_jones_array(
        &self,
        azels: &[AzEl],
        freq_hz: f64,
        tile_index: Option<usize>,
        latitude_rad: f64,
        time: Option<BeamTime>,
    ) -> Result<Vec<Jones<f64>>, BeamError> {
        let mut results = vec![Jones::default(); azels.len()];
        self.calc_jones_array_inner(
            azels,
            freq_hz,
            tile_index,
            latitude_rad,
            time,
            &mut results,
        )?;
        Ok(results)
    }

    fn calc_jones_array_inner(
        &self,
        azels: &[AzEl],
        freq_hz: f64,
        tile_index: Option<usize>,
        latitude_rad: f64,
        time: Option<BeamTime>,
        results: &mut [Jones<f64>],
    ) -> Result<(), BeamError> {
        // EveryBeam stations are not "ideal"; if no tile is specified, use the
        // first station.
        let station = tile_index.unwrap_or(0);
        let results = ArrayViewMut3::from_shape((1, 1, azels.len()), results)
            .expect("results has the same length as azels");
        self.inner.calc_jones_azel(
            azels,
            &[freq_hz],
            &[station],
            latitude_rad,
            time,
            results,
        )
    }

    fn calc_jones_tiles_freqs(
        &self,
        azels: &[AzEl],
        freqs_hz: &[f64],
        tile_indices: &[usize],
        latitude_rad: f64,
        time: Option<BeamTime>,
        results: ArrayViewMut3<Jones<f64>>,
    ) -> Result<(), BeamError> {
        self.inner
            .calc_jones_azel(azels, freqs_hz, tile_indices, latitude_rad, time, results)
    }

    fn find_closest_freq(&self, desired_freq_hz: f64) -> f64 {
        desired_freq_hz
    }

    fn empty_coeff_cache(&self) {}

    fn get_unique_tiles(
        &self,
        total_num_tiles: usize,
        flagged_tiles: &std::collections::HashSet<usize>,
    ) -> (Vec<usize>, Vec<usize>) {
        // Every station is assumed to have a unique response.
        let mut unique_tiles = Vec::with_capacity(total_num_tiles);
        let mut map = Vec::with_capacity(total_num_tiles);
        for i_tile in 0..total_num_tiles {
            if flagged_tiles.contains(&i_tile) {
                map.push(0);
            } else {
                map.push(unique_tiles.len());
                unique_tiles.push(i_tile);
            }
        }
        (unique_tiles, map)
    }

    #[cfg(any(feature = "cuda", feature = "hip"))]
    fn prepare_gpu_beam(&self, freqs_hz: &[u32]) -> Result<Box<dyn BeamGpu>, BeamError> {
        // All stations and frequencies are treated as unique.
        let tile_map: Vec<i32> = (0..self.inner.num_stations as i32).collect();
        let freq_map: Vec<i32> = (0..freqs_hz.len() as i32).collect();
        Ok(Box::new(EveryBeamGpu {
            beam: Arc::clone(&self.inner),
            freqs_hz: freqs_hz.iter().map(|&f| f as f64).collect(),
            tile_map: DevicePointer::copy_to_device(&tile_map)?,
            freq_map: DevicePointer::copy_to_device(&freq_map)?,
        }))
    }
}

/// An EveryBeam "GPU" beam. EveryBeam can't run on a GPU, so beam responses
/// are calculated on the CPU and copied to the device.
#[cfg(any(feature = "cuda", feature = "hip"))]
struct EveryBeamGpu {
    beam: Arc<EveryBeamInner>,
    freqs_hz: Vec<f64>,
    tile_map: DevicePointer<i32>,
    freq_map: DevicePointer<i32>,
}

#[cfg(any(feature = "cuda", feature = "hip"))]
impl BeamGpu for EveryBeamGpu {
    unsafe fn calc_jones_pair(
        &self,
        az_rad: &[GpuFloat],
        za_rad: &[GpuFloat],
        latitude_rad: f64,
        time: Option<BeamTime>,
        d_jones: *mut std::ffi::c_void,
    ) -> Result<(), BeamError> {
        #[cfg(feature = "cuda")]
        use cuda_runtime_sys::{
            cudaMemcpy as gpuMemcpy,
            cudaMemcpyKind::cudaMemcpyHostToDevice as gpuMemcpyHostToDevice,
        };
        #[cfg(feature = "hip")]
        use hip_sys::hiprt::{
            hipMemcpy as gpuMemcpy, hipMemcpyKind::hipMemcpyHostToDevice as gpuMemcpyHostToDevice,
        };

        let beam = &self.beam;
        let azels: Vec<AzEl> = az_rad
            .iter()
            .zip(za_rad.iter())
            .map(|(&az, &za)| {
                AzEl::from_radians(az as f64, std::f64::consts::FRAC_PI_2 - za as f64)
            })
            .collect();
        let stations: Vec<usize> = (0..beam.num_stations).collect();
        let mut results = Array3::default((stations.len(), self.freqs_hz.len(), azels.len()));
        beam.calc_jones_azel(
            &azels,
            &self.freqs_hz,
            &stations,
            latitude_rad,
            time,
            results.view_mut(),
        )?;

        // The results are ordered tile, frequency, direction, slowest to
        // fastest, matching what the GPU code expects.
        let results: Vec<Jones<GpuFloat>> = results
            .into_iter()
            .map(|j| {
                Jones::from([
                    num_complex::Complex::new(j[0].re as GpuFloat, j[0].im as GpuFloat),
                    num_complex::Complex::new(j[1].re as GpuFloat, j[1].im as GpuFloat),
                    num_complex::Complex::new(j[2].re as GpuFloat, j[2].im as GpuFloat),
                    num_complex::Complex::new(j[3].re as GpuFloat, j[3].im as GpuFloat),
                ])
            })
            .collect();
        gpuMemcpy(
            d_jones,
            results.as_ptr().cast(),
            results.len() * std::mem::size_of::<Jones<GpuFloat>>(),
            gpuMemcpyHostToDevice,
        );
        crate::gpu::check_for_errors(crate::gpu::GpuCall::CopyToDevice)?;
        Ok(())
    }

    fn get_beam_type(&self) -> BeamType {
        BeamType::EveryBeam
    }

    fn get_tile_map(&self) -> *const i32 {
        self.tile_map.get()
    }

    fn get_freq_map(&self) -> *const i32 {
        self.freq_map.get()
    }

    fn get_num_unique_tiles(&self) -> i32 {
        self.beam.num_stations as i32
    }

    fn get_num_unique_freqs(&self) -> i32 {
        self.freqs_hz.len() as i32
    }
}
