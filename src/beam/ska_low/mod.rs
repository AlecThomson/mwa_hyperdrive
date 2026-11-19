// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! SKA-Low station beams, calculated natively (on the CPU or GPU) with
//! hyperbeam's SKA-Low element models.
//!
//! | hyperdrive beam type | EveryBeam `element_response_model` |
//! |----------------------|------------------------------------|
//! | `ska-low-analytic`   | `oskar_dipole_cos`                 |
//! | `ska-low-eep-grid`   | `gridded_eep`                      |
//! | `ska-low-eep-sph`    | `skalow_feko`                      |
//!
//! As with EveryBeam, the telescope is described by a measurement set: each
//! station's element layout, element flags and rigid rotation come from the
//! `PHASED_ARRAY` table, and the stations are pointed at the `DELAY_DIR` of
//! the `FIELD` table (field 0).
//!
//! # Conventions
//!
//! The responses are in the same convention as the MWA FEE beam in hyperdrive
//! (hyperbeam's `calc_jones` with the parallactic-angle correction and
//! non-IAU order): the columns are the sky's (East, North) components and the
//! rows are each station's own X and Y feeds, in the station's (rotated)
//! frame. This is the same convention as hyperdrive's EveryBeam beam with its
//! default ("none") normalisation, and the responses are likewise
//! unnormalised (e.g. the EEP models are in their files' units); a normalised
//! copy (scaled by the element response at the zenith) is used for vetoing
//! sources.
//!
//! # Station frames
//!
//! Each station's response is calculated in its own local frame (from its
//! `COORDINATE_AXES`), so stations far from the array centre (whose horizons
//! differ from the centre's) are handled as EveryBeam handles them.
//! Directions given to the beam are relative to the array position derived
//! from the measurement set (its OBSERVATION table's ARRAY_CENTER, if
//! present, otherwise the mean ANTENNA position), as hyperdrive uses by
//! default.

#[cfg(test)]
mod tests;

use std::{
    collections::HashMap,
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
};

use log::debug;
use marlu::{rubbl_casatables, AzEl, Jones, RADec, XyzGeocentric};
use mwa_hyperbeam::ska::{
    rotation_from_p_axis, SkaBeamError, SkaElementResponse, SkaLowBeam as HbSkaLowBeam, SkaStation,
};
use ndarray::prelude::*;
use rubbl_casatables::{Table, TableOpenMode};

use super::{beam_knots_from_file_freqs, beam_knots_on_grid, Beam, BeamError, BeamTime, BeamType};

#[cfg(any(feature = "cuda", feature = "hip"))]
use super::{BeamGpu, DevicePointer, GpuFloat};

/// The number of guard knots used for the analytic SKA-Low model, which is
/// defined at any frequency (see [`Beam::get_beam_knots`]).
const ANALYTIC_NUM_GUARD_KNOTS: usize = 2;

/// hyperbeam's runtime-selected SKA-Low element model, shared between
/// stations.
#[derive(Clone)]
struct Model(Arc<HbSkaLowBeam>);

impl SkaElementResponse for Model {
    fn element_response(
        &self,
        element_id: usize,
        freq_hz: f64,
        theta_rad: f64,
        phi_rad: f64,
    ) -> Result<Jones<f64>, SkaBeamError> {
        self.0
            .element_response(element_id, freq_hz, theta_rad, phi_rad)
    }

    fn num_elements(&self) -> Option<usize> {
        self.0.num_elements()
    }
}

/// The geometry of a station, from a measurement set.
#[derive(Debug, Clone)]
pub(crate) struct StationGeometry {
    /// Element positions in the station's local East, North, Up frame
    /// \[metres\].
    pub(crate) positions_enu: Vec<[f64; 3]>,
    /// Weights of the X and Y elements (0 if flagged, otherwise 1).
    pub(crate) weights_x: Vec<f64>,
    pub(crate) weights_y: Vec<f64>,
    /// The angle of the station's X dipoles from its local East towards North
    /// \[radians\].
    pub(crate) rotation_rad: f64,
    /// The rows are the station's local East, North and Up unit vectors in
    /// the array's East, North, Up frame.
    pub(crate) local_frame: [[f64; 3]; 3],
}

/// The telescope layout needed for SKA-Low station beams.
#[derive(Debug, Clone)]
pub(crate) struct SkaLowLayout {
    pub(crate) stations: Vec<StationGeometry>,
    /// The direction the stations are pointed at (J2000).
    pub(crate) delay_dir: RADec,
}

type V3 = [f64; 3];

fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

fn normalise(a: V3) -> V3 {
    let n = norm(a);
    [a[0] / n, a[1] / n, a[2] / n]
}

/// The (East, North, Up) unit vectors (in ITRF) of a local frame whose
/// vertical is `up` (East is horizontal).
fn enu_from_up(up: V3) -> [V3; 3] {
    let up = normalise(up);
    let east = normalise([-up[1], up[0], 0.0]);
    let north = cross(up, east);
    [east, north, up]
}

/// The geodetic (East, North, Up) unit vectors (in ITRF) at a position.
fn geodetic_enu(xyz: XyzGeocentric) -> [V3; 3] {
    let pos = xyz.to_earth_wgs84();
    let (s_lat, c_lat) = pos.latitude_rad.sin_cos();
    let (s_lon, c_lon) = pos.longitude_rad.sin_cos();
    enu_from_up([c_lat * c_lon, c_lat * s_lon, s_lat])
}

fn ms_err(ms: &Path, what: impl std::fmt::Display) -> BeamError {
    BeamError::SkaLowMs(format!("{}: {what}", ms.display()))
}

fn open_table(ms: &Path, name: &str) -> Result<Table, BeamError> {
    Table::open(ms.join(name), TableOpenMode::Read)
        .map_err(|e| ms_err(ms, format!("couldn't open the {name} table: {e}")))
}

fn get_vec<T: rubbl_casatables::CasaScalarData>(
    ms: &Path,
    table: &mut Table,
    name: &str,
    column: &str,
    row: u64,
) -> Result<Vec<T>, BeamError> {
    table
        .get_cell_as_vec(column, row)
        .map_err(|e| ms_err(ms, format!("couldn't read {name} {column} row {row}: {e}")))
}

/// Read the station layout and pointing from a measurement set.
pub(crate) fn read_layout(ms: &Path) -> Result<SkaLowLayout, BeamError> {
    // The array reference position, as hyperdrive's measurement-set reader
    // determines it by default.
    let mut antenna = open_table(ms, "ANTENNA")?;
    let num_antennas = antenna.n_rows();
    let mut antenna_xyzs = Vec::with_capacity(num_antennas as usize);
    for row in 0..num_antennas {
        let p: Vec<f64> = get_vec(ms, &mut antenna, "ANTENNA", "POSITION", row)?;
        if p.len() != 3 {
            return Err(ms_err(ms, "ANTENNA POSITION doesn't have 3 values"));
        }
        antenna_xyzs.push([p[0], p[1], p[2]]);
    }
    let mut observation = open_table(ms, "OBSERVATION")?;
    let array_xyz = match observation.get_cell_as_vec::<f64>("ARRAY_CENTER", 0) {
        Ok(v) if v.len() >= 3 => [v[0], v[1], v[2]],
        _ => {
            if antenna_xyzs.is_empty() {
                return Err(ms_err(ms, "the ANTENNA table is empty"));
            }
            let n = antenna_xyzs.len() as f64;
            let mut mean = [0.0; 3];
            for p in &antenna_xyzs {
                for k in 0..3 {
                    mean[k] += p[k] / n;
                }
            }
            mean
        }
    };
    let array_enu = geodetic_enu(XyzGeocentric {
        x: array_xyz[0],
        y: array_xyz[1],
        z: array_xyz[2],
    });

    let mut field = open_table(ms, "FIELD")?;
    let delay_dir: Vec<f64> = get_vec(ms, &mut field, "FIELD", "DELAY_DIR", 0)?;
    if delay_dir.len() < 2 {
        return Err(ms_err(ms, "FIELD DELAY_DIR doesn't have 2 values"));
    }
    let delay_dir = RADec::from_radians(delay_dir[0], delay_dir[1]);

    let mut pa = open_table(ms, "PHASED_ARRAY").map_err(|e| {
        BeamError::SkaLowMs(format!(
            "{e}; SKA-Low beams need a PHASED_ARRAY table describing each station's elements"
        ))
    })?;
    let num_stations = pa.n_rows();
    if num_stations != num_antennas {
        return Err(ms_err(
            ms,
            format!("PHASED_ARRAY has {num_stations} rows, but ANTENNA has {num_antennas}"),
        ));
    }
    let mut stations = Vec::with_capacity(num_stations as usize);
    for row in 0..num_stations {
        let position: Vec<f64> = get_vec(ms, &mut pa, "PHASED_ARRAY", "POSITION", row)?;
        let axes: Vec<f64> = get_vec(ms, &mut pa, "PHASED_ARRAY", "COORDINATE_AXES", row)?;
        let offsets: Vec<f64> = get_vec(ms, &mut pa, "PHASED_ARRAY", "ELEMENT_OFFSET", row)?;
        if position.len() != 3
            || axes.len() != 9
            || !offsets.len().is_multiple_of(3)
            || offsets.is_empty()
        {
            return Err(ms_err(
                ms,
                format!("PHASED_ARRAY row {row} has unexpected POSITION, COORDINATE_AXES or ELEMENT_OFFSET shapes"),
            ));
        }
        let num_elements = offsets.len() / 3;
        let flags: Vec<bool> = match pa.get_cell_as_vec::<bool>("ELEMENT_FLAG", row) {
            Ok(f) if f.len() == 2 * num_elements => f,
            Ok(f) => {
                return Err(ms_err(
                    ms,
                    format!(
                    "PHASED_ARRAY row {row} has {} ELEMENT_FLAG values for {num_elements} elements",
                    f.len()
                ),
                ))
            }
            Err(_) => vec![false; 2 * num_elements],
        };

        // casacore arrays are column major, so the flattened COORDINATE_AXES
        // are the p, q and r axes in turn (each in ITRF), the element offsets
        // are (x, y, z) per element, and the flags are (X, Y) per element.
        let p = [axes[0], axes[1], axes[2]];
        let r = [axes[6], axes[7], axes[8]];
        let station_enu = if norm(r) > 0.5 {
            enu_from_up(r)
        } else {
            geodetic_enu(XyzGeocentric {
                x: position[0],
                y: position[1],
                z: position[2],
            })
        };
        let to_station = |v: V3| {
            [
                dot(v, station_enu[0]),
                dot(v, station_enu[1]),
                dot(v, station_enu[2]),
            ]
        };
        let positions_enu = offsets
            .as_chunks::<3>()
            .0
            .iter()
            .map(|&o| to_station(o))
            .collect();
        let (weights_x, weights_y) = flags
            .as_chunks::<2>()
            .0
            .iter()
            .map(|f| (if f[0] { 0.0 } else { 1.0 }, if f[1] { 0.0 } else { 1.0 }))
            .unzip();
        let mut local_frame = [[0.0; 3]; 3];
        for (i, s) in station_enu.iter().enumerate() {
            for (j, a) in array_enu.iter().enumerate() {
                local_frame[i][j] = dot(*s, *a);
            }
        }
        stations.push(StationGeometry {
            positions_enu,
            weights_x,
            weights_y,
            rotation_rad: rotation_from_p_axis(to_station(p)),
            local_frame,
        });
    }
    debug!(
        "Read {} SKA-Low stations from {}; delay direction {delay_dir}",
        stations.len(),
        ms.display()
    );

    Ok(SkaLowLayout {
        stations,
        delay_dir,
    })
}

/// One station: its beam with the X element weights, and (if the Y elements
/// have different weights) its beam with the Y element weights.
struct Station {
    x: SkaStation<Model>,
    y: Option<SkaStation<Model>>,
}

struct Inner {
    beam_type: BeamType,
    model: Model,
    layout: SkaLowLayout,
    stations: Vec<Station>,
    ms: PathBuf,
    /// Normalisation factors (see [`Inner::zenith_norm`]), keyed by the bits
    /// of the frequency.
    norms: Mutex<HashMap<u64, f64>>,
}

/// SKA-Low station beams from hyperbeam, implementing the [`Beam`] trait.
pub(crate) struct SkaLowBeam {
    inner: Arc<Inner>,
    /// Whether responses are normalised by the element response at the
    /// zenith.
    normalised: bool,
    /// A normalised copy, used for vetoing (if this beam isn't normalised).
    veto: Option<Box<SkaLowBeam>>,
}

impl SkaLowBeam {
    /// Set up an SKA-Low beam of `beam_type` for the telescope described by
    /// the measurement set `ms`. `data_path` is the model's data (the EEP grid
    /// HDF5 file for `ska-low-eep-grid`, or the directory of FEKO coefficient
    /// files for `ska-low-eep-sph`); if not given, hyperbeam looks for it in
    /// the environment (`SKA_LOW_EEP_GRID_FILE`, `SKA_LOW_EEP_SPH_PATH` or
    /// `EVERYBEAM_DATADIR`). If `num_tiles` is given, it must match the number
    /// of stations.
    pub(crate) fn new(
        beam_type: BeamType,
        ms: &Path,
        data_path: Option<&Path>,
        num_tiles: Option<usize>,
    ) -> Result<SkaLowBeam, BeamError> {
        let model = match beam_type {
            BeamType::SkaLowAnalytic => HbSkaLowBeam::new_analytic(true),
            BeamType::SkaLowEepGrid => HbSkaLowBeam::new_eep_grid(data_path, true)?,
            BeamType::SkaLowEepSph => HbSkaLowBeam::new_eep_sph(data_path, None)?,
            _ => panic!("{beam_type} is not an SKA-Low beam type"),
        };
        let layout = read_layout(ms)?;
        Self::from_layout(beam_type, model, layout, ms, num_tiles)
    }

    pub(crate) fn from_layout(
        beam_type: BeamType,
        model: HbSkaLowBeam,
        layout: SkaLowLayout,
        ms: &Path,
        num_tiles: Option<usize>,
    ) -> Result<SkaLowBeam, BeamError> {
        if let Some(num_tiles) = num_tiles {
            if num_tiles != layout.stations.len() {
                return Err(BeamError::StationCountMismatch {
                    ms: layout.stations.len(),
                    tiles: num_tiles,
                });
            }
        }
        let model = Model(Arc::new(model));
        let stations = layout
            .stations
            .iter()
            .map(|g| {
                let make = |weights: &[f64]| -> Result<SkaStation<Model>, SkaBeamError> {
                    Ok(SkaStation::new(
                        model.clone(),
                        g.positions_enu.clone(),
                        Some(weights.to_vec()),
                    )?
                    .with_rotation(g.rotation_rad)
                    .with_local_frame(g.local_frame))
                };
                let x = make(&g.weights_x)?;
                let y = if g.weights_y != g.weights_x {
                    Some(make(&g.weights_y)?)
                } else {
                    None
                };
                Ok(Station { x, y })
            })
            .collect::<Result<Vec<_>, SkaBeamError>>()?;
        let inner = Arc::new(Inner {
            beam_type,
            model,
            layout,
            stations,
            ms: ms.to_path_buf(),
            norms: Mutex::new(HashMap::new()),
        });
        Ok(SkaLowBeam {
            inner: Arc::clone(&inner),
            normalised: false,
            veto: Some(Box::new(SkaLowBeam {
                inner,
                normalised: true,
                veto: None,
            })),
        })
    }
}

impl Inner {
    /// The pointing direction at a time, in the same frame as hyperdrive's
    /// [`AzEl`]s (i.e. from the same LST and latitude).
    fn pointing(&self, latitude_rad: f64, time: BeamTime) -> AzEl {
        self.layout
            .delay_dir
            .to_hadec(time.lst_rad)
            .to_azel(latitude_rad)
    }

    /// The scalar amplitude of the element response at the zenith (averaged
    /// over elements, for per-element models), used to normalise responses.
    fn zenith_norm(&self, freq_hz: f64) -> Result<f64, BeamError> {
        let key = freq_hz.to_bits();
        if let Some(n) = self.norms.lock().expect("not poisoned").get(&key) {
            return Ok(*n);
        }
        let num_elements = self.model.num_elements().unwrap_or(1);
        let mut sum = 0.0;
        for e in 0..num_elements {
            let j = self.model.element_response(e, freq_hz, 0.0, 0.0)?;
            sum += j.iter().map(|c| c.norm_sqr()).sum::<f64>() / 2.0;
        }
        let n = (sum / num_elements as f64).sqrt();
        self.norms.lock().expect("not poisoned").insert(key, n);
        Ok(n)
    }

    fn station(&self, tile_index: Option<usize>) -> Result<&Station, BeamError> {
        let i = tile_index.unwrap_or(0);
        self.stations.get(i).ok_or(BeamError::BadTileIndex {
            got: i,
            num_tiles: self.stations.len(),
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn calc(
        &self,
        azels: &[AzEl],
        freq_hz: f64,
        tile_index: Option<usize>,
        latitude_rad: f64,
        time: Option<BeamTime>,
        normalised: bool,
        results: &mut [Jones<f64>],
    ) -> Result<(), BeamError> {
        assert_eq!(azels.len(), results.len());
        let time = time.ok_or(BeamError::NeedsTime)?;
        let station = self.station(tile_index)?;
        let pointing = self.pointing(latitude_rad, time);
        let (az, za): (Vec<f64>, Vec<f64>) = azels.iter().map(|a| (a.az, a.za())).unzip();
        let calc = |s: &SkaStation<Model>| {
            s.calc_jones_array_pair(
                &az,
                &za,
                freq_hz,
                pointing.az,
                pointing.za(),
                Some(latitude_rad),
                false,
            )
        };
        let x = calc(&station.x)?;
        match &station.y {
            None => results.copy_from_slice(&x),
            Some(y) => {
                // The rows are the X and Y feeds (the parallactic-angle
                // correction only mixes the columns).
                let y = calc(y)?;
                for ((r, x), y) in results.iter_mut().zip(x).zip(y) {
                    *r = Jones::from([x[0], x[1], y[2], y[3]]);
                }
            }
        }
        if normalised {
            let scale = 1.0 / self.zenith_norm(freq_hz)?;
            for r in results.iter_mut() {
                *r *= scale;
            }
        }
        Ok(())
    }
}

impl Beam for SkaLowBeam {
    fn get_beam_type(&self) -> BeamType {
        self.inner.beam_type
    }

    fn get_num_tiles(&self) -> usize {
        self.inner.stations.len()
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

    fn veto_beam(&self) -> Option<&dyn Beam> {
        self.veto.as_deref().map(|b| b as &dyn Beam)
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
        self.calc_jones_array_inner(azels, freq_hz, tile_index, latitude_rad, time, &mut results)?;
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
        self.inner.calc(
            azels,
            freq_hz,
            tile_index,
            latitude_rad,
            time,
            self.normalised,
            results,
        )
    }

    fn find_closest_freq(&self, desired_freq_hz: f64) -> f64 {
        self.inner.model.0.find_closest_freq(desired_freq_hz)
    }

    fn get_beam_knots(&self, freqs_hz: &[f64]) -> Vec<f64> {
        match self.inner.model.0.get_freqs() {
            // Defined at any frequency: knots on the coarse-channel grid.
            None => beam_knots_on_grid(freqs_hz, ANALYTIC_NUM_GUARD_KNOTS),
            // Defined at the frequencies of its data.
            Some(beam_freqs) => beam_knots_from_file_freqs(&beam_freqs, freqs_hz),
        }
    }

    fn empty_coeff_cache(&self) {
        self.inner.model.0.empty_cache();
    }

    fn get_unique_tiles(
        &self,
        total_num_tiles: usize,
        flagged_tiles: &std::collections::HashSet<usize>,
    ) -> (Vec<usize>, Vec<usize>) {
        // SKA-Low stations are rotated differently, so every station is
        // treated as unique.
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
        Ok(Box::new(SkaLowBeamGpu::new(&self.inner, freqs_hz)?))
    }
}

/// SKA-Low station beams on a GPU, using hyperbeam's GPU code. All stations
/// and frequencies are treated as unique.
#[cfg(any(feature = "cuda", feature = "hip"))]
struct SkaLowBeamGpu {
    inner: Arc<Inner>,
    /// The stations with their X element weights.
    gpu_x: mwa_hyperbeam::ska::SkaLowBeamGpu,
    /// The stations with their Y element weights, if any differ from the X
    /// weights.
    gpu_y: Option<mwa_hyperbeam::ska::SkaLowBeamGpu>,
    num_freqs: usize,
    tile_map: DevicePointer<i32>,
    freq_map: DevicePointer<i32>,
}

#[cfg(any(feature = "cuda", feature = "hip"))]
impl SkaLowBeamGpu {
    fn new(inner: &Arc<Inner>, freqs_hz: &[u32]) -> Result<SkaLowBeamGpu, BeamError> {
        use mwa_hyperbeam::ska::SkaGpuTiles;

        let stations = &inner.layout.stations;
        let num_stations = stations.len();
        // Stations with fewer elements are padded with zero-weight elements.
        let max_elements = stations
            .iter()
            .map(|s| s.positions_enu.len())
            .max()
            .unwrap_or(0);
        let mut positions = Array3::zeros((num_stations, max_elements, 3));
        let mut weights_x = Array2::zeros((num_stations, max_elements));
        let mut weights_y = Array2::zeros((num_stations, max_elements));
        for (i, s) in stations.iter().enumerate() {
            for (e, p) in s.positions_enu.iter().enumerate() {
                for k in 0..3 {
                    positions[(i, e, k)] = p[k];
                }
                weights_x[(i, e)] = s.weights_x[e];
                weights_y[(i, e)] = s.weights_y[e];
            }
        }
        let rotations: Vec<f64> = stations.iter().map(|s| s.rotation_rad).collect();
        let frames: Vec<[[f64; 3]; 3]> = stations.iter().map(|s| s.local_frame).collect();
        let freqs: Vec<f64> = freqs_hz.iter().map(|&f| f as f64).collect();
        let prepare = |weights: &Array2<f64>| unsafe {
            inner.model.0.gpu_prepare(
                &freqs,
                SkaGpuTiles::Stations {
                    positions_enu: positions.view(),
                    weights: Some(weights.view()),
                    rotations_rad: Some(&rotations),
                    local_frames: Some(&frames),
                },
            )
        };
        let gpu_x = prepare(&weights_x)?;
        let gpu_y = if weights_y != weights_x {
            Some(prepare(&weights_y)?)
        } else {
            None
        };
        let tile_map: Vec<i32> = (0..num_stations as i32).collect();
        let freq_map: Vec<i32> = (0..freqs.len() as i32).collect();
        Ok(SkaLowBeamGpu {
            inner: Arc::clone(inner),
            gpu_x,
            gpu_y,
            num_freqs: freqs.len(),
            tile_map: DevicePointer::copy_to_device(&tile_map)?,
            freq_map: DevicePointer::copy_to_device(&freq_map)?,
        })
    }
}

#[cfg(any(feature = "cuda", feature = "hip"))]
impl BeamGpu for SkaLowBeamGpu {
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
            cudaMemcpy2D as gpuMemcpy2D,
            cudaMemcpyKind::cudaMemcpyDeviceToDevice as gpuMemcpyDeviceToDevice,
        };
        #[cfg(feature = "hip")]
        use hip_sys::hiprt::{
            hipMemcpy2D as gpuMemcpy2D,
            hipMemcpyKind::hipMemcpyDeviceToDevice as gpuMemcpyDeviceToDevice,
        };

        let time = time.ok_or(BeamError::NeedsTime)?;
        let pointing = self.inner.pointing(latitude_rad, time);
        let num_directions: i32 = az_rad.len().try_into().expect("not bigger than i32::MAX");
        let d_az = DevicePointer::copy_to_device(az_rad)?;
        let d_za = DevicePointer::copy_to_device(za_rad)?;
        let d_latitude = DevicePointer::copy_to_device(&[latitude_rad as GpuFloat])?;
        let calc = |gpu: &mwa_hyperbeam::ska::SkaLowBeamGpu, d_out: *mut std::ffi::c_void| {
            gpu.calc_jones_device_pair_inner(
                d_az.get(),
                d_za.get(),
                num_directions,
                pointing.az as GpuFloat,
                pointing.za() as GpuFloat,
                d_latitude.get(),
                false,
                d_out,
            )
        };
        calc(&self.gpu_x, d_jones)?;
        if let Some(gpu_y) = &self.gpu_y {
            // Calculate the responses with the Y weights, and copy their Y
            // rows (the second half of each Jones matrix) into the results.
            let n = self.inner.stations.len() * self.num_freqs * az_rad.len();
            let jones_size = std::mem::size_of::<Jones<GpuFloat>>();
            let mut d_y = DevicePointer::<Jones<GpuFloat>>::malloc(n * jones_size)?;
            calc(gpu_y, d_y.get_mut().cast())?;
            let half = jones_size / 2;
            gpuMemcpy2D(
                d_jones.cast::<u8>().add(half).cast(),
                jones_size,
                d_y.get().cast::<u8>().add(half).cast(),
                jones_size,
                half,
                n,
                gpuMemcpyDeviceToDevice,
            );
            crate::gpu::check_for_errors(crate::gpu::GpuCall::CopyToDevice)?;
        }
        Ok(())
    }

    fn get_beam_type(&self) -> BeamType {
        self.inner.beam_type
    }

    fn get_tile_map(&self) -> *const i32 {
        self.tile_map.get()
    }

    fn get_freq_map(&self) -> *const i32 {
        self.freq_map.get()
    }

    fn get_num_unique_tiles(&self) -> i32 {
        self.inner.stations.len() as i32
    }

    fn get_num_unique_freqs(&self) -> i32 {
        self.num_freqs as i32
    }
}
