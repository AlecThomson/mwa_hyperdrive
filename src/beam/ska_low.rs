// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Code for SKA-Low station beam calculations.

use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

use log::debug;
use marlu::{AzEl, Jones};
use ndarray::prelude::*;
use num_complex::Complex64 as c64;

use super::{Beam, BeamError, BeamType};
use crate::context::{PhasedArray, PolConvention};

#[cfg(any(feature = "cuda", feature = "hip"))]
use super::BeamGpu;

/// Where the per-element patterns come from.
#[derive(Debug, Clone)]
pub(crate) enum SkaLowSource {
    /// FEKO spherical-wave coefficients, one HDF5 file per frequency.
    Swe { use_ticra_convention: bool },

    /// Gridded embedded element patterns.
    Grid {
        filebase: String,
        format: mwa_hyperbeam::ska_low::GridFormat,
        normalise: bool,
    },
}

/// A wrapper of hyperbeam's `SkaLowBeam` that implements the [`Beam`] trait.
///
/// hyperbeam's object holds one station's element layout, so we keep one per
/// *distinct* layout and map tiles onto them. Each object owns a coefficient
/// cache, so sharing them matters when every station reads the same files.
pub(crate) struct SkaLowBeam {
    stations: Vec<mwa_hyperbeam::ska_low::SkaLowBeam>,

    /// Maps a tile index to an index into `stations`.
    station_map: Vec<usize>,

    /// Per-element excitations, one entry per tile. Dead elements are zero.
    weights: Vec<Vec<c64>>,

    num_tiles: usize,
    dir: PathBuf,
}

impl SkaLowBeam {
    pub(crate) fn new(
        dir: Option<&Path>,
        num_tiles: usize,
        phased_array: &PhasedArray,
        source: &SkaLowSource,
    ) -> Result<SkaLowBeam, BeamError> {
        let dir = match dir {
            Some(d) => d.to_path_buf(),
            None => PathBuf::from(
                std::env::var("SKA_LOW_BEAM_DIR").map_err(|_| BeamError::NoSkaLowBeamDir)?,
            ),
        };
        debug!("SKA-Low beam directory: {}", dir.display());

        let num_stations = phased_array.element_offsets.len();
        if num_stations != num_tiles {
            return Err(BeamError::InconsistentStations {
                num_stations,
                num_tiles,
            });
        }

        let mut stations = Vec::new();
        let mut station_map = Vec::with_capacity(num_tiles);
        let mut weights = Vec::with_capacity(num_tiles);
        // ponytail: dedupe on an exact bit-for-bit layout match. Stations that
        // differ only by rounding get their own object and their own cache.
        let mut seen: HashMap<Vec<u64>, usize> = HashMap::new();

        for i in 0..num_tiles {
            let offsets = &phased_array.element_offsets[i];
            let axes = phased_array.coordinate_axes.as_ref().map(|a| a[i]);

            let mut key: Vec<u64> = offsets.iter().map(|v| v.to_bits()).collect();
            if let Some(a) = axes.as_ref() {
                key.extend(a.iter().flatten().map(|v| v.to_bits()));
            }

            let index = match seen.get(&key) {
                Some(&index) => index,
                None => {
                    let beam = match source {
                        SkaLowSource::Swe {
                            use_ticra_convention,
                        } => mwa_hyperbeam::ska_low::SkaLowBeam::new(
                            &dir,
                            offsets.view(),
                            axes,
                            *use_ticra_convention,
                        )?,
                        SkaLowSource::Grid {
                            filebase,
                            format,
                            normalise,
                        } => mwa_hyperbeam::ska_low::SkaLowBeam::new_from_grid(
                            &dir,
                            filebase,
                            *format,
                            offsets.view(),
                            axes,
                            *normalise,
                        )?,
                    };
                    stations.push(beam);
                    seen.insert(key, stations.len() - 1);
                    stations.len() - 1
                }
            };
            station_map.push(index);
            weights.push(element_weights(phased_array, i));
        }

        debug!(
            "{num_tiles} stations use {} distinct element layout(s)",
            stations.len()
        );

        Ok(SkaLowBeam {
            stations,
            station_map,
            weights,
            num_tiles,
            dir,
        })
    }

    /// Resolve a tile index to its hyperbeam object and element weights.
    /// Without an index we use the first station, mirroring how the FEE beam
    /// falls back to ideal delays.
    fn station(
        &self,
        tile_index: Option<usize>,
    ) -> Result<(&mwa_hyperbeam::ska_low::SkaLowBeam, &[c64]), BeamError> {
        let i = tile_index.unwrap_or(0);
        if i >= self.num_tiles {
            return Err(BeamError::BadTileIndex {
                got: i,
                max: self.num_tiles.saturating_sub(1),
            });
        }
        Ok((&self.stations[self.station_map[i]], &self.weights[i]))
    }
}

/// Fold a station's element flags into hyperbeam's complex weights. hyperbeam
/// has no flag concept: a dead element is one with a weight of zero.
fn element_weights(phased_array: &PhasedArray, tile_index: usize) -> Vec<c64> {
    let num_elements = phased_array.element_offsets[tile_index].len_of(Axis(1));
    match &phased_array.element_flags {
        None => vec![c64::new(1.0, 0.0); num_elements],
        Some(flags) => {
            let f = &flags[tile_index];
            (0..num_elements)
                .map(|e| {
                    // The MS flags X and Y separately, but hyperbeam takes one
                    // weight per element, so either one killing it is fatal.
                    if f[[0, e]] || f[[1, e]] {
                        c64::new(0.0, 0.0)
                    } else {
                        c64::new(1.0, 0.0)
                    }
                })
                .collect()
        }
    }
}

impl Beam for SkaLowBeam {
    fn get_beam_type(&self) -> BeamType {
        BeamType::SkaLow
    }

    fn get_pol_convention(&self) -> Option<PolConvention> {
        // We always ask hyperbeam for IAU-ordered Jones matrices below; SKA-Low
        // is an IAU instrument.
        Some(PolConvention::Iau)
    }

    fn get_num_tiles(&self) -> usize {
        self.num_tiles
    }

    fn get_ideal_dipole_delays(&self) -> Option<[u32; 16]> {
        None
    }

    fn get_dipole_delays(&self) -> Option<ArcArray<u32, Dim<[usize; 2]>>> {
        None
    }

    fn get_dipole_gains(&self) -> Option<ArcArray<f64, Dim<[usize; 2]>>> {
        None
    }

    fn get_beam_file(&self) -> Option<&Path> {
        Some(&self.dir)
    }

    fn calc_jones(
        &self,
        azel: AzEl,
        freq_hz: f64,
        tile_index: Option<usize>,
        latitude_rad: f64,
    ) -> Result<Jones<f64>, BeamError> {
        let (station, weights) = self.station(tile_index)?;
        Ok(station.calc_jones(azel, freq_hz, weights, Some(latitude_rad), true)?)
    }

    fn calc_jones_array(
        &self,
        azels: &[AzEl],
        freq_hz: f64,
        tile_index: Option<usize>,
        latitude_rad: f64,
    ) -> Result<Vec<Jones<f64>>, BeamError> {
        let (station, weights) = self.station(tile_index)?;
        Ok(station.calc_jones_array(azels, freq_hz, weights, Some(latitude_rad), true)?)
    }

    fn calc_jones_array_inner(
        &self,
        azels: &[AzEl],
        freq_hz: f64,
        tile_index: Option<usize>,
        latitude_rad: f64,
        results: &mut [Jones<f64>],
    ) -> Result<(), BeamError> {
        let (station, weights) = self.station(tile_index)?;
        station.calc_jones_array_inner(
            azels,
            freq_hz,
            weights,
            Some(latitude_rad),
            true,
            results,
        )?;
        Ok(())
    }

    fn find_closest_freq(&self, desired_freq_hz: f64) -> f64 {
        self.stations[0].find_closest_freq(desired_freq_hz)
    }

    fn empty_coeff_cache(&self) {
        self.stations.iter().for_each(|s| s.empty_cache());
    }

    #[cfg(any(feature = "cuda", feature = "hip"))]
    fn prepare_gpu_beam(&self, _freqs_hz: &[u32]) -> Result<Box<dyn BeamGpu>, BeamError> {
        Err(BeamError::SkaLowNoGpu)
    }
}
