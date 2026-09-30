// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Tests for EveryBeam beam responses. These use a small synthetic SKA-Low
//! measurement set (see `test_files/everybeam/make_skalow_ms.py`), whose
//! stations are rigidly rotated by 0, 30 and 75 degrees, and whose third
//! station has some flagged X elements.

use std::collections::HashSet;

use approx::assert_abs_diff_eq;
use hifitime::Epoch;
use marlu::{c64, AzEl, Jones, RADec};
use ndarray::prelude::*;

use super::{reference::*, *};

const MS: &str = "test_files/everybeam/skalow_mini.ms";

fn get_beam(normalisation: &str) -> EveryBeam {
    EveryBeam::new(
        Path::new(MS),
        Some(3),
        EveryBeamOptions {
            beam_normalisation_mode: Some(normalisation.to_string()),
            ..Default::default()
        },
    )
    .unwrap()
}

fn ref_epoch() -> Epoch {
    Epoch::from_mjd_utc(REF_TIME_MJD_S / 86400.0)
}

fn ref_radecs() -> Vec<RADec> {
    REF_RADECS
        .iter()
        .map(|&(ra, dec)| RADec::from_radians(ra, dec))
        .collect()
}

/// Convert a reference Jones matrix (EveryBeam's convention) to hyperdrive's
/// convention (swap the columns).
fn ref_to_hyperdrive(r: &[[f64; 2]; 4]) -> Jones<f64> {
    let c = |i: usize| c64::new(r[i][0], r[i][1]);
    Jones::from([c(1), c(0), c(3), c(2)])
}

fn check_against_reference(normalisation: &str, reference: &[[[[[f64; 2]; 4]; 3]; 2]; 3]) {
    let beam = get_beam(normalisation);
    let radecs = ref_radecs();
    let stations = [0, 1, 2];
    let mut results = Array3::default((stations.len(), REF_FREQS.len(), radecs.len()));
    beam.inner
        .calc_jones_radec(
            &radecs,
            &REF_FREQS,
            &stations,
            BeamTime {
                epoch: ref_epoch(),
                lst_rad: 0.0,
            },
            results.view_mut(),
        )
        .unwrap();

    for (i_station, station_ref) in reference.iter().enumerate() {
        for (i_freq, freq_ref) in station_ref.iter().enumerate() {
            for (i_dir, dir_ref) in freq_ref.iter().enumerate() {
                let expected = ref_to_hyperdrive(dir_ref);
                let got = results[(i_station, i_freq, i_dir)];
                // EveryBeam calculates in single precision.
                assert_abs_diff_eq!(got, expected, epsilon = 1e-6);
            }
        }
    }
}

#[test]
fn test_load_and_station_count() {
    let beam = get_beam("amplitude");
    assert_eq!(beam.get_num_tiles(), 3);
    assert!(matches!(beam.get_beam_type(), BeamType::EveryBeam));
    assert_eq!(beam.get_beam_file(), Some(Path::new(MS)));
    assert!(beam.get_dipole_delays().is_none());
    assert!(beam.get_dipole_gains().is_none());
    assert_eq!(beam.find_closest_freq(123.456e6), 123.456e6);
}

#[test]
fn test_station_count_mismatch() {
    let result = EveryBeam::new(Path::new(MS), Some(4), EveryBeamOptions::default());
    assert!(matches!(
        result,
        Err(BeamError::StationCountMismatch { ms: 3, tiles: 4 })
    ));
}

#[test]
fn test_bad_ms() {
    let result = EveryBeam::new(
        Path::new("test_files/everybeam/does_not_exist.ms"),
        None,
        EveryBeamOptions::default(),
    );
    assert!(matches!(result, Err(BeamError::EveryBeam(_))));
}

#[test]
fn test_bad_element_model() {
    let result = EveryBeam::new(
        Path::new(MS),
        None,
        EveryBeamOptions {
            element_response_model: Some("not_a_model".to_string()),
            ..Default::default()
        },
    );
    assert!(matches!(result, Err(BeamError::EveryBeam(_))));
}

#[test]
fn test_against_everybeam_python_no_normalisation() {
    check_against_reference("none", &REF_NONE);
}

#[test]
fn test_against_everybeam_python_amplitude_normalisation() {
    check_against_reference("amplitude", &REF_AMPLITUDE);
}

#[test]
fn test_rotated_stations() {
    // At the beam centre (zenith), a station's X dipole responds to the sky's
    // East and North components with cos(rotation) and sin(rotation), because
    // the stations are rigidly rotated. The rows are in each station's own
    // frame, so they aren't de-rotated.
    let beam = get_beam("amplitude");
    let radecs = ref_radecs();
    let mut results = Array3::default((3, 1, 1));
    beam.inner
        .calc_jones_radec(
            &radecs[..1],
            &[REF_FREQS[0]],
            &[0, 1, 2],
            BeamTime {
                epoch: ref_epoch(),
                lst_rad: 0.0,
            },
            results.view_mut(),
        )
        .unwrap();
    for (i_station, rotation_deg) in [0.0_f64, 30.0, 75.0].into_iter().enumerate() {
        let j = results[(i_station, 0, 0)];
        let (s, c) = rotation_deg.to_radians().sin_cos();
        // hyperdrive's columns are (East, North). Allow for the zenith not
        // being exactly at the beam centre after 10 seconds, and for the sign
        // convention of the X dipole.
        assert_abs_diff_eq!(j[0].norm(), c, epsilon = 5e-3);
        assert_abs_diff_eq!(j[1].norm(), s, epsilon = 5e-3);
        assert_abs_diff_eq!(j[2].norm(), s, epsilon = 5e-3);
        assert_abs_diff_eq!(j[3].norm(), c, epsilon = 5e-3);
    }
}

#[test]
fn test_azel_matches_radec() {
    let beam = get_beam("amplitude");
    let radecs = ref_radecs();
    let epoch = ref_epoch();
    // Any LST and latitude will do, as long as the same values are used to
    // generate the AzEls and supplied as the beam time.
    let lst_rad = 1.2;
    let latitude_rad = -0.4682;
    let azels: Vec<AzEl> = radecs
        .iter()
        .map(|r| r.to_hadec(lst_rad).to_azel(latitude_rad))
        .collect();
    let time = BeamTime { epoch, lst_rad };

    let mut expected = Array3::default((1, 1, radecs.len()));
    beam.inner
        .calc_jones_radec(&radecs, &[110e6], &[1], time, expected.view_mut())
        .unwrap();
    let got = beam
        .calc_jones_array(&azels, 110e6, Some(1), latitude_rad, Some(time))
        .unwrap();
    for (g, e) in got.iter().zip(expected.iter()) {
        assert_abs_diff_eq!(*g, *e, epsilon = 1e-6);
    }

    let single = beam
        .calc_jones(azels[1], 110e6, Some(1), latitude_rad, Some(time))
        .unwrap();
    assert_abs_diff_eq!(single, expected[(0, 0, 1)], epsilon = 1e-6);

    // No tile index means the first station.
    let first = beam
        .calc_jones(azels[1], 110e6, None, latitude_rad, Some(time))
        .unwrap();
    let first_explicit = beam
        .calc_jones(azels[1], 110e6, Some(0), latitude_rad, Some(time))
        .unwrap();
    assert_abs_diff_eq!(first, first_explicit);
}

#[test]
fn test_needs_time() {
    let beam = get_beam("amplitude");
    let result = beam.calc_jones(AzEl::from_degrees(0.0, 80.0), 110e6, Some(0), -0.4682, None);
    assert!(matches!(result, Err(BeamError::NeedsTime)));
}

#[test]
fn test_bad_station_index() {
    let beam = get_beam("amplitude");
    let time = BeamTime {
        epoch: ref_epoch(),
        lst_rad: 0.0,
    };
    let result = beam.calc_jones(
        AzEl::from_degrees(0.0, 80.0),
        110e6,
        Some(3),
        -0.4682,
        Some(time),
    );
    assert!(matches!(
        result,
        Err(BeamError::BadTileIndex { got: 3, max: 2 })
    ));
}

#[test]
fn test_many_directions_tiles_freqs() {
    // Use more directions than a single chunk, so that the multi-threaded code
    // is exercised, and check against direction-by-direction results.
    let beam = get_beam("amplitude");
    let epoch = ref_epoch();
    let lst_rad = 1.2;
    let latitude_rad = -0.4682;
    let time = BeamTime { epoch, lst_rad };
    let azels: Vec<AzEl> = (0..3 * DIRECTION_CHUNK_SIZE + 7)
        .map(|i| AzEl::from_degrees(i as f64 * 1.7 % 360.0, 30.0 + (i % 60) as f64))
        .collect();
    let freqs = [100e6, 150e6];
    let tiles = [2, 0];

    let mut results = Array3::default((tiles.len(), freqs.len(), azels.len()));
    beam.calc_jones_tiles_freqs(
        &azels,
        &freqs,
        &tiles,
        latitude_rad,
        Some(time),
        results.view_mut(),
    )
    .unwrap();

    for (i_tile, &tile) in tiles.iter().enumerate() {
        for (i_freq, &freq) in freqs.iter().enumerate() {
            for (i_dir, &azel) in azels.iter().enumerate().step_by(13) {
                let expected = beam
                    .calc_jones(azel, freq, Some(tile), latitude_rad, Some(time))
                    .unwrap();
                assert_abs_diff_eq!(results[(i_tile, i_freq, i_dir)], expected);
            }
        }
    }
}

#[test]
fn test_unique_tiles() {
    let beam = get_beam("amplitude");
    // Every station is unique; flagged tiles aren't included.
    let flagged = HashSet::from([1]);
    let (unique_tiles, map) = beam.get_unique_tiles(3, &flagged);
    assert_eq!(unique_tiles, vec![0, 2]);
    assert_eq!(map, vec![0, 0, 1]);

    let (unique_tiles, map) = beam.get_unique_tiles(3, &HashSet::new());
    assert_eq!(unique_tiles, vec![0, 1, 2]);
    assert_eq!(map, vec![0, 1, 2]);
}
