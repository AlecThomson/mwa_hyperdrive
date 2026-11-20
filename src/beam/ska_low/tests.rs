// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Tests for the SKA-Low beams. These use the same small synthetic SKA-Low
//! measurement set as the EveryBeam tests (see
//! `test_files/everybeam/make_mini_ms.py`), and compare against EveryBeam's
//! responses for it (which don't need EveryBeam to be compiled in).

use approx::assert_abs_diff_eq;
use hifitime::{Duration, Epoch};
use marlu::{c64, precession::precess_time, AzEl, Jones, RADec, XyzGeocentric};

use super::*;

// EveryBeam's reference responses for the test measurement set. (They
// don't need EveryBeam itself, so they're also used without it.)
#[cfg(not(feature = "everybeam"))]
#[path = "../everybeam/reference.rs"]
#[allow(dead_code)]
mod reference;
#[cfg(feature = "everybeam")]
use crate::beam::everybeam::reference;
use reference::*;

const MS: &str = "test_files/everybeam/mini.ms";

fn get_beam() -> SkaLowBeam {
    SkaLowBeam::new(BeamType::SkaLowAnalytic, Path::new(MS), None, Some(8)).unwrap()
}

/// The array position hyperdrive would use for the measurement set.
fn array_position() -> marlu::LatLngHeight {
    let mut t = Table::open(Path::new(MS).join("ANTENNA"), TableOpenMode::Read).unwrap();
    let n = t.n_rows();
    let mut mean = XyzGeocentric::default();
    for row in 0..n {
        let p: Vec<f64> = t.get_cell_as_vec("POSITION", row).unwrap();
        mean.x += p[0] / n as f64;
        mean.y += p[1] / n as f64;
        mean.z += p[2] / n as f64;
    }
    mean.to_earth_wgs84()
}

/// EveryBeam's reference Jones matrix in hyperdrive's convention (its columns
/// are the (North, East) sky basis; hyperdrive's are (East, North)).
fn ref_to_hyperdrive(r: &[[f64; 2]; 4]) -> Jones<f64> {
    let c = |i: usize| c64::new(r[i][0], r[i][1]);
    Jones::from([c(1), c(0), c(3), c(2)])
}

/// Directions, LST and latitude as hyperdrive's modellers make them (with
/// precession, the default).
fn reference_azels() -> (Vec<AzEl>, BeamTime, f64) {
    let pos = array_position();
    let epoch = Epoch::from_mjd_utc(REF_TIME_MJD_S / 86400.0);
    let radecs: Vec<RADec> = REF_RADECS
        .iter()
        .map(|&(ra, dec)| RADec::from_radians(ra, dec))
        .collect();
    let p = precess_time(
        pos.longitude_rad,
        pos.latitude_rad,
        radecs[0],
        epoch,
        Duration::default(),
    );
    let azels = radecs
        .iter()
        .map(|r| r.to_hadec(p.lmst_j2000).to_azel(p.array_latitude_j2000))
        .collect();
    (
        azels,
        BeamTime {
            epoch,
            lst_rad: p.lmst_j2000,
        },
        p.array_latitude_j2000,
    )
}

#[test]
fn layout_from_mini_ms() {
    let layout = read_layout(Path::new(MS)).unwrap();
    assert_eq!(layout.stations.len(), 8);
    let rotations_deg = [0.0, 30.0, 75.0, 15.0, 60.0, 105.0, 140.0, 45.0];
    for (s, r) in layout.stations.iter().zip(rotations_deg) {
        assert_eq!(s.positions_enu.len(), 16);
        // The stations are rotated rigidly: the element grid and the dipoles.
        assert_abs_diff_eq!(s.rotation_rad.to_degrees(), r, epsilon = 1e-6);
        // Element offsets are on a 1.5 m grid in the station's rotated frame.
        let (sr, cr) = s.rotation_rad.sin_cos();
        let p = s.positions_enu[0];
        let (u, v) = (cr * p[0] + sr * p[1], -sr * p[0] + cr * p[1]);
        assert_abs_diff_eq!(u, -2.25, epsilon = 1e-6);
        assert_abs_diff_eq!(v, -2.25, epsilon = 1e-6);
        assert_abs_diff_eq!(p[2], 0.0, epsilon = 1e-6);
        // These stations are close to the array centre, so their frames are
        // close to the array's.
        for (i, row) in s.local_frame.iter().enumerate() {
            assert_abs_diff_eq!(row[i], 1.0, epsilon = 1e-6);
        }
    }
    // Station 2 has X elements 0, 5 and 10 flagged.
    let s2 = &layout.stations[2];
    for e in 0..16 {
        let expected_x = if [0, 5, 10].contains(&e) { 0.0 } else { 1.0 };
        assert_eq!(s2.weights_x[e], expected_x);
        assert_eq!(s2.weights_y[e], 1.0);
    }
    assert!(layout.stations[0].weights_x.iter().all(|&w| w == 1.0));
}

/// The norms of a Jones matrix's rows. A rotation of the sky's polarisation
/// basis (which mixes the columns) doesn't change them.
fn row_norms(j: Jones<f64>) -> [f64; 2] {
    [
        (j[0].norm_sqr() + j[1].norm_sqr()).sqrt(),
        (j[2].norm_sqr() + j[3].norm_sqr()).sqrt(),
    ]
}

/// Compare with EveryBeam (`oskar_dipole_cos`, its default for SKA-Low, with
/// no normalisation) for three stations with different rotations and
/// flagged elements, two frequencies and three directions.
///
/// EveryBeam refers its responses to the J2000 celestial pole exactly, whereas
/// hyperdrive's directions (like those given to all of its beams) are in a
/// horizon frame precessed to J2000, whose azimuth reference is rotated from
/// the true one by about the precession of the pole since J2000 (0.13 deg for
/// this 2023 measurement set). That rotates the sky's polarisation basis
/// (mixing the Jones columns) by about 2e-3 radians; everything else (the
/// array factor, element pattern, rotations and flags) is compared with the
/// row norms, which this rotation doesn't change.
#[test]
fn analytic_matches_everybeam() {
    let beam = get_beam();
    let (azels, time, latitude) = reference_azels();
    let (mut max_diff, mut max_norm_diff): (f64, f64) = (0.0, 0.0);
    for (i_station, station_ref) in REF_NONE.iter().enumerate() {
        for (i_freq, freq_ref) in station_ref.iter().enumerate() {
            let got = beam
                .calc_jones_array(
                    &azels,
                    REF_FREQS[i_freq],
                    Some(i_station),
                    latitude,
                    Some(time),
                )
                .unwrap();
            for (g, dir_ref) in got.iter().zip(freq_ref) {
                let expected = ref_to_hyperdrive(dir_ref);
                for k in 0..4 {
                    max_diff = max_diff.max((g[k] - expected[k]).norm());
                }
                for (a, b) in row_norms(*g).into_iter().zip(row_norms(expected)) {
                    max_norm_diff = max_norm_diff.max((a - b).abs());
                }
            }
        }
    }
    assert!(max_diff < 3e-3, "max diff {max_diff:e}");
    assert!(max_norm_diff < 1e-3, "max row-norm diff {max_norm_diff:e}");
}

fn copy_dir(src: &Path, dst: &Path) {
    std::fs::create_dir_all(dst).unwrap();
    for entry in std::fs::read_dir(src).unwrap() {
        let entry = entry.unwrap();
        let to = dst.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &to);
        } else {
            std::fs::copy(entry.path(), &to).unwrap();
        }
    }
}

/// Some measurement sets store ELEMENT_OFFSET as 3 x N (and ELEMENT_FLAG as
/// 2 x N) rather than N x 3 (and N x 2); both must give the same stations.
#[test]
fn transposed_phased_array_cells() {
    use ndarray::Array2;

    let dir = tempfile::tempdir().unwrap();
    let ms = dir.path().join("transposed.ms");
    copy_dir(Path::new(MS), &ms);
    {
        let mut pa = Table::open(ms.join("PHASED_ARRAY"), TableOpenMode::ReadWrite).unwrap();
        for row in 0..pa.n_rows() {
            let offsets: Vec<f64> = pa.get_cell_as_vec("ELEMENT_OFFSET", row).unwrap();
            let offsets = Array2::from_shape_vec((16, 3), offsets).unwrap();
            pa.put_cell(
                "ELEMENT_OFFSET",
                row,
                &offsets.reversed_axes().as_standard_layout().into_owned(),
            )
            .unwrap();
            let flags: Vec<bool> = pa.get_cell_as_vec("ELEMENT_FLAG", row).unwrap();
            let flags = Array2::from_shape_vec((16, 2), flags).unwrap();
            pa.put_cell(
                "ELEMENT_FLAG",
                row,
                &flags.reversed_axes().as_standard_layout().into_owned(),
            )
            .unwrap();
        }
        // The cells really are transposed now.
        let offsets: Vec<f64> = pa.get_cell_as_vec("ELEMENT_OFFSET", 0).unwrap();
        let normal: Vec<f64> = Table::open(Path::new(MS).join("PHASED_ARRAY"), TableOpenMode::Read)
            .unwrap()
            .get_cell_as_vec("ELEMENT_OFFSET", 0)
            .unwrap();
        assert_eq!(offsets[1], normal[3]);
    }
    let normal = read_layout(Path::new(MS)).unwrap();
    let transposed = read_layout(&ms).unwrap();
    for (a, b) in normal.stations.iter().zip(&transposed.stations) {
        assert_eq!(a.positions_enu, b.positions_enu);
        assert_eq!(a.weights_x, b.weights_x);
        assert_eq!(a.weights_y, b.weights_y);
        assert_eq!(a.rotation_rad, b.rotation_rad);
    }
    // Station 2's flagged X elements survive the transpose.
    assert_eq!(transposed.stations[2].weights_x[5], 0.0);
}

/// How well does splining from knots on a frequency grid reproduce an SKA-Low
/// station beam (256 elements across 38 m) evaluated at every fine channel?
/// The station's sidelobes change with frequency much faster than an MWA
/// tile's, so this measures the error at increasing angles from the
/// pointing centre (zenith), for several knot spacings.
#[test]
fn ska_low_station_knot_spline_error() {
    use mwa_hyperbeam::ska::SkaLowAnalytic;

    use crate::beam::{beam_knots_on_grid_with_spacing, BeamFreqSpline};

    // A 256-element Vogel spiral of radius 19 m, like an SKA-Low station.
    let positions: Vec<[f64; 3]> = (0..256)
        .map(|i| {
            let r = 19.0 * ((i as f64 + 0.5) / 256.0).sqrt();
            let a = i as f64 * 2.399963229728653;
            [r * a.cos(), r * a.sin(), 0.0]
        })
        .collect();
    let station = SkaStation::new(SkaLowAnalytic::new(), positions, None).unwrap();
    // 10 kHz fine channels over 20 MHz.
    let freqs: Vec<f64> = (0..2000).map(|i| 150e6 + 10e3 * i as f64).collect();
    let lat = -0.4660608448386394;
    let peak = station
        .calc_jones_pair(0.0, 0.0, 150e6, 0.0, 0.0, Some(lat), false)
        .unwrap()
        .iter()
        .map(|c| c.norm())
        .fold(0.0, f64::max);
    for (paz, pza) in [(0.0, 0.0), (1.0, 0.6)] {
        let respond = |f: f64, za: f64| {
            station
                .calc_jones_pair(paz + 0.3, za, f, paz, pza, Some(lat), false)
                .unwrap()
        };
        for (spacing_hz, max_allowed) in [(1.28e6, 1e-4), (0.64e6, 1e-5), (0.32e6, 1e-6)] {
            let knots = beam_knots_on_grid_with_spacing(&freqs, 2, spacing_hz);
            let spline = BeamFreqSpline::new(knots.clone(), &freqs);
            let s = spline.second_derivs.as_ref().unwrap();
            let mut worst: f64 = 0.0;
            let mut line = format!(
            "pointing za {pza} rad, knots every {} kHz: max |spline - exact| / peak by zenith angle:",
            spacing_hz / 1e3
        );
            for za_deg in [5.0f64, 15.0, 30.0, 45.0, 60.0, 75.0] {
                let za = za_deg.to_radians();
                let at_knots: Vec<Jones<f64>> = knots.iter().map(|&f| respond(f, za)).collect();
                let mut err: f64 = 0.0;
                for (&f, &(lo, hi, [a, b, c, d])) in freqs.iter().zip(&spline.interp) {
                    let m = |i: usize| {
                        let mut j = Jones::default();
                        for (k, &w) in s.row(i).iter().enumerate() {
                            j += at_knots[k] * w;
                        }
                        j
                    };
                    let interp = at_knots[lo] * a + at_knots[hi] * b + m(lo) * c + m(hi) * d;
                    let exact = respond(f, za);
                    for k in 0..4 {
                        err = err.max((interp[k] - exact[k]).norm() / peak);
                    }
                }
                line += &format!(" {za_deg}°: {err:.1e}");
                worst = worst.max(err);
            }
            eprintln!("{line}");
            assert!(worst < max_allowed, "{line}");
        }
    }
}
