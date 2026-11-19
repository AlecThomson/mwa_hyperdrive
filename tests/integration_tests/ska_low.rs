// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! End-to-end tests with the SKA-Low beams (calculated natively with
//! hyperbeam), using the same small synthetic SKA-Low measurement set as the
//! EveryBeam tests (see `test_files/everybeam/make_mini_ms.py`).

use std::path::Path;

use approx::assert_abs_diff_eq;
use marlu::{c64, Jones};
use ndarray::prelude::*;
use tempfile::TempDir;

use crate::{get_cmd_output, hyperdrive};
use mwa_hyperdrive::CalibrationSolutions;

const MS: &str = "test_files/everybeam/mini.ms";
const SRCLIST: &str = "test_files/everybeam/mini_srclist.yaml";
const NUM_STATIONS: usize = 8;
const NUM_CHANS: usize = 3;
const BEAM: &str = "ska-low-analytic";

/// The synthetic measurement set has no visibilities. Fill a copy with a sky
/// model generated with the SKA-Low beam by "subtracting" a negated sky
/// model.
fn simulate(tmp_dir: &Path) -> std::path::PathBuf {
    let srclist = std::fs::read_to_string(SRCLIST).unwrap();
    let negated = srclist.replace("        i: ", "        i: -");
    assert_ne!(srclist, negated);
    let neg_srclist = tmp_dir.join("negated.yaml");
    std::fs::write(&neg_srclist, negated).unwrap();

    let sim = tmp_dir.join("sim.ms");
    let cmd = hyperdrive()
        .args(["vis-subtract", "--data", MS, "--source-list"])
        .arg(&neg_srclist)
        .args(["--named-sources", "zenith_source", "src_n", "src_s"])
        .args(["src_e", "src_w", "src_ne"])
        .args(["--beam-type", BEAM])
        .arg("--outputs")
        .arg(&sim)
        .arg("--no-progress-bars")
        .ok();
    assert!(cmd.is_ok(), "{:?}", get_cmd_output(cmd));
    sim
}

/// Simulate visibilities with the SKA-Low beam, corrupt them with known
/// direction-independent gains, and check that calibrating with the SKA-Low
/// beam recovers them.
#[test]
fn test_ska_low_simulate_and_calibrate() {
    let tmp_dir = TempDir::new().unwrap();
    let sim = simulate(tmp_dir.path());

    // Corrupt the visibilities. Applying solutions G gives
    // G^-1 V (G^-1)^H, so calibration should find G^-1.
    let mut gains = Array3::from_elem((1, NUM_STATIONS, NUM_CHANS), Jones::identity());
    for (i_station, mut station_gains) in gains.axis_iter_mut(Axis(1)).enumerate() {
        let g = 1.0 + 0.05 * i_station as f64;
        let phase = 0.05 * i_station as f64;
        station_gains.fill(Jones::from([
            c64::from_polar(g, phase),
            c64::new(0.0, 0.0),
            c64::new(0.0, 0.0),
            c64::from_polar(1.0 / g, -phase),
        ]));
    }
    let corruption = tmp_dir.path().join("corruption.fits");
    CalibrationSolutions {
        di_jones: gains.clone(),
        ..Default::default()
    }
    .write_solutions_from_ext::<&Path>(&corruption)
    .unwrap();
    let corrupted = tmp_dir.path().join("corrupted.ms");
    let cmd = hyperdrive()
        .arg("solutions-apply")
        .arg("--data")
        .arg(&sim)
        .arg("--solutions")
        .arg(&corruption)
        .arg("--outputs")
        .arg(&corrupted)
        .arg("--no-progress-bars")
        .ok();
    assert!(cmd.is_ok(), "{:?}", get_cmd_output(cmd));

    // Calibrate. hyperdrive-written measurement sets keep the input's
    // telescope description (including PHASED_ARRAY), so the SKA-Low beam can
    // use them directly.
    let solutions = tmp_dir.path().join("sols.fits");
    let cmd = hyperdrive()
        .arg("di-calibrate")
        .arg("--data")
        .arg(&corrupted)
        .args(["--source-list", SRCLIST])
        .args(["--beam-type", BEAM])
        .arg("--outputs")
        .arg(&solutions)
        .arg("--no-progress-bars")
        .ok();
    assert!(cmd.is_ok(), "{:?}", get_cmd_output(cmd));

    let sols =
        CalibrationSolutions::read_solutions_from_ext::<&Path, &Path>(&solutions, None).unwrap();
    assert_eq!(sols.di_jones.dim(), (1, NUM_STATIONS, NUM_CHANS));
    assert!(!sols.di_jones.iter().any(|j| j.any_nan()));

    // Corrupting with G and applying gives G^-1 V (G^-1)^H. With an
    // unpolarised sky, full-Jones calibration only determines the solutions up
    // to J_i = G_i^-1 B_i W B_i^-1, where B_i is station i's beam (which
    // differ, because the stations are rigidly rotated) and W is common to
    // all stations. G_i J_i is therefore similar to W, so its trace and
    // determinant must be the same for every station.
    let trace = |j: Jones<f64>| j[0] + j[3];
    let det = |j: Jones<f64>| j[0] * j[3] - j[1] * j[2];
    for i_chan in 0..NUM_CHANS {
        let gj = |i_station: usize| {
            gains[(0, i_station, i_chan)] * sols.di_jones[(0, i_station, i_chan)]
        };
        let gj0 = gj(0);
        for i_station in 1..NUM_STATIONS {
            assert_abs_diff_eq!(trace(gj(i_station)), trace(gj0), epsilon = 1e-3);
            assert_abs_diff_eq!(det(gj(i_station)), det(gj0), epsilon = 1e-3);
        }
    }
}

/// uvfits files don't describe SKA-Low stations, so the SKA-Low beam needs
/// --beam-ms; check that a helpful error is given otherwise, and that
/// --beam-ms works.
#[test]
fn test_ska_low_needs_telescope_ms() {
    let tmp_dir = TempDir::new().unwrap();
    let sim = tmp_dir.path().join("sim.uvfits");
    let cmd = hyperdrive()
        .args(["vis-convert", "--data", MS, "--outputs"])
        .arg(&sim)
        .arg("--no-progress-bars")
        .ok();
    assert!(cmd.is_ok(), "{:?}", get_cmd_output(cmd));
    let cmd = hyperdrive()
        .arg("di-calibrate")
        .arg("--data")
        .arg(&sim)
        .args(["--source-list", SRCLIST])
        .args(["--beam-type", BEAM])
        .arg("--outputs")
        .arg(tmp_dir.path().join("sols.fits"))
        .arg("--no-progress-bars")
        .ok();
    assert!(cmd.is_err());
    let (_, stderr) = get_cmd_output(cmd);
    assert!(stderr.contains("--beam-ms"), "{stderr}");

    let cmd = hyperdrive()
        .arg("di-calibrate")
        .arg("--data")
        .arg(&sim)
        .args(["--source-list", SRCLIST])
        .args(["--beam-type", BEAM])
        .args(["--beam-ms", MS])
        .arg("--outputs")
        .arg(tmp_dir.path().join("sols.fits"))
        .arg("--no-progress-bars")
        .ok();
    assert!(cmd.is_ok(), "{:?}", get_cmd_output(cmd));
}

/// Model the test measurement set's sky with the SKA-Low beam on the CPU
/// and on the GPU, and check that they agree. Station 2 has flagged X (but
/// not Y) elements, which takes a separate path on the GPU.
#[cfg(any(feature = "cuda", feature = "hip"))]
#[test]
fn test_ska_low_cpu_gpu_match() {
    use marlu::{
        c32,
        rubbl_casatables::{Table, TableOpenMode},
    };

    let tmp_dir = TempDir::new().unwrap();
    let srclist = std::fs::read_to_string(SRCLIST).unwrap();
    // "Subtracting" a negated sky model from the (zero) visibilities gives
    // the sky model.
    let neg_srclist = tmp_dir.path().join("negated.yaml");
    std::fs::write(&neg_srclist, srclist.replace("        i: ", "        i: -")).unwrap();
    let model = |name: &str, cpu: bool| -> Vec<c32> {
        let out = tmp_dir.path().join(name);
        let mut cmd = hyperdrive();
        cmd.args(["vis-subtract", "--data", MS, "--source-list"])
            .arg(&neg_srclist)
            .args(["--beam-type", BEAM])
            .arg("--outputs")
            .arg(&out)
            .arg("--no-progress-bars");
        if cpu {
            cmd.arg("--cpu");
        }
        let result = cmd.ok();
        assert!(result.is_ok(), "{:?}", get_cmd_output(result));
        let mut t = Table::open(&out, TableOpenMode::Read).unwrap();
        let mut vis = vec![];
        for row in 0..t.n_rows() {
            vis.extend(t.get_cell_as_vec::<c32>("DATA", row).unwrap());
        }
        vis
    };
    let cpu = model("cpu.ms", true);
    let gpu = model("gpu.ms", false);
    assert_eq!(cpu.len(), gpu.len());
    let max = cpu.iter().map(|v| v.norm()).fold(0.0, f32::max);
    assert!(max > 0.0);
    let max_diff = cpu
        .iter()
        .zip(&gpu)
        .map(|(c, g)| (c - g).norm())
        .fold(0.0, f32::max);
    // The visibilities are written in single precision.
    #[cfg(not(feature = "gpu-single"))]
    let tolerance = 1e-5;
    #[cfg(feature = "gpu-single")]
    let tolerance = 1e-3;
    eprintln!(
        "SKA-Low CPU vs GPU: max |diff| / max |vis| = {:e}",
        max_diff / max
    );
    assert!(max_diff / max < tolerance, "{max_diff} vs max {max}");
}
