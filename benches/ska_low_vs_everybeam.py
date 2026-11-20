# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "python-casacore"]
# ///

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""Head-to-head of hyperdrive's SKA-Low beams (hyperbeam) and its EveryBeam
beam, on real data.

For each SKA-Low element model, this runs the same `hyperdrive di-calibrate`
command with:

- hyperbeam on the GPU (`--beam-type ska-low-*`);
- hyperbeam on the CPU (the same, with `--cpu`);
- EveryBeam (`--beam-type everybeam --everybeam-element-model ...`). EveryBeam
  only runs on the CPU; in a GPU build of hyperdrive, its responses are
  copied to the GPU for modelling.

and reports the time spent reading data and sky modelling (the stage that
the beam dominates; it's taken from hyperdrive's log), the total run time and
the calibration convergence. The model visibilities (`--model-filenames`)
are compared against EveryBeam's.

| hyperbeam          | EveryBeam          |
|--------------------|--------------------|
| `ska-low-analytic` | `oskar_dipole_cos` |
| `ska-low-eep-sph`  | `skalow_feko`      |
| `ska-low-eep-grid` | `gridded_eep`      |

hyperdrive must have been built with GPU support and EveryBeam, e.g.

    cargo build --release --locked --features cuda,hdf5-static,cfitsio-static,everybeam-vendored

Expected differences between the model visibilities:

- hyperdrive's directions are in a J2000-precessed horizon frame (for all of
  its beams), whereas EveryBeam refers responses to the J2000 pole exactly.
  This rotates the sky's polarisation basis by roughly the precession of the
  pole since J2000 (~0.15 deg in 2026), i.e. differences of a few 1e-3.
- The vendored EveryBeam (0.9.0) predates two `gridded_eep` fixes on
  EveryBeam's master branch (the correlation order and the SKA-Low
  conjugation), which hyperbeam follows; `ska-low-eep-grid` won't match it,
  but its timings are still comparable.

Usage:

    uv run benches/ska_low_vs_everybeam.py --hyperdrive target/release/hyperdrive \\
        --ms 3C063.ms --source-list 3C063.gleam_model.txt \\
        --feko-dir /path/to/skalow --grid-file /path/to/eep_grid.h5 \\
        -- --uvw-min 500m --freq-average 32

Arguments after `--` are passed to every `hyperdrive di-calibrate` run.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

EVERYBEAM_MODELS = {
    "ska-low-analytic": "oskar_dipole_cos",
    "ska-low-eep-sph": "skalow_feko",
    "ska-low-eep-grid": "gridded_eep",
}

LOG_TIME = re.compile(r"^\[(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d)Z")


def log_time(line):
    m = LOG_TIME.match(line)
    return datetime.fromisoformat(m.group(1)) if m else None


def run(cmd, log_path):
    """Run a hyperdrive command, returning (wall seconds, modelling seconds,
    converged string, log text)."""
    start = time.perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    wall = time.perf_counter() - start
    log_path.write_text(proc.stdout)
    if proc.returncode != 0:
        tail = "\n".join(proc.stdout.splitlines()[-15:])
        raise RuntimeError(f"{' '.join(cmd)} failed (exit {proc.returncode}); see {log_path}:\n{tail}")
    t0 = t1 = None
    converged = "?"
    for line in proc.stdout.splitlines():
        if "Reading input data and sky modelling" in line:
            t0 = log_time(line)
        elif "Finished reading input data and sky modelling" in line:
            t1 = log_time(line)
        elif "chanblocks converged" in line:
            converged = line.split("]", 1)[1].strip()
    modelling = (t1 - t0).total_seconds() if t0 and t1 else float("nan")
    return wall, modelling, converged, proc.stdout


def read_model(ms):
    import casacore.tables as t

    with t.table(str(ms), ack=False) as tab:
        data = tab.getcol("DATA")
        flags = tab.getcol("FLAG")
    return np.where(flags, 0, data)


def main():
    argv = sys.argv[1:]
    extra = []
    if "--" in argv:
        i = argv.index("--")
        argv, extra = argv[:i], argv[i + 1 :]
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hyperdrive", required=True, help="hyperdrive binary built with GPU support and EveryBeam")
    p.add_argument("--ms", required=True, help="SKA-Low measurement set")
    p.add_argument("--source-list", required=True)
    p.add_argument("--feko-dir", help="FEKO coefficient directory (ska-low-eep-sph / skalow_feko); skipped without it")
    p.add_argument("--grid-file", help="Gridded EEP HDF5 file (ska-low-eep-grid / gridded_eep); skipped without it")
    p.add_argument("--models", nargs="+", default=list(EVERYBEAM_MODELS), choices=list(EVERYBEAM_MODELS))
    p.add_argument("--skip-cpu", action="store_true", help="Don't run hyperbeam on the CPU")
    p.add_argument("--out-dir", default="ska_low_vs_everybeam", help="Where to write solutions, models and logs")
    args = p.parse_args(argv)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    version = subprocess.run([args.hyperdrive, "--version"], capture_output=True, text=True).stdout.strip()
    print(f"# SKA-Low beams in hyperdrive: hyperbeam vs EveryBeam\n")
    print(f"- hyperdrive: {version} ({args.hyperdrive})")
    print(f"- data: {args.ms}; sky model: {args.source_list}")
    print(f"- extra di-calibrate arguments: {' '.join(extra) or '(none)'}")
    gpu = shutil.which("nvidia-smi")
    if gpu:
        name = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"], capture_output=True, text=True).stdout.strip()
        print(f"- GPU: {name}")
    print(f"- CPUs: {os.cpu_count()}\n")

    rows = []
    comparisons = []
    for model in args.models:
        data = {"ska-low-eep-sph": args.feko_dir, "ska-low-eep-grid": args.grid_file}.get(model, "")
        if data is None:
            print(f"Skipping {model} (no data given)", file=sys.stderr)
            continue
        runs = [("hyperbeam GPU", ["--beam-type", model] + (["--beam-file", data] if data else []))]
        if not args.skip_cpu:
            runs.append(("hyperbeam CPU", runs[0][1] + ["--cpu"]))
        runs.append(
            (
                "EveryBeam",
                ["--beam-type", "everybeam", "--everybeam-element-model", EVERYBEAM_MODELS[model]]
                + (["--everybeam-coeff-path", data] if data else []),
            )
        )
        models = {}
        for label, beam_args in runs:
            tag = f"{model}_{label.replace(' ', '_')}"
            sols, model_ms = out / f"sols_{tag}.fits", out / f"model_{tag}.ms"
            if model_ms.exists():
                shutil.rmtree(model_ms)
            cmd = [args.hyperdrive, "di-calibrate", "-d", args.ms, "-s", args.source_list, *beam_args, *extra]
            cmd += ["-o", str(sols), "--model-filenames", str(model_ms), "--no-progress-bars"]
            print(f"Running {model} with {label}...", file=sys.stderr)
            wall, modelling, converged, _ = run(cmd, out / f"log_{tag}.txt")
            rows.append((model, label, modelling, wall, converged))
            models[label] = model_ms
        ref = read_model(models["EveryBeam"])
        scale = np.abs(ref).max()
        for label in models:
            if label == "EveryBeam":
                continue
            m = read_model(models[label])
            comparisons.append((model, label, np.abs(m - ref).max() / scale, np.abs(m).max() / scale))

    print("## Timings\n")
    print("| model | beam | reading + sky modelling [s] | total [s] | convergence |")
    print("|---|---|---|---|---|")
    for model, label, modelling, wall, converged in rows:
        print(f"| {model} | {label} | {modelling:.0f} | {wall:.1f} | {converged} |")
    print("\nThe reading and sky-modelling time is from hyperdrive's log (1 s resolution).\n")
    print("## Model visibilities vs EveryBeam\n")
    print("| model | beam | max abs. diff / max abs. EveryBeam vis. | max abs. vis. / max abs. EveryBeam vis. |")
    print("|---|---|---|---|")
    for model, label, diff, ratio in comparisons:
        print(f"| {model} | {label} | {diff:.2e} | {ratio:.4f} |")
    print(
        "\nExpected: a few 1e-3 from the polarisation basis (hyperdrive's J2000-precessed frame vs EveryBeam's"
        " exact J2000), except for ska-low-eep-grid against EveryBeam 0.9.0 (see the script's docstring)."
    )


if __name__ == "__main__":
    main()
