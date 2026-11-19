#!/bin/bash
# Minimal repro for hyperdrive di-calibrate (EveryBeam build) on this machine.
# Usage: test_files/everybeam/repro_di_calibrate.sh /path/to/hyperdrive
# Mirrors the release smoke test in .github/workflows/releases.yml.
set -uo pipefail

HYP=$(readlink -f "${1:-$(command -v hyperdrive)}")
WORK=$(mktemp -d)
echo "hyperdrive: $HYP"
echo "work dir:   $WORK"

echo "== System"
grep PRETTY_NAME /etc/os-release
ldd --version | head -1

echo "== Shared libraries (missing or casacore/EveryBeam/HDF5)"
ldd "$HYP" | grep -Ei 'not found|casa|everybeam|hdf5' || echo "(none)"
echo "RUNPATH: $(readelf -d "$HYP" 2>/dev/null | grep -Eo 'R(UN)?PATH.*' || echo none)"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}"

echo "== --version"
"$HYP" --version || { echo "FAILED to start (see libraries above)"; exit 1; }

echo "== Fetching test data"
git clone -q --depth 1 --filter=blob:none --sparse \
    https://github.com/AlecThomson/mwa_hyperdrive.git "$WORK/repo"
git -C "$WORK/repo" sparse-checkout set test_files/everybeam
TF="$WORK/repo/test_files/everybeam"
cd "$WORK"

echo "== Simulate data (vis-subtract with negated fluxes)"
sed 's/        i: /        i: -/' "$TF/mini_srclist.yaml" > neg.yaml
"$HYP" vis-subtract -d "$TF/mini.ms" -s neg.yaml \
    --beam-type everybeam --everybeam-element-model skala40_wave \
    --named-sources zenith_source src_n src_s src_e src_w src_ne -o sim.ms \
    || { echo "vis-subtract FAILED"; exit 1; }

echo "== di-calibrate"
"$HYP" di-calibrate -d sim.ms -s "$TF/mini_srclist.yaml" \
    --beam-type everybeam --everybeam-element-model skala40_wave \
    --beam-ms "$TF/mini.ms" -o sols.fits \
    || { echo "di-calibrate FAILED (exit $?)"; exit 1; }

echo "== OK: $WORK/sols.fits"
