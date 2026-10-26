This is a release of `mwa_hyperdrive` with
[EveryBeam](https://git.astron.nl/RD/EveryBeam) support, obtained from the
[GitHub releases page](https://github.com/MWATelescope/mwa_hyperdrive/releases).
EveryBeam provides beam models for telescopes other than the MWA, e.g. SKA-Low
and LOFAR (use `--beam-type everybeam`).

Documentation on `hyperdrive` can be found
[here](https://mwatelescope.github.io/mwa_hyperdrive/index.html); EveryBeam
usage is described
[here](https://mwatelescope.github.io/mwa_hyperdrive/defs/beam.html#everybeam).

Run `bin/hyperdrive`. The libraries it needs are in `lib/`, and casacore's
measures data are in `share/casacore/data` (these are used automatically,
unless casacore's measures data are configured elsewhere, e.g. in `~/.casarc`).
These binaries should work on any x86-64 Linux with glibc 2.28 or newer.

The MWA FEE beam also works, but needs the MWA FEE beam HDF5 file:

  `wget http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5`
  `export MWA_BEAM_FILE=/path/to/mwa_full_embedded_element_pattern.h5`

# Licensing

`hyperdrive` is licensed under the [Mozilla Public License 2.0 (MPL
2.0)](https://www.mozilla.org/en-US/MPL/2.0/) (LICENSE-hyperdrive). This binary
includes EveryBeam and aocommon, which are licensed under the GNU General Public
License 3.0 (LICENSE-EveryBeam, LICENSE-aocommon), and casacore, which is
licensed under the GNU Lesser General Public License (LICENSE-casacore). As a
combined work, this binary is therefore distributed under the terms of the
GPL-3.0. The source code of all of these is publicly available:

- hyperdrive: https://github.com/MWATelescope/mwa_hyperdrive
- EveryBeam: https://git.astron.nl/RD/EveryBeam
- aocommon: https://gitlab.com/aroffringa/aocommon
- casacore: https://github.com/casacore/casacore

Shared libraries in `lib/` (e.g. HDF5, FFTW, BLAS/LAPACK) are included under
their own licenses, from the AlmaLinux 8 and EPEL repositories.
