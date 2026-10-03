# Installing with EveryBeam support

`hyperdrive` can optionally use
[EveryBeam](https://git.astron.nl/RD/EveryBeam) for beam responses. This
allows beam models other than the MWA FEE beam to be used, e.g. SKA-Low,
LOFAR and OSKAR-simulated arrays; see [EveryBeam beam
responses](../defs/beam.md#everybeam) for usage.

There are three ways to get EveryBeam support:

1. [Pre-compiled binaries](#pre-compiled-binaries) (Linux x86-64; easiest);
2. [Building from source with a vendored EveryBeam](#from-source-vendored-everybeam)
   (recommended if building from source); or
3. [Building from source against an installed EveryBeam](#from-source-installed-everybeam).

~~~admonish warning title="Licensing"
EveryBeam is licensed under the GPL-3.0, whereas `hyperdrive` is licensed under
the MPL-2.0. `hyperdrive` binaries that include EveryBeam are therefore subject
to the GPL-3.0. Binaries without EveryBeam (the default) are unaffected.
~~~

## Pre-compiled binaries

Each `hyperdrive` [release](https://github.com/MWATelescope/mwa_hyperdrive/releases)
includes `...-everybeam.tar.gz` tarballs. These contain a `hyperdrive`
binary with EveryBeam and casacore statically linked, the other libraries it
needs (in `lib/`), and casacore's measures data (in `share/`). They run on any
x86-64 Linux with glibc 2.28 or newer (e.g. RHEL/Rocky 8+, Ubuntu 20.04+).
Extract the tarball and run `bin/hyperdrive`; no other setup is required.

## From source: vendored EveryBeam

With the `everybeam-vendored` feature, the build downloads EveryBeam and
casacore (and some header-only dependencies), verifies them, builds them and
links them statically. Nothing needs to be installed except some common system
packages:

| Distribution | Packages |
| ------------ | -------- |
| Debian/Ubuntu | `build-essential cmake curl git gfortran flex bison libboost-dev libhdf5-dev libfftw3-dev libgsl-dev libblas-dev liblapack-dev casacore-data` |
| Fedora/RHEL | `gcc-c++ gcc-gfortran cmake curl git flex bison boost-devel hdf5-devel fftw-devel gsl-devel blas-devel lapack-devel` |

(plus `hyperdrive`'s usual [dependencies](from_source.md)). On RHEL-like
systems, `hdf5-devel` comes from EPEL and the BLAS/LAPACK packages from the
PowerTools/CRB repository, and a newer GCC (e.g. `gcc-toolset-13`) is needed on
RHEL 8, because EveryBeam needs a C++20 compiler (GCC 10 or later). Then:

```shell
# From a clone of the hyperdrive repo:
cargo install --path . --locked --features everybeam-vendored
```

~~~admonish info title="Build time"
The first build compiles casacore and EveryBeam, which adds a few minutes (on
the order of 5-10 minutes on 4 cores; less with more cores). These are built
inside cargo's target directory (using about 2.5 GB of disk space), so later
`cargo build`s don't rebuild them.
To avoid rebuilding them for every `cargo install`:
- use a persistent target directory, e.g.
  `CARGO_TARGET_DIR=~/.cache/hyperdrive-target cargo install ...`; and/or
- use a compiler cache: `export CMAKE_C_COMPILER_LAUNCHER=sccache
  CMAKE_CXX_COMPILER_LAUNCHER=sccache` (or `ccache`) before building.
~~~

For offline builds (or mirrors), any of the downloaded sources can instead be
supplied as an extracted directory with `EVERYBEAM_SYS_<NAME>_SRC`, where
`<NAME>` is one of `CASACORE`, `EVERYBEAM`, `AOCOMMON`, `SKA_SDP_FUNC`, `XTL`,
`XTENSOR` or `EIGEN`; see `crates/everybeam-sys/build.rs` for the versions.

EveryBeam's element-response coefficients (e.g. for SKA-Low and LOFAR element
models) are embedded in the binary, and are used automatically.

## From source: installed EveryBeam

If EveryBeam (>= 0.9) and casacore (>= 3.6) are already installed, use the
`everybeam` feature:

```shell
# Only needed if they aren't in $CONDA_PREFIX, /opt/everybeam, /opt/casacore,
# /usr/local or /usr:
export EVERYBEAM_DIR=/path/to/everybeam/prefix
export CASACORE_DIR=/path/to/casacore/prefix

cargo install --path . --locked --features everybeam
```

The binary records where the libraries are (an rpath), so `LD_LIBRARY_PATH`
isn't needed. Other build-time environment variables are `EVERYBEAM_CXXFLAGS`
(extra C++ compiler flags) and `EVERYBEAM_LIBS` (extra libraries to link).

<details>
<summary>Installing casacore and EveryBeam from source</summary>

At the time of writing, Ubuntu and Debian only package casacore 3.5, so
casacore must be built from source. EveryBeam needs a C++20 compiler (e.g.
GCC 10 or later).

```shell
sudo apt install -y build-essential cmake git wget gfortran flex bison \
    libboost-all-dev libhdf5-dev libfftw3-dev libblas-dev liblapack-dev \
    libcfitsio-dev wcslib-dev libgsl-dev casacore-data

git clone --depth 1 --branch v3.8.2 https://github.com/casacore/casacore.git
mkdir casacore/build && cd casacore/build
cmake .. -DCMAKE_INSTALL_PREFIX=/opt/casacore -DCMAKE_BUILD_TYPE=Release \
    -DMODULE=ms -DBUILD_PYTHON=OFF -DBUILD_PYTHON3=OFF -DBUILD_TESTING=OFF \
    -DBUILD_SISCO=OFF -DDATA_DIR=/usr/share/casacore/data
make -j"$(nproc)" install
cd ../..

git clone --recursive --branch v0.9.0 https://git.astron.nl/RD/EveryBeam.git
mkdir EveryBeam/build && cd EveryBeam/build
cmake .. -DCMAKE_INSTALL_PREFIX=/opt/everybeam -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH=/opt/casacore
make -j"$(nproc)" install
```

See [EveryBeam's build
instructions](https://everybeam.readthedocs.io/en/latest/build-instructions.html)
for more options. Don't build against the `everybeam` Python wheels on PyPI;
they contain no headers and use a renamed casacore namespace.
</details>

## Runtime: casacore's measures data

casacore needs its "measures" data (leap seconds, Earth-orientation tables) to
convert coordinates. The pre-compiled binaries include it. Otherwise, install it
(e.g. `casacore-data` on Debian/Ubuntu) and, if casacore can't find it (errors
mentioning `TAI_UTC` or `IERS`), tell casacore where it is:

```shell
echo "measures.directory: /usr/share/casacore/data" >> ~/.casarc
```

The data can also be downloaded from ASTRON
(`ftp://ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar`).

## Check that it worked

```shell
hyperdrive di-calibrate --help | grep -i everybeam
```

should list `everybeam` as a beam type, along with the `BEAM (EVERYBEAM)`
options. From a clone of the repo, the EveryBeam tests can be run with e.g.

```shell
cargo test --features everybeam-vendored --lib beam::everybeam
cargo test --features everybeam-vendored --test integration_tests everybeam
```
