# Installing with EveryBeam support

`hyperdrive` can optionally use
[EveryBeam](https://git.astron.nl/RD/EveryBeam) for beam responses. This
allows beam models other than the MWA FEE beam to be used, e.g. SKA-Low,
LOFAR, and OSKAR-simulated arrays. EveryBeam support is enabled with the
`everybeam` cargo feature, and is not included in the pre-compiled binaries.

EveryBeam is a C++ library, and it depends on
[casacore](https://github.com/casacore/casacore). The instructions below build
both from source and install them under `/opt`; adjust the prefixes to suit.
They were written for Ubuntu 24.04, but other distributions work similarly.

~~~admonish warning title="Versions"
- EveryBeam **0.9 or later** is required.
- EveryBeam >= 0.6.2 requires casacore **3.6 or later**. At the time of
  writing, Ubuntu and Debian only package casacore 3.5, so casacore must be
  built from source.
- A C++20 compiler is required (e.g. GCC 10 or later).
~~~

## 1. System dependencies

```shell
sudo apt install -y build-essential cmake git wget gfortran flex bison \
    libboost-all-dev libhdf5-dev libfftw3-dev libblas-dev liblapack-dev \
    libcfitsio-dev wcslib-dev libgsl-dev libreadline-dev \
    casacore-data
```

`casacore-data` provides the measures data (leap seconds, Earth orientation
tables, etc.) that casacore needs at runtime; see [step 5](#5-runtime-setup).

## 2. Build casacore (>= 3.6)

Only the modules that EveryBeam needs are built here.

```shell
git clone --depth 1 --branch v3.6.1 https://github.com/casacore/casacore.git
mkdir casacore/build && cd casacore/build
cmake .. \
    -DCMAKE_INSTALL_PREFIX=/opt/casacore \
    -DCMAKE_BUILD_TYPE=Release \
    -DMODULE=ms \
    -DBUILD_PYTHON=OFF -DBUILD_PYTHON3=OFF \
    -DBUILD_TESTING=OFF \
    -DUSE_OPENMP=ON \
    -DDATA_DIR=/usr/share/casacore/data
make -j"$(nproc)" install
cd ../..
```

## 3. Build EveryBeam (>= 0.9)

EveryBeam has git submodules, so it must be cloned recursively.

```shell
git clone --recursive --branch v0.9.0 https://git.astron.nl/RD/EveryBeam.git
mkdir EveryBeam/build && cd EveryBeam/build
cmake .. \
    -DCMAKE_INSTALL_PREFIX=/opt/everybeam \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH=/opt/casacore
make -j"$(nproc)" install
cd ../..
```

If the binaries will run on a different CPU to the one they're compiled on,
add `-DPORTABLE=ON`. For LOFAR's LOBEs element model, add
`-DDOWNLOAD_LOBES=ON`. See [EveryBeam's build
instructions](https://everybeam.readthedocs.io/en/latest/build-instructions.html)
for more options.

This installs the headers in `/opt/everybeam/include/EveryBeam`, the libraries
in `/opt/everybeam/lib` and the element-response coefficients in
`/opt/everybeam/share/everybeam`.

## 4. Build `hyperdrive`

Tell the build where EveryBeam and casacore are, and enable the `everybeam`
feature (it can be combined with other features, e.g. `cuda`):

```shell
export EVERYBEAM_DIR=/opt/everybeam
export CASACORE_DIR=/opt/casacore

# From a clone of the hyperdrive repo:
cargo install --path . --locked --features everybeam
```

To avoid needing `LD_LIBRARY_PATH` at runtime (see below), the library paths
can be baked into the binary:

```shell
RUSTFLAGS="-C link-args=-Wl,-rpath,/opt/everybeam/lib:/opt/casacore/lib" \
    cargo install --path . --locked --features everybeam
```

The build is controlled by these environment variables:

| Variable             | Description                                                                                              | Default            |
| -------------------- | -------------------------------------------------------------------------------------------------------- | ------------------ |
| `EVERYBEAM_DIR`      | EveryBeam install prefix(es), `:`-separated. Headers in `include` and `include/EveryBeam`, libs in `lib` or `lib64`. | `/usr/local:/usr` |
| `CASACORE_DIR`       | casacore install prefix(es), `:`-separated.                                                              | `/usr/local:/usr` |
| `EVERYBEAM_CXXFLAGS` | Extra flags for the C++ compiler, whitespace separated.                                                  |                    |
| `EVERYBEAM_LIBS`     | Extra libraries to link, whitespace separated.                                                           |                    |

~~~admonish tip title="Don't link against the EveryBeam Python wheel"
The `everybeam` wheels on PyPI contain a copy of `libeverybeam.so` but no
headers. They also rename casacore's namespace. Building against them is not
supported; use a proper installation as described above.
~~~

## 5. Runtime setup

- The EveryBeam and casacore libraries must be found at runtime, unless an
  rpath was used above:

  ```shell
  export LD_LIBRARY_PATH=/opt/everybeam/lib:/opt/casacore/lib:$LD_LIBRARY_PATH
  ```

- EveryBeam looks for its coefficient files in the data directory set at
  compile time (`/opt/everybeam/share/everybeam` above). If EveryBeam has been
  moved, set `EVERYBEAM_DATADIR`.

- casacore needs its measures data for coordinate conversions. If you see
  errors about `TAI_UTC` or `IERS` tables, make sure `casacore-data` is
  installed and point casacore at it, e.g.

  ```shell
  echo "measures.directory: /usr/share/casacore/data" >> ~/.casarc
  ```

## 6. Check that it worked

```shell
hyperdrive di-calibrate --help | grep -i everybeam
```

should list `everybeam` as a beam type, along with the `BEAM (EVERYBEAM)`
options. See [EveryBeam beam responses](../defs/beam.md#everybeam) for usage.

To run `hyperdrive`'s EveryBeam tests from a clone of the repo:

```shell
cargo test --features everybeam --lib beam::everybeam
```
