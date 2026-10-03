# Beam responses

Beam responses are given by
[`mwa_hyperbeam`](https://github.com/MWATelescope/mwa_hyperbeam) (the default)
or, if `hyperdrive` was compiled with the `everybeam` feature,
[EveryBeam](#everybeam).

To function, MWA beam code needs a few things:

- The [dipole delays](mwa/delays.md);
- The [dipole gains](mwa/dead_dipoles.md) (usually dead dipoles are 0, others are 1);
- The direction we want the beam response as an Azimuth-Elevation coordinate; and
- A frequency.

In addition, the FEE beam code needs an HDF5 file to function. See the
[post-installation instructions](../installation/post.md) for information on
getting that set up.

## Errors

Beam code usually does not error, but if it does it's likely because:

1. There aren't exactly 16 dipole delays;
2. There aren't exactly 16 or 32 dipole gains per tile; or
3. There's something wrong with the FEE HDF5 file. The official file is well
   tested.
## EveryBeam

If `hyperdrive` was compiled with the `everybeam` feature (see [Installing with
EveryBeam support](../installation/everybeam.md)), `--beam-type everybeam`
uses EveryBeam for beam responses. This supports any telescope that EveryBeam
supports (e.g. SKA-Low, LOFAR, OSKAR-simulated arrays). EveryBeam reads the
telescope's description (station positions, element layouts and orientations,
pointing) from a measurement set:

- If the input visibilities are a measurement set, it is used by default.
- Otherwise (or to use a different one), supply it with `--beam-ms`.

The number of stations in the measurement set must match the number of tiles
in the input data, and the station order is assumed to be the same.

~~~admonish warning title="Measurement sets written by hyperdrive"
Measurement sets written by `hyperdrive` (e.g. by `vis-subtract` or
`solutions-apply`) don't describe non-MWA telescopes, so EveryBeam can't load
the telescope from them. When processing such outputs, point `--beam-ms` at the
original measurement set.
~~~

```shell
hyperdrive di-calibrate -d obs.ms -s srclist.yaml --beam-type everybeam
```

### Options

| Option | Description |
| ------ | ----------- |
| `--beam-ms` | The measurement set describing the telescope. |
| `--everybeam-element-model` | The element response model, e.g. `default`, `oskar_dipole_cos`, `oskar_dipole`, `skala40_wave`, `skalow_feko`, `hamaker`, `lobes`. The default depends on the telescope; for SKA-Low/OSKAR it is `oskar_dipole_cos`, which is the best match to real SKA-Low stations. |
| `--everybeam-mode` | `full` (default), `array_factor` or `element`. |
| `--everybeam-normalisation` | `none` (default), `amplitude`, `full`, `preapplied` or `preapplied_or_full`. See below. |
| `--everybeam-coeff-path` | Path to element-response coefficients (telescope dependent). |
| `--everybeam-field-id` | The measurement set `FIELD` used for the beam pointing (default 0). |
| `--everybeam-data-column` | The data column used to check for a pre-applied beam (LOFAR only). |
| `--everybeam-subband-frequency` | Use the subband frequency for the station beamformer (LOFAR only). |
| `--everybeam-frequency-interpolation` | Interpolate the beam over frequency (MWA only). |

### Polarisation conventions and normalisation

`hyperdrive`'s sky model uses the MWA [polarisation](pols.md) convention (X is
East-West, Y is North-South), whereas the IAU convention has X as North-South.
Only the sky side (the columns of the beam Jones matrices) is converted to
`hyperdrive`'s convention. The feed side (the rows) is left in the basis of the
telescope's feeds, so that it matches the visibilities.

- `none` (the default) uses each station's unnormalised response, so
  calibration solutions don't contain the beam's spectral response. `amplitude`
  instead scales each station's response, separately at each time and
  frequency, to unit amplitude in the direction of the measurement set's
  `FIELD` `REFERENCE_DIR`.
- With `none` or `amplitude` normalisation, the rows are each
  station's own feeds. For example, SKA-Low stations are rigidly rotated with
  respect to each other (as described in the measurement set's `PHASED_ARRAY`
  table), so each station's X and Y feeds are generally not East-West and
  North-South. This is appropriate for data that have not had a beam
  correction applied.
- With `full` (or `preapplied`) normalisation, EveryBeam multiplies the
  response by the inverse of the response in the `REFERENCE_DIR` direction. The rows are
  then in the (North, East) sky basis, i.e. the IAU order. This is only
  appropriate for data that have already had the beam at the phase centre
  corrected (e.g. with DP3's `applybeam`).

### Performance

EveryBeam only runs on the CPU. When GPU modelling is used, EveryBeam beam
responses are calculated on the CPU and copied to the GPU, which is slower than
the MWA FEE beam's GPU code.
