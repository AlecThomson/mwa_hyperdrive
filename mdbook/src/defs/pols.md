# Instrumental polarisations

In `hyperdrive` (and [`mwalib`](https://github.com/MWATelescope/mwalib) and
[`hyperbeam`](https://github.com/MWATelescope/mwa_hyperbeam)), the X
polarisation refers to the East-West dipoles and the Y refers to North-South.
Note that this contrasts with the IAU definition of X and Y, which is opposite
to this. However, this is consistent within the MWA.

MWA visibilities in raw data products are ordered XX, XY, YX, YY where X is
East-West and Y is North-South. `Birli` and `cotter` also write pre-processed
visibilities this way.

`wsclean` expects its input measurement sets to be in the IAU order, meaning
that, by default, `hyperdrive` outputs are (somewhat) inappropriate for usage
with `wsclean`. We are discussing how to move forward given the history of MWA
data processing and expectations in the community. See
[Polarisation conventions](#polarisation-conventions) below for how to model in
the IAU convention.

We expect that any input data contains 4 cross-correlation polarisations (XX XY
YX YY), but `hyperdrive` is able to read the following combinations out of the
supported [input data types](./vis_formats_read.md):
- XX
- YY
- XX YY
- XX XY YY

In addition, uvfits files need not have a weight associated with each
polarisation.

# Stokes polarisations

By default (`--pol-convention mwa`), `hyperdrive` converts a sky-model
component's Stokes flux densities into instrumental polarisations as:
- \\( \text{XX} = \text{I} - \text{Q} \\)
- \\( \text{XY} = \text{U} - i\text{V} \\)
- \\( \text{YX} = \text{U} + i\text{V} \\)
- \\( \text{YY} = \text{I} + \text{Q} \\)

where \\( \text{I} \\), \\( \text{Q} \\), \\( \text{U} \\), \\( \text{V} \\) are
Stokes polarisations and \\( i \\) is the imaginary unit.

# Polarisation conventions

`--pol-convention` declares which physical dipole the label "X" refers to in the
visibility data being modelled. It is available on `di-calibrate`,
`vis-simulate`, `vis-subtract` and `peel`.

| Value | X | Y | Used by |
| ----- | - | - | ------- |
| `mwa` (default) | East-West | North-South | MWA, `mwalib`, `hyperbeam` |
| `iau` | North-South | East-West | IAU, TMS, SKA |

With `--pol-convention iau` the conversion above becomes:
- \\( \text{XX} = \text{I} + \text{Q} \\)
- \\( \text{XY} = \text{U} + i\text{V} \\)
- \\( \text{YX} = \text{U} - i\text{V} \\)
- \\( \text{YY} = \text{I} - \text{Q} \\)

Note that only the *feed labelling* differs; `hyperdrive`'s sky-frame Stokes
convention is IAU either way (\\( \text{Q} > 0 \\) means North-South linear
polarisation). Equivalently, the two conventions are the same brightness matrix
with its rows and columns swapped.

The MWA FEE beam produces Jones matrices with X as East-West, so it cannot be
combined with `--pol-convention iau`; `hyperdrive` rejects that combination.
Non-MWA instruments should therefore use `--beam-type none` for now.
