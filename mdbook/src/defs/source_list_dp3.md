# The DP3 (makesourcedb/BBS) source list format

This is the comma-separated text format used by DP3, makesourcedb and BBS, and
written by e.g. WSClean (`-save-source-list`). `hyperdrive` can read it (it is
detected automatically, or use the `dp3` type), but not write it; use
`srclist-convert` to convert it to another format.

The first line describes the columns, e.g.

```plaintext
format = Name, Type, Patch, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='183984261.97', MajorAxis, MinorAxis, Orientation
J131139-221640,POINT,,13:11:39.3334008,-22.16.41.315016,43.94590718123622,[],true,183984261.97,,,
,,cluster,13:00:00.0,-25.00.00.0
cluster_g,GAUSSIAN,cluster,13:00:00.0,-25.00.00.0,2.5,[-0.8,0.05],true,,120.0,60.0,30.0
```

- Column names are case insensitive, and can be in any order. The
  `# (Name, Type, ...) = format` form of the first line is also accepted.
- A value in quotes after a column name (e.g. `ReferenceFrequency='150e6'`) is
  that column's default, used when a line's value is empty or missing.
- Used columns: `Name`, `Type` (`POINT` or `GAUSSIAN`), `Ra`, `Dec`, `I`, `Q`,
  `U`, `V` (0 if absent), `ReferenceFrequency` \[Hz\], `SpectralIndex`,
  `LogarithmicSI` (default `true`), `MajorAxis` and `MinorAxis` (FWHM, arcsec),
  `Orientation` (degrees) and `Patch`. Other columns are ignored.
- `Ra` is `hh:mm:ss.s` (or `XXhYYmZZs`); `Dec` is `dd.mm.ss.s` or `dd:mm:ss.s`
  (or `XXdYYmZZs`). Either can instead be a number with a `deg` or `rad` unit.
- Lines with an empty name and type define patches. Components with the same
  `Patch` become one `hyperdrive` source, named after the patch; components
  without a patch are each their own source.

## Spectra

DP3's spectral shape is applied to all Stokes parameters.

| DP3 spectrum | `hyperdrive` [flux-density type](fd_types.md) |
| ------------ | -------------------------------------------- |
| `SpectralIndex` `[]` | power law with a spectral index of 0 |
| `LogarithmicSI=true`, `[a]` | power law, \\( \alpha = a \\) |
| `LogarithmicSI=true`, `[a, b]` | curved power law, \\( \alpha = a \\), \\( q = b / \ln 10 \\) |
| more terms, or `LogarithmicSI=false` | list, sampled between 50 and 350 MHz (a warning is printed) |

With `LogarithmicSI=true`, DP3's spectrum is
\\( \log_{10} S = \log_{10} S_0 + \sum_i c_i \left[\log_{10}(\nu/\nu_0)\right]^{i+1} \\);
with `LogarithmicSI=false`, it is
\\( S = S_0 + \sum_i c_i (\nu/\nu_0 - 1)^{i+1} \\).
