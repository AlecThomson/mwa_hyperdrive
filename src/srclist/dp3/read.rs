// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Parsing of DP3/makesourcedb ("BBS") source lists.
//!
//! These are comma-separated text files whose first line describes the columns,
//! e.g.
//!
//! ```text
//! format = Name, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='150e6', MajorAxis, MinorAxis, Orientation
//! J131139-221640,POINT,13:11:39.33,-22.16.41.31,43.9,[-0.8],true,183984261.97,,,
//! ```
//!
//! A value in quotes after a column name is that column's default. Lines with
//! an empty name and type define a patch; components with a patch are combined
//! into a single source named after the patch.

use std::f64::consts::LN_10;

use indexmap::IndexMap;
use marlu::{sexagesimal::sexagesimal_hms_string_to_degrees, RADec};
use vec1::Vec1;

use crate::{
    cli::Warn,
    srclist::{
        error::{ReadSourceListCommonError, ReadSourceListDp3Error, ReadSourceListError},
        ComponentType, FluxDensity, FluxDensityType, Source, SourceComponent, SourceList,
    },
};

/// When a spectrum can't be represented by a (curved) power law, it is sampled
/// at these frequencies \[Hz\] into a list.
const LIST_FREQS_HZ: (f64, f64, usize) = (50e6, 350e6, 61);

/// The columns of a DP3 source list.
struct Columns {
    /// The lower-case name of each column.
    names: Vec<String>,
    /// The default value of each column (empty if there's no default).
    defaults: Vec<String>,
}

impl Columns {
    fn index(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }
}

/// Split `s` on commas, except those inside brackets or quotes.
fn split_fields(s: &str) -> Vec<&str> {
    let mut fields = vec![];
    let mut depth = 0;
    let mut in_quotes = false;
    let mut start = 0;
    for (i, c) in s.char_indices() {
        match c {
            '\'' | '"' => in_quotes = !in_quotes,
            '[' if !in_quotes => depth += 1,
            ']' if !in_quotes => depth -= 1,
            ',' if !in_quotes && depth == 0 => {
                fields.push(s[start..i].trim());
                start = i + 1;
            }
            _ => (),
        }
    }
    fields.push(s[start..].trim());
    fields
}

/// If `line` describes the columns, return the column descriptions. Both
/// "format = Name, Type, ..." and "# (Name, Type, ...) = format" are accepted.
fn parse_format_line(line: &str) -> Option<Columns> {
    let line = line.trim();
    let spec = if let Some(rest) = line.strip_prefix('#') {
        let rest = rest.trim();
        let (lhs, rhs) = rest.rsplit_once('=')?;
        if !rhs.trim().eq_ignore_ascii_case("format") {
            return None;
        }
        lhs.trim().strip_prefix('(')?.strip_suffix(')')?
    } else {
        let (lhs, rhs) = line.split_once('=')?;
        if !lhs.trim().eq_ignore_ascii_case("format") {
            return None;
        }
        rhs
    };

    let mut names = vec![];
    let mut defaults = vec![];
    for field in split_fields(spec) {
        let (name, default) = match field.split_once('=') {
            Some((name, default)) => (name, default.trim().trim_matches(['\'', '"'])),
            None => (field, ""),
        };
        // Strip any type annotations (e.g. "Name:string").
        let name = name.split(':').next().unwrap_or(name).trim();
        if name.is_empty() {
            return None;
        }
        names.push(name.to_ascii_lowercase());
        defaults.push(default.to_string());
    }
    Some(Columns { names, defaults })
}

/// Parse an angle with an explicit unit suffix ("deg" or "rad") into degrees.
fn parse_angle_with_unit(s: &str) -> Option<f64> {
    if let Some(v) = s.strip_suffix("deg") {
        v.trim().parse().ok()
    } else if let Some(v) = s.strip_suffix("rad") {
        v.trim().parse::<f64>().ok().map(f64::to_degrees)
    } else {
        None
    }
}

/// Parse a sexagesimal "d m s" triple (as strings) with a sign taken from the
/// string (so that e.g. "-00" is negative).
fn dms_to_degrees(negative: bool, d: &str, m: &str, s: &str) -> Option<f64> {
    let d: f64 = d.trim_start_matches(['-', '+']).parse().ok()?;
    let m: f64 = m.parse().ok()?;
    let s: f64 = s.parse().ok()?;
    if !(0.0..60.0).contains(&m) || !(0.0..60.0).contains(&s) {
        return None;
    }
    let v = d + m / 60.0 + s / 3600.0;
    Some(if negative { -v } else { v })
}

/// Parse a DP3 RA ("hh:mm:ss.s", "XXhYYmZZs", or with a "deg"/"rad" unit) into
/// degrees.
fn parse_ra(s: &str) -> Option<f64> {
    let s = s.trim();
    if let Some(v) = parse_angle_with_unit(s) {
        return Some(v);
    }
    if s.contains('h') {
        return sexagesimal_hms_string_to_degrees(s).ok();
    }
    let parts: Vec<&str> = s.split(':').collect();
    match parts.as_slice() {
        [h, m, sec] => dms_to_degrees(s.starts_with('-'), h, m, sec).map(|h| h * 15.0),
        _ => None,
    }
}

/// Parse a DP3 Dec ("dd.mm.ss.s", "dd:mm:ss.s", "XXdYYmZZs", or with a
/// "deg"/"rad" unit) into degrees.
fn parse_dec(s: &str) -> Option<f64> {
    let s = s.trim();
    if let Some(v) = parse_angle_with_unit(s) {
        return Some(v);
    }
    let negative = s.starts_with('-');
    if s.contains('d') {
        let (d, rest) = s.split_once('d')?;
        let (m, rest) = rest.split_once('m')?;
        return dms_to_degrees(negative, d, m, rest.trim_end_matches('s'));
    }
    let parts: Vec<&str> = if s.contains(':') {
        s.split(':').collect()
    } else {
        s.split('.').collect()
    };
    match parts.as_slice() {
        [d, m, sec] => dms_to_degrees(negative, d, m, sec),
        // "dd.mm.ss.sss"; the last part is the fractional seconds.
        [d, m, sec, frac] if !s.contains(':') => {
            dms_to_degrees(negative, d, m, &format!("{sec}.{frac}"))
        }
        _ => None,
    }
}

/// Parse a spectral index list, e.g. "[-0.7, 0.01]" or "[]".
fn parse_spectral_index(s: &str) -> Option<Vec<f64>> {
    let inner = s.trim().strip_prefix('[')?.strip_suffix(']')?.trim();
    if inner.is_empty() {
        return Some(vec![]);
    }
    inner.split(',').map(|v| v.trim().parse().ok()).collect()
}

fn parse_bool(s: &str) -> Option<bool> {
    match s.to_ascii_lowercase().as_str() {
        "true" | "t" | "1" | "yes" => Some(true),
        "false" | "f" | "0" | "no" => Some(false),
        _ => None,
    }
}

/// Convert a DP3 spectrum into a hyperdrive flux-density type. `fd` holds the
/// Stokes flux densities at the reference frequency; DP3 applies the same
/// spectral shape to all of them.
fn spectrum(
    fd: FluxDensity,
    terms: &[f64],
    logarithmic: bool,
) -> (FluxDensityType, bool /* sampled */) {
    match (logarithmic, terms) {
        (_, []) => (FluxDensityType::PowerLaw { si: 0.0, fd }, false),
        (true, [si]) => (FluxDensityType::PowerLaw { si: *si, fd }, false),
        // DP3: log10(S) = log10(S0) + c0*x + c1*x^2, x = log10(ν/ν0).
        // hyperdrive: S = S0 (ν/ν0)^si exp(q ln²(ν/ν0)), so q = c1 / ln(10).
        (true, [si, c1]) => (
            FluxDensityType::CurvedPowerLaw {
                si: *si,
                fd,
                q: c1 / LN_10,
            },
            false,
        ),
        _ => {
            let (start, end, n) = LIST_FREQS_HZ;
            let fds = (0..n)
                .map(|i| {
                    let freq = start + (end - start) * i as f64 / (n - 1) as f64;
                    let ratio = if logarithmic {
                        let x = (freq / fd.freq).log10();
                        let exponent: f64 = terms
                            .iter()
                            .enumerate()
                            .map(|(i, c)| c * x.powi(i as i32 + 1))
                            .sum();
                        10_f64.powf(exponent)
                    } else {
                        // DP3: S = S0 + sum_i c_i (ν/ν0 - 1)^(i+1).
                        let x = freq / fd.freq - 1.0;
                        let sum: f64 = terms
                            .iter()
                            .enumerate()
                            .map(|(i, c)| c * x.powi(i as i32 + 1))
                            .sum();
                        1.0 + sum / fd.i
                    };
                    FluxDensity {
                        freq,
                        i: fd.i * ratio,
                        q: fd.q * ratio,
                        u: fd.u * ratio,
                        v: fd.v * ratio,
                    }
                })
                .collect();
            (
                FluxDensityType::List(Vec1::try_from_vec(fds).expect("not empty")),
                true,
            )
        }
    }
}

/// Parse a buffer containing a DP3 source list into a [`SourceList`].
pub(crate) fn parse_source_list<T: std::io::BufRead>(
    buf: &mut T,
) -> Result<SourceList, ReadSourceListError> {
    let mut line = String::new();
    let mut line_num: u32 = 0;
    let mut columns: Option<Columns> = None;
    // Components, grouped by source (patch) name.
    let mut sources: IndexMap<String, Vec<SourceComponent>> = IndexMap::new();
    let mut num_sampled = 0;

    while buf.read_line(&mut line)? > 0 {
        line_num += 1;
        let this_line = line.trim();
        if this_line.is_empty() {
            line.clear();
            continue;
        }

        // The first line that isn't blank must describe the columns (it may
        // be a "# (...) = format" comment).
        let cols = match columns.as_ref() {
            Some(c) => c,
            None => {
                match parse_format_line(this_line) {
                    Some(c) => columns = Some(c),
                    None if this_line.starts_with('#') => (),
                    None => return Err(ReadSourceListDp3Error::NoFormatLine(line_num).into()),
                }
                line.clear();
                continue;
            }
        };
        if this_line.starts_with('#') {
            line.clear();
            continue;
        }

        let fields = split_fields(this_line);
        if fields.len() > cols.names.len() {
            return Err(ReadSourceListDp3Error::TooManyFields {
                line_num,
                got: fields.len(),
                expected: cols.names.len(),
            }
            .into());
        }
        // A field's value, or the column's default if it's empty or absent.
        let get = |name: &str| -> Option<&str> {
            let i = cols.index(name)?;
            let v = fields.get(i).copied().unwrap_or("");
            let v = if v.is_empty() {
                cols.defaults[i].as_str()
            } else {
                v
            };
            let v = v.trim().trim_matches(['\'', '"']).trim();
            (!v.is_empty()).then_some(v)
        };
        let require = |name: &'static str| -> Result<&str, ReadSourceListError> {
            get(name).ok_or_else(|| {
                ReadSourceListDp3Error::MissingValue {
                    line_num,
                    column: name,
                }
                .into()
            })
        };
        let parse_float = |name: &'static str, s: &str| -> Result<f64, ReadSourceListError> {
            s.parse().map_err(|_| {
                ReadSourceListCommonError::ParseFloatError {
                    line_num,
                    string: format!("{s} ({name})"),
                }
                .into()
            })
        };
        let optional_float = |name: &'static str| -> Result<f64, ReadSourceListError> {
            get(name).map_or(Ok(0.0), |s| parse_float(name, s))
        };

        // Lines with no name and type define patches; we only need the
        // patch names of the components.
        let name = get("name");
        let comp_type = get("type");
        if name.is_none() && comp_type.is_none() {
            line.clear();
            continue;
        }
        let name = require("name")?;
        let comp_type = match require("type")?.to_ascii_uppercase().as_str() {
            "POINT" => ComponentType::Point,
            "GAUSSIAN" => ComponentType::Gaussian {
                maj: optional_float("majoraxis")?.to_radians() / 3600.0,
                min: optional_float("minoraxis")?.to_radians() / 3600.0,
                pa: optional_float("orientation")?.to_radians(),
            },
            other => {
                return Err(ReadSourceListDp3Error::UnsupportedType {
                    line_num,
                    comp_type: other.to_string(),
                }
                .into())
            }
        };

        let ra_str = require("ra")?;
        let mut ra = parse_ra(ra_str).ok_or_else(|| ReadSourceListDp3Error::InvalidRa {
            line_num,
            value: ra_str.to_string(),
        })?;
        if ra < 0.0 {
            ra += 360.0;
        }
        let dec_str = require("dec")?;
        let dec = parse_dec(dec_str).ok_or_else(|| ReadSourceListDp3Error::InvalidDec {
            line_num,
            value: dec_str.to_string(),
        })?;
        if !(0.0..=360.0).contains(&ra) {
            return Err(ReadSourceListError::InvalidRa(ra));
        }
        if !(-90.0..=90.0).contains(&dec) {
            return Err(ReadSourceListError::InvalidDec(dec));
        }

        let freq = parse_float("ReferenceFrequency", require("referencefrequency")?)?;
        let fd = FluxDensity {
            freq,
            i: parse_float("I", require("i")?)?,
            q: optional_float("q")?,
            u: optional_float("u")?,
            v: optional_float("v")?,
        };
        let terms = match get("spectralindex") {
            None => vec![],
            Some(s) => {
                parse_spectral_index(s).ok_or_else(|| ReadSourceListDp3Error::InvalidValue {
                    line_num,
                    column: "SpectralIndex",
                    value: s.to_string(),
                })?
            }
        };
        let logarithmic = match get("logarithmicsi") {
            None => true,
            Some(s) => parse_bool(s).ok_or_else(|| ReadSourceListDp3Error::InvalidValue {
                line_num,
                column: "LogarithmicSI",
                value: s.to_string(),
            })?,
        };
        let (flux_type, sampled) = spectrum(fd, &terms, logarithmic);
        if sampled {
            num_sampled += 1;
        }

        let source_name = get("patch").unwrap_or(name).to_string();
        sources
            .entry(source_name)
            .or_default()
            .push(SourceComponent {
                radec: RADec::from_degrees(ra, dec),
                comp_type,
                flux_type,
            });

        line.clear();
    }

    if columns.is_none() {
        return Err(ReadSourceListDp3Error::NoFormatLine(line_num).into());
    }
    if sources.is_empty() {
        return Err(ReadSourceListCommonError::NoSources(line_num).into());
    }
    if num_sampled > 0 {
        let (start, end, n) = LIST_FREQS_HZ;
        format!(
            "{num_sampled} DP3 source-list components have spectra that can't be represented as (curved) power laws; they have been sampled at {n} frequencies between {} and {} MHz",
            start / 1e6,
            end / 1e6
        )
        .warn();
    }

    let mut source_list = SourceList::new();
    for (name, components) in sources {
        source_list.insert(
            name,
            Source {
                components: components.into_boxed_slice(),
            },
        );
    }
    Ok(source_list)
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use approx::assert_abs_diff_eq;

    use super::*;

    fn parse(s: &str) -> Result<SourceList, ReadSourceListError> {
        parse_source_list(&mut Cursor::new(s))
    }

    #[test]
    fn test_user_example() {
        let sl = parse(
            "format = Name, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='183984261.97193283', MajorAxis, MinorAxis, Orientation
J131139-221640,POINT,13:11:39.3334008,-22.16.41.315016,43.94590718123622,[],true,183984261.97193283,,,
",
        )
        .unwrap();
        assert_eq!(sl.len(), 1);
        let comp = &sl["J131139-221640"].components[0];
        assert_abs_diff_eq!(comp.radec.ra.to_degrees(), 197.91388917, epsilon = 1e-8);
        assert_abs_diff_eq!(comp.radec.dec.to_degrees(), -22.27814306, epsilon = 1e-8);
        assert!(matches!(comp.comp_type, ComponentType::Point));
        match &comp.flux_type {
            FluxDensityType::PowerLaw { si, fd } => {
                assert_abs_diff_eq!(*si, 0.0);
                assert_abs_diff_eq!(fd.freq, 183984261.97193283);
                assert_abs_diff_eq!(fd.i, 43.94590718123622);
            }
            _ => panic!("expected a power law"),
        }
    }

    #[test]
    fn test_coordinates() {
        assert_abs_diff_eq!(
            parse_ra("13:11:39.3334").unwrap(),
            197.913889,
            epsilon = 1e-6
        );
        assert_abs_diff_eq!(
            parse_ra("13h11m39.3334s").unwrap(),
            197.913889,
            epsilon = 1e-6
        );
        assert_abs_diff_eq!(parse_ra("197.5deg").unwrap(), 197.5);
        assert_abs_diff_eq!(parse_ra("1rad").unwrap(), 1.0_f64.to_degrees());
        assert!(parse_ra("197.5").is_none());

        assert_abs_diff_eq!(
            parse_dec("-22.16.41.315").unwrap(),
            -22.278143,
            epsilon = 1e-6
        );
        assert_abs_diff_eq!(
            parse_dec("-22:16:41.315").unwrap(),
            -22.278143,
            epsilon = 1e-6
        );
        assert_abs_diff_eq!(parse_dec("+22.16.41").unwrap(), 22.278056, epsilon = 1e-6);
        assert_abs_diff_eq!(
            parse_dec("-22d16m41.315s").unwrap(),
            -22.278143,
            epsilon = 1e-6
        );
        // The sign of "-00" must not be lost.
        assert_abs_diff_eq!(parse_dec("-00.30.00").unwrap(), -0.5);
        assert_abs_diff_eq!(parse_dec("-00:30:00").unwrap(), -0.5);
        assert_abs_diff_eq!(parse_dec("-45.5deg").unwrap(), -45.5);
        assert!(parse_dec("-22.5").is_none());
        assert!(parse_dec("-22.75.00").is_none());
    }

    #[test]
    fn test_spectra() {
        let sl = parse(
            "FORMAT = Name, Type, Ra, Dec, I, Q, U, V, SpectralIndex='[]', LogarithmicSI='true', ReferenceFrequency='100e6'
pl, POINT, 00:00:00, +00.00.00, 2.0, 0.2, 0.1, 0.0, [-0.7]
cpl, POINT, 01:00:00, +10.00.00, 2.0, , , , [-0.7, 0.1]
cubic, POINT, 02:00:00, +20.00.00, 2.0, , , , [-0.7, 0.1, 0.01]
lin, POINT, 03:00:00, +30.00.00, 2.0, , , , [0.5, -0.1], false
flat, POINT, 04:00:00, +40.00.00, 2.0
",
        )
        .unwrap();
        assert_eq!(sl.len(), 5);
        let dp3_log = |terms: &[f64], freq: f64| {
            let x = (freq / 100e6).log10();
            2.0 * 10_f64.powf(
                terms
                    .iter()
                    .enumerate()
                    .map(|(i, c)| c * x.powi(i as i32 + 1))
                    .sum(),
            )
        };
        let dp3_lin = |terms: &[f64], freq: f64| {
            let x = freq / 100e6 - 1.0;
            2.0 + terms
                .iter()
                .enumerate()
                .map(|(i, c)| c * x.powi(i as i32 + 1))
                .sum::<f64>()
        };
        // Frequencies in the sampled list's grid are exact; check the
        // analytic ones elsewhere too.
        for freq in [80e6, 155e6, 200e6] {
            let at = |name: &str| sl[name].components[0].flux_type.estimate_at_freq(freq);
            assert_abs_diff_eq!(at("pl").i, dp3_log(&[-0.7], freq), epsilon = 1e-10);
            assert_abs_diff_eq!(at("pl").q, 0.1 * dp3_log(&[-0.7], freq), epsilon = 1e-10);
            assert_abs_diff_eq!(at("cpl").i, dp3_log(&[-0.7, 0.1], freq), epsilon = 1e-10);
            assert_abs_diff_eq!(
                at("cubic").i,
                dp3_log(&[-0.7, 0.1, 0.01], freq),
                epsilon = 1e-10
            );
            assert_abs_diff_eq!(at("lin").i, dp3_lin(&[0.5, -0.1], freq), epsilon = 1e-10);
            assert_abs_diff_eq!(at("flat").i, 2.0);
        }
        assert!(matches!(
            sl["cpl"].components[0].flux_type,
            FluxDensityType::CurvedPowerLaw { .. }
        ));
        assert!(matches!(
            sl["cubic"].components[0].flux_type,
            FluxDensityType::List(_)
        ));
    }

    #[test]
    fn test_gaussians_and_patches() {
        let sl = parse(
            "# (Name, Type, Patch, Ra, Dec, I, ReferenceFrequency='150e6', MajorAxis, MinorAxis, Orientation) = format

# A patch definition, then its components.
, , patch_a, 12:00:00, -30.00.00
a1, GAUSSIAN, patch_a, 12:00:00, -30.00.00, 1.0, , 120.0, 60.0, 45.0
a2, POINT, patch_a, 12:00:10, -30.00.10, 0.5
b, POINT, , 13:00:00, -31.00.00, 3.0
",
        )
        .unwrap();
        assert_eq!(sl.len(), 2);
        assert_eq!(sl["patch_a"].components.len(), 2);
        assert_eq!(sl["b"].components.len(), 1);
        match sl["patch_a"].components[0].comp_type {
            ComponentType::Gaussian { maj, min, pa } => {
                assert_abs_diff_eq!(maj.to_degrees() * 3600.0, 120.0, epsilon = 1e-10);
                assert_abs_diff_eq!(min.to_degrees() * 3600.0, 60.0, epsilon = 1e-10);
                assert_abs_diff_eq!(pa.to_degrees(), 45.0, epsilon = 1e-10);
            }
            _ => panic!("expected a Gaussian"),
        }
    }

    #[test]
    fn test_errors() {
        // Other formats are rejected immediately (important for
        // auto-detection).
        assert!(matches!(
            parse("skymodel fileformat 1.1\nsource {\n"),
            Err(ReadSourceListError::Dp3(
                ReadSourceListDp3Error::NoFormatLine(1)
            ))
        ));
        assert!(matches!(
            parse("SOURCE a P 1.0 2.0\n"),
            Err(ReadSourceListError::Dp3(
                ReadSourceListDp3Error::NoFormatLine(1)
            ))
        ));
        assert!(matches!(
            parse("format = Name, Type, Ra, Dec, I\n"),
            Err(ReadSourceListError::Common(
                ReadSourceListCommonError::NoSources(1)
            ))
        ));
        let header = "format = Name, Type, Ra, Dec, I, ReferenceFrequency\n";
        assert!(matches!(
            parse(&format!(
                "{header}a, SHAPELET, 00:00:00, +00.00.00, 1.0, 1e8\n"
            )),
            Err(ReadSourceListError::Dp3(
                ReadSourceListDp3Error::UnsupportedType { .. }
            ))
        ));
        assert!(matches!(
            parse(&format!("{header}a, POINT, 00:00:00, +00.00.00, 1.0\n")),
            Err(ReadSourceListError::Dp3(
                ReadSourceListDp3Error::MissingValue {
                    column: "referencefrequency",
                    ..
                }
            ))
        ));
        assert!(matches!(
            parse(&format!("{header}a, POINT, 0.0, +00.00.00, 1.0, 1e8\n")),
            Err(ReadSourceListError::Dp3(
                ReadSourceListDp3Error::InvalidRa { .. }
            ))
        ));
        assert!(matches!(
            parse(&format!("{header}a, POINT, 00:00:00, +00.00.00, x, 1e8\n")),
            Err(ReadSourceListError::Common(
                ReadSourceListCommonError::ParseFloatError { line_num: 2, .. }
            ))
        ));
        assert!(matches!(
            parse(&format!(
                "{header}a, POINT, 00:00:00, +00.00.00, 1.0, 1e8, 7\n"
            )),
            Err(ReadSourceListError::Dp3(
                ReadSourceListDp3Error::TooManyFields { .. }
            ))
        ));
    }
}

#[cfg(test)]
mod file_tests {
    use std::path::Path;

    use crate::srclist::{read::read_source_list_file, SourceListType};

    #[test]
    fn test_auto_detect() {
        let (sl, sl_type) =
            read_source_list_file(Path::new("test_files/dp3_srclist.txt"), None).unwrap();
        assert_eq!(sl_type, SourceListType::Dp3);
        assert_eq!(sl.len(), 2);
        assert_eq!(sl["cluster"].components.len(), 2);

        // Other formats are still detected as themselves.
        for (file, expected) in [
            ("test_files/ao_cluster_srclist.txt", SourceListType::AO),
            (
                "test_files/srclist_1099334672_100_rts.txt",
                SourceListType::Rts,
            ),
            (
                "test_files/srclist_1099334672_100.yaml",
                SourceListType::Hyperdrive,
            ),
        ] {
            let (_, sl_type) = read_source_list_file(Path::new(file), None).unwrap();
            assert_eq!(sl_type, expected, "{file}");
        }
    }
}
