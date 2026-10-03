// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Code to handle DP3/makesourcedb ("BBS") source list files, as written by
//! e.g. DP3 and WSClean.
//!
//! See for more info:
//! <https://mwatelescope.github.io/mwa_hyperdrive/defs/source_list_dp3.html>

mod read;

// Re-exports.
pub(crate) use read::parse_source_list;
