// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Raw bindings to the C shim around EveryBeam (see `shim.h`).

use std::ffi::{c_char, c_double, c_int};

/// An opaque handle to a loaded EveryBeam telescope.
#[repr(C)]
pub(super) struct eb_telescope {
    _private: [u8; 0],
}

#[repr(C)]
pub(super) struct eb_options {
    pub(super) element_response_model: *const c_char,
    pub(super) beam_mode: *const c_char,
    pub(super) beam_normalisation_mode: *const c_char,
    pub(super) coeff_path: *const c_char,
    pub(super) data_column_name: *const c_char,
    pub(super) use_channel_frequency: c_int,
    pub(super) frequency_interpolation: c_int,
}

extern "C" {
    pub(super) fn eb_load(
        ms_path: *const c_char,
        options: *const eb_options,
        err: *mut c_char,
        err_len: usize,
    ) -> *mut eb_telescope;

    pub(super) fn eb_free(telescope: *mut eb_telescope);

    pub(super) fn eb_num_stations(telescope: *const eb_telescope) -> usize;

    #[allow(clippy::too_many_arguments)]
    pub(super) fn eb_point_responses(
        telescope: *const eb_telescope,
        time_mjd_s: c_double,
        field_id: usize,
        ras: *const c_double,
        decs: *const c_double,
        num_directions: usize,
        freqs_hz: *const c_double,
        num_freqs: usize,
        stations: *const usize,
        num_stations: usize,
        out: *mut c_double,
        station_stride: usize,
        freq_stride: usize,
        err: *mut c_char,
        err_len: usize,
    ) -> c_int;
}
