// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Errors associated with beam calculations.

use thiserror::Error;

use super::{BeamType, BEAM_TYPES_COMMA_SEPARATED};
use crate::context::{PolConvention, POL_CONVENTIONS_COMMA_SEPARATED};

#[derive(Error, Debug)]
pub enum BeamError {
    #[error(
        "Unrecognised beam model '{_0}'; supported beam models are: {}",
        *BEAM_TYPES_COMMA_SEPARATED
    )]
    Unrecognised(String),

    #[error(
        "Unrecognised polarisation convention '{_0}'; supported conventions are: {}",
        *POL_CONVENTIONS_COMMA_SEPARATED
    )]
    UnrecognisedPolConvention(String),

    #[error("The '{beam_type}' beam produces Jones matrices in the '{beam}' polarisation convention, but '{requested}' was requested. Use '--beam-type none' for non-MWA instruments")]
    PolConventionMismatch {
        beam_type: BeamType,
        beam: PolConvention,
        requested: PolConvention,
    },

    #[error(
        "Tried to set up a '{0}' beam, which requires MWA dipole delays, but none are available"
    )]
    NoDelays(&'static str),

    #[error(
        "The specified MWA dipole delays aren't valid; there should be 16 values between 0 and 32"
    )]
    BadDelays,

    #[error("There are dipole delays specified for {num_rows} tiles, but when creating the beam object, {num_tiles} was specified as the number of tiles; refusing to continue")]
    InconsistentDelays { num_rows: usize, num_tiles: usize },

    #[error("The number of delays per tile ({delays}) didn't match the number of gains per tile ({gains})")]
    DelayGainsDimensionMismatch { delays: usize, gains: usize },

    #[error("Got tile index {got}, but the biggest tile index is {max}")]
    BadTileIndex { got: usize, max: usize },

    #[error("hyperbeam error: {0}")]
    Hyperbeam(#[from] mwa_hyperbeam::fee::FEEBeamError),

    #[error("hyperbeam init error: {0}")]
    HyperbeamInit(#[from] mwa_hyperbeam::fee::InitFEEBeamError),

    #[error("hyperbeam SKA-Low error: {0}")]
    HyperbeamSkaLow(#[from] mwa_hyperbeam::ska_low::SkaLowBeamError),

    #[error("hyperbeam SKA-Low init error: {0}")]
    HyperbeamSkaLowInit(#[from] mwa_hyperbeam::ska_low::InitSkaLowBeamError),

    #[error("Tried to set up a SKA-Low beam, but the input data has no PHASED_ARRAY table; only measurement sets can supply station element positions")]
    NoPhasedArray,

    #[error("No SKA-Low beam directory was given; use --ska-low-beam-dir or set SKA_LOW_BEAM_DIR")]
    NoSkaLowBeamDir,

    #[error("The PHASED_ARRAY table describes {num_stations} stations, but there are {num_tiles} tiles; refusing to continue")]
    InconsistentStations {
        num_stations: usize,
        num_tiles: usize,
    },

    #[error("Unrecognised SKA-Low gridded-EEP format '{0}'; supported formats are: npy, npz")]
    UnrecognisedGridFormat(String),

    #[error("--ska-low-grid-filebase is required alongside --ska-low-grid-format")]
    NoGridFilebase,

    #[error("The SKA-Low beam has no GPU support; run without GPU acceleration")]
    SkaLowNoGpu,

    #[cfg(any(feature = "cuda", feature = "hip"))]
    #[error(transparent)]
    Gpu(#[from] crate::gpu::GpuError),
}
