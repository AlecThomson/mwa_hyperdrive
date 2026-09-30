// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Errors associated with beam calculations.

use thiserror::Error;

use super::BEAM_TYPES_COMMA_SEPARATED;

#[derive(Error, Debug)]
pub enum BeamError {
    #[error(
        "Unrecognised beam model '{0}'; supported beam models are: {}",
        *BEAM_TYPES_COMMA_SEPARATED
    )]
    Unrecognised(String),

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

    #[error("EveryBeam error: {0}")]
    EveryBeam(String),

    #[error("The EveryBeam beam requires a measurement set to describe the telescope, but the input data isn't a measurement set; please supply one with --beam-ms")]
    NeedsBeamMs,

    #[error("This beam requires the time of the beam-response calculation, but it wasn't supplied; this is a hyperdrive bug or an unsupported beam for this feature")]
    NeedsTime,

    #[error("The beam's measurement set describes {ms} stations, but the input data has {tiles} tiles")]
    StationCountMismatch { ms: usize, tiles: usize },

    #[cfg(any(feature = "cuda", feature = "hip"))]
    #[error(transparent)]
    Gpu(#[from] crate::gpu::GpuError),
}
