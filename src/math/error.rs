// src/math/error.rs
use thiserror::Error;

#[derive(Error, Debug)]
pub enum MathError {
    #[error("Insufficient points for operation: expected at least {expected}, got {actual}")]
    InsufficientPoints { expected: usize, actual: usize },

    #[error("Invalid configuration: {message}")]
    InvalidConfiguration { message: String },

    #[error("Triangulation failed: {reason}")]
    TriangulationFailed { reason: String },

    #[error("Geometric calculation failed: {operation}")]
    GeometricFailure { operation: String },

    #[error("Conversion error: {from} to {to}")]
    ConversionError { from: String, to: String },

    #[error("No Voronoi cells found.")]
    EmptyVoronoiCells,
}

pub type MathResult<T> = Result<T, MathError>;
