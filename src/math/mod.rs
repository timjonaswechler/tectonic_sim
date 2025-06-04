pub mod error;
pub mod geometry;
pub mod probability;
pub mod types;
pub mod utils;

// Re-exports für einfache Verwendung
pub use error::{MathError, MathResult};
pub use types::*;

// Öffentliche API
pub mod prelude {
    pub use super::{
        error::{MathError, MathResult},
        geometry::{
            polygon::operations::*,
            sphere::{coordinates::*, projection::*, sampling::*},
            voronoi::{builder::VoronoiGenerator, config::VoronoiConfig},
        },
        probability::{/*distributions::*, sampling::**/ noise::*},
        types::*,
    };
}
