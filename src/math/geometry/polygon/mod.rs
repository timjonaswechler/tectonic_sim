pub mod core;
pub mod operations;
pub mod smoothing;
pub mod transformations;
pub mod utils;

// Re-exports für einfache Verwendung
pub use core::polygon::Polygon;
pub use core::properties::PolygonProperties;
pub use smoothing::*;
pub use transformations::*;
pub use utils::builders::PolygonBuilder;
