// src/math/geometry/metalballs/mod.rs

pub mod builder; // MetaballsBuilder
pub mod components; // GridConfig
pub mod field; // MetaballField
pub mod influence; // Trait FieldInfluence und Implementierungen wie Composite, Line, Point
pub mod marching_squares; // MarchingSquares Algorithmus und Contour
pub mod metalballs;

// Wichtige Re-Exporte für die einfache Nutzung des Metaball-Systems
pub use self::builder::MetaballsBuilder;
pub use self::components::GridConfig;
pub use self::field::MetaballField;
pub use self::influence::{
    CombinationMode, CompositeInfluence, FieldInfluence, InfluenceUtils, LineInfluence,
    MetaballFalloff, MetaballSource, PointInfluence, PolygonFalloff, PolygonInfluence,
    PolygonInfluenceBuilder,
};
pub use self::marching_squares::{Contour, MarchingSquares, MarchingSquaresIterator};
pub use self::metalballs::Metaballs;
