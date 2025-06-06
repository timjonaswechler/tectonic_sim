// src/math/geometry/polygon/mod.rs

// Deklaration der Untermodule für Polygon-spezifische Funktionalität
pub mod builder;
pub mod core; // Enthält die Polygon-Struktur selbst
pub mod properties; // Enthält den PolygonProperties-Trait
pub mod validation; // Enthält den PolygonValidator und verwandte Typen // Enthält den PolygonBuilder (optional, falls benötigt)

// Optional: Transformations- und Operationsmodule, falls sie sehr Polygon-spezifisch sind
// Alternativ könnten diese später unter math/algorithms landen, wenn sie generischer sind.
// Fürs Erste nehmen wir an, sie sind Polygon-spezifisch:
pub mod operations;
pub mod transformations; // Enthält affine, distortion, projection für Polygone // Enthält boolean, clipping, convex_hull, triangulation für Polygone

// Re-Exporte für den einfachen Zugriff auf die wichtigsten Polygon-Elemente
pub use self::core::Polygon; // Die Polygon-Struktur
pub use self::properties::{Orientation, PolygonProperties}; // Der Eigenschaften-Trait und Enum
pub use self::validation::{
    PolygonValidator, QuickValidation, ValidationError, ValidationLevel, ValidationReport,
    ValidationSuggestion, ValidationWarning,
}; // Validator und Ergebnis-Typen

// Re-Exporte für Transformationen (wenn sie hier bleiben)
pub use self::transformations::{
    affine::{AffineTransform, AffineTransformable, TransformBuilder as AffineTransformBuilder}, // Umbenannt, um Klarheit zu schaffen
    distortion::{Distortion, DistortionEffects, DistortionFn, DistortionType, SineAxis},
    projection::{
        Complex as ProjectionComplex, ConformalType, Polygon2DProjector, Polygon3DProjector,
        Projection2DType, ProjectionBatch, ProjectionEffects as PolygonProjectionEffects,
    }, // Umbenannt
};

// Re-Exporte für Operationen (wenn sie hier bleiben)
// Diese könnten auch in ein allgemeineres `algorithms`-Modul verschoben werden,
// wenn sie nicht strikt an die `Polygon`-Struktur gebunden sind.
// Fürs Erste belassen wir sie hier, da sie oft stark mit der Polygon-Definition interagieren.
pub use self::operations::{
    boolean::{BooleanAnalysis, BooleanOperation, BooleanOperations, BooleanStats, PolygonBoolean},
    clipping::{
        ClippingAlgorithm, ClippingAnalysis, ClippingOperations, ClippingStats, PolygonClipper,
    },
    convex_hull::{ConvexHullAlgorithm, ConvexHullComputer, ConvexityDefect, ConvexityUtils},
    triangulation::{PolygonTriangulator, Triangle, TriangulationAlgorithm, TriangulationUtils},
};

// --- Builders ---
pub use self::builder::{PolygonBuilder, ShapeGenerators}; // Der PolygonBuilder und sein Fehler-Typ

// Dann hier re-exportieren:
// pub use self::builder::PolygonBuilder;
