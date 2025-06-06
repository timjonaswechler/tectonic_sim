// src/math/algorithms/smoothing/mod.rs

// Deklaration der verschiedenen Smoothing-Algorithmus-Module
pub mod bezier; // Bézier-Kurven-basierte Glättung
pub mod catmull_rom;
pub mod chaikin; // Chaikin-Algorithmus
pub mod traits; // Enthält den Smoothing-Trait und verwandte Hilfstypen // Catmull-Rom Spline-basierte Glättung
// Hier könnten später weitere Smoothing-Algorithmen hinzukommen,
// z.B. Laplacian Smoothing, etc.

// Re-Exporte der wichtigsten Elemente für die einfache Nutzung des Smoothing-Moduls
pub use self::bezier::{
    BezierSmoother, /*, BezierCurve, BezierType // Falls diese intern bleiben sollen */
};
pub use self::catmull_rom::{
    CatmullRomSmoother, /*, CatmullRomVariants, AdaptiveCatmullRomSmoother // Falls diese Varianten implementiert werden */
};
pub use self::chaikin::ChaikinSmoother;
pub use self::traits::{
    AlgorithmParameters, // Konfigurationsstrukturen
    // SmoothingOutput, SmoothingQualityMetrics, // Ergebnis- und Qualitätsstrukturen (falls QualityControlled implementiert wird)
    Smoothing, // Der Haupt-Trait
    // AdaptiveSmoothing, QualityControlledSmoothing, // Falls diese Traits weiter ausdefiniert werden
    SmoothingParameters,
};

// Optional: Ein Enum, um verschiedene Smoother-Typen zur Laufzeit auszuwählen
/*
#[derive(Debug, Clone)]
pub enum AnySmoother {
    Chaikin(ChaikinSmoother),
    Bezier(BezierSmoother),
    CatmullRom(CatmullRomSmoother),
}

impl Smoothing for AnySmoother {
    fn smooth_points(&self, points: &[bevy::math::Vec2], is_closed: bool) -> crate::math::error::MathResult<Vec<bevy::math::Vec2>> {
        match self {
            AnySmoother::Chaikin(s) => s.smooth_points(points, is_closed),
            AnySmoother::Bezier(s) => s.smooth_points(points, is_closed),
            AnySmoother::CatmullRom(s) => s.smooth_points(points, is_closed),
        }
    }
    // supports_iterations und estimate_iterations müssten ebenfalls implementiert werden
}
*/
