// src/math/algorithms/smoothing/mod.rs

//! Modul für verschiedene Algorithmen zur Glättung von Linienzügen und Polygonen.
//!
//! Dieses Modul stellt Implementierungen von Glättungsalgorithmen sowie die notwendigen
//! Traits und Hilfsstrukturen bereit.
//!
//! # Verfügbare Algorithmen:
//!
//! - [`BezierSmoother`]: Glättung basierend auf kubischen Bézier-Kurven.
//! - [`CatmullRomSmoother`]: Glättung basierend auf Catmull-Rom-Splines.
//! - [`ChaikinSmoother`]: Glättung basierend auf dem iterativen Chaikin-Algorithmus.
//!
//! # Haupt-Trait:
//!
//! - [`Smoothing`]: Definiert die Basisschnittstelle für alle Glättungsalgorithmen.
//!
//! # Verwendung
//!
//! Die Smoother können instanziiert und konfiguriert werden. Anschließend wird die
//! Methode `smooth_points` (aus dem `Smoothing`-Trait) aufgerufen, um einen Linienzug
//! zu glätten, oder `smooth_polygon` für ganze Polygone.
//!
//! ```rust,ignore
//! use bevy::math::Vec2;
//! use your_crate::math::algorithms::smoothing::{ChaikinSmoother, Smoothing};
//! use your_crate::math::geometry::polygon::Polygon; // Angenommen, Polygon ist verfügbar
//!
//! let points = vec![Vec2::new(0.0, 0.0), Vec2::new(10.0, 10.0), Vec2::new(20.0, 0.0)];
//! let mut polygon = Polygon::new(points.clone()).unwrap(); // Beispiel
//!
//! // Chaikin-Glättung
//! let chaikin_smoother = ChaikinSmoother::new(2).with_cut_ratio(0.25);
//! let smoothed_points_chaikin = chaikin_smoother.smooth_points(&points, false).unwrap();
//! let smoothed_polygon_chaikin = chaikin_smoother.smooth_polygon(&polygon).unwrap();
//!
//! // Weitere Smoother könnten ähnlich verwendet werden...
//! ```

// Deklaration der verschiedenen Smoothing-Algorithmus-Module
pub mod bezier;
pub mod catmull_rom;
pub mod chaikin;
pub mod traits;

// Re-Exporte der wichtigsten Elemente für die einfache Nutzung des Smoothing-Moduls
pub use self::bezier::BezierSmoother;
pub use self::catmull_rom::CatmullRomSmoother;
pub use self::chaikin::ChaikinSmoother;
pub use self::traits::{
    AlgorithmParameters, // Konfigurationsstrukturen für spezifische Algorithmen
    Smoothing,           // Der Haupt-Trait für Glättungsoperationen
    SmoothingParameters, // Allgemeine Parameter für Glättungsoperationen
};

#[derive(Debug, Clone)]
pub enum AnySmoother {
    Chaikin(ChaikinSmoother),
    Bezier(BezierSmoother),
    CatmullRom(CatmullRomSmoother),
}

impl Smoothing for AnySmoother {
    fn smooth_points(
        &self,
        points: &[bevy::math::Vec2],
        is_closed: bool,
    ) -> crate::math::error::MathResult<Vec<bevy::math::Vec2>> {
        match self {
            AnySmoother::Chaikin(s) => s.smooth_points(points, is_closed),
            AnySmoother::Bezier(s) => s.smooth_points(points, is_closed),
            AnySmoother::CatmullRom(s) => s.smooth_points(points, is_closed),
        }
    }

    fn supports_iterations(&self) -> bool {
        match self {
            AnySmoother::Chaikin(s) => s.supports_iterations(),
            AnySmoother::Bezier(s) => s.supports_iterations(),
            AnySmoother::CatmullRom(s) => s.supports_iterations(),
        }
    }

    fn estimate_iterations_for_polygon(
        &self,
        polygon: &crate::math::geometry::polygon::core::Polygon,
    ) -> usize {
        match self {
            AnySmoother::Chaikin(s) => s.estimate_iterations_for_polygon(polygon),
            AnySmoother::Bezier(s) => s.estimate_iterations_for_polygon(polygon),
            AnySmoother::CatmullRom(s) => s.estimate_iterations_for_polygon(polygon),
        }
    }
}
