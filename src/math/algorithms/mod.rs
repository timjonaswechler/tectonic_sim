// src/math/algorithms/mod.rs

// Deklaration der verschiedenen Algorithmus-Kategorien
pub mod boolean;
pub mod clipping;
pub mod convex_hull;
pub mod smoothing;
pub mod triangulation; // Für Boolean-Operationen auf Polygon-Punktlisten

// Platzhalter für zukünftige Algorithmen-Module
pub mod point_relaxation;
// pub mod pathfinding;
// pub mod optimization;

// Re-Exporte für den direkten Zugriff auf die Algorithmus-Module
// oder spezifische Algorithmen und deren Konfigurationen.

// Smoothing Algorithmen und Traits
pub use self::smoothing::{
    AlgorithmParameters as SmoothingAlgorithmParameters, // Spezifische Parameter für Smoother
    BezierSmoother,
    CatmullRomSmoother,
    ChaikinSmoother,
    Smoothing,           // Haupt-Trait
    SmoothingParameters, // Allgemeine Parameter (falls noch verwendet)
};

// Konvexe Hülle
pub use self::convex_hull::{
    ConvexHullAlgorithm,
    // ConvexityUtils und ConvexityDefect wurden entfernt, da sie eher Polygon-Analyse sind
    ConvexHullComputer,
};

// Triangulation
pub use self::triangulation::{
    PolygonTriangulator,
    Triangle, // Die Dreiecksstruktur, die von Triangulatoren verwendet wird
    TriangulationAlgorithm,
    TriangulationUtils,
};

// Clipping
pub use self::clipping::{
    ClippingAlgorithm,
    ClippingOperations, // Utility-Funktionen
    // ClippingStats, ClippingAnalysis, // Falls benötigt und implementiert
    PolygonClipper,
};

// Boolean-Operationen
pub use self::boolean::{
    BooleanOpType, // Enum für die Operationsart
    // BooleanOperations (Utility-Struct) wurde entfernt, da PolygonBooleanOps dies übernimmt
    // BooleanAnalysis, BooleanStats, // Falls benötigt und implementiert
    PolygonBooleanOps, // Umbenannt von PolygonBoolean zur Unterscheidung von Polygon-Struktur
};

pub use self::point_relaxation::{
    BoundaryHandling, // Enum für Grenzbehandlung
                      // LloydStatistics, // Wird in lloyd.rs definiert und hier re-exportiert
    LloydConfig, // Die Konfiguration für Lloyd
    LloydRelaxation,
};
