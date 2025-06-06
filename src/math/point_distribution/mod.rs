// src/math/point_distribution/mod.rs

// Deklaration der verschiedenen Punktverteilungs-Methoden/Module
pub mod voronoi;
// Hier könnten später andere Punktverteilungsalgorithmen hinzukommen,
// z.B. Poisson Disk Sampling (als eigenständiger Algorithmus, nicht nur für Kugel),
// Grid-basierte Verteilungen etc.

// Re-Exporte der wichtigsten Elemente aus den Untermodulen
pub use self::voronoi::{
    VoronoiBuilder,
    // LloydConfig, BoundaryHandling, // Bereits durch VoronoiConfig indirekt abgedeckt oder spezifisch im Builder verwendet
    VoronoiCell,
    VoronoiConfig, // Konfiguration für den gesamten Voronoi-Prozess
    VoronoiDiagram,
    // PolygonMerger, MergeStrategy, // Sind eher interne Komponenten des Builders
    VoronoiEdge,
};

// Wenn LloydRelaxation als allgemeiner Algorithmus betrachtet wird,
// wird er über `math::algorithms::prelude` oder direkt `math::algorithms::point_relaxation` importiert.
// Wenn er spezifisch für die Punktverteilung hier relevant ist, könnte man ihn auch hier re-exportieren,
// aber das würde die Trennung zu `algorithms` etwas aufweichen.
