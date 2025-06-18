// src/math/point_distribution/voronoi/mod.rs

// Deklaration der Untermodule für Voronoi-spezifische Funktionalität
pub mod builder;
pub mod config; // Enthält VoronoiConfig (und LloydConfig, falls es dort bleibt)
pub mod error;
pub mod merger;
pub mod voronoi_diagram; // Enthält VoronoiCell, VoronoiEdge, VoronoiDiagram, VoronoiExtractor

// Re-Exporte für den einfachen Zugriff auf die wichtigsten Voronoi-Elemente
pub use self::builder::VoronoiBuilder;
pub use self::config::{BoundaryHandling, LloydConfig, VoronoiConfig}; // LloydConfig und BoundaryHandling hier, da eng mit Voronoi-Prozess verbunden
pub use self::merger::{MergeStrategy, PolygonMerger};
pub use self::voronoi_diagram::{VoronoiCell, VoronoiDiagram, VoronoiEdge, VoronoiExtractor};
