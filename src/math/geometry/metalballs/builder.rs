// src/math/geometry/metalballs/builder.rs

use crate::math::geometry::metalballs::{
    components::GridConfig,
    influence::FieldInfluence, // Der Trait
    influence::MetaballSource, // Die konkrete Metaball-Struktur
    metalballs::Metaballs,     // Die zu bauende Struktur
};
use std::sync::Arc;

/// Builder zum komfortablen Erstellen und Konfigurieren von `Metaballs`-Instanzen.
#[derive(Debug, Clone)] // Clone ist möglich, da Vec<Arc<...>> cloneable ist
pub struct MetaballsBuilder {
    pub config: GridConfig,
    sources: Vec<Arc<dyn FieldInfluence>>,
}

impl Default for MetaballsBuilder {
    fn default() -> Self {
        Self {
            config: GridConfig::default(),
            sources: Vec::new(),
        }
    }
}

impl MetaballsBuilder {
    /// Erstellt einen neuen `MetaballsBuilder` mit Standardkonfiguration.
    pub fn new() -> Self {
        Self::default()
    }

    /// Setzt die Breite und Höhe des Grids in Zellen.
    pub fn grid_dimensions(mut self, width: usize, height: usize) -> Self {
        self.config.width = width.max(1); // Mindestens 1x1 Grid
        self.config.height = height.max(1);
        self
    }

    /// Setzt die Größe einer einzelnen Zelle in Weltkoordinaten.
    pub fn cell_size(mut self, size: f32) -> Self {
        self.config.cell_size = size.max(crate::math::utils::constants::EPSILON); // Muss positiv sein
        self
    }

    /// Setzt den Schwellenwert für die Iso-Kontur-Extraktion.
    pub fn threshold(mut self, threshold_value: f32) -> Self {
        self.config.threshold = threshold_value;
        self
    }

    /// Setzt die gesamte Grid-Konfiguration.
    pub fn with_config(mut self, config: GridConfig) -> Self {
        self.config = config;
        self
    }

    /// Fügt eine beliebige Feld-Einflussquelle hinzu.
    /// Die Quelle muss das `FieldInfluence`-Trait implementieren und `'static` sein,
    /// um in einem `Arc<dyn ...>` gespeichert werden zu können.
    pub fn add_source<I: FieldInfluence + 'static>(mut self, source: I) -> Self {
        self.sources.push(Arc::new(source));
        self
    }

    /// Fügt eine bereits in `Arc` gepackte Feld-Einflussquelle hinzu.
    pub fn add_source_arc(mut self, source_arc: Arc<dyn FieldInfluence>) -> Self {
        self.sources.push(source_arc);
        self
    }

    /// (Beibehalten für API-Kompatibilität, obwohl `add_source` generischer ist)
    /// Fügt einen spezifischen `Metaball` als Einflussquelle hinzu.
    pub fn add_ball(self, ball: MetaballSource) -> Self {
        self.add_source(ball)
    }

    /// Fügt mehrere Quellen auf einmal hinzu.
    pub fn add_sources<I>(mut self, sources_iter: I) -> Self
    where
        I: IntoIterator<Item = Arc<dyn FieldInfluence>>,
    {
        self.sources.extend(sources_iter);
        self
    }

    /// Erstellt die `Metaballs`-Instanz mit der aktuellen Konfiguration und den hinzugefügten Quellen.
    pub fn build(self) -> Metaballs {
        Metaballs::from_builder_data(self.config, self.sources)
    }
}
