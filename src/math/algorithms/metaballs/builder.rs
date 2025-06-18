// src/math/algorithms/metaballs/builder.rs

use crate::math::algorithms::metaballs::{
    influence::{FieldInfluence, MetaballSource}, // MetaballSource nur für add_ball
    metaballs::Metaballs,
};
use crate::math::utils::constants; // Für EPSILON
use std::sync::Arc;

/// Builder zum komfortablen Erstellen und Konfigurieren von `Metaballs`-Instanzen.
#[derive(Debug, Clone)]
pub struct MetaballsBuilder {
    // Felder von GridConfig direkt hier (außer threshold, der geht in Metaballs struct oder bleibt hier für den Builder)
    // oder wir behalten config hier für Einfachheit und der Builder kopiert es in Metaballs
    width: usize,
    height: usize,
    cell_size: f32,
    threshold: f32, // Threshold ist nun ein direktes Konfigurationsfeld des Builders
    sources: Vec<Arc<dyn FieldInfluence>>,
}

impl Default for MetaballsBuilder {
    fn default() -> Self {
        Self {
            width: 100,
            height: 100,
            cell_size: 1.0,
            threshold: 0.5,
            sources: Vec::new(),
        }
    }
}

impl MetaballsBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn grid_dimensions(mut self, width: usize, height: usize) -> Self {
        self.width = width.max(1);
        self.height = height.max(1);
        self
    }

    pub fn cell_size(mut self, size: f32) -> Self {
        self.cell_size = size.max(constants::EPSILON);
        self
    }

    pub fn threshold(mut self, threshold_value: f32) -> Self {
        self.threshold = threshold_value;
        self
    }

    // with_config wird nicht mehr benötigt, da es keine GridConfig-Struct mehr gibt

    pub fn add_source<I: FieldInfluence + 'static>(mut self, source: I) -> Self {
        self.sources.push(Arc::new(source));
        self
    }

    pub fn add_source_arc(mut self, source_arc: Arc<dyn FieldInfluence>) -> Self {
        self.sources.push(source_arc);
        self
    }

    pub fn add_ball(self, ball: MetaballSource) -> Self {
        self.add_source(ball)
    }

    pub fn add_sources<I>(mut self, sources_iter: I) -> Self
    where
        I: IntoIterator<Item = Arc<dyn FieldInfluence>>,
    {
        self.sources.extend(sources_iter);
        self
    }

    pub fn build(self) -> Metaballs {
        // Übergebe die Felder direkt an Metaballs::new oder eine from_builder_data
        Metaballs::from_builder_data(
            self.width,
            self.height,
            self.cell_size,
            self.threshold,
            self.sources,
        )
    }
}
