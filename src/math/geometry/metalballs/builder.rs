use super::{GridConfig, influence::FieldInfluence};

#[derive(Debug, Clone)]
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
    pub fn new() -> Self {
        Self::default()
    }

    pub fn grid_size(mut self, w: usize, h: usize) -> Self {
        self.config.width = w;
        self.config.height = h;
        self
    }

    pub fn cell_size(mut self, size: f32) -> Self {
        self.config.cell_size = size;
        self
    }

    pub fn threshold(mut self, t: f32) -> Self {
        self.config.threshold = t;
        self
    }

    /// Fügt eine beliebige Feld-Quelle hinzu.
    pub fn add_source<I: FieldInfluence + 'static>(mut self, s: I) -> Self {
        self.sources.push(Arc::new(s));
        self
    }

    /// Komfort-Alias speziell für Kreise – bricht alte Nutzer-API nicht.
    pub fn add_ball(mut self, ball: crate::metaball::Metaball) -> Self {
        self.sources.push(Arc::new(ball));
        self
    }

    pub fn build(self) -> crate::metaballs::Metaballs {
        crate::metaballs::Metaballs::with_config(self.config, self.sources)
    }
}
