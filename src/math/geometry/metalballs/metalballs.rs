use super::{FieldInfluence, GridConfig, MetaballField};
use bevy::prelude::*;

#[derive(Debug, Clone)]
pub struct Metaballs {
    pub config: GridConfig,
    pub sources: Vec<Arc<dyn FieldInfluence>>,
    pub field: MetaballField,
}

impl Metaballs {
    pub fn with_config(config: GridConfig, sources: Vec<Arc<dyn FieldInfluence>>) -> Self {
        let field = MetaballField::new(config.width, config.height, config.cell_size);
        Self {
            config,
            sources,
            field,
        }
    }

    pub fn recompute_field(&mut self) {
        let (w, h) = (self.config.width, self.config.height);

        for y in 0..h {
            for x in 0..w {
                let p = vec2(
                    x as f32 * self.config.cell_size,
                    y as f32 * self.config.cell_size,
                );
                let val = self.sources.iter().map(|s| s.influence_at(p)).sum();
                self.field.set(x, y, val);
            }
        }
    }
}
