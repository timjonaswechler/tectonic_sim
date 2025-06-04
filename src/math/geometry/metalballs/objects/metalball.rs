use bevy::prelude::*;

#[derive(Debug, Clone)]
pub struct Metaball {
    pub position: Vec2,
    pub radius: f32,
    pub strength: f32,
}

impl Metaball {
    pub fn new<P: Into<Vec2>>(position: P, radius: f32, strength: f32) -> Self {
        Self {
            position: position.into(),
            radius,
            strength,
        }
    }

    /// Einflussfunktion ∝ strength / (dist² + 1)
    pub fn influence_at(&self, p: Vec2) -> f32 {
        let d = self.position.distance(p);
        self.strength / (d * d + 1.0)
    }
}
