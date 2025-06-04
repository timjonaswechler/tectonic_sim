use bevy::prelude::*;

/// Basis-Metaball Implementierung
#[derive(Debug, Clone)]
pub struct Metaball {
    pub position: Vec2,
    pub radius: f32,
    pub strength: f32,
    pub falloff_type: MetaballFalloff,
}

/// Verschiedene Falloff-Funktionen für Metaballs
#[derive(Debug, Clone, Copy)]
pub enum MetaballFalloff {
    /// 1/r² Falloff
    InverseSquare,
    /// Exponentieller Falloff
    Exponential,
    /// Gaussian Falloff
    Gaussian,
    /// Linearer Falloff
    Linear,
    /// Soft-Step Falloff
    SmoothStep,
    /// Blobby Falloff (Wyvill)
    Wyvill,
}

impl Metaball {
    pub fn new(position: Vec2, radius: f32, strength: f32) -> Self {
        Self {
            position,
            radius,
            strength,
            falloff_type: MetaballFalloff::InverseSquare,
        }
    }
    /// Einflussfunktion ∝ strength / (dist² + 1)
    pub fn influence_at(&self, p: Vec2) -> f32 {
        let d = self.position.distance(p);
        self.strength / (d * d + 1.0)
    }
    pub fn with_falloff(mut self, falloff: MetaballFalloff) -> Self {
        self.falloff_type = falloff;
        self
    }

    /// Berechnet den Falloff basierend auf dem Typ
    fn calculate_falloff(&self, distance: f32) -> f32 {
        if distance >= self.radius {
            return 0.0;
        }

        let normalized_distance = distance / self.radius;

        match self.falloff_type {
            MetaballFalloff::InverseSquare => self.strength / (distance * distance + 1.0),
            MetaballFalloff::Exponential => self.strength * (-normalized_distance * 4.0).exp(),
            MetaballFalloff::Gaussian => {
                self.strength * (-(normalized_distance * normalized_distance) * 2.0).exp()
            }
            MetaballFalloff::Linear => self.strength * (1.0 - normalized_distance),
            MetaballFalloff::SmoothStep => {
                let t = 1.0 - normalized_distance;
                self.strength * t * t * (3.0 - 2.0 * t)
            }
            MetaballFalloff::Wyvill => {
                // Wyvill falloff function: (1 - (r/R)^6)^2
                let t = normalized_distance;
                let t6 = t.powi(6);
                let falloff = (1.0 - t6).max(0.0);
                self.strength * falloff * falloff
            }
        }
    }
}

impl FieldInfluence for Metaball {
    fn influence_at(&self, point: Vec2) -> f32 {
        let distance = self.position.distance(point);
        self.calculate_falloff(distance)
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        let min = Point2D::new(self.position.x - self.radius, self.position.y - self.radius);
        let max = Point2D::new(self.position.x + self.radius, self.position.y + self.radius);
        Some(Bounds2D::from_points(min, max))
    }

    fn max_influence_distance(&self) -> Option<f32> {
        Some(self.radius)
    }

    fn influence_type(&self) -> &'static str {
        "Metaball"
    }
}
