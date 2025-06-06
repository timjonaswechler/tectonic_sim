// src/math/geometry/metalballs/influence/point.rs

use super::FieldInfluence;
use crate::math::types::Bounds2D;
use crate::math::utils::constants;
use bevy::math::Vec2;

/// Erzeugt einen einfachen punktbasierten Einfluss mit linearem Falloff.
#[derive(Debug, Clone)]
pub struct PointInfluence {
    pub position: Vec2,
    pub strength: f32,
    pub radius: f32, // Radius, innerhalb dessen der Einfluss wirkt
    cached_bounds: Option<Bounds2D>,
}

impl PointInfluence {
    pub fn new(position: Vec2, strength: f32, radius: f32) -> Self {
        let r = radius.max(constants::EPSILON);
        let min_coord = position - Vec2::splat(r);
        let max_coord = position + Vec2::splat(r);
        Self {
            position,
            strength,
            radius: r,
            cached_bounds: Bounds2D::new(min_coord, max_coord).ok(),
        }
    }
}

impl FieldInfluence for PointInfluence {
    fn influence_at(&self, point: Vec2) -> f32 {
        if let Some(bounds) = self.cached_bounds {
            if !bounds.contains_point(point) {
                return 0.0;
            }
        }

        let distance_sq = self.position.distance_squared(point);
        let radius_sq = self.radius * self.radius;

        if distance_sq >= radius_sq {
            return 0.0;
        }

        // Einfacher linearer Falloff
        let distance = distance_sq.sqrt();
        self.strength * (1.0 - distance / self.radius).max(0.0)
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        self.cached_bounds
    }

    fn max_influence_distance(&self) -> Option<f32> {
        Some(self.radius)
    }

    fn influence_type(&self) -> &'static str {
        "PointInfluence"
    }
}
