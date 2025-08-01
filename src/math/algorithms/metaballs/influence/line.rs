use super::FieldInfluence;
use crate::math::types::Bounds2D;
use crate::math::utils::constants; // Für EPSILON und EPSILON_SQUARED
use bevy::math::Vec2;

/// Erzeugt einen linearen Einflussbereich zwischen zwei Punkten.
#[derive(Debug, Clone)]
pub struct LineInfluence {
    pub start: Vec2,
    pub end: Vec2,
    pub strength: f32,
    /// Breite der Linie, innerhalb derer die volle Stärke (oder ein Teil davon) wirkt.
    pub width: f32,
    /// Distanz vom Rand der `width`, über die der Einfluss auf Null abfällt.
    pub falloff_distance: f32,
    cached_bounds_expanded: Option<Bounds2D>,
}

impl LineInfluence {
    pub fn new(start: Vec2, end: Vec2, strength: f32, width: f32) -> Self {
        let mut new_self = Self {
            start,
            end,
            strength,
            width: width.max(0.0),
            // Falloff sollte mind. so groß wie Breite sein, oder einen sinnvollen Mindestwert haben
            falloff_distance: (width * 0.5).max(constants::EPSILON), // Angepasst: Falloff kann kleiner als Breite sein
            cached_bounds_expanded: None,
        };
        new_self.recalculate_cached_bounds();
        new_self
    }

    pub fn with_falloff_distance(mut self, distance: f32) -> Self {
        // Erlaube beliebige positive Falloff-Distanz
        self.falloff_distance = distance.max(constants::EPSILON);
        self.recalculate_cached_bounds();
        self
    }

    fn recalculate_cached_bounds(&mut self) {
        let total_reach = self.width * 0.5 + self.falloff_distance;
        let min_x = self.start.x.min(self.end.x) - total_reach;
        let max_x = self.start.x.max(self.end.x) + total_reach;
        let min_y = self.start.y.min(self.end.y) - total_reach;
        let max_y = self.start.y.max(self.end.y) + total_reach;

        self.cached_bounds_expanded =
            Bounds2D::new(Vec2::new(min_x, min_y), Vec2::new(max_x, max_y)).ok();
    }

    /// Berechnet den kürzesten Abstand von einem Punkt zum Liniensegment.
    fn distance_to_line_segment(&self, point: Vec2) -> f32 {
        let line_vec = self.end - self.start;
        let point_vec = point - self.start;

        let line_length_sq = line_vec.length_squared();

        if line_length_sq < constants::EPSILON_SQUARED {
            // Verwende EPSILON_SQUARED
            return point.distance(self.start);
        }

        let t = (point_vec.dot(line_vec) / line_length_sq).clamp(0.0, 1.0);
        let projection_on_segment = self.start + line_vec * t;

        point.distance(projection_on_segment)
    }
}

impl FieldInfluence for LineInfluence {
    fn influence_at(&self, point: Vec2) -> f32 {
        if let Some(bounds) = self.cached_bounds_expanded {
            if !bounds.contains_point(point) {
                return 0.0;
            }
        }

        let distance = self.distance_to_line_segment(point);
        let half_width = self.width * 0.5;

        if distance <= half_width {
            self.strength
        } else if distance <= half_width + self.falloff_distance {
            // Stelle sicher, dass falloff_distance positiv ist, um Division durch Null zu vermeiden
            if self.falloff_distance < constants::EPSILON {
                // Wenn kein Falloff-Bereich, aber außerhalb der Kernbreite -> 0
                return 0.0;
            }
            let dist_into_falloff = distance - half_width;
            let normalized_falloff_dist = dist_into_falloff / self.falloff_distance;
            self.strength * (1.0 - normalized_falloff_dist).max(0.0)
        } else {
            0.0
        }
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        self.cached_bounds_expanded
    }

    fn max_influence_distance(&self) -> Option<f32> {
        Some(self.width * 0.5 + self.falloff_distance)
    }

    fn influence_type(&self) -> &'static str {
        "LineInfluence"
    }
}
