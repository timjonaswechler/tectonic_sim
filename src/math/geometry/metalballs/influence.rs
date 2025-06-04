// src/math/geometry/metalballs/influence.rs

use crate::math::{types::*, utils::*};
use bevy::prelude::*;
use std::sync::Arc;

/// Trait für alle Objekte, die ein Skalarfeld beeinflussen können
pub trait FieldInfluence: Send + Sync + std::fmt::Debug {
    /// Berechnet den Einfluss an einem gegebenen Punkt
    fn influence_at(&self, point: Vec2) -> f32;

    /// Optionale Bounding Box für Optimierungen
    fn bounding_box(&self) -> Option<Bounds2D> {
        None
    }

    /// Prüft ob ein Punkt möglicherweise beeinflusst wird (grobe Schätzung)
    fn might_influence(&self, point: Vec2) -> bool {
        if let Some(bounds) = self.bounding_box() {
            bounds.contains_point(Point2D::new(point.x, point.y))
        } else {
            true // Konservative Schätzung
        }
    }

    /// Gibt die maximale Reichweite des Einflusses zurück
    fn max_influence_distance(&self) -> Option<f32> {
        None
    }

    /// Gibt den Namen/Typ der Influence zurück für Debugging
    fn influence_type(&self) -> &'static str {
        "Unknown"
    }
}

/// Linien-basierter Einfluss
#[derive(Debug, Clone)]
pub struct LineInfluence {
    pub start: Vec2,
    pub end: Vec2,
    pub strength: f32,
    pub width: f32,
    pub falloff_distance: f32,
}

impl LineInfluence {
    pub fn new(start: Vec2, end: Vec2, strength: f32, width: f32) -> Self {
        Self {
            start,
            end,
            strength,
            width,
            falloff_distance: width * 2.0,
        }
    }

    pub fn with_falloff_distance(mut self, distance: f32) -> Self {
        self.falloff_distance = distance;
        self
    }

    /// Berechnet den Abstand zur Linie
    fn distance_to_line(&self, point: Vec2) -> f32 {
        let line_vec = self.end - self.start;
        let point_vec = point - self.start;

        let line_length_sq = line_vec.length_squared();

        if line_length_sq < constants::EPSILON {
            return point.distance(self.start);
        }

        let t = (point_vec.dot(line_vec) / line_length_sq).clamp(0.0, 1.0);
        let projection = self.start + line_vec * t;

        point.distance(projection)
    }
}

impl FieldInfluence for LineInfluence {
    fn influence_at(&self, point: Vec2) -> f32 {
        let distance = self.distance_to_line(point);

        if distance >= self.falloff_distance {
            return 0.0;
        }

        // Lineare Falloff innerhalb der Breite, dann exponentieller Falloff
        if distance <= self.width {
            self.strength
        } else {
            let normalized = (distance - self.width) / (self.falloff_distance - self.width);
            self.strength * (-normalized * 3.0).exp()
        }
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        let min_x = self.start.x.min(self.end.x) - self.falloff_distance;
        let max_x = self.start.x.max(self.end.x) + self.falloff_distance;
        let min_y = self.start.y.min(self.end.y) - self.falloff_distance;
        let max_y = self.start.y.max(self.end.y) + self.falloff_distance;

        Some(Bounds2D::from_points(
            Point2D::new(min_x, min_y),
            Point2D::new(max_x, max_y),
        ))
    }

    fn max_influence_distance(&self) -> Option<f32> {
        Some(self.falloff_distance)
    }

    fn influence_type(&self) -> &'static str {
        "Line"
    }
}

/// Punkt-basierter Einfluss (vereinfachter Metaball)
#[derive(Debug, Clone)]
pub struct PointInfluence {
    pub position: Vec2,
    pub strength: f32,
    pub radius: f32,
}

impl PointInfluence {
    pub fn new(position: Vec2, strength: f32, radius: f32) -> Self {
        Self {
            position,
            strength,
            radius,
        }
    }
}

impl FieldInfluence for PointInfluence {
    fn influence_at(&self, point: Vec2) -> f32 {
        let distance = self.position.distance(point);

        if distance >= self.radius {
            return 0.0;
        }

        // Einfacher linearer Falloff
        self.strength * (1.0 - distance / self.radius)
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        Some(Bounds2D::from_points(
            Point2D::new(self.position.x - self.radius, self.position.y - self.radius),
            Point2D::new(self.position.x + self.radius, self.position.y + self.radius),
        ))
    }

    fn max_influence_distance(&self) -> Option<f32> {
        Some(self.radius)
    }

    fn influence_type(&self) -> &'static str {
        "Point"
    }
}

/// Kombinierter Einfluss aus mehreren Quellen
#[derive(Debug)]
pub struct CompositeInfluence {
    influences: Vec<Arc<dyn FieldInfluence>>,
    combination_mode: CombinationMode,
}

/// Verschiedene Modi zum Kombinieren von Einflüssen
#[derive(Debug, Clone, Copy)]
pub enum CombinationMode {
    /// Addiere alle Einflüsse
    Add,
    /// Multipliziere alle Einflüsse
    Multiply,
    /// Nimm das Maximum
    Maximum,
    /// Nimm das Minimum
    Minimum,
    /// Gewichtetes Mittel
    Average,
}

impl CompositeInfluence {
    pub fn new(mode: CombinationMode) -> Self {
        Self {
            influences: Vec::new(),
            combination_mode: mode,
        }
    }

    pub fn add_influence<I: FieldInfluence + 'static>(&mut self, influence: I) {
        self.influences.push(Arc::new(influence));
    }

    pub fn with_influence<I: FieldInfluence + 'static>(mut self, influence: I) -> Self {
        self.add_influence(influence);
        self
    }
}

impl FieldInfluence for CompositeInfluence {
    fn influence_at(&self, point: Vec2) -> f32 {
        if self.influences.is_empty() {
            return 0.0;
        }

        let values: Vec<f32> = self
            .influences
            .iter()
            .map(|influence| influence.influence_at(point))
            .collect();

        match self.combination_mode {
            CombinationMode::Add => values.iter().sum(),
            CombinationMode::Multiply => values.iter().product(),
            CombinationMode::Maximum => values.iter().fold(0.0, |acc, &x| acc.max(x)),
            CombinationMode::Minimum => values.iter().fold(f32::INFINITY, |acc, &x| acc.min(x)),
            CombinationMode::Average => values.iter().sum::<f32>() / values.len() as f32,
        }
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        let mut combined_bounds: Option<Bounds2D> = None;

        for influence in &self.influences {
            if let Some(bounds) = influence.bounding_box() {
                if let Some(existing) = combined_bounds {
                    combined_bounds = Some(existing.union(&bounds));
                } else {
                    combined_bounds = Some(bounds);
                }
            }
        }

        combined_bounds
    }

    fn max_influence_distance(&self) -> Option<f32> {
        self.influences
            .iter()
            .filter_map(|i| i.max_influence_distance())
            .fold(None, |acc, dist| {
                Some(match acc {
                    None => dist,
                    Some(current) => current.max(dist),
                })
            })
    }

    fn influence_type(&self) -> &'static str {
        "Composite"
    }
}

/// Hilfsfunktionen für Field Influence
pub struct InfluenceUtils;

impl InfluenceUtils {
    /// Erstellt einen optimierten Composite-Einfluss basierend auf räumlicher Nähe
    pub fn create_spatial_composite(
        influences: Vec<Arc<dyn FieldInfluence>>,
        mode: CombinationMode,
    ) -> CompositeInfluence {
        // TODO: Hier könnte eine räumliche Datenstruktur (Quadtree) implementiert werden
        // für bessere Performance bei vielen Einflüssen

        let mut composite = CompositeInfluence::new(mode);
        for influence in influences {
            composite.influences.push(influence);
        }
        composite
    }

    /// Optimiert eine Liste von Einflüssen durch Removal von überlappenden oder schwachen Einflüssen
    pub fn optimize_influences(
        influences: Vec<Arc<dyn FieldInfluence>>,
        min_strength_threshold: f32,
    ) -> Vec<Arc<dyn FieldInfluence>> {
        influences
            .into_iter()
            .filter(|influence| {
                // Einfache Heuristik: teste Einfluss am eigenen Zentrum
                if let Some(bounds) = influence.bounding_box() {
                    let center = bounds.center();
                    let center_vec2 = Vec2::new(center.x, center.y);
                    influence.influence_at(center_vec2) >= min_strength_threshold
                } else {
                    true // Behalte Einflüsse ohne Bounds
                }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_metaball_influence() {
        let metaball = Metaball::new(Vec2::ZERO, 5.0, 1.0);

        // Am Zentrum sollte maximaler Einfluss sein
        let center_influence = metaball.influence_at(Vec2::ZERO);
        assert!(center_influence > 0.0);

        // Außerhalb des Radius sollte kein Einfluss sein
        let far_influence = metaball.influence_at(Vec2::new(10.0, 0.0));
        assert_eq!(far_influence, 0.0);

        // Bounding box sollte korrekt sein
        let bounds = metaball.bounding_box().unwrap();
        assert_eq!(bounds.min.x, -5.0);
        assert_eq!(bounds.max.x, 5.0);
    }

    #[test]
    fn test_line_influence() {
        let line = LineInfluence::new(Vec2::new(-5.0, 0.0), Vec2::new(5.0, 0.0), 1.0, 1.0);

        // Auf der Linie sollte voller Einfluss sein
        let on_line = line.influence_at(Vec2::ZERO);
        assert_eq!(on_line, 1.0);

        // Weit von der Linie entfernt sollte kein Einfluss sein
        let far_from_line = line.influence_at(Vec2::new(0.0, 10.0));
        assert_eq!(far_from_line, 0.0);
    }

    #[test]
    fn test_composite_influence() {
        let mut composite = CompositeInfluence::new(CombinationMode::Add);
        composite.add_influence(Metaball::new(Vec2::new(-1.0, 0.0), 2.0, 0.5));
        composite.add_influence(Metaball::new(Vec2::new(1.0, 0.0), 2.0, 0.5));

        // Am Punkt zwischen beiden sollte die Summe der Einflüsse sein
        let combined_influence = composite.influence_at(Vec2::ZERO);
        assert!(combined_influence > 0.0);
    }

    #[test]
    fn test_influence_optimization() {
        let influences: Vec<Arc<dyn FieldInfluence>> = vec![
            Arc::new(Metaball::new(Vec2::ZERO, 1.0, 0.1)), // Schwach
            Arc::new(Metaball::new(Vec2::new(5.0, 0.0), 1.0, 1.0)), // Stark
        ];

        let optimized = InfluenceUtils::optimize_influences(influences, 0.5);
        assert_eq!(optimized.len(), 1); // Nur der starke Einfluss sollte bleiben
    }
}
