// src/math/geometry/metalballs/objects/polygon.rs

use crate::math::{
    error::*,
    geometry::polygon::{Polygon, PolygonProperties},
    types::*,
    utils::*,
};
use bevy::prelude::*;

/// Verschiedene Falloff-Funktionen für Polygon-Einfluss
#[derive(Debug, Clone, Copy)]
pub enum FalloffType {
    /// Linearer Falloff
    Linear,
    /// Quadratischer Falloff (1/r²)
    Quadratic,
    /// Exponentieller Falloff
    Exponential,
    /// Gaussian Falloff
    Gaussian,
    /// Konstanter Einfluss (Schwellwert)
    Step,
    /// Sigmoid-Funktion für sanften Übergang
    Sigmoid,
}

/// Eigenschaften für Polygon-Einfluss
#[derive(Debug, Clone)]
pub struct PolygonProperties {
    /// Basis-Stärke des Einflusses
    pub strength: f32,
    /// Falloff-Distanz
    pub falloff_distance: f32,
    /// Art des Falloffs
    pub falloff_type: FalloffType,
    /// Ob der Einfluss nur innerhalb des Polygons wirkt
    pub interior_only: bool,
    /// Zusätzlicher Offset für den Einfluss
    pub offset: f32,
    /// Skalierungsfaktor für nicht-lineare Effekte
    pub scale_factor: f32,
}

impl Default for PolygonProperties {
    fn default() -> Self {
        Self {
            strength: 1.0,
            falloff_distance: 10.0,
            falloff_type: FalloffType::Quadratic,
            interior_only: false,
            offset: 0.0,
            scale_factor: 1.0,
        }
    }
}

/// Wrapper für Polygon + Einfluss-Eigenschaften
#[derive(Clone)]
pub struct PolygonInfluence {
    pub polygon: Polygon,
    pub properties: PolygonProperties,
    /// Cache für berechnete Bounding Box
    cached_bounds: Option<Bounds2D>,
}

impl PolygonInfluence {
    /// Erstellt eine neue PolygonInfluence
    pub fn new(polygon: Polygon, properties: PolygonProperties) -> Self {
        Self {
            polygon,
            properties,
            cached_bounds: None,
        }
    }

    /// Erstellt eine PolygonInfluence mit Standard-Eigenschaften
    pub fn with_strength(polygon: Polygon, strength: f32) -> Self {
        let properties = PolygonProperties {
            strength,
            ..Default::default()
        };
        Self::new(polygon, properties)
    }

    /// Setzt die Falloff-Eigenschaften
    pub fn with_falloff(mut self, distance: f32, falloff_type: FalloffType) -> Self {
        self.properties.falloff_distance = distance;
        self.properties.falloff_type = falloff_type;
        self
    }

    /// Aktiviert Interior-Only Modus
    pub fn interior_only(mut self) -> Self {
        self.properties.interior_only = true;
        self
    }

    /// Berechnet die Bounding Box (mit Caching)
    fn get_bounds(&mut self) -> Bounds2D {
        if let Some(bounds) = self.cached_bounds {
            return bounds;
        }

        if let Some((min, max)) = self.polygon.bounds() {
            let bounds = Bounds2D::from_points(min, max);
            // Erweitere um Falloff-Distanz
            let expanded = bounds.expand(self.properties.falloff_distance);
            self.cached_bounds = Some(expanded);
            expanded
        } else {
            // Fallback für leere Polygone
            Bounds2D::from_points(Point2D::ZERO, Point2D::ZERO)
        }
    }

    /// Berechnet den minimalen Abstand von einem Punkt zum Polygon
    pub fn distance_to_polygon(&self, point: Vec2) -> f32 {
        let point_2d = Point2D::new(point.x, point.y);

        // Prüfe zuerst ob der Punkt innerhalb des Polygons liegt
        if self.polygon.contains_point(point_2d) {
            if self.properties.interior_only {
                return 0.0; // Innerhalb = kein Abstand
            } else {
                // Berechne Abstand zur nächsten Kante (negativ für innen)
                return -self.distance_to_boundary(point_2d);
            }
        }

        // Punkt ist außerhalb - berechne Abstand zur nächsten Kante
        self.distance_to_boundary(point_2d)
    }

    /// Berechnet den Abstand zur Polygon-Grenze
    fn distance_to_boundary(&self, point: Point2D) -> f32 {
        let vertices = self.polygon.vertices();
        let mut min_distance = f32::INFINITY;

        for i in 0..vertices.len() {
            let j = (i + 1) % vertices.len();
            let edge_start = vertices[i];
            let edge_end = vertices[j];

            let distance = self.point_to_line_segment_distance(point, edge_start, edge_end);
            min_distance = min_distance.min(distance);
        }

        min_distance
    }

    /// Berechnet den Abstand von einem Punkt zu einem Liniensegment
    fn point_to_line_segment_distance(
        &self,
        point: Point2D,
        line_start: Point2D,
        line_end: Point2D,
    ) -> f32 {
        let line_vec = line_end - line_start;
        let point_vec = point - line_start;

        let line_length_sq = line_vec.length_squared();

        if line_length_sq < constants::EPSILON {
            // Degenerate line segment
            return point.distance(line_start);
        }

        // Projiziere Punkt auf Linie
        let t = (point_vec.dot(line_vec) / line_length_sq).clamp(0.0, 1.0);
        let projection = line_start + line_vec * t;

        point.distance(projection)
    }

    /// Wendet die Falloff-Funktion an
    fn apply_falloff(&self, distance: f32) -> f32 {
        if distance >= self.properties.falloff_distance {
            return 0.0;
        }

        let normalized_distance = distance / self.properties.falloff_distance;
        let base_strength = self.properties.strength;

        match self.properties.falloff_type {
            FalloffType::Linear => base_strength * (1.0 - normalized_distance),
            FalloffType::Quadratic => base_strength / (distance * distance + 1.0),
            FalloffType::Exponential => {
                base_strength * (-normalized_distance * self.properties.scale_factor).exp()
            }
            FalloffType::Gaussian => {
                let variance = self.properties.falloff_distance * 0.5;
                base_strength * (-(distance * distance) / (2.0 * variance * variance)).exp()
            }
            FalloffType::Step => {
                if distance <= self.properties.falloff_distance * 0.5 {
                    base_strength
                } else {
                    0.0
                }
            }
            FalloffType::Sigmoid => {
                let steepness = self.properties.scale_factor.max(0.1);
                let midpoint = self.properties.falloff_distance * 0.5;
                base_strength / (1.0 + ((distance - midpoint) * steepness).exp())
            }
        }
    }
}

/// Implementation des FieldInfluence Traits
impl crate::math::geometry::metalballs::FieldInfluence for PolygonInfluence {
    fn influence_at(&self, point: Vec2) -> f32 {
        // Frühe Ausstiege für Performance
        let point_2d = Point2D::new(point.x, point.y);

        // Prüfe Bounding Box
        let mut bounds = self.cached_bounds.unwrap_or_else(|| {
            if let Some((min, max)) = self.polygon.bounds() {
                Bounds2D::from_points(min, max).expand(self.properties.falloff_distance)
            } else {
                return 0.0;
            }
        });

        if !bounds.contains_point(point_2d) {
            return 0.0;
        }

        // Berechne Abstand zum Polygon
        let distance = self.distance_to_polygon(point);

        // Wenn Interior-Only und Punkt außerhalb, kein Einfluss
        if self.properties.interior_only && distance > 0.0 {
            return 0.0;
        }

        // Wende Falloff an
        let mut influence = self.apply_falloff(distance.abs());

        // Addiere Offset
        influence += self.properties.offset;

        // Stelle sicher, dass der Wert nicht negativ wird
        influence.max(0.0)
    }
}

/// Spezielle Polygon-Influence Typen
pub struct PolygonShapes;

impl PolygonShapes {
    /// Erstellt einen Kreis als Polygon-Influence
    pub fn circle(center: Vec2, radius: f32, strength: f32, segments: usize) -> PolygonInfluence {
        let mut vertices = Vec::with_capacity(segments);

        for i in 0..segments {
            let angle = (i as f32 / segments as f32) * std::f32::consts::TAU;
            let x = center.x + radius * angle.cos();
            let y = center.y + radius * angle.sin();
            vertices.push(Point2D::new(x, y));
        }

        let polygon = Polygon::closed(vertices).unwrap_or_else(|_| {
            // Fallback: einfaches Dreieck
            Polygon::closed(vec![
                Point2D::new(center.x - radius, center.y),
                Point2D::new(center.x + radius, center.y),
                Point2D::new(center.x, center.y + radius),
            ])
            .unwrap()
        });

        PolygonInfluence::with_strength(polygon, strength)
    }

    /// Erstellt ein Rechteck als Polygon-Influence
    pub fn rectangle(center: Vec2, width: f32, height: f32, strength: f32) -> PolygonInfluence {
        let half_w = width * 0.5;
        let half_h = height * 0.5;

        let vertices = vec![
            Point2D::new(center.x - half_w, center.y - half_h),
            Point2D::new(center.x + half_w, center.y - half_h),
            Point2D::new(center.x + half_w, center.y + half_h),
            Point2D::new(center.x - half_w, center.y + half_h),
        ];

        let polygon = Polygon::closed(vertices).unwrap();
        PolygonInfluence::with_strength(polygon, strength)
    }

    /// Erstellt einen Stern als Polygon-Influence
    pub fn star(
        center: Vec2,
        inner_radius: f32,
        outer_radius: f32,
        points: usize,
        strength: f32,
    ) -> PolygonInfluence {
        let mut vertices = Vec::with_capacity(points * 2);

        for i in 0..(points * 2) {
            let angle = (i as f32 / (points * 2) as f32) * std::f32::consts::TAU;
            let radius = if i % 2 == 0 {
                outer_radius
            } else {
                inner_radius
            };

            let x = center.x + radius * angle.cos();
            let y = center.y + radius * angle.sin();
            vertices.push(Point2D::new(x, y));
        }

        let polygon = Polygon::closed(vertices).unwrap();
        PolygonInfluence::with_strength(polygon, strength)
    }

    /// Erstellt ein regelmäßiges n-Eck
    pub fn regular_polygon(
        center: Vec2,
        radius: f32,
        sides: usize,
        strength: f32,
    ) -> PolygonInfluence {
        if sides < 3 {
            return Self::circle(center, radius, strength, 8);
        }

        let mut vertices = Vec::with_capacity(sides);

        for i in 0..sides {
            let angle = (i as f32 / sides as f32) * std::f32::consts::TAU;
            let x = center.x + radius * angle.cos();
            let y = center.y + radius * angle.sin();
            vertices.push(Point2D::new(x, y));
        }

        let polygon = Polygon::closed(vertices).unwrap();
        PolygonInfluence::with_strength(polygon, strength)
    }

    /// Erstellt eine elliptische Form
    pub fn ellipse(
        center: Vec2,
        width: f32,
        height: f32,
        strength: f32,
        segments: usize,
    ) -> PolygonInfluence {
        let mut vertices = Vec::with_capacity(segments);
        let half_w = width * 0.5;
        let half_h = height * 0.5;

        for i in 0..segments {
            let angle = (i as f32 / segments as f32) * std::f32::consts::TAU;
            let x = center.x + half_w * angle.cos();
            let y = center.y + half_h * angle.sin();
            vertices.push(Point2D::new(x, y));
        }

        let polygon = Polygon::closed(vertices).unwrap();
        PolygonInfluence::with_strength(polygon, strength)
    }
}

/// Builder für komplexe Polygon-Influences
pub struct PolygonInfluenceBuilder {
    polygon: Option<Polygon>,
    properties: PolygonProperties,
}

impl PolygonInfluenceBuilder {
    pub fn new() -> Self {
        Self {
            polygon: None,
            properties: PolygonProperties::default(),
        }
    }

    pub fn polygon(mut self, polygon: Polygon) -> Self {
        self.polygon = Some(polygon);
        self
    }

    pub fn vertices(self, vertices: Vec<Point2D>) -> Self {
        if let Ok(polygon) = Polygon::closed(vertices) {
            self.polygon(polygon)
        } else {
            self
        }
    }

    pub fn strength(mut self, strength: f32) -> Self {
        self.properties.strength = strength;
        self
    }

    pub fn falloff_distance(mut self, distance: f32) -> Self {
        self.properties.falloff_distance = distance;
        self
    }

    pub fn falloff_type(mut self, falloff_type: FalloffType) -> Self {
        self.properties.falloff_type = falloff_type;
        self
    }

    pub fn interior_only(mut self, interior_only: bool) -> Self {
        self.properties.interior_only = interior_only;
        self
    }

    pub fn offset(mut self, offset: f32) -> Self {
        self.properties.offset = offset;
        self
    }

    pub fn scale_factor(mut self, scale_factor: f32) -> Self {
        self.properties.scale_factor = scale_factor;
        self
    }

    pub fn build(self) -> Option<PolygonInfluence> {
        if let Some(polygon) = self.polygon {
            Some(PolygonInfluence::new(polygon, self.properties))
        } else {
            None
        }
    }
}

impl Default for PolygonInfluenceBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polygon_influence_basic() {
        let square = PolygonShapes::rectangle(Vec2::new(0.0, 0.0), 2.0, 2.0, 1.0);

        // Punkt im Zentrum sollte vollen Einfluss haben
        let center_influence = square.influence_at(Vec2::new(0.0, 0.0));
        assert!(center_influence > 0.0);

        // Punkt weit entfernt sollte keinen Einfluss haben
        let far_influence = square.influence_at(Vec2::new(100.0, 100.0));
        assert_eq!(far_influence, 0.0);
    }

    #[test]
    fn test_falloff_types() {
        let circle =
            PolygonShapes::circle(Vec2::ZERO, 1.0, 1.0, 8).with_falloff(5.0, FalloffType::Linear);

        let linear_influence = circle.influence_at(Vec2::new(2.0, 0.0));

        let circle_exp = PolygonShapes::circle(Vec2::ZERO, 1.0, 1.0, 8)
            .with_falloff(5.0, FalloffType::Exponential);

        let exp_influence = circle_exp.influence_at(Vec2::new(2.0, 0.0));

        // Beide sollten Einfluss haben, aber unterschiedlich stark
        assert!(linear_influence > 0.0);
        assert!(exp_influence > 0.0);
        assert_ne!(linear_influence, exp_influence);
    }

    #[test]
    fn test_interior_only_mode() {
        let square = PolygonShapes::rectangle(Vec2::ZERO, 2.0, 2.0, 1.0).interior_only();

        // Punkt innerhalb sollte Einfluss haben
        let inside_influence = square.influence_at(Vec2::new(0.5, 0.5));
        assert!(inside_influence > 0.0);

        // Punkt außerhalb sollte keinen Einfluss haben
        let outside_influence = square.influence_at(Vec2::new(2.0, 2.0));
        assert_eq!(outside_influence, 0.0);
    }

    #[test]
    fn test_polygon_shapes() {
        let circle = PolygonShapes::circle(Vec2::ZERO, 1.0, 1.0, 12);
        let star = PolygonShapes::star(Vec2::ZERO, 0.5, 1.0, 5, 1.0);
        let hexagon = PolygonShapes::regular_polygon(Vec2::ZERO, 1.0, 6, 1.0);

        // Alle sollten am Zentrum Einfluss haben
        assert!(circle.influence_at(Vec2::ZERO) > 0.0);
        assert!(star.influence_at(Vec2::ZERO) > 0.0);
        assert!(hexagon.influence_at(Vec2::ZERO) > 0.0);
    }

    #[test]
    fn test_builder_pattern() {
        let polygon_influence = PolygonInfluenceBuilder::new()
            .vertices(vec![
                Point2D::new(-1.0, -1.0),
                Point2D::new(1.0, -1.0),
                Point2D::new(1.0, 1.0),
                Point2D::new(-1.0, 1.0),
            ])
            .strength(2.0)
            .falloff_distance(5.0)
            .falloff_type(FalloffType::Gaussian)
            .build()
            .unwrap();

        assert_eq!(polygon_influence.properties.strength, 2.0);
        assert_eq!(polygon_influence.properties.falloff_distance, 5.0);
        assert!(matches!(
            polygon_influence.properties.falloff_type,
            FalloffType::Gaussian
        ));
    }
}
