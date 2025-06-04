// src/math/geometry/metalballs/objects/polygon.rs

use crate::math::{
    geometry::polygon::{Polygon /*, PolygonProperties // Entfernt */},
    types::*,
    utils::*, // Für Bounds2D, Point2D, constants::EPSILON
};
use bevy::prelude::*; // Für Vec2

/// Falloff-Typen für Polygon-Einfluss
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PolygonFalloff {
    Linear,
    Quadratic,
    Exponential,
    Gaussian,
    Step,
    Sigmoid,
}

impl Default for PolygonFalloff {
    fn default() -> Self {
        PolygonFalloff::Quadratic
    }
}

/// Polygon-basierter Einfluss
#[derive(Clone)]
pub struct PolygonInfluence {
    pub polygon: Polygon,
    pub strength: f32,
    pub falloff_distance: f32,
    pub falloff_type: PolygonFalloff,
    pub interior_only: bool,
    pub offset: f32,       // Hinzugefügt basierend auf alter influence_at Logik
    pub scale_factor: f32, // Hinzugefügt basierend auf alter apply_falloff Logik
    cached_bounds: Option<Bounds2D>,
}

impl PolygonInfluence {
    /// Privater Konstruktor, der vom Builder verwendet wird.
    fn new(
        polygon: Polygon,
        strength: f32,
        falloff_distance: f32,
        falloff_type: PolygonFalloff,
        interior_only: bool,
        offset: f32,
        scale_factor: f32,
    ) -> Self {
        Self {
            polygon,
            strength,
            falloff_distance,
            falloff_type,
            interior_only,
            offset,
            scale_factor,
            cached_bounds: None, // Cache wird bei Bedarf gefüllt
        }
    }

    /// Berechnet die Bounding Box (mit Caching).
    /// Diese Methode benötigt &mut self, um den Cache zu aktualisieren.
    pub fn get_bounds_mut(&mut self) -> Bounds2D {
        if let Some(bounds) = self.cached_bounds {
            return bounds;
        }

        if let Some((min, max)) = self.polygon.bounds() {
            let bounds = Bounds2D::from_points(min, max);
            // Erweitere um Falloff-Distanz
            let expanded = bounds.expand(self.falloff_distance);
            self.cached_bounds = Some(expanded);
            expanded
        } else {
            // Fallback für leere Polygone
            let fallback_bounds = Bounds2D::from_points(Point2D::ZERO, Point2D::ZERO);
            self.cached_bounds = Some(fallback_bounds);
            fallback_bounds
        }
    }

    /// Berechnet den minimalen Abstand von einem Punkt zum Polygon.
    /// Positiv für außen, negativ für innen (wenn nicht interior_only), 0.0 wenn interior_only und innen.
    fn distance_to_polygon(&self, point: Vec2) -> f32 {
        let point_2d = Point2D::new(point.x, point.y);

        if self.polygon.contains_point(point_2d) {
            if self.interior_only {
                0.0
            } else {
                // Negativer Abstand zur nächsten Kante für Punkte innerhalb
                -self.distance_to_boundary(point_2d)
            }
        } else {
            // Positiver Abstand zur nächsten Kante für Punkte außerhalb
            self.distance_to_boundary(point_2d)
        }
    }

    /// Berechnet den Abstand zur Polygon-Grenze (immer positiv).
    fn distance_to_boundary(&self, point: Point2D) -> f32 {
        let vertices = self.polygon.vertices();
        if vertices.is_empty() {
            return f32::INFINITY;
        }
        if vertices.len() == 1 {
            return point.distance(vertices[0]);
        }

        let mut min_distance_sq = f32::INFINITY;

        for i in 0..vertices.len() {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % vertices.len()];
            min_distance_sq =
                min_distance_sq.min(self.point_to_line_segment_distance_sq(point, p1, p2));
        }
        min_distance_sq.sqrt()
    }

    /// Berechnet den quadratischen Abstand von einem Punkt zu einem Liniensegment.
    /// Effizienter, da sqrt erst am Ende in distance_to_boundary aufgerufen wird.
    fn point_to_line_segment_distance_sq(
        &self,
        point: Point2D,
        line_start: Point2D,
        line_end: Point2D,
    ) -> f32 {
        let line_vec = line_end - line_start;
        let point_vec = point - line_start;

        let line_length_sq = line_vec.length_squared();

        if line_length_sq < constants::EPSILON * constants::EPSILON {
            // Vergleiche mit quadriertem Epsilon
            // Degeneriertes Liniensegment (Start und Ende sind quasi gleich)
            return point.distance_squared(line_start);
        }

        // Projiziere Punkt auf die Linie, die das Segment enthält
        // t ist der Parameter der Projektion entlang des Segments
        // t = dot(point_vec, line_vec) / |line_vec|^2
        let t = (point_vec.dot(line_vec) / line_length_sq).clamp(0.0, 1.0);

        // Der nächstgelegene Punkt auf dem Segment ist line_start + t * line_vec
        let projection = line_start + line_vec * t;

        point.distance_squared(projection)
    }

    /// Wendet die Falloff-Funktion an.
    fn apply_falloff(&self, distance: f32) -> f32 {
        // distance hier ist immer positiv (abs() wurde in influence_at aufgerufen)
        if distance >= self.falloff_distance {
            return 0.0;
        }
        if self.falloff_distance < constants::EPSILON {
            // Vermeide Division durch Null
            return if distance < constants::EPSILON {
                self.strength
            } else {
                0.0
            };
        }

        let normalized_distance = distance / self.falloff_distance;

        match self.falloff_type {
            PolygonFalloff::Linear => self.strength * (1.0 - normalized_distance),
            PolygonFalloff::Quadratic => {
                // Vermeide Division durch fast 0, wenn distance sehr klein ist
                // Original war: self.strength / (distance * distance + 1.0)
                // Dies skaliert mit der Distanz, nicht normalisiert.
                // Eine normalisierte Variante könnte sein:
                self.strength * (1.0 - normalized_distance * normalized_distance);
                // Oder, wenn die Idee ist, dass es bei distance=0 strength ist und bei falloff_distance = 0:
                // self.strength / (normalized_distance * normalized_distance * (self.falloff_distance * self.falloff_distance) + 1.0)
                // Die gegebene Formel ist aber häufig: strength / (d^2/k + 1) oder strength * (1 - (d/D)^2)
                // Behalten wir die Variante, die weniger stark abfällt als Linear:
                // self.strength * (1.0 - normalized_distance).powi(2) // Quadratischer Abfall zu 0
                self.strength / (distance * distance + 1.0) // Diese Variante ist üblich für Metaballs, aber nicht normalisiert
            }
            PolygonFalloff::Exponential => {
                // scale_factor steuert die Stärke des Abfalls
                self.strength * (-normalized_distance * self.scale_factor.max(0.1)).exp()
            }
            PolygonFalloff::Gaussian => {
                // variance könnte mit scale_factor verbunden sein, oder ein fester Teil der falloff_distance
                // Original: variance = self.falloff_distance * 0.5
                // Ein typischerer Ansatz ist, dass bei falloff_distance der Wert einen bestimmten Bruchteil erreicht.
                // Hier wird scale_factor als eine Art "Breite" der Glockenkurve interpretiert.
                let variance_param = (self.falloff_distance / self.scale_factor.max(0.1)).powi(2); // scale_factor > 0
                self.strength * (-(distance * distance) / (2.0 * variance_param)).exp()
            }
            PolygonFalloff::Step => {
                if distance <= self.falloff_distance * 0.5 {
                    // Harter Übergang in der Mitte
                    self.strength
                } else {
                    0.0
                }
            }
            PolygonFalloff::Sigmoid => {
                // steepness wird hier durch scale_factor gesteuert
                // midpoint ist die Mitte der falloff_distance
                let steepness = self.scale_factor.max(0.1); // Sichert > 0
                let midpoint = self.falloff_distance * 0.5;
                self.strength
                    / (1.0 + ((distance - midpoint) * steepness / self.falloff_distance).exp()) // Normalisierte Steilheit
            }
        }
    }
}

impl std::fmt::Debug for PolygonInfluence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PolygonInfluence")
            .field("polygon_vertices", &self.polygon.vertices().len()) // Zeige nur Anzahl für Kürze
            .field("strength", &self.strength)
            .field("falloff_distance", &self.falloff_distance)
            .field("falloff_type", &self.falloff_type)
            .field("interior_only", &self.interior_only)
            .field("offset", &self.offset)
            .field("scale_factor", &self.scale_factor)
            .field("cached_bounds", &self.cached_bounds.is_some())
            .finish()
    }
}

impl crate::math::geometry::metalballs::FieldInfluence for PolygonInfluence {
    fn influence_at(&self, point: Vec2) -> f32 {
        // Früher Ausstieg, wenn Punkt außerhalb der erweiterten Bounding Box liegt.
        // Die Bounding Box wird hier bei Bedarf berechnet, aber nicht im self gecacht (da &self).
        // Für Performance bei wiederholten Aufrufen an vielen Punkten sollte get_bounds_mut() einmal vorher aufgerufen werden.
        let bounds = if let Some(cached) = self.cached_bounds {
            cached
        } else if let Some((min, max)) = self.polygon.bounds() {
            Bounds2D::from_points(min, max).expand(self.falloff_distance)
        } else {
            // Leeres Polygon oder Polygon ohne Bounds hat keinen Einflussbereich
            return 0.0;
        };

        let point_2d = Point2D::new(point.x, point.y);
        if !bounds.contains_point(point_2d) {
            return 0.0;
        }

        let distance = self.distance_to_polygon(point); // Kann negativ sein

        // Wenn interior_only und Punkt außerhalb (distance > 0), kein Einfluss
        if self.interior_only && distance > constants::EPSILON {
            // kleine Toleranz für Fließkomma
            return 0.0;
        }

        let mut influence_val = self.apply_falloff(distance.abs());

        // Addiere Offset nur, wenn ein Einfluss vorhanden ist oder interior_only und wir sind drinnen
        if influence_val > constants::EPSILON
            || (self.interior_only && distance <= constants::EPSILON)
        {
            influence_val += self.offset;
        }

        influence_val.max(0.0) // Stelle sicher, dass der Wert nicht negativ wird
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        if let Some(cached) = self.cached_bounds {
            return Some(cached);
        }
        // Wenn nicht gecacht, berechne es frisch (ohne zu cachen, da &self)
        if let Some((min, max)) = self.polygon.bounds() {
            Some(Bounds2D::from_points(min, max).expand(self.falloff_distance))
        } else {
            None // Leeres Polygon hat keine Bounds
        }
    }

    fn max_influence_distance(&self) -> Option<f32> {
        Some(self.falloff_distance)
    }

    fn influence_type(&self) -> &'static str {
        "Polygon"
    }
}

/// Builder für komplexe Polygon-Influences
pub struct PolygonInfluenceBuilder {
    polygon: Option<Polygon>,
    strength: f32,
    falloff_distance: f32,
    falloff_type: PolygonFalloff,
    interior_only: bool,
    offset: f32,
    scale_factor: f32,
}

impl PolygonInfluenceBuilder {
    pub fn new() -> Self {
        Self {
            polygon: None,
            strength: 1.0,
            falloff_distance: 10.0,                  // Standardwert
            falloff_type: PolygonFalloff::default(), // Quadratic
            interior_only: false,
            offset: 0.0,
            scale_factor: 1.0, // Neutraler Standard für die meisten Falloffs (z.B. Exp, Gaussian)
                               // Für Sigmoid könnte ein anderer Default besser sein, z.B. 5.0-10.0
        }
    }

    pub fn polygon(mut self, polygon: Polygon) -> Self {
        self.polygon = Some(polygon);
        self
    }

    pub fn vertices(mut self, vertices: Vec<Point2D>) -> Self {
        match Polygon::closed(vertices) {
            Ok(polygon) => self.polygon = Some(polygon),
            Err(_) => {
                // Fehlerbehandlung: Loggen oder Panik, oder Builder in invaliden Zustand?
                // Fürs Erste: Polygon bleibt None, build() wird None zurückgeben.
                eprintln!("Failed to create polygon from vertices for PolygonInfluenceBuilder");
            }
        }
        self
    }

    pub fn strength(mut self, strength: f32) -> Self {
        self.strength = strength;
        self
    }

    pub fn falloff_distance(mut self, distance: f32) -> Self {
        self.falloff_distance = distance.max(0.0); // Sicherstellen, dass nicht negativ
        self
    }

    pub fn falloff_type(mut self, falloff_type: PolygonFalloff) -> Self {
        self.falloff_type = falloff_type;
        self
    }

    pub fn interior_only(mut self, interior_only: bool) -> Self {
        self.interior_only = interior_only;
        self
    }

    pub fn offset(mut self, offset: f32) -> Self {
        self.offset = offset;
        self
    }

    pub fn scale_factor(mut self, scale_factor: f32) -> Self {
        self.scale_factor = scale_factor;
        self
    }

    pub fn build(self) -> Option<PolygonInfluence> {
        if let Some(polygon) = self.polygon {
            Some(PolygonInfluence::new(
                polygon,
                self.strength,
                self.falloff_distance,
                self.falloff_type,
                self.interior_only,
                self.offset,
                self.scale_factor,
            ))
        } else {
            None // Kein Polygon definiert
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
    // Annahme: FieldInfluence ist hier im Scope für Tests
    use crate::math::geometry::metalballs::FieldInfluence;

    /// Hilfsfunktionen zum Erstellen von Polygon-Shapes für Tests
    struct TestShapes;

    impl TestShapes {
        fn circle_vertices(center: Vec2, radius: f32, segments: usize) -> Vec<Point2D> {
            let mut vertices = Vec::with_capacity(segments);
            for i in 0..segments {
                let angle = (i as f32 / segments as f32) * std::f32::consts::TAU;
                let x = center.x + radius * angle.cos();
                let y = center.y + radius * angle.sin();
                vertices.push(Point2D::new(x, y));
            }
            vertices
        }

        pub fn circle_influence_builder(
            center: Vec2,
            radius: f32,
            strength: f32,
            segments: usize,
        ) -> PolygonInfluenceBuilder {
            let vertices = Self::circle_vertices(center, radius, segments);
            PolygonInfluenceBuilder::new()
                .vertices(vertices)
                .strength(strength)
        }

        pub fn rectangle_influence_builder(
            center: Vec2,
            width: f32,
            height: f32,
            strength: f32,
        ) -> PolygonInfluenceBuilder {
            let half_w = width * 0.5;
            let half_h = height * 0.5;
            let vertices = vec![
                Point2D::new(center.x - half_w, center.y - half_h),
                Point2D::new(center.x + half_w, center.y - half_h),
                Point2D::new(center.x + half_w, center.y + half_h),
                Point2D::new(center.x - half_w, center.y + half_h),
            ];
            PolygonInfluenceBuilder::new()
                .vertices(vertices)
                .strength(strength)
        }
        pub fn star_influence_builder(
            center: Vec2,
            inner_radius: f32,
            outer_radius: f32,
            points: usize,
            strength: f32,
        ) -> PolygonInfluenceBuilder {
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
            PolygonInfluenceBuilder::new()
                .vertices(vertices)
                .strength(strength)
        }

        pub fn regular_polygon_influence_builder(
            center: Vec2,
            radius: f32,
            sides: usize,
            strength: f32,
        ) -> PolygonInfluenceBuilder {
            let mut builder = PolygonInfluenceBuilder::new().strength(strength);
            if sides < 3 {
                builder = builder.vertices(Self::circle_vertices(center, radius, 8)); // Fallback
            } else {
                let mut vertices = Vec::with_capacity(sides);
                for i in 0..sides {
                    let angle = (i as f32 / sides as f32) * std::f32::consts::TAU;
                    let x = center.x + radius * angle.cos();
                    let y = center.y + radius * angle.sin();
                    vertices.push(Point2D::new(x, y));
                }
                builder = builder.vertices(vertices);
            }
            builder
        }
    }

    #[test]
    fn test_polygon_influence_basic() {
        let mut square = TestShapes::rectangle_influence_builder(Vec2::ZERO, 2.0, 2.0, 1.0)
            .falloff_distance(5.0) // Explizit setzen für Test
            .build()
            .unwrap();

        // Cache muss manuell initialisiert werden, wenn er vor influence_at() genutzt werden soll
        // oder wenn bounding_box() den gecachten Wert zurückgeben soll.
        // influence_at() berechnet es intern, wenn nicht vorhanden, aber cacht nicht.
        square.get_bounds_mut();

        // Punkt im Zentrum: Abstand sollte ca. -1 sein (Hälfte der Dicke),
        // abs(Abstand) = 1. Falloff wird angewendet.
        let center_influence = square.influence_at(Vec2::new(0.0, 0.0));
        // Da strength=1 und falloff_distance=5, und distance_to_boundary(0,0) = 1
        // d_abs = 1. normalized_d = 1/5 = 0.2
        // Quadratic: 1.0 / (1*1 + 1) = 0.5
        assert!(
            center_influence > 0.4 && center_influence < 0.6,
            "Center influence: {}",
            center_influence
        );

        // Punkt weit entfernt sollte keinen Einfluss haben
        let far_influence = square.influence_at(Vec2::new(100.0, 100.0));
        assert_eq!(far_influence, 0.0, "Far influence should be 0");
    }

    #[test]
    fn test_falloff_types() {
        let base_builder =
            TestShapes::circle_influence_builder(Vec2::ZERO, 1.0, 1.0, 16).falloff_distance(5.0); // Distanz vom Rand des Kreises

        let linear_influence_obj = base_builder
            .clone() // Builder ist Cloneable wenn Polygon Option<Arc<Polygon>> wäre oder Polygon Clone ist
            .falloff_type(PolygonFalloff::Linear)
            .build()
            .unwrap();

        let exp_influence_obj = base_builder
            .clone()
            .falloff_type(PolygonFalloff::Exponential)
            .scale_factor(3.0) // für Exponential relevant
            .build()
            .unwrap();

        // Punkt außerhalb des Kreises (Radius 1), aber innerhalb Falloff (5)
        // Abstand zum Rand des Kreises ist 1.0 (Punkt (2,0) zu Kreisrand (1,0))
        let test_point = Vec2::new(2.0, 0.0);
        // distance_to_polygon für (2,0) bei Kreis (0,0) Radius 1 ist 1.0

        let linear_influence = linear_influence_obj.influence_at(test_point); // d=1, D=5, norm_d=0.2. Str*(1-0.2)=0.8
        let exp_influence = exp_influence_obj.influence_at(test_point); // Str*exp(-0.2 * 3) = exp(-0.6) approx 0.548

        assert!(linear_influence > 0.0, "Linear influence should be > 0");
        assert!(
            linear_influence > 0.79 && linear_influence < 0.81,
            "Linear: {}",
            linear_influence
        );

        assert!(exp_influence > 0.0, "Exponential influence should be > 0");
        assert!(
            exp_influence > 0.54 && exp_influence < 0.56,
            "Exp: {}",
            exp_influence
        );

        assert_ne!(
            linear_influence.to_bits(),
            exp_influence.to_bits(),
            "Influences should differ"
        );
    }

    #[test]
    fn test_interior_only_mode() {
        let square = TestShapes::rectangle_influence_builder(Vec2::ZERO, 2.0, 2.0, 1.0)
            .interior_only(true)
            .falloff_distance(1.0) // Falloff ist hier weniger relevant, da außen 0
            .build()
            .unwrap();

        // Punkt innerhalb sollte Einfluss haben
        // distance_to_polygon ist 0.0. apply_falloff(0) sollte strength (1.0) geben.
        let inside_influence = square.influence_at(Vec2::new(0.5, 0.5));
        assert!(
            inside_influence > 0.99 && inside_influence < 1.01,
            "Inside influence: {}",
            inside_influence
        );

        // Punkt außerhalb sollte keinen Einfluss haben
        let outside_influence = square.influence_at(Vec2::new(2.0, 2.0));
        assert_eq!(
            outside_influence, 0.0,
            "Outside influence should be 0 for interior_only"
        );
    }

    #[test]
    fn test_polygon_shapes() {
        let circle = TestShapes::circle_influence_builder(Vec2::ZERO, 1.0, 1.0, 12)
            .build()
            .unwrap();
        let star = TestShapes::star_influence_builder(Vec2::ZERO, 0.5, 1.0, 5, 1.0)
            .build()
            .unwrap();
        let hexagon = TestShapes::regular_polygon_influence_builder(Vec2::ZERO, 1.0, 6, 1.0)
            .build()
            .unwrap();

        // Alle sollten am Zentrum einen Wert > 0 haben (Details hängen von Falloff ab)
        // Standard Falloff ist Quadratic. strength / (d*d + 1)
        // Für Punkt im Zentrum (0,0) ist d = distance_to_boundary.
        // Für Kreis(R=1): d=1 -> 1/(1*1+1)=0.5
        // Für Star(inner=0.5): d=0.5 -> 1/(0.5*0.5+1) = 1/1.25 = 0.8
        // Für Hexagon(R=1): d=cos(PI/6) = sqrt(3)/2 ~0.866 -> 1/(0.866^2+1) = 1/(0.75+1) = 1/1.75 ~ 0.57
        assert!(circle.influence_at(Vec2::ZERO) > 0.4 && circle.influence_at(Vec2::ZERO) < 0.6);
        assert!(star.influence_at(Vec2::ZERO) > 0.7 && star.influence_at(Vec2::ZERO) < 0.9);
        assert!(hexagon.influence_at(Vec2::ZERO) > 0.5 && hexagon.influence_at(Vec2::ZERO) < 0.7);
    }

    #[test]
    fn test_builder_pattern_full() {
        let polygon_influence = PolygonInfluenceBuilder::new()
            .vertices(vec![
                Point2D::new(-1.0, -1.0),
                Point2D::new(1.0, -1.0),
                Point2D::new(1.0, 1.0),
                Point2D::new(-1.0, 1.0),
            ])
            .strength(2.0)
            .falloff_distance(5.0)
            .falloff_type(PolygonFalloff::Gaussian)
            .interior_only(true)
            .offset(0.1)
            .scale_factor(0.8)
            .build()
            .unwrap();

        assert_eq!(polygon_influence.strength, 2.0);
        assert_eq!(polygon_influence.falloff_distance, 5.0);
        assert_eq!(polygon_influence.falloff_type, PolygonFalloff::Gaussian);
        assert_eq!(polygon_influence.interior_only, true);
        assert_eq!(polygon_influence.offset, 0.1);
        assert_eq!(polygon_influence.scale_factor, 0.8);
    }
    #[test]
    fn test_empty_polygon_influence() {
        let influence = PolygonInfluenceBuilder::new()
            .vertices(vec![]) // Leeres Polygon
            .build();
        assert!(
            influence.is_none(),
            "Building with empty vertices should fail or return None"
        );

        // Test mit Polygon, das keine Bounds hat (z.B. nur ein Punkt, obwohl Polygon::closed das verhindern sollte)
        // Angenommen Polygon::closed würde ein Polygon mit einem Punkt erstellen:
        // let single_point_poly = Polygon::closed(vec![Point2D::new(0.0,0.0)]).unwrap();
        // let influence_single = PolygonInfluenceBuilder::new().polygon(single_point_poly).build().unwrap();
        // assert_eq!(influence_single.influence_at(Vec2::new(10.0, 10.0)), 0.0);
        // assert!(influence_single.bounding_box().is_some()); // Sollte kleine Box um den Punkt sein
    }

    #[test]
    fn test_distance_to_boundary_precision() {
        // Quadrat von (-1,-1) bis (1,1)
        let square_verts = vec![
            Point2D::new(-1.0, -1.0),
            Point2D::new(1.0, -1.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(-1.0, 1.0),
        ];
        let polygon = Polygon::closed(square_verts).unwrap();
        let influence = PolygonInfluenceBuilder::new()
            .polygon(polygon)
            .strength(1.0)
            .build()
            .unwrap();

        let p_center = Point2D::new(0.0, 0.0);
        let dist_center = influence.distance_to_boundary(p_center);
        assert!(
            (dist_center - 1.0).abs() < constants::EPSILON,
            "Center distance: {}",
            dist_center
        );

        let p_edge_mid = Point2D::new(0.5, -1.0); // Direkt auf der unteren Kante
        let dist_edge_mid = influence.distance_to_boundary(p_edge_mid);
        assert!(
            dist_edge_mid < constants::EPSILON,
            "Edge mid distance: {}",
            dist_edge_mid
        );

        let p_corner = Point2D::new(1.0, 1.0); // Direkt auf Ecke
        let dist_corner = influence.distance_to_boundary(p_corner);
        assert!(
            dist_corner < constants::EPSILON,
            "Corner distance: {}",
            dist_corner
        );

        let p_outside_edge = Point2D::new(0.0, -2.0); // Außerhalb, mittig zur unteren Kante
        let dist_outside_edge = influence.distance_to_boundary(p_outside_edge);
        assert!(
            (dist_outside_edge - 1.0).abs() < constants::EPSILON,
            "Outside edge distance: {}",
            dist_outside_edge
        );

        let p_outside_corner = Point2D::new(2.0, 2.0); // Außerhalb, diagonal zur Ecke (1,1)
        // Distanz zu (1,1) ist sqrt((2-1)^2 + (2-1)^2) = sqrt(1+1) = sqrt(2)
        let dist_outside_corner = influence.distance_to_boundary(p_outside_corner);
        assert!(
            (dist_outside_corner - std::f32::consts::SQRT_2).abs() < constants::EPSILON,
            "Outside corner distance: {}",
            dist_outside_corner
        );
    }
}
