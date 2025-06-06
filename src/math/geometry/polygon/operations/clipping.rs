// src/math/geometry/polygon/operations/clipping.rs

use super::super::Polygon;
use crate::math::{error::*, types::*, utils::*};
use bevy::math::Vec2;

/// Verschiedene Clipping-Algorithmen
#[derive(Debug, Clone, Copy)]
pub enum ClippingAlgorithm {
    /// Sutherland-Hodgman Clipping
    SutherlandHodgman,
    /// Cohen-Sutherland (für rechteckige Clipper)
    CohenSutherland,
    /// Liang-Barsky (für rechteckige Clipper)
    LiangBarsky,
    /// Weiler-Atherton (für beliebige Polygone)
    WeilerAtherton,
}

/// Edge-Klassifikation für Clipping
#[derive(Debug, Clone, Copy, PartialEq)]
enum EdgeRelation {
    Inside,
    Outside,
    Entering,
    Leaving,
}

/// Clipping-Engine für Polygone
pub struct PolygonClipper {
    algorithm: ClippingAlgorithm,
    tolerance: f32,
}

impl PolygonClipper {
    /// Erstellt einen neuen Clipper
    pub fn new(algorithm: ClippingAlgorithm) -> Self {
        Self {
            algorithm,
            tolerance: constants::EPSILON * 1000.0,
        }
    }

    /// Setzt die Toleranz
    pub fn with_tolerance(mut self, tolerance: f32) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Clippt ein Polygon gegen ein Clipping-Polygon
    pub fn clip_polygon(&self, subject: &Polygon, clipper: &Polygon) -> MathResult<Vec<Polygon>> {
        match self.algorithm {
            ClippingAlgorithm::SutherlandHodgman => self.sutherland_hodgman_clip(subject, clipper),
            ClippingAlgorithm::WeilerAtherton => self.weiler_atherton_clip(subject, clipper),
            _ => Err(MathError::InvalidConfiguration {
                message: "Algorithm not supported for polygon-polygon clipping".to_string(),
            }),
        }
    }

    /// Clippt ein Polygon gegen ein Rechteck
    pub fn clip_against_rectangle(
        &self,
        subject: &Polygon,
        bounds: Bounds2D,
    ) -> MathResult<Vec<Polygon>> {
        match self.algorithm {
            ClippingAlgorithm::SutherlandHodgman => {
                self.clip_against_rect_sutherland_hodgman(subject, bounds)
            }
            ClippingAlgorithm::CohenSutherland => {
                self.clip_against_rect_cohen_sutherland(subject, bounds)
            }
            ClippingAlgorithm::LiangBarsky => self.clip_against_rect_liang_barsky(subject, bounds),
            _ => Err(MathError::InvalidConfiguration {
                message: "Algorithm not supported for rectangle clipping".to_string(),
            }),
        }
    }

    // === Sutherland-Hodgman Algorithmus ===

    fn sutherland_hodgman_clip(
        &self,
        subject: &Polygon,
        clipper: &Polygon,
    ) -> MathResult<Vec<Polygon>> {
        let mut output_vertices = subject.vertices().to_vec();
        let clipper_vertices = clipper.vertices();
        let n_clipper =
            if clipper.is_closed() && clipper_vertices.first() == clipper_vertices.last() {
                clipper_vertices.len() - 1
            } else {
                clipper_vertices.len()
            };

        // Für jede Kante des Clipping-Polygons
        for i in 0..n_clipper {
            if output_vertices.is_empty() {
                break;
            }

            let clip_vertex1 = clipper_vertices[i];
            let clip_vertex2 = clipper_vertices[(i + 1) % n_clipper];

            let input_vertices = output_vertices.clone();
            output_vertices.clear();

            if input_vertices.is_empty() {
                continue;
            }

            let mut s = input_vertices[input_vertices.len() - 1];

            for &e in &input_vertices {
                if self.is_inside(e, clip_vertex1, clip_vertex2) {
                    if !self.is_inside(s, clip_vertex1, clip_vertex2) {
                        // Entering
                        if let Some(intersection) =
                            self.line_intersection(s, e, clip_vertex1, clip_vertex2)
                        {
                            output_vertices.push(intersection);
                        }
                    }
                    output_vertices.push(e);
                } else if self.is_inside(s, clip_vertex1, clip_vertex2) {
                    // Leaving
                    if let Some(intersection) =
                        self.line_intersection(s, e, clip_vertex1, clip_vertex2)
                    {
                        output_vertices.push(intersection);
                    }
                }
                s = e;
            }
        }

        if output_vertices.len() >= 3 {
            Ok(vec![Polygon::closed(output_vertices)?])
        } else {
            Ok(vec![])
        }
    }

    fn clip_against_rect_sutherland_hodgman(
        &self,
        subject: &Polygon,
        bounds: Bounds2D,
    ) -> MathResult<Vec<Polygon>> {
        let mut output_vertices = subject.vertices().to_vec();

        // Clippe gegen jede Seite des Rechtecks
        let clip_edges = [
            (bounds.min, Vec2::new(bounds.max.x, bounds.min.y)), // Bottom
            (Vec2::new(bounds.max.x, bounds.min.y), bounds.max), // Right
            (bounds.max, Vec2::new(bounds.min.x, bounds.max.y)), // Top
            (Vec2::new(bounds.min.x, bounds.max.y), bounds.min), // Left
        ];

        for (edge_start, edge_end) in &clip_edges {
            if output_vertices.is_empty() {
                break;
            }

            let input_vertices = output_vertices.clone();
            output_vertices.clear();

            if input_vertices.is_empty() {
                continue;
            }

            let mut s = input_vertices[input_vertices.len() - 1];

            for &e in &input_vertices {
                if self.is_inside(e, *edge_start, *edge_end) {
                    if !self.is_inside(s, *edge_start, *edge_end) {
                        if let Some(intersection) =
                            self.line_intersection(s, e, *edge_start, *edge_end)
                        {
                            output_vertices.push(intersection);
                        }
                    }
                    output_vertices.push(e);
                } else if self.is_inside(s, *edge_start, *edge_end) {
                    if let Some(intersection) = self.line_intersection(s, e, *edge_start, *edge_end)
                    {
                        output_vertices.push(intersection);
                    }
                }
                s = e;
            }
        }

        if output_vertices.len() >= 3 {
            Ok(vec![Polygon::closed(output_vertices)?])
        } else {
            Ok(vec![])
        }
    }

    // === Cohen-Sutherland für Rechteck-Clipping ===

    fn clip_against_rect_cohen_sutherland(
        &self,
        subject: &Polygon,
        bounds: Bounds2D,
    ) -> MathResult<Vec<Polygon>> {
        let vertices = subject.vertices();
        let mut clipped_edges = Vec::new();

        for i in 0..vertices.len() {
            let j = (i + 1) % vertices.len();
            let start = vertices[i];
            let end = vertices[j];

            if let Some((clipped_start, clipped_end)) =
                self.cohen_sutherland_line_clip(start, end, bounds)
            {
                clipped_edges.push((clipped_start, clipped_end));
            }
        }

        // Verbinde die geclippten Kanten zu Polygonen
        self.reconstruct_polygons_from_edges(&clipped_edges)
    }

    fn cohen_sutherland_line_clip(
        &self,
        mut p1: Vec2,
        mut p2: Vec2,
        bounds: Bounds2D,
    ) -> Option<(Vec2, Vec2)> {
        let mut outcode1 = self.compute_outcode(p1, bounds);
        let mut outcode2 = self.compute_outcode(p2, bounds);

        loop {
            if (outcode1 | outcode2) == 0 {
                // Beide Punkte sind inside
                return Some((p1, p2));
            } else if (outcode1 & outcode2) != 0 {
                // Beide Punkte sind outside auf derselben Seite
                return None;
            } else {
                // Mindestens ein Punkt ist outside
                let outcode_out = if outcode1 != 0 { outcode1 } else { outcode2 };

                let (x, y) = if (outcode_out & 8) != 0 {
                    // Top
                    let x = p1.x + (p2.x - p1.x) * (bounds.max.y - p1.y) / (p2.y - p1.y);
                    (x, bounds.max.y)
                } else if (outcode_out & 4) != 0 {
                    // Bottom
                    let x = p1.x + (p2.x - p1.x) * (bounds.min.y - p1.y) / (p2.y - p1.y);
                    (x, bounds.min.y)
                } else if (outcode_out & 2) != 0 {
                    // Right
                    let y = p1.y + (p2.y - p1.y) * (bounds.max.x - p1.x) / (p2.x - p1.x);
                    (bounds.max.x, y)
                } else {
                    // Left
                    let y = p1.y + (p2.y - p1.y) * (bounds.min.x - p1.x) / (p2.x - p1.x);
                    (bounds.min.x, y)
                };

                if outcode_out == outcode1 {
                    p1 = Vec2::new(x, y);
                    outcode1 = self.compute_outcode(p1, bounds);
                } else {
                    p2 = Vec2::new(x, y);
                    outcode2 = self.compute_outcode(p2, bounds);
                }
            }
        }
    }

    fn compute_outcode(&self, point: Vec2, bounds: Bounds2D) -> u8 {
        let mut code = 0;

        if point.x < bounds.min.x {
            code |= 1; // Left
        }
        if point.x > bounds.max.x {
            code |= 2; // Right
        }
        if point.y < bounds.min.y {
            code |= 4; // Bottom
        }
        if point.y > bounds.max.y {
            code |= 8; // Top
        }

        code
    }

    // === Liang-Barsky für Rechteck-Clipping ===

    fn clip_against_rect_liang_barsky(
        &self,
        subject: &Polygon,
        bounds: Bounds2D,
    ) -> MathResult<Vec<Polygon>> {
        let vertices = subject.vertices();
        let mut clipped_edges = Vec::new();

        for i in 0..vertices.len() {
            let j = (i + 1) % vertices.len();
            let start = vertices[i];
            let end = vertices[j];

            if let Some((clipped_start, clipped_end)) =
                self.liang_barsky_line_clip(start, end, bounds)
            {
                clipped_edges.push((clipped_start, clipped_end));
            }
        }

        self.reconstruct_polygons_from_edges(&clipped_edges)
    }

    fn liang_barsky_line_clip(&self, p1: Vec2, p2: Vec2, bounds: Bounds2D) -> Option<(Vec2, Vec2)> {
        let dx = p2.x - p1.x;
        let dy = p2.y - p1.y;

        let p = [-dx, dx, -dy, dy];
        let q = [
            p1.x - bounds.min.x,
            bounds.max.x - p1.x,
            p1.y - bounds.min.y,
            bounds.max.y - p1.y,
        ];

        let mut u1 = 0.0_f32;
        let mut u2 = 1.0_f32;

        for i in 0..4 {
            if p[i] == 0.0 {
                if q[i] < 0.0 {
                    return None; // Linie ist parallel und außerhalb
                }
            } else {
                let t = q[i] / p[i];
                if p[i] < 0.0 {
                    u1 = u1.max(t);
                } else {
                    u2 = u2.min(t);
                }
            }
        }

        if u1 <= u2 {
            let clipped_start = Vec2::new(p1.x + u1 * dx, p1.y + u1 * dy);
            let clipped_end = Vec2::new(p1.x + u2 * dx, p1.y + u2 * dy);
            Some((clipped_start, clipped_end))
        } else {
            None
        }
    }

    // === Weiler-Atherton für Polygon-Polygon Clipping ===

    fn weiler_atherton_clip(
        &self,
        subject: &Polygon,
        clipper: &Polygon,
    ) -> MathResult<Vec<Polygon>> {
        // Vereinfachte Implementierung von Weiler-Atherton
        // Für vollständige Implementierung wäre eine komplexere Datenstruktur nötig

        let subject_vertices = subject.vertices();
        let clipper_vertices = clipper.vertices();

        // Finde alle Schnittpunkte
        let mut intersections = Vec::new();

        for i in 0..subject_vertices.len() {
            let j = (i + 1) % subject_vertices.len();
            let subject_edge = (subject_vertices[i], subject_vertices[j]);

            for k in 0..clipper_vertices.len() {
                let l = (k + 1) % clipper_vertices.len();
                let clipper_edge = (clipper_vertices[k], clipper_vertices[l]);

                if let Some(intersection) = self.line_intersection(
                    subject_edge.0,
                    subject_edge.1,
                    clipper_edge.0,
                    clipper_edge.1,
                ) {
                    intersections.push(intersection);
                }
            }
        }

        // Vereinfachte Behandlung: verwende Sutherland-Hodgman als Fallback
        self.sutherland_hodgman_clip(subject, clipper)
    }

    // === Hilfsfunktionen ===

    fn is_inside(&self, point: Vec2, edge_start: Vec2, edge_end: Vec2) -> bool {
        // Punkt ist "inside" wenn er auf der linken Seite der gerichteten Kante liegt
        let cross = (edge_end.x - edge_start.x) * (point.y - edge_start.y)
            - (edge_end.y - edge_start.y) * (point.x - edge_start.x);
        cross >= -self.tolerance
    }

    fn line_intersection(&self, p1: Vec2, p2: Vec2, p3: Vec2, p4: Vec2) -> Option<Vec2> {
        let d1 = p2 - p1;
        let d2 = p4 - p3;

        let denominator = d1.x * d2.y - d1.y * d2.x;

        if denominator.abs() < self.tolerance {
            return None; // Parallele Linien
        }

        let t = ((p3.x - p1.x) * d2.y - (p3.y - p1.y) * d2.x) / denominator;

        Some(p1 + d1 * t)
    }

    fn reconstruct_polygons_from_edges(&self, edges: &[(Vec2, Vec2)]) -> MathResult<Vec<Polygon>> {
        if edges.is_empty() {
            return Ok(vec![]);
        }

        // Vereinfachte Rekonstruktion: verbinde alle Kanten zu einem Polygon
        let mut vertices = Vec::new();

        for (start, end) in edges {
            if vertices.is_empty() || vertices.last() != Some(start) {
                vertices.push(*start);
            }
            vertices.push(*end);
        }

        if vertices.len() >= 3 {
            Ok(vec![Polygon::new(vertices)?])
        } else {
            Ok(vec![])
        }
    }
}

/// Spezielle Clipping-Operationen
pub struct ClippingOperations;

impl ClippingOperations {
    /// Clippt ein Polygon gegen ein Rechteck (vereinfachte API)
    pub fn clip_to_rectangle(polygon: &Polygon, bounds: Bounds2D) -> MathResult<Vec<Polygon>> {
        let clipper = PolygonClipper::new(ClippingAlgorithm::SutherlandHodgman);
        clipper.clip_against_rectangle(polygon, bounds)
    }

    /// Clippt ein Polygon gegen einen Kreis
    pub fn clip_to_circle(
        polygon: &Polygon,
        center: Vec2,
        radius: f32,
    ) -> MathResult<Vec<Polygon>> {
        let vertices = polygon.vertices();
        let mut clipped_vertices = Vec::new();

        for &vertex in vertices {
            let distance = vertex.distance(center);
            if distance <= radius {
                clipped_vertices.push(vertex);
            } else {
                // Projiziere auf Kreis-Rand
                let direction = (vertex - center).normalize();
                clipped_vertices.push(center + direction * radius);
            }
        }

        if clipped_vertices.len() >= 3 {
            Ok(vec![if polygon.is_closed() {
                Polygon::closed(clipped_vertices)?
            } else {
                Polygon::new(clipped_vertices)?
            }])
        } else {
            Ok(vec![])
        }
    }

    /// Multi-Level Clipping (hierarchisches Clipping)
    pub fn multi_level_clip(polygon: &Polygon, clippers: &[Bounds2D]) -> MathResult<Vec<Polygon>> {
        let mut current_polygons = vec![polygon.clone()];

        for &bounds in clippers {
            let mut next_polygons = Vec::new();

            for current_polygon in &current_polygons {
                let clipped = Self::clip_to_rectangle(current_polygon, bounds)?;
                next_polygons.extend(clipped);
            }

            current_polygons = next_polygons;

            if current_polygons.is_empty() {
                break;
            }
        }

        Ok(current_polygons)
    }

    /// Adaptive Clipping (wählt Algorithmus basierend auf Polygon-Eigenschaften)
    pub fn adaptive_clip(subject: &Polygon, clipper: &Polygon) -> MathResult<Vec<Polygon>> {
        // Wähle Algorithmus basierend auf Komplexität
        let algorithm = if subject.len() < 10 && clipper.len() < 10 {
            ClippingAlgorithm::SutherlandHodgman
        } else {
            ClippingAlgorithm::WeilerAtherton
        };

        let clipper_engine = PolygonClipper::new(algorithm);
        clipper_engine.clip_polygon(subject, clipper)
    }
}

/// Clipping-Analyse und Debugging
pub struct ClippingAnalysis;

impl ClippingAnalysis {
    /// Analysiert warum ein Clipping-Vorgang fehlgeschlagen ist
    pub fn analyze_clipping_failure(subject: &Polygon, clipper: &Polygon) -> String {
        use crate::math::geometry::polygon::PolygonProperties;

        let mut issues = Vec::new();

        if subject.len() < 3 {
            issues.push("Subject polygon has fewer than 3 vertices".to_string());
        }

        if clipper.len() < 3 {
            issues.push("Clipper polygon has fewer than 3 vertices".to_string());
        }

        if subject.area() < constants::EPSILON {
            issues.push("Subject polygon has zero area".to_string());
        }

        if clipper.area() < constants::EPSILON {
            issues.push("Clipper polygon has zero area".to_string());
        }

        // Prüfe auf Überschneidung
        if let (Some(subject_bounds), Some(clipper_bounds)) = (subject.bounds(), clipper.bounds()) {
            let subject_bounds = Bounds2D::from_points(subject_bounds.min, subject_bounds.max);
            let clipper_bounds = Bounds2D::from_points(clipper_bounds.min, clipper_bounds.max);

            if !subject_bounds.intersects(&clipper_bounds) {
                issues.push("Polygons do not overlap (bounding boxes don't intersect)".to_string());
            }
        }

        if issues.is_empty() {
            "No obvious issues detected".to_string()
        } else {
            issues.join("; ")
        }
    }

    /// Berechnet Clipping-Statistiken
    pub fn clipping_statistics(
        subject: &Polygon,
        clipper: &Polygon,
        result: &[Polygon],
    ) -> ClippingStats {
        use crate::math::geometry::polygon::PolygonProperties;

        let original_area = subject.area();
        let clipped_area: f32 = result.iter().map(|p| p.area()).sum();
        let area_ratio = if original_area > 0.0 {
            clipped_area / original_area
        } else {
            0.0
        };

        ClippingStats {
            original_vertex_count: subject.len(),
            clipper_vertex_count: clipper.len(),
            result_polygon_count: result.len(),
            total_result_vertices: result.iter().map(|p| p.len()).sum(),
            original_area,
            clipped_area,
            area_retention_ratio: area_ratio,
        }
    }
}

/// Statistiken über einen Clipping-Vorgang
#[derive(Debug, Clone)]
pub struct ClippingStats {
    pub original_vertex_count: usize,
    pub clipper_vertex_count: usize,
    pub result_polygon_count: usize,
    pub total_result_vertices: usize,
    pub original_area: f32,
    pub clipped_area: f32,
    pub area_retention_ratio: f32,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::geometry::polygon::Polygon;

    #[test]
    fn test_rectangle_clipping() {
        let square = Polygon::closed(vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(2.0, 0.0),
            Vec2::new(2.0, 2.0),
            Vec2::new(0.0, 2.0),
        ])
        .unwrap();

        let clip_bounds = Bounds2D::from_points(Vec2::new(1.0, 1.0), Vec2::new(3.0, 3.0));

        let clipped = ClippingOperations::clip_to_rectangle(&square, clip_bounds).unwrap();

        assert_eq!(clipped.len(), 1);
        assert!(clipped[0].len() >= 3);
    }

    #[test]
    fn test_sutherland_hodgman() {
        let subject = Polygon::closed(vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(2.0, 0.0),
            Vec2::new(2.0, 2.0),
            Vec2::new(0.0, 2.0),
        ])
        .unwrap();

        let clipper = Polygon::closed(vec![
            Vec2::new(1.0, 1.0),
            Vec2::new(3.0, 1.0),
            Vec2::new(3.0, 3.0),
            Vec2::new(1.0, 3.0),
        ])
        .unwrap();

        let clipper_engine = PolygonClipper::new(ClippingAlgorithm::SutherlandHodgman);
        let result = clipper_engine.clip_polygon(&subject, &clipper).unwrap();

        assert!(!result.is_empty());
    }

    #[test]
    fn test_circle_clipping() {
        let square = Polygon::closed(vec![
            Vec2::new(-1.0, -1.0),
            Vec2::new(1.0, -1.0),
            Vec2::new(1.0, 1.0),
            Vec2::new(-1.0, 1.0),
        ])
        .unwrap();

        let clipped = ClippingOperations::clip_to_circle(&square, Vec2::ZERO, 1.5).unwrap();

        assert_eq!(clipped.len(), 1);
        // All vertices should be within or on the circle
        for vertex in clipped[0].vertices() {
            assert!(vertex.distance(Vec2::ZERO) <= 1.5 + constants::EPSILON);
        }
    }

    #[test]
    fn test_cohen_sutherland_outcode() {
        let clipper = PolygonClipper::new(ClippingAlgorithm::CohenSutherland);
        let bounds = Bounds2D::from_points(Vec2::new(0.0, 0.0), Vec2::new(10.0, 10.0));

        assert_eq!(clipper.compute_outcode(Vec2::new(5.0, 5.0), bounds), 0); // Inside
        assert_eq!(clipper.compute_outcode(Vec2::new(-1.0, 5.0), bounds), 1); // Left
        assert_eq!(clipper.compute_outcode(Vec2::new(11.0, 5.0), bounds), 2); // Right
        assert_eq!(clipper.compute_outcode(Vec2::new(5.0, -1.0), bounds), 4); // Bottom
        assert_eq!(clipper.compute_outcode(Vec2::new(5.0, 11.0), bounds), 8); // Top
    }
}
