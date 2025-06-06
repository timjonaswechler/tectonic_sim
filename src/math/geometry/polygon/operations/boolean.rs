// src/math/geometry/polygon/operations/boolean.rs

use super::super::Polygon;
use super::clipping::{ClippingAlgorithm, PolygonClipper};
use crate::math::{error::*, types::*, utils::*};
use bevy::math::Vec2;

/// Boolean-Operationstypen
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BooleanOperation {
    /// Union (A ∪ B)
    Union,
    /// Intersection (A ∩ B)
    Intersection,
    /// Difference (A - B)
    Difference,
    /// Symmetric Difference (A ⊕ B)
    SymmetricDifference,
    /// Exclusive Or (gleich wie Symmetric Difference)
    Xor,
}

/// Boolean-Engine für Polygon-Operationen
pub struct PolygonBoolean {
    tolerance: f32,
    simplify_result: bool,
}

impl PolygonBoolean {
    /// Erstellt eine neue Boolean-Engine
    pub fn new() -> Self {
        Self {
            tolerance: constants::EPSILON * 1000.0,
            simplify_result: true,
        }
    }

    /// Setzt die Toleranz
    pub fn with_tolerance(mut self, tolerance: f32) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Aktiviert/deaktiviert Ergebnis-Vereinfachung
    pub fn simplify_result(mut self, simplify: bool) -> Self {
        self.simplify_result = simplify;
        self
    }

    /// Führt eine Boolean-Operation durch
    pub fn execute(
        &self,
        polygon_a: &Polygon,
        polygon_b: &Polygon,
        operation: BooleanOperation,
    ) -> MathResult<Vec<Polygon>> {
        match operation {
            BooleanOperation::Union => self.union(polygon_a, polygon_b),
            BooleanOperation::Intersection => self.intersection(polygon_a, polygon_b),
            BooleanOperation::Difference => self.difference(polygon_a, polygon_b),
            BooleanOperation::SymmetricDifference | BooleanOperation::Xor => {
                self.symmetric_difference(polygon_a, polygon_b)
            }
        }
    }

    // === Boolean-Operationen ===

    /// Union (Vereinigung)
    fn union(&self, polygon_a: &Polygon, polygon_b: &Polygon) -> MathResult<Vec<Polygon>> {
        // Prüfe triviale Fälle
        if !self.polygons_overlap(polygon_a, polygon_b)? {
            return Ok(vec![polygon_a.clone(), polygon_b.clone()]);
        }

        if self.polygon_contains_polygon(polygon_a, polygon_b)? {
            return Ok(vec![polygon_a.clone()]);
        }

        if self.polygon_contains_polygon(polygon_b, polygon_a)? {
            return Ok(vec![polygon_b.clone()]);
        }

        // Komplexe Union mit Sutherland-Hodgman und post-processing
        let mut result = self.complex_union(polygon_a, polygon_b)?;

        if self.simplify_result {
            result = self.simplify_polygons(result)?;
        }

        Ok(result)
    }

    /// Intersection (Schnittmenge)
    fn intersection(&self, polygon_a: &Polygon, polygon_b: &Polygon) -> MathResult<Vec<Polygon>> {
        // Verwende Clipping für Intersection
        let clipper = PolygonClipper::new(ClippingAlgorithm::SutherlandHodgman)
            .with_tolerance(self.tolerance);

        let result = clipper.clip_polygon(polygon_a, polygon_b)?;

        if self.simplify_result {
            self.simplify_polygons(result)
        } else {
            Ok(result)
        }
    }

    /// Difference (Differenz)
    fn difference(&self, polygon_a: &Polygon, polygon_b: &Polygon) -> MathResult<Vec<Polygon>> {
        // A - B = A ∩ (complement of B)
        // Verwende inverse Clipping
        let result = self.difference_operation(polygon_a, polygon_b)?;

        if self.simplify_result {
            self.simplify_polygons(result)
        } else {
            Ok(result)
        }
    }

    /// Symmetric Difference (Symmetrische Differenz)
    fn symmetric_difference(
        &self,
        polygon_a: &Polygon,
        polygon_b: &Polygon,
    ) -> MathResult<Vec<Polygon>> {
        // A ⊕ B = (A - B) ∪ (B - A)
        let diff_a_b = self.difference(polygon_a, polygon_b)?;
        let diff_b_a = self.difference(polygon_b, polygon_a)?;

        let mut result = diff_a_b;
        result.extend(diff_b_a);

        if self.simplify_result {
            self.simplify_polygons(result)
        } else {
            Ok(result)
        }
    }

    // === Komplexe Operationen ===

    fn complex_union(&self, polygon_a: &Polygon, polygon_b: &Polygon) -> MathResult<Vec<Polygon>> {
        // Vereinfachte Union-Implementierung
        // Für vollständige Implementierung wäre ein robuster Boolean-Algorithmus wie Martinez-Rueda nötig

        // 1. Finde alle Schnittpunkte
        let intersections = self.find_intersections(polygon_a, polygon_b)?;

        if intersections.is_empty() {
            // Keine Schnittpunkte: entweder getrennt oder eines enthält das andere
            return Ok(vec![polygon_a.clone(), polygon_b.clone()]);
        }

        // 2. Erstelle vereinigtes Polygon durch Konvexe Hülle als Approximation
        let mut all_vertices = polygon_a.vertices().to_vec();
        all_vertices.extend_from_slice(polygon_b.vertices());

        use super::convex_hull::{ConvexHullAlgorithm, ConvexHullComputer};
        let hull_computer = ConvexHullComputer::new(ConvexHullAlgorithm::AndrewMonotone);
        let hull = hull_computer.compute_hull(&all_vertices)?;

        Ok(vec![hull])
    }

    fn difference_operation(
        &self,
        polygon_a: &Polygon,
        polygon_b: &Polygon,
    ) -> MathResult<Vec<Polygon>> {
        // Vereinfachte Differenz-Implementierung
        let intersections = self.find_intersections(polygon_a, polygon_b)?;

        if intersections.is_empty() {
            // Keine Überschneidung
            if self.polygon_contains_polygon(polygon_b, polygon_a)? {
                return Ok(vec![]); // A ist komplett in B
            } else {
                return Ok(vec![polygon_a.clone()]); // Keine Überschneidung
            }
        }

        // Verwende inverse Clipping
        let vertices_a = polygon_a.vertices();
        let mut result_vertices = Vec::new();

        for &vertex in vertices_a {
            if !self.point_in_polygon(vertex, polygon_b)? {
                result_vertices.push(vertex);
            }
        }

        // Füge Schnittpunkte hinzu die auf dem Rand von A liegen
        for intersection in &intersections {
            if self.point_on_polygon_boundary(*intersection, polygon_a)? {
                result_vertices.push(*intersection);
            }
        }

        if result_vertices.len() >= 3 {
            Ok(vec![Polygon::new(result_vertices)?])
        } else {
            Ok(vec![])
        }
    }

    // === Hilfsfunktionen ===

    fn polygons_overlap(&self, polygon_a: &Polygon, polygon_b: &Polygon) -> MathResult<bool> {
        // Prüfe Bounding Box Überschneidung
        if let (Some(bounds_a), Some(bounds_b)) = (polygon_a.bounds(), polygon_b.bounds()) {
            let bounds_a = Bounds2D::from_points(bounds_a.min, bounds_a.max);
            let bounds_b = Bounds2D::from_points(bounds_b.min, bounds_b.max);

            if !bounds_a.intersects(&bounds_b) {
                return Ok(false);
            }
        }

        // Prüfe auf echte Überschneidung
        let intersections = self.find_intersections(polygon_a, polygon_b)?;
        Ok(!intersections.is_empty())
    }

    fn polygon_contains_polygon(&self, outer: &Polygon, inner: &Polygon) -> MathResult<bool> {
        use crate::math::geometry::polygon::PolygonProperties;

        // Prüfe ob alle Vertices von inner in outer liegen
        for &vertex in inner.vertices() {
            if !outer.contains_point(vertex) {
                return Ok(false);
            }
        }

        Ok(true)
    }

    fn point_in_polygon(&self, point: Vec2, polygon: &Polygon) -> MathResult<bool> {
        use crate::math::geometry::polygon::PolygonProperties;
        Ok(polygon.contains_point(point))
    }

    fn point_on_polygon_boundary(&self, point: Vec2, polygon: &Polygon) -> MathResult<bool> {
        let vertices = polygon.vertices();

        for i in 0..vertices.len() {
            let j = (i + 1) % vertices.len();
            let distance = self.point_line_distance(point, vertices[i], vertices[j]);

            if distance < self.tolerance {
                // Prüfe ob Punkt auf dem Liniensegment liegt
                let on_segment = self.point_on_line_segment(point, vertices[i], vertices[j]);
                if on_segment {
                    return Ok(true);
                }
            }
        }

        Ok(false)
    }

    fn point_line_distance(&self, point: Vec2, line_start: Vec2, line_end: Vec2) -> f32 {
        let line_vec = line_end - line_start;
        let point_vec = point - line_start;

        let line_length = line_vec.length();
        if line_length < self.tolerance {
            return point_vec.length();
        }

        let cross = line_vec.x * point_vec.y - line_vec.y * point_vec.x;
        cross.abs() / line_length
    }

    fn point_on_line_segment(&self, point: Vec2, segment_start: Vec2, segment_end: Vec2) -> bool {
        let segment_vec = segment_end - segment_start;
        let point_vec = point - segment_start;

        let dot = segment_vec.x * point_vec.x + segment_vec.y * point_vec.y;
        let segment_length_sq = segment_vec.x * segment_vec.x + segment_vec.y * segment_vec.y;

        if segment_length_sq < self.tolerance {
            return point.distance(segment_start) < self.tolerance;
        }

        let t = dot / segment_length_sq;
        t >= 0.0 && t <= 1.0
    }

    fn find_intersections(
        &self,
        polygon_a: &Polygon,
        polygon_b: &Polygon,
    ) -> MathResult<Vec<Vec2>> {
        let mut intersections = Vec::new();

        let vertices_a = polygon_a.vertices();
        let vertices_b = polygon_b.vertices();

        for i in 0..vertices_a.len() {
            let j = (i + 1) % vertices_a.len();
            let edge_a = (vertices_a[i], vertices_a[j]);

            for k in 0..vertices_b.len() {
                let l = (k + 1) % vertices_b.len();
                let edge_b = (vertices_b[k], vertices_b[l]);

                if let Some(intersection) = self.line_segment_intersection(edge_a, edge_b) {
                    intersections.push(intersection);
                }
            }
        }

        // Entferne Duplikate
        intersections.sort_by(|a, b| {
            a.x.partial_cmp(&b.x)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
        });

        intersections.dedup_by(|a, b| {
            (a.x - b.x).abs() < self.tolerance && (a.y - b.y).abs() < self.tolerance
        });

        Ok(intersections)
    }

    fn line_segment_intersection(
        &self,
        edge_a: (Vec2, Vec2),
        edge_b: (Vec2, Vec2),
    ) -> Option<Vec2> {
        let (p1, p2) = edge_a;
        let (p3, p4) = edge_b;

        let d1 = p2 - p1;
        let d2 = p4 - p3;

        let denominator = d1.x * d2.y - d1.y * d2.x;

        if denominator.abs() < self.tolerance {
            return None; // Parallele Linien
        }

        let t = ((p3.x - p1.x) * d2.y - (p3.y - p1.y) * d2.x) / denominator;
        let u = ((p3.x - p1.x) * d1.y - (p3.y - p1.y) * d1.x) / denominator;

        // Prüfe ob Schnittpunkt auf beiden Segmenten liegt
        if t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0 {
            Some(p1 + d1 * t)
        } else {
            None
        }
    }

    fn simplify_polygons(&self, polygons: Vec<Polygon>) -> MathResult<Vec<Polygon>> {
        // Entferne sehr kleine Polygone
        let min_area = self.tolerance * self.tolerance;
        let mut simplified = Vec::new();

        for polygon in polygons {
            use crate::math::geometry::polygon::PolygonProperties;

            if polygon.area() >= min_area {
                // Optional: weitere Vereinfachungen hier
                simplified.push(polygon);
            }
        }

        Ok(simplified)
    }
}

impl Default for PolygonBoolean {
    fn default() -> Self {
        Self::new()
    }
}

/// Vereinfachte API für häufige Boolean-Operationen
pub struct BooleanOperations;

impl BooleanOperations {
    /// Union von zwei Polygonen
    pub fn union(polygon_a: &Polygon, polygon_b: &Polygon) -> MathResult<Vec<Polygon>> {
        let boolean_engine = PolygonBoolean::new();
        boolean_engine.union(polygon_a, polygon_b)
    }

    /// Intersection von zwei Polygonen
    pub fn intersection(polygon_a: &Polygon, polygon_b: &Polygon) -> MathResult<Vec<Polygon>> {
        let boolean_engine = PolygonBoolean::new();
        boolean_engine.intersection(polygon_a, polygon_b)
    }

    /// Difference von zwei Polygonen
    pub fn difference(polygon_a: &Polygon, polygon_b: &Polygon) -> MathResult<Vec<Polygon>> {
        let boolean_engine = PolygonBoolean::new();
        boolean_engine.difference(polygon_a, polygon_b)
    }

    /// Symmetric Difference von zwei Polygonen
    pub fn symmetric_difference(
        polygon_a: &Polygon,
        polygon_b: &Polygon,
    ) -> MathResult<Vec<Polygon>> {
        let boolean_engine = PolygonBoolean::new();
        boolean_engine.symmetric_difference(polygon_a, polygon_b)
    }

    /// Batch Boolean-Operationen
    pub fn batch_operation(
        polygons: &[Polygon],
        operation: BooleanOperation,
    ) -> MathResult<Vec<Polygon>> {
        if polygons.is_empty() {
            return Ok(vec![]);
        }

        if polygons.len() == 1 {
            return Ok(vec![polygons[0].clone()]);
        }

        let boolean_engine = PolygonBoolean::new();
        let mut result = vec![polygons[0].clone()];

        for polygon in &polygons[1..] {
            let mut new_result = Vec::new();

            for existing in &result {
                let operation_result = boolean_engine.execute(existing, polygon, operation)?;
                new_result.extend(operation_result);
            }

            result = new_result;
        }

        Ok(result)
    }

    /// Optimierte Union für viele Polygone
    pub fn union_many(polygons: &[Polygon]) -> MathResult<Vec<Polygon>> {
        Self::batch_operation(polygons, BooleanOperation::Union)
    }
}

/// Statistiken für Boolean-Operationen
#[derive(Debug)]
pub struct BooleanStats {
    pub operation: BooleanOperation,
    pub input_area_a: f32,
    pub input_area_b: f32,
    pub result_area: f32,
    pub expected_area: f32,
    pub result_polygon_count: usize,
    pub area_efficiency: f32,
}

/// Boolean-Statistiken und Analyse
pub struct BooleanAnalysis;

impl BooleanAnalysis {
    /// Berechnet Boolean-Operationsstatistiken
    pub fn operation_statistics(
        polygon_a: &Polygon,
        polygon_b: &Polygon,
        result: &[Polygon],
        operation: BooleanOperation,
    ) -> BooleanStats {
        use crate::math::geometry::polygon::PolygonProperties;

        let area_a = polygon_a.area();
        let area_b = polygon_b.area();
        let result_area: f32 = result.iter().map(|p| p.area()).sum();

        let expected_area = match operation {
            BooleanOperation::Union => area_a + area_b, // Approximation
            BooleanOperation::Intersection => area_a.min(area_b), // Approximation
            BooleanOperation::Difference => area_a,     // Maximum possible
            BooleanOperation::SymmetricDifference | BooleanOperation::Xor => area_a + area_b,
        };

        BooleanStats {
            operation,
            input_area_a: area_a,
            input_area_b: area_b,
            result_area,
            expected_area,
            result_polygon_count: result.len(),
            area_efficiency: if expected_area > 0.0 {
                result_area / expected_area
            } else {
                0.0
            },
        }
    }

    /// Validiert Boolean-Operationsergebnis
    pub fn validate_result(
        polygon_a: &Polygon,
        polygon_b: &Polygon,
        result: &[Polygon],
        operation: BooleanOperation,
    ) -> Vec<String> {
        let mut issues = Vec::new();

        // Grundlegende Validierungen
        for (i, polygon) in result.iter().enumerate() {
            if polygon.len() < 3 {
                issues.push(format!("Result polygon {} has fewer than 3 vertices", i));
            }

            use crate::math::geometry::polygon::PolygonProperties;
            if polygon.area() < constants::EPSILON {
                issues.push(format!("Result polygon {} has zero or negative area", i));
            }
        }

        // Operationsspezifische Validierungen
        match operation {
            BooleanOperation::Intersection => {
                // Intersection sollte höchstens so groß sein wie die kleinere Eingabe
                use crate::math::geometry::polygon::PolygonProperties;
                let area_a = polygon_a.area();
                let area_b = polygon_b.area();
                for (i, polygon) in result.iter().enumerate() {
                    let result_area = polygon.area();
                    if result_area > area_a {
                        issues.push(format!("Result polygon {} has larger area than input A", i));
                    }
                    if result_area > area_b {
                        issues.push(format!("Result polygon {} has larger area than input B", i));
                    }
                }
            }
            BooleanOperation::Union => {
                // Union sollte mindestens so groß sein wie die größere Eingabe
                use crate::math::geometry::polygon::PolygonProperties;
                let area_a = polygon_a.area();
                let area_b = polygon_b.area();
                for (i, polygon) in result.iter().enumerate() {
                    let result_area = polygon.area();
                    if result_area < area_a && result_area < area_b {
                        issues.push(format!(
                            "Result polygon {} has smaller area than both inputs",
                            i
                        ));
                    }
                }
            }
            BooleanOperation::Difference => {
                // Difference sollte höchstens so groß sein wie A
                use crate::math::geometry::polygon::PolygonProperties;
                let area_a = polygon_a.area();
                for (i, polygon) in result.iter().enumerate() {
                    let result_area = polygon.area();
                    if result_area > area_a {
                        issues.push(format!("Result polygon {} has larger area than input A", i));
                    }
                }
            }
            BooleanOperation::SymmetricDifference | BooleanOperation::Xor => {
                // Symmetric Difference sollte kleiner als Summe der Eingaben sein
                use crate::math::geometry::polygon::PolygonProperties;
                let area_a = polygon_a.area();
                let area_b = polygon_b.area();
                let expected_area = area_a + area_b;

                for (i, polygon) in result.iter().enumerate() {
                    let result_area = polygon.area();
                    if result_area > expected_area {
                        issues.push(format!(
                            "Result polygon {} has larger area than expected",
                            i
                        ));
                    }
                }
            }
        }
        // Weitere Validierungen können hier hinzugefügt werden
        issues
    }
}
