// src/math/geometry/polygon/operations/triangulation.rs
use super::super::Polygon;
use crate::math::{error::*, utils::*};
use bevy::math::Vec2;
use std::collections::HashMap;

/// Verschiedene Triangulations-Algorithmen
#[derive(Debug, Clone, Copy)]
pub enum TriangulationAlgorithm {
    /// Ear Clipping (O(n²))
    EarClipping,
    /// Delaunay Triangulation (O(n log n))
    Delaunay,
    /// Monotone Polygon Triangulation (O(n))
    Monotone,
    /// Fan Triangulation (O(n), nur für konvexe Polygone)
    Fan,
}

/// Triangle representation
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Triangle {
    pub a: Vec2,
    pub b: Vec2,
    pub c: Vec2,
}

impl Triangle {
    pub fn new(a: Vec2, b: Vec2, c: Vec2) -> Self {
        Self { a, b, c }
    }

    /// Berechnet die Fläche des Dreiecks
    pub fn area(&self) -> f32 {
        0.5 * ((self.b.x - self.a.x) * (self.c.y - self.a.y)
            - (self.c.x - self.a.x) * (self.b.y - self.a.y))
            .abs()
    }

    /// Berechnet den Umkreis-Mittelpunkt
    pub fn circumcenter(&self) -> Option<Vec2> {
        let ax = self.a.x;
        let ay = self.a.y;
        let bx = self.b.x;
        let by = self.b.y;
        let cx = self.c.x;
        let cy = self.c.y;

        let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));

        if d.abs() < constants::EPSILON {
            return None; // Kollineare Punkte
        }

        let ux = ((ax * ax + ay * ay) * (by - cy)
            + (bx * bx + by * by) * (cy - ay)
            + (cx * cx + cy * cy) * (ay - by))
            / d;
        let uy = ((ax * ax + ay * ay) * (cx - bx)
            + (bx * bx + by * by) * (ax - cx)
            + (cx * cx + cy * cy) * (bx - ax))
            / d;

        Some(Vec2::new(ux, uy))
    }

    /// Prüft ob ein Punkt im Umkreis liegt
    pub fn point_in_circumcircle(&self, point: Vec2) -> bool {
        if let Some(center) = self.circumcenter() {
            let radius_sq = center.distance_squared(self.a);
            center.distance_squared(point) < radius_sq
        } else {
            false
        }
    }

    /// Prüft ob ein Punkt im Dreieck liegt
    pub fn contains_point(&self, point: Vec2) -> bool {
        let sign = |p1: Vec2, p2: Vec2, p3: Vec2| -> f32 {
            (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)
        };

        let d1 = sign(point, self.a, self.b);
        let d2 = sign(point, self.b, self.c);
        let d3 = sign(point, self.c, self.a);

        let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
        let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);

        !(has_neg && has_pos)
    }
}

/// Triangulations-Engine
pub struct PolygonTriangulator {
    algorithm: TriangulationAlgorithm,
    tolerance: f32,
}

impl PolygonTriangulator {
    /// Erstellt einen neuen Triangulator
    pub fn new(algorithm: TriangulationAlgorithm) -> Self {
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

    /// Trianguliert ein Polygon
    pub fn triangulate(&self, polygon: &Polygon) -> MathResult<Vec<Triangle>> {
        let vertices = polygon.vertices();

        if vertices.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: vertices.len(),
            });
        }

        // Entferne duplizierten letzten Punkt bei geschlossenen Polygonen
        let clean_vertices = if polygon.is_closed() && vertices.first() == vertices.last() {
            &vertices[..vertices.len() - 1]
        } else {
            vertices
        };

        if clean_vertices.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: clean_vertices.len(),
            });
        }

        match self.algorithm {
            TriangulationAlgorithm::EarClipping => self.ear_clipping(clean_vertices),
            TriangulationAlgorithm::Fan => self.fan_triangulation(clean_vertices),
            TriangulationAlgorithm::Delaunay => self.delaunay_triangulation(clean_vertices),
            TriangulationAlgorithm::Monotone => Err(MathError::InvalidConfiguration {
                message: "Monotone triangulation not yet implemented".to_string(),
            }),
        }
    }

    // === Algorithmus-Implementierungen ===

    /// Ear Clipping Algorithmus
    fn ear_clipping(&self, vertices: &[Vec2]) -> MathResult<Vec<Triangle>> {
        let mut triangles = Vec::new();
        let mut remaining_vertices = vertices.to_vec();

        while remaining_vertices.len() > 3 {
            let mut ear_found = false;

            for i in 0..remaining_vertices.len() {
                let prev = remaining_vertices
                    [(i + remaining_vertices.len() - 1) % remaining_vertices.len()];
                let curr = remaining_vertices[i];
                let next = remaining_vertices[(i + 1) % remaining_vertices.len()];

                if self.is_ear(prev, curr, next, &remaining_vertices) {
                    // Erstelle Dreieck
                    triangles.push(Triangle::new(prev, curr, next));

                    // Entferne den Ear-Vertex
                    remaining_vertices.remove(i);
                    ear_found = true;
                    break;
                }
            }

            if !ear_found {
                return Err(MathError::GeometricFailure {
                    operation: "No ear found in polygon (possibly self-intersecting)".to_string(),
                });
            }
        }

        // Letztes Dreieck
        if remaining_vertices.len() == 3 {
            triangles.push(Triangle::new(
                remaining_vertices[0],
                remaining_vertices[1],
                remaining_vertices[2],
            ));
        }

        Ok(triangles)
    }

    /// Prüft ob drei aufeinanderfolgende Vertices ein "Ear" bilden
    fn is_ear(&self, prev: Vec2, curr: Vec2, next: Vec2, vertices: &[Vec2]) -> bool {
        // 1. Prüfe ob das Dreieck konvex ist (Links-Kurve)
        let cross = (curr.x - prev.x) * (next.y - prev.y) - (curr.y - prev.y) * (next.x - prev.x);
        if cross <= 0.0 {
            return false; // Reflex vertex
        }

        // 2. Prüfe ob andere Vertices im Dreieck liegen
        let triangle = Triangle::new(prev, curr, next);

        for &vertex in vertices {
            if vertex == prev || vertex == curr || vertex == next {
                continue;
            }

            if triangle.contains_point(vertex) {
                return false;
            }
        }

        true
    }

    /// Fan Triangulation (für konvexe Polygone)
    fn fan_triangulation(&self, vertices: &[Vec2]) -> MathResult<Vec<Triangle>> {
        if vertices.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: vertices.len(),
            });
        }

        let mut triangles = Vec::new();
        let first_vertex = vertices[0];

        for i in 1..(vertices.len() - 1) {
            triangles.push(Triangle::new(first_vertex, vertices[i], vertices[i + 1]));
        }

        Ok(triangles)
    }

    /// Vereinfachte Delaunay Triangulation
    fn delaunay_triangulation(&self, vertices: &[Vec2]) -> MathResult<Vec<Triangle>> {
        if vertices.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: vertices.len(),
            });
        }

        // Für kleine Anzahl von Punkten: brute force
        if vertices.len() <= 4 {
            return self.fan_triangulation(vertices);
        }

        // Bowyer-Watson Algorithmus (vereinfacht)
        let mut triangles = Vec::new();

        // Super-Triangle das alle Punkte umschließt
        let bounds = self.compute_bounding_box(vertices);
        let super_triangle = self.create_super_triangle(bounds);
        triangles.push(super_triangle);

        // Füge Punkte einen nach dem anderen hinzu
        for &point in vertices {
            let mut bad_triangles = Vec::new();

            // Finde Dreiecke deren Umkreis den Punkt enthält
            for (i, triangle) in triangles.iter().enumerate() {
                if triangle.point_in_circumcircle(point) {
                    bad_triangles.push(i);
                }
            }

            // Sammle die Kanten der "bad triangles"
            let mut polygon_edges = Vec::new();
            for &i in &bad_triangles {
                let triangle = triangles[i];
                polygon_edges.push((triangle.a, triangle.b));
                polygon_edges.push((triangle.b, triangle.c));
                polygon_edges.push((triangle.c, triangle.a));
            }

            // Entferne doppelte Kanten
            polygon_edges.sort_by(|a, b| {
                a.0.x
                    .partial_cmp(&b.0.x)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| {
                        a.0.y
                            .partial_cmp(&b.0.y)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .then_with(|| {
                        a.1.x
                            .partial_cmp(&b.1.x)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .then_with(|| {
                        a.1.y
                            .partial_cmp(&b.1.y)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
            });

            let mut unique_edges = Vec::new();
            let mut i = 0;
            while i < polygon_edges.len() {
                let current_edge = polygon_edges[i];
                let reverse_edge = (current_edge.1, current_edge.0);

                if i + 1 < polygon_edges.len()
                    && (polygon_edges[i + 1] == current_edge
                        || polygon_edges[i + 1] == reverse_edge)
                {
                    i += 2; // Skip beide Kanten (doppelt)
                } else {
                    unique_edges.push(current_edge);
                    i += 1;
                }
            }

            // Entferne bad triangles (rückwärts iterieren für korrekte Indizes)
            bad_triangles.sort_by(|a, b| b.cmp(a));
            for &i in &bad_triangles {
                triangles.remove(i);
            }

            // Erstelle neue Dreiecke
            for edge in unique_edges {
                triangles.push(Triangle::new(edge.0, edge.1, point));
            }
        }

        // Entferne Dreiecke die Vertices des Super-Triangles enthalten
        let super_vertices = [super_triangle.a, super_triangle.b, super_triangle.c];
        triangles.retain(|triangle| {
            !super_vertices.contains(&triangle.a)
                && !super_vertices.contains(&triangle.b)
                && !super_vertices.contains(&triangle.c)
        });

        Ok(triangles)
    }

    fn compute_bounding_box(&self, vertices: &[Vec2]) -> (Vec2, Vec2) {
        let mut min = vertices[0];
        let mut max = vertices[0];

        for &vertex in vertices {
            min.x = min.x.min(vertex.x);
            min.y = min.y.min(vertex.y);
            max.x = max.x.max(vertex.x);
            max.y = max.y.max(vertex.y);
        }

        (min, max)
    }

    fn create_super_triangle(&self, bounds: (Vec2, Vec2)) -> Triangle {
        let (min, max) = bounds;
        let width = max.x - min.x;
        let height = max.y - min.y;
        let margin = (width + height) * 2.0;

        Triangle::new(
            Vec2::new(min.x - margin, min.y - margin),
            Vec2::new(max.x + margin, min.y - margin),
            Vec2::new((min.x + max.x) * 0.5, max.y + margin),
        )
    }
}

/// Triangulations-Utilities
pub struct TriangulationUtils;

impl TriangulationUtils {
    /// Berechnet die Gesamtfläche einer Triangulation
    pub fn total_area(triangles: &[Triangle]) -> f32 {
        triangles.iter().map(|t| t.area()).sum()
    }

    /// Findet das Dreieck das einen Punkt enthält
    pub fn find_containing_triangle(triangles: &[Triangle], point: Vec2) -> Option<usize> {
        triangles
            .iter()
            .position(|triangle| triangle.contains_point(point))
    }

    /// Konvertiert Triangulation zu Index-basierter Repräsentation
    pub fn to_indexed_triangulation(triangles: &[Triangle]) -> (Vec<Vec2>, Vec<[usize; 3]>) {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();
        let mut vertex_map = HashMap::new();

        for triangle in triangles {
            let mut triangle_indices = [0; 3];

            for (i, &vertex) in [triangle.a, triangle.b, triangle.c].iter().enumerate() {
                let key = ((vertex.x * 1e6) as i64, (vertex.y * 1e6) as i64);

                let index = *vertex_map.entry(key).or_insert_with(|| {
                    vertices.push(vertex);
                    vertices.len() - 1
                });

                triangle_indices[i] = index;
            }

            indices.push(triangle_indices);
        }

        (vertices, indices)
    }

    /// Erstellt Wireframe-Kanten aus Triangulation
    pub fn extract_edges(triangles: &[Triangle]) -> Vec<(Vec2, Vec2)> {
        let mut edges = Vec::new();

        for triangle in triangles {
            edges.push((triangle.a, triangle.b));
            edges.push((triangle.b, triangle.c));
            edges.push((triangle.c, triangle.a));
        }

        // Entferne Duplikate
        edges.sort_by(|a, b| {
            a.0.x
                .partial_cmp(&b.0.x)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    a.0.y
                        .partial_cmp(&b.0.y)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| {
                    a.1.x
                        .partial_cmp(&b.1.x)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| {
                    a.1.y
                        .partial_cmp(&b.1.y)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });
        edges.dedup();

        edges
    }

    /// Validiert eine Triangulation
    pub fn validate_triangulation(triangles: &[Triangle], original_area: f32) -> bool {
        let triangulation_area = Self::total_area(triangles);
        (triangulation_area - original_area).abs() < constants::EPSILON * 100.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::geometry::polygon::{Polygon, PolygonProperties};

    #[test]
    fn test_triangle_basic_operations() {
        let triangle = Triangle::new(
            Vec2::new(0.0, 0.0),
            Vec2::new(1.0, 0.0),
            Vec2::new(0.5, 1.0),
        );

        assert!(comparison::nearly_equal(triangle.area(), 0.5));
        assert!(triangle.contains_point(Vec2::new(0.5, 0.3)));
        assert!(!triangle.contains_point(Vec2::new(0.0, 1.0)));
    }

    #[test]
    fn test_ear_clipping_square() {
        let square = Polygon::closed(vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(1.0, 0.0),
            Vec2::new(1.0, 1.0),
            Vec2::new(0.0, 1.0),
        ])
        .unwrap();

        let triangulator = PolygonTriangulator::new(TriangulationAlgorithm::EarClipping);
        let triangles = triangulator.triangulate(&square).unwrap();

        assert_eq!(triangles.len(), 2); // Square should give 2 triangles

        let total_area = TriangulationUtils::total_area(&triangles);
        assert!(comparison::nearly_equal(total_area, 1.0));
    }

    #[test]
    fn test_fan_triangulation() {
        let pentagon = Polygon::closed(vec![
            Vec2::new(1.0, 0.0),
            Vec2::new(0.31, 0.95),
            Vec2::new(-0.81, 0.59),
            Vec2::new(-0.81, -0.59),
            Vec2::new(0.31, -0.95),
        ])
        .unwrap();

        let triangulator = PolygonTriangulator::new(TriangulationAlgorithm::Fan);
        let triangles = triangulator.triangulate(&pentagon).unwrap();

        assert_eq!(triangles.len(), 3); // Pentagon should give 3 triangles

        // Total area should be close to original polygon area
        let original_area = pentagon.area();
        let _triangulation_area = TriangulationUtils::total_area(&triangles);
        assert!(TriangulationUtils::validate_triangulation(
            &triangles,
            original_area
        ));
    }

    #[test]
    fn test_delaunay_triangulation() {
        let triangle = Polygon::closed(vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(2.0, 0.0),
            Vec2::new(1.0, 2.0),
        ])
        .unwrap();

        let triangulator = PolygonTriangulator::new(TriangulationAlgorithm::Delaunay);
        let triangles = triangulator.triangulate(&triangle).unwrap();

        assert_eq!(triangles.len(), 1); // Triangle should give 1 triangle
        assert!(comparison::nearly_equal(
            triangles[0].area(),
            triangle.area()
        ));
    }

    #[test]
    fn test_triangulation_utilities() {
        let triangles = vec![
            Triangle::new(
                Vec2::new(0.0, 0.0),
                Vec2::new(1.0, 0.0),
                Vec2::new(0.5, 1.0),
            ),
            Triangle::new(
                Vec2::new(1.0, 0.0),
                Vec2::new(2.0, 0.0),
                Vec2::new(1.5, 1.0),
            ),
        ];

        let total_area = TriangulationUtils::total_area(&triangles);
        assert!(comparison::nearly_equal(total_area, 1.0));

        let containing =
            TriangulationUtils::find_containing_triangle(&triangles, Vec2::new(0.5, 0.3));
        assert_eq!(containing, Some(0));

        let (vertices, indices) = TriangulationUtils::to_indexed_triangulation(&triangles);
        assert_eq!(indices.len(), 2);
        assert!(vertices.len() <= 6); // Maximum 6 unique vertices
    }
}
