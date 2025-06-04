// src/math/geometry/polygon/smoothing/chaikin.rs

use super::super::Polygon;
use crate::math::{error::MathResult, types::*};

/// Chaikin-Subdivision Algorithmus für Polygon-Glättung
pub struct ChaikinSmoother {
    /// Anzahl der Iterationen
    pub iterations: usize,
    /// Schnitt-Verhältnis (normalerweise 0.25)
    pub cut_ratio: f32,
    /// Ob das Polygon nach dem Glätten vereinfacht werden soll
    pub simplify: bool,
}

impl Default for ChaikinSmoother {
    fn default() -> Self {
        Self {
            iterations: 1,
            cut_ratio: 0.25,
            simplify: false,
        }
    }
}

impl Smoothing for ChaikinSmoother {
    fn smooth(&self, polygon: &Polygon) -> MathResult<Polygon> {
        self.smooth(polygon)
    }
}

impl ChaikinSmoother {
    pub fn new(iterations: usize) -> Self {
        Self {
            iterations,
            ..Default::default()
        }
    }

    pub fn with_cut_ratio(mut self, ratio: f32) -> Self {
        self.cut_ratio = ratio.clamp(0.0, 0.5);
        self
    }

    pub fn with_simplification(mut self, simplify: bool) -> Self {
        self.simplify = simplify;
        self
    }

    /// Glättet ein Polygon mit Chaikin-Subdivision
    pub fn smooth(&self, polygon: &Polygon) -> MathResult<Polygon> {
        if polygon.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: polygon.len(),
            });
        }

        let mut current_vertices = polygon.vertices().to_vec();

        // Für geschlossene Polygone: letzten Punkt temporär entfernen
        let was_closed = polygon.is_closed();
        if was_closed && current_vertices.len() > 1 {
            if current_vertices.first() == current_vertices.last() {
                current_vertices.pop();
            }
        }

        // Iterative Glättung
        for _ in 0..self.iterations {
            current_vertices = self.chaikin_iteration(&current_vertices, was_closed);
        }

        // Polygon wieder zusammensetzen
        let result = if was_closed {
            Polygon::closed(current_vertices)?
        } else {
            Polygon::new(current_vertices)?
        };

        // Optional: Vereinfachung
        if self.simplify {
            self.simplify_polygon(result)
        } else {
            Ok(result)
        }
    }

    fn chaikin_iteration(&self, vertices: &[Point2D], is_closed: bool) -> Vec<Point2D> {
        if vertices.len() < 2 {
            return vertices.to_vec();
        }

        let mut smoothed = Vec::new();
        let n = vertices.len();

        for i in 0..n {
            let current = vertices[i];
            let next = if is_closed {
                vertices[(i + 1) % n]
            } else if i == n - 1 {
                break; // Letzter Punkt bei offenen Polygonen
            } else {
                vertices[i + 1]
            };

            // Zwei neue Punkte zwischen current und next
            let p1 = current + (next - current) * self.cut_ratio;
            let p2 = current + (next - current) * (1.0 - self.cut_ratio);

            smoothed.push(p1);
            smoothed.push(p2);
        }

        // Bei offenen Polygonen: ersten und letzten Original-Punkt beibehalten
        if !is_closed && !vertices.is_empty() {
            smoothed.insert(0, vertices[0]);
            smoothed.push(vertices[vertices.len() - 1]);
        }

        smoothed
    }

    fn simplify_polygon(&self, polygon: Polygon) -> MathResult<Polygon> {
        // Douglas-Peucker ähnliche Vereinfachung
        let simplified = self.douglas_peucker(polygon.vertices(), 0.1);

        if polygon.is_closed() {
            Polygon::closed(simplified)
        } else {
            Polygon::new(simplified)
        }
    }

    fn douglas_peucker(&self, points: &[Point2D], epsilon: f32) -> Vec<Point2D> {
        if points.len() <= 2 {
            return points.to_vec();
        }

        // Finde den Punkt mit dem größten Abstand zur Linie
        let mut max_distance = 0.0;
        let mut max_index = 0;

        let start = points[0];
        let end = points[points.len() - 1];

        for (i, point) in points.iter().enumerate().skip(1).take(points.len() - 2) {
            let distance = self.point_line_distance(*point, start, end);
            if distance > max_distance {
                max_distance = distance;
                max_index = i;
            }
        }

        // Wenn der maximale Abstand größer als epsilon ist, rekursiv unterteilen
        if max_distance > epsilon {
            let mut result1 = self.douglas_peucker(&points[0..=max_index], epsilon);
            let result2 = self.douglas_peucker(&points[max_index..], epsilon);

            result1.pop(); // Überschneidungspunkt entfernen
            result1.extend(result2);
            result1
        } else {
            vec![start, end]
        }
    }

    fn point_line_distance(&self, point: Point2D, line_start: Point2D, line_end: Point2D) -> f32 {
        let line_vec = line_end - line_start;
        let point_vec = point - line_start;

        let line_len = line_vec.length();
        if line_len < 1e-6 {
            return point_vec.length();
        }

        let line_unit = line_vec / line_len;
        let projection_length = point_vec.dot(line_unit);
        let projection = line_start + line_unit * projection_length;

        point.distance(projection)
    }
}
