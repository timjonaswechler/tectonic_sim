// src/math/geometry/polygon/smoothing/bezier.rs

use super::{
    super::Polygon,
    traits::{Smoothing, SmoothingParameters, SmoothingQuality, SmoothingResult, SmoothingUtils},
};
use crate::math::{error::*, types::*, utils::*};

/// Bézier-basierte Polygon-Glättung
#[derive(Debug, Clone)]
pub struct BezierSmoother {
    /// Verhältnis für Kontrollpunkt-Abstand (0.1 bis 0.9)
    pub control_point_ratio: f32,
    /// Anzahl der Samples pro Kurve
    pub samples_per_curve: usize,
    /// Spannungs-Parameter für Kontrollpunkte
    pub tension: f32,
    /// Ob Endpunkte bei offenen Polygonen fixiert werden sollen
    pub preserve_endpoints: bool,
}

impl Default for BezierSmoother {
    fn default() -> Self {
        Self {
            control_point_ratio: 0.3,
            samples_per_curve: 10,
            tension: 0.0,
            preserve_endpoints: true,
        }
    }
}

impl BezierSmoother {
    /// Erstellt einen neuen Bézier-Smoother
    pub fn new(control_point_ratio: f32) -> Self {
        Self {
            control_point_ratio: control_point_ratio.clamp(0.1, 0.9),
            ..Default::default()
        }
    }

    /// Setzt die Anzahl der Samples pro Kurve
    pub fn with_samples(mut self, samples: usize) -> Self {
        self.samples_per_curve = samples.max(3);
        self
    }

    /// Setzt den Spannungs-Parameter
    pub fn with_tension(mut self, tension: f32) -> Self {
        self.tension = tension.clamp(-1.0, 1.0);
        self
    }

    /// Aktiviert/deaktiviert Endpunkt-Erhaltung
    pub fn preserve_endpoints(mut self, preserve: bool) -> Self {
        self.preserve_endpoints = preserve;
        self
    }

    /// Berechnet Kontrollpunkte für einen Vertex basierend auf seinen Nachbarn
    fn calculate_control_points(
        &self,
        prev: Point2D,
        current: Point2D,
        next: Point2D,
    ) -> (Point2D, Point2D) {
        // Tangenten-Vektoren
        let incoming = current - prev;
        let outgoing = next - current;

        // Normalisierte Tangenten
        let incoming_len = incoming.length();
        let outgoing_len = outgoing.length();

        let incoming_normalized = if incoming_len > constants::EPSILON {
            incoming / incoming_len
        } else {
            Point2D::new(1.0, 0.0)
        };

        let outgoing_normalized = if outgoing_len > constants::EPSILON {
            outgoing / outgoing_len
        } else {
            Point2D::new(1.0, 0.0)
        };

        // Durchschnitts-Tangente für Smoothness
        let avg_tangent = (incoming_normalized + outgoing_normalized).normalize();

        // Skalierung basierend auf Segmentlängen
        let scale_in = incoming_len * self.control_point_ratio;
        let scale_out = outgoing_len * self.control_point_ratio;

        // Anwendung von Tension
        let tension_factor = 1.0 - self.tension.abs();

        let control_in = if self.tension >= 0.0 {
            current - avg_tangent * scale_in * tension_factor
        } else {
            current - incoming_normalized * scale_in * tension_factor
        };

        let control_out = if self.tension >= 0.0 {
            current + avg_tangent * scale_out * tension_factor
        } else {
            current + outgoing_normalized * scale_out * tension_factor
        };

        (control_in, control_out)
    }

    /// Erstellt Bézier-Kurven für alle Segmente
    fn create_bezier_curves(&self, polygon: &Polygon) -> MathResult<Vec<BezierCurve>> {
        let vertices = polygon.vertices();
        if vertices.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: vertices.len(),
            });
        }

        let mut curves = Vec::new();
        let n = vertices.len();
        let effective_n = if polygon.is_closed() && vertices.first() == vertices.last() {
            n - 1
        } else {
            n
        };

        for i in 0..effective_n {
            let current = vertices[i];
            let next = vertices[(i + 1) % effective_n];

            // Berechne Kontrollpunkte
            let prev = vertices[(i + effective_n - 1) % effective_n];
            let next_next = vertices[(i + 2) % effective_n];

            let (_, control_out) = self.calculate_control_points(prev, current, next);
            let (control_in, _) = self.calculate_control_points(current, next, next_next);

            curves.push(BezierCurve::cubic(current, control_out, control_in, next));
        }

        Ok(curves)
    }

    /// Sampelt alle Bézier-Kurven zu einem neuen Polygon
    fn sample_curves(&self, curves: &[BezierCurve]) -> Vec<Point2D> {
        let mut sampled_points = Vec::new();

        for (i, curve) in curves.iter().enumerate() {
            // Für die letzte Kurve bei geschlossenen Polygonen: überspringe den letzten Punkt
            let end_t = if i == curves.len() - 1 {
                1.0 - (1.0 / self.samples_per_curve as f32)
            } else {
                1.0
            };

            for j in 0..self.samples_per_curve {
                let t = j as f32 / self.samples_per_curve as f32;
                if t <= end_t {
                    sampled_points.push(curve.evaluate(t));
                }
            }
        }

        sampled_points
    }
}

impl Smoothing for BezierSmoother {
    fn smooth(&self, polygon: &Polygon) -> MathResult<Polygon> {
        if polygon.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: polygon.len(),
            });
        }

        // Spezialbehandlung für offene Polygone mit Endpunkt-Erhaltung
        if !polygon.is_closed() && self.preserve_endpoints {
            return self.smooth_open_polygon_preserving_endpoints(polygon);
        }

        let curves = self.create_bezier_curves(polygon)?;
        let smoothed_points = self.sample_curves(&curves);

        if polygon.is_closed() {
            Polygon::closed(smoothed_points)
        } else {
            Polygon::new(smoothed_points)
        }
    }

    fn supports_iterations(&self) -> bool {
        true
    }

    fn estimate_iterations(&self, polygon: &Polygon) -> usize {
        let smoothness = SmoothingUtils::calculate_smoothness(polygon);
        if smoothness < 0.3 {
            3
        } else if smoothness < 0.7 {
            2
        } else {
            1
        }
    }
}

impl BezierSmoother {
    /// Spezielle Behandlung für offene Polygone mit Endpunkt-Erhaltung
    fn smooth_open_polygon_preserving_endpoints(&self, polygon: &Polygon) -> MathResult<Polygon> {
        let vertices = polygon.vertices();
        if vertices.len() < 3 {
            return Ok(polygon.clone());
        }

        let mut smoothed_points = Vec::new();

        // Ersten Punkt beibehalten
        smoothed_points.push(vertices[0]);

        // Mittlere Segmente mit Bézier glätten
        if vertices.len() > 2 {
            let mut middle_vertices = vertices[1..vertices.len() - 1].to_vec();

            // Erstelle temporäres geschlossenes Polygon für die Mitte
            let temp_polygon = Polygon::new(middle_vertices)?;
            let curves = self.create_bezier_curves(&temp_polygon)?;
            let middle_smoothed = self.sample_curves(&curves);

            smoothed_points.extend(middle_smoothed);
        }

        // Letzten Punkt beibehalten
        smoothed_points.push(vertices[vertices.len() - 1]);

        Polygon::new(smoothed_points)
    }
}

/// Bézier-Kurve Repräsentation
#[derive(Debug, Clone)]
pub struct BezierCurve {
    control_points: Vec<Point2D>,
    curve_type: BezierType,
}

#[derive(Debug, Clone, Copy)]
pub enum BezierType {
    Quadratic,
    Cubic,
}

impl BezierCurve {
    /// Erstellt eine quadratische Bézier-Kurve
    pub fn quadratic(p0: Point2D, p1: Point2D, p2: Point2D) -> Self {
        Self {
            control_points: vec![p0, p1, p2],
            curve_type: BezierType::Quadratic,
        }
    }

    /// Erstellt eine kubische Bézier-Kurve
    pub fn cubic(p0: Point2D, p1: Point2D, p2: Point2D, p3: Point2D) -> Self {
        Self {
            control_points: vec![p0, p1, p2, p3],
            curve_type: BezierType::Cubic,
        }
    }

    /// Evaluiert die Kurve an Parameter t (0.0 bis 1.0)
    pub fn evaluate(&self, t: f32) -> Point2D {
        let t = t.clamp(0.0, 1.0);

        match self.curve_type {
            BezierType::Quadratic => self.evaluate_quadratic(t),
            BezierType::Cubic => self.evaluate_cubic(t),
        }
    }

    /// Evaluiert quadratische Bézier-Kurve
    fn evaluate_quadratic(&self, t: f32) -> Point2D {
        if self.control_points.len() != 3 {
            return Point2D::ZERO;
        }

        let u = 1.0 - t;
        let p0 = self.control_points[0];
        let p1 = self.control_points[1];
        let p2 = self.control_points[2];

        p0 * (u * u) + p1 * (2.0 * u * t) + p2 * (t * t)
    }

    /// Evaluiert kubische Bézier-Kurve
    fn evaluate_cubic(&self, t: f32) -> Point2D {
        if self.control_points.len() != 4 {
            return Point2D::ZERO;
        }

        let u = 1.0 - t;
        let tt = t * t;
        let uu = u * u;
        let uuu = uu * u;
        let ttt = tt * t;

        let p0 = self.control_points[0];
        let p1 = self.control_points[1];
        let p2 = self.control_points[2];
        let p3 = self.control_points[3];

        p0 * uuu + p1 * (3.0 * uu * t) + p2 * (3.0 * u * tt) + p3 * ttt
    }

    /// Berechnet die Tangente an Parameter t
    pub fn tangent(&self, t: f32) -> Point2D {
        let t = t.clamp(0.0, 1.0);

        match self.curve_type {
            BezierType::Quadratic => self.tangent_quadratic(t),
            BezierType::Cubic => self.tangent_cubic(t),
        }
    }

    fn tangent_quadratic(&self, t: f32) -> Point2D {
        if self.control_points.len() != 3 {
            return Point2D::new(1.0, 0.0);
        }

        let p0 = self.control_points[0];
        let p1 = self.control_points[1];
        let p2 = self.control_points[2];

        (p1 - p0) * (2.0 * (1.0 - t)) + (p2 - p1) * (2.0 * t)
    }

    fn tangent_cubic(&self, t: f32) -> Point2D {
        if self.control_points.len() != 4 {
            return Point2D::new(1.0, 0.0);
        }

        let u = 1.0 - t;
        let p0 = self.control_points[0];
        let p1 = self.control_points[1];
        let p2 = self.control_points[2];
        let p3 = self.control_points[3];

        (p1 - p0) * (3.0 * u * u) + (p2 - p1) * (6.0 * u * t) + (p3 - p2) * (3.0 * t * t)
    }

    /// Berechnet die Länge der Kurve (Näherung durch Sampling)
    pub fn length(&self, samples: usize) -> f32 {
        if samples < 2 {
            return 0.0;
        }

        let mut length = 0.0;
        let mut prev_point = self.evaluate(0.0);

        for i in 1..=samples {
            let t = i as f32 / samples as f32;
            let current_point = self.evaluate(t);
            length += prev_point.distance(current_point);
            prev_point = current_point;
        }

        length
    }
}

/// Erweiterte Bézier-Smoothing Varianten
pub struct BezierVariants;

impl BezierVariants {
    /// Sanftes Bézier-Smoothing
    pub fn gentle() -> BezierSmoother {
        BezierSmoother::new(0.2).with_tension(0.3)
    }

    /// Standard Bézier-Smoothing
    pub fn standard() -> BezierSmoother {
        BezierSmoother::default()
    }

    /// Aggressives Bézier-Smoothing
    pub fn aggressive() -> BezierSmoother {
        BezierSmoother::new(0.5).with_tension(-0.2)
    }

    /// Hochauflösendes Bézier-Smoothing
    pub fn high_resolution() -> BezierSmoother {
        BezierSmoother::default().with_samples(20)
    }

    /// Adaptives Bézier-Smoothing basierend auf Polygon-Eigenschaften
    pub fn adaptive(polygon: &Polygon) -> BezierSmoother {
        use crate::math::geometry::polygon::PolygonProperties;

        let vertex_count = polygon.len();
        let area = polygon.area();
        let perimeter = polygon.perimeter();
        let complexity = if perimeter > 0.0 {
            area / (perimeter * perimeter)
        } else {
            0.0
        };

        let control_ratio = if vertex_count < 10 {
            0.4 // Mehr Glättung für einfache Polygone
        } else if complexity < 0.01 {
            0.2 // Weniger Glättung für komplexe/elongierte Polygone
        } else {
            0.3 // Standard
        };

        let samples = if vertex_count > 50 {
            15 // Mehr Samples für komplexe Polygone
        } else {
            10
        };

        BezierSmoother::new(control_ratio).with_samples(samples)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bezier_curve_evaluation() {
        let curve = BezierCurve::cubic(
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(2.0, 1.0),
            Point2D::new(3.0, 0.0),
        );

        let start = curve.evaluate(0.0);
        let end = curve.evaluate(1.0);
        let mid = curve.evaluate(0.5);

        assert!(comparison::nearly_equal(start.x, 0.0));
        assert!(comparison::nearly_equal(end.x, 3.0));
        assert!(mid.x > 0.0 && mid.x < 3.0);
    }

    #[test]
    fn test_bezier_smoothing_basic() {
        let square = Polygon::closed(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(0.0, 1.0),
        ])
        .unwrap();

        let smoother = BezierSmoother::new(0.3);
        let smoothed = smoother.smooth(&square).unwrap();

        // Sollte mehr Vertices haben als das Original
        assert!(smoothed.len() > square.len());
        assert!(smoothed.is_closed());
    }

    #[test]
    fn test_bezier_variants() {
        let triangle = Polygon::closed(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(0.5, 1.0),
        ])
        .unwrap();

        let gentle = BezierVariants::gentle().smooth(&triangle).unwrap();
        let aggressive = BezierVariants::aggressive().smooth(&triangle).unwrap();
        let adaptive = BezierVariants::adaptive(&triangle)
            .smooth(&triangle)
            .unwrap();

        // Alle sollten gültig sein
        assert!(gentle.len() >= 3);
        assert!(aggressive.len() >= 3);
        assert!(adaptive.len() >= 3);
    }

    #[test]
    fn test_bezier_curve_length() {
        let line = BezierCurve::cubic(
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(2.0, 0.0),
            Point2D::new(3.0, 0.0),
        );

        let length = line.length(100);
        assert!((length - 3.0).abs() < 0.1); // Sollte ungefähr 3.0 sein
    }
}
