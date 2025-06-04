// src/math/geometry/polygon/smoothing/catmull_rom.rs

use super::super::Polygon;
use super::Smoothing;
use crate::math::{error::*, types::*, utils::*};

/// Catmull-Rom Spline Smoother für Polygone
pub struct CatmullRomSmoother {
    /// Anzahl der Interpolationspunkte zwischen jedem Vertex-Paar
    pub segments_per_edge: usize,
    /// Tension-Parameter (0.0 = sehr smooth, 1.0 = weniger smooth)
    pub tension: f32,
    /// Ob die Endpunkte bei offenen Polygonen fixiert bleiben sollen
    pub preserve_endpoints: bool,
    /// Alpha-Parameter für centripetal/chordal parameterization
    pub alpha: f32,
}

impl Default for CatmullRomSmoother {
    fn default() -> Self {
        Self {
            segments_per_edge: 10,
            tension: 0.0,
            preserve_endpoints: true,
            alpha: 0.5, // Centripetal parameterization
        }
    }
}

impl CatmullRomSmoother {
    /// Erstellt einen neuen Catmull-Rom Smoother
    pub fn new(segments_per_edge: usize) -> Self {
        Self {
            segments_per_edge,
            ..Default::default()
        }
    }

    /// Setzt den Tension-Parameter
    pub fn with_tension(mut self, tension: f32) -> Self {
        self.tension = tension.clamp(0.0, 1.0);
        self
    }

    /// Setzt die Alpha-Parametrisierung
    pub fn with_alpha(mut self, alpha: f32) -> Self {
        self.alpha = alpha.clamp(0.0, 1.0);
        self
    }

    /// Aktiviert/deaktiviert Endpunkt-Erhaltung
    pub fn preserve_endpoints(mut self, preserve: bool) -> Self {
        self.preserve_endpoints = preserve;
        self
    }

    /// Interpoliert zwischen vier Kontrollpunkten
    fn catmull_rom_interpolate(
        &self,
        p0: Point2D,
        p1: Point2D,
        p2: Point2D,
        p3: Point2D,
        t: f32,
    ) -> Point2D {
        // Standard Catmull-Rom mit Tension
        let t2 = t * t;
        let t3 = t2 * t;

        // Basis-Funktionen für Catmull-Rom
        let b0 = -self.tension * t + 2.0 * self.tension * t2 - self.tension * t3;
        let b1 = 1.0 + (self.tension - 3.0) * t2 + (2.0 - self.tension) * t3;
        let b2 = self.tension * t + (3.0 - 2.0 * self.tension) * t2 + (self.tension - 2.0) * t3;
        let b3 = -self.tension * t2 + self.tension * t3;

        p0 * b0 + p1 * b1 + p2 * b2 + p3 * b3
    }

    /// Berechnet Parameterisierungs-Werte für centripetal/chordal parameterization
    fn calculate_parameterization(
        &self,
        p0: Point2D,
        p1: Point2D,
        p2: Point2D,
        p3: Point2D,
    ) -> (f32, f32, f32, f32) {
        if self.alpha == 0.0 {
            // Uniform parameterization
            return (0.0, 1.0, 2.0, 3.0);
        }

        let d01 = p0.distance(p1).powf(self.alpha);
        let d12 = p1.distance(p2).powf(self.alpha);
        let d23 = p2.distance(p3).powf(self.alpha);

        (0.0, d01, d01 + d12, d01 + d12 + d23)
    }

    /// Centripetal/Chordal Catmull-Rom Interpolation
    fn centripetal_catmull_rom(
        &self,
        p0: Point2D,
        p1: Point2D,
        p2: Point2D,
        p3: Point2D,
        t: f32,
    ) -> Point2D {
        let (t0, t1, t2, t3) = self.calculate_parameterization(p0, p1, p2, p3);

        if (t2 - t1).abs() < constants::EPSILON {
            return p1; // Degenerate case
        }

        // Map t from [0,1] to [t1, t2]
        let t_mapped = t1 + t * (t2 - t1);

        // Lagrange interpolation
        let a1 = p0 * ((t1 - t_mapped) / (t1 - t0)) + p1 * ((t_mapped - t0) / (t1 - t0));
        let a2 = p1 * ((t2 - t_mapped) / (t2 - t1)) + p2 * ((t_mapped - t1) / (t2 - t1));
        let a3 = p2 * ((t3 - t_mapped) / (t3 - t2)) + p3 * ((t_mapped - t2) / (t3 - t2));

        let b1 = a1 * ((t2 - t_mapped) / (t2 - t0)) + a2 * ((t_mapped - t0) / (t2 - t0));
        let b2 = a2 * ((t3 - t_mapped) / (t3 - t1)) + a3 * ((t_mapped - t1) / (t3 - t1));

        b1 * ((t2 - t_mapped) / (t2 - t1)) + b2 * ((t_mapped - t1) / (t2 - t1))
    }

    /// Erzeugt Kontrollpunkte für offene Splines
    fn generate_open_control_points(&self, vertices: &[Point2D]) -> Vec<Point2D> {
        if vertices.len() < 2 {
            return vertices.to_vec();
        }

        let mut control_points = Vec::with_capacity(vertices.len() + 2);

        // Erster Kontrollpunkt (extrapoliert)
        if vertices.len() >= 2 {
            let extrapolated = vertices[0] - (vertices[1] - vertices[0]);
            control_points.push(extrapolated);
        } else {
            control_points.push(vertices[0]);
        }

        // Originale Vertices
        control_points.extend_from_slice(vertices);

        // Letzter Kontrollpunkt (extrapoliert)
        if vertices.len() >= 2 {
            let n = vertices.len();
            let extrapolated = vertices[n - 1] + (vertices[n - 1] - vertices[n - 2]);
            control_points.push(extrapolated);
        } else {
            control_points.push(vertices[vertices.len() - 1]);
        }

        control_points
    }

    /// Erzeugt Kontrollpunkte für geschlossene Splines
    fn generate_closed_control_points(&self, vertices: &[Point2D]) -> Vec<Point2D> {
        if vertices.len() < 3 {
            return vertices.to_vec();
        }

        let mut control_points = Vec::with_capacity(vertices.len() + 2);

        // Für geschlossene Splines: verwende die letzten/ersten Punkte als Kontrollpunkte
        let n = vertices.len();
        control_points.push(vertices[n - 1]); // Letzter Punkt als erster Kontrollpunkt
        control_points.extend_from_slice(vertices);
        control_points.push(vertices[0]); // Erster Punkt als letzter Kontrollpunkt

        control_points
    }
}

impl Smoothing for CatmullRomSmoother {
    fn smooth(&self, polygon: &Polygon) -> MathResult<Polygon> {
        let vertices = polygon.vertices();

        if vertices.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: vertices.len(),
            });
        }

        let mut smoothed_vertices = Vec::new();

        if polygon.is_closed() {
            self.smooth_closed_polygon(vertices, &mut smoothed_vertices)?;
        } else {
            self.smooth_open_polygon(vertices, &mut smoothed_vertices)?;
        }

        if polygon.is_closed() {
            Polygon::closed(smoothed_vertices)
        } else {
            Polygon::new(smoothed_vertices)
        }
    }
}

impl CatmullRomSmoother {
    fn smooth_closed_polygon(
        &self,
        vertices: &[Point2D],
        output: &mut Vec<Point2D>,
    ) -> MathResult<()> {
        let n = if vertices.first() == vertices.last() {
            vertices.len() - 1 // Remove duplicate last vertex
        } else {
            vertices.len()
        };

        if n < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: n,
            });
        }

        let control_points = self.generate_closed_control_points(&vertices[..n]);

        for i in 0..n {
            let p0 = control_points[i];
            let p1 = control_points[i + 1];
            let p2 = control_points[i + 2];
            let p3 = control_points[(i + 3) % (control_points.len())];

            // Interpoliere zwischen p1 und p2
            for j in 0..self.segments_per_edge {
                let t = j as f32 / self.segments_per_edge as f32;

                let interpolated = if self.alpha > 0.0 {
                    self.centripetal_catmull_rom(p0, p1, p2, p3, t)
                } else {
                    self.catmull_rom_interpolate(p0, p1, p2, p3, t)
                };

                output.push(interpolated);
            }
        }

        Ok(())
    }

    fn smooth_open_polygon(
        &self,
        vertices: &[Point2D],
        output: &mut Vec<Point2D>,
    ) -> MathResult<()> {
        if vertices.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: vertices.len(),
            });
        }

        let control_points = self.generate_open_control_points(vertices);
        let n = vertices.len();

        // Erster Punkt (optional: preserve endpoint)
        if self.preserve_endpoints {
            output.push(vertices[0]);
        }

        // Interpoliere zwischen allen Vertex-Paaren
        for i in 0..(n - 1) {
            let p0 = control_points[i];
            let p1 = control_points[i + 1];
            let p2 = control_points[i + 2];
            let p3 = control_points[i + 3];

            // Interpoliere zwischen p1 und p2
            let start_j = if self.preserve_endpoints && i == 0 {
                1
            } else {
                0
            };
            let end_j = if self.preserve_endpoints && i == n - 2 {
                self.segments_per_edge
            } else {
                self.segments_per_edge + 1
            };

            for j in start_j..end_j {
                let t = j as f32 / self.segments_per_edge as f32;

                let interpolated = if self.alpha > 0.0 {
                    self.centripetal_catmull_rom(p0, p1, p2, p3, t)
                } else {
                    self.catmull_rom_interpolate(p0, p1, p2, p3, t)
                };

                output.push(interpolated);
            }
        }

        // Letzter Punkt (optional: preserve endpoint)
        if self.preserve_endpoints {
            output.push(vertices[n - 1]);
        }

        Ok(())
    }
}

/// Erweiterte Catmull-Rom Varianten
pub struct CatmullRomVariants;

impl CatmullRomVariants {
    /// Uniform Catmull-Rom (Standard)
    pub fn uniform(segments_per_edge: usize) -> CatmullRomSmoother {
        CatmullRomSmoother::new(segments_per_edge).with_alpha(0.0)
    }

    /// Centripetal Catmull-Rom (weniger Loops und Überschwinger)
    pub fn centripetal(segments_per_edge: usize) -> CatmullRomSmoother {
        CatmullRomSmoother::new(segments_per_edge).with_alpha(0.5)
    }

    /// Chordal Catmull-Rom (noch weniger Artifacts)
    pub fn chordal(segments_per_edge: usize) -> CatmullRomSmoother {
        CatmullRomSmoother::new(segments_per_edge).with_alpha(1.0)
    }

    /// Tense Catmull-Rom (weniger smooth, mehr originalgetreu)
    pub fn tense(segments_per_edge: usize, tension: f32) -> CatmullRomSmoother {
        CatmullRomSmoother::new(segments_per_edge)
            .with_alpha(0.5)
            .with_tension(tension)
    }

    /// Adaptive Catmull-Rom (passt Segmente an Krümmung an)
    pub fn adaptive(base_segments: usize, max_segments: usize) -> AdaptiveCatmullRomSmoother {
        AdaptiveCatmullRomSmoother::new(base_segments, max_segments)
    }
}

/// Adaptive Catmull-Rom Smoother (variiert Segmentanzahl basierend auf Krümmung)
pub struct AdaptiveCatmullRomSmoother {
    base_segments: usize,
    max_segments: usize,
    curvature_threshold: f32,
    alpha: f32,
}

impl AdaptiveCatmullRomSmoother {
    pub fn new(base_segments: usize, max_segments: usize) -> Self {
        Self {
            base_segments,
            max_segments,
            curvature_threshold: 0.1,
            alpha: 0.5,
        }
    }

    pub fn with_curvature_threshold(mut self, threshold: f32) -> Self {
        self.curvature_threshold = threshold;
        self
    }

    /// Berechnet die Krümmung an einem Punkt
    fn calculate_curvature(&self, p0: Point2D, p1: Point2D, p2: Point2D) -> f32 {
        let v1 = p1 - p0;
        let v2 = p2 - p1;

        let v1_len = v1.length();
        let v2_len = v2.length();

        if v1_len < constants::EPSILON || v2_len < constants::EPSILON {
            return 0.0;
        }

        let cross = v1.x * v2.y - v1.y * v2.x;
        let dot = v1.x * v2.x + v1.y * v2.y;

        let angle = cross.atan2(dot).abs();
        angle / (v1_len + v2_len) * 2.0
    }

    /// Bestimmt die Anzahl der Segmente basierend auf der Krümmung
    fn adaptive_segment_count(&self, p0: Point2D, p1: Point2D, p2: Point2D, p3: Point2D) -> usize {
        let curvature1 = self.calculate_curvature(p0, p1, p2);
        let curvature2 = self.calculate_curvature(p1, p2, p3);
        let max_curvature = curvature1.max(curvature2);

        if max_curvature < self.curvature_threshold {
            self.base_segments
        } else {
            let factor = (max_curvature / self.curvature_threshold).min(3.0);
            ((self.base_segments as f32 * factor) as usize).min(self.max_segments)
        }
    }
}

impl Smoothing for AdaptiveCatmullRomSmoother {
    fn smooth(&self, polygon: &Polygon) -> MathResult<Polygon> {
        let vertices = polygon.vertices();

        if vertices.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: vertices.len(),
            });
        }

        let mut smoothed_vertices = Vec::new();
        let n = if polygon.is_closed() && vertices.first() == vertices.last() {
            vertices.len() - 1
        } else {
            vertices.len()
        };

        // Für geschlossene Polygone
        if polygon.is_closed() {
            for i in 0..n {
                let p0 = vertices[(i + n - 1) % n];
                let p1 = vertices[i];
                let p2 = vertices[(i + 1) % n];
                let p3 = vertices[(i + 2) % n];

                let segment_count = self.adaptive_segment_count(p0, p1, p2, p3);
                let smoother = CatmullRomSmoother::new(segment_count).with_alpha(self.alpha);

                for j in 0..segment_count {
                    let t = j as f32 / segment_count as f32;
                    let interpolated = smoother.centripetal_catmull_rom(p0, p1, p2, p3, t);
                    smoothed_vertices.push(interpolated);
                }
            }

            Polygon::closed(smoothed_vertices)
        } else {
            return Err(MathError::InvalidConfiguration {
                message: "Adaptive smoothing not yet implemented for open polygons".to_string(),
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::geometry::polygon::Polygon;

    #[test]
    fn test_catmull_rom_basic() {
        let square = Polygon::closed(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(0.0, 1.0),
        ])
        .unwrap();

        let smoother = CatmullRomSmoother::new(5);
        let smoothed = smoother.smooth(&square).unwrap();

        // Should have more vertices than original
        assert!(smoothed.len() > square.len());
        assert!(smoothed.is_closed());
    }

    #[test]
    fn test_catmull_rom_variants() {
        let triangle = Polygon::closed(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(0.5, 1.0),
        ])
        .unwrap();

        let uniform = CatmullRomVariants::uniform(8).smooth(&triangle).unwrap();
        let centripetal = CatmullRomVariants::centripetal(8)
            .smooth(&triangle)
            .unwrap();
        let chordal = CatmullRomVariants::chordal(8).smooth(&triangle).unwrap();

        // All should be valid and have similar vertex counts
        assert!(uniform.len() > triangle.len());
        assert!(centripetal.len() > triangle.len());
        assert!(chordal.len() > triangle.len());
    }

    #[test]
    fn test_adaptive_smoothing() {
        let jagged_polygon = Polygon::closed(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.1),
            Point2D::new(2.0, 0.0),
            Point2D::new(1.8, 1.0),
            Point2D::new(1.0, 0.9),
            Point2D::new(0.2, 1.0),
        ])
        .unwrap();

        let adaptive_smoother = CatmullRomVariants::adaptive(3, 12);
        let smoothed = adaptive_smoother.smooth(&jagged_polygon).unwrap();

        assert!(smoothed.len() > jagged_polygon.len());
        assert!(smoothed.is_closed());
    }
}
