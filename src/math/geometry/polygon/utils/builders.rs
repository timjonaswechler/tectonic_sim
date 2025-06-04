// src/math/geometry/polygon/utils/builders.rs

use super::super::Polygon;
use crate::math::{error::MathResult, types::*};
use std::f32::consts::TAU;

/// Builder für komplexe Polygon-Erstellung
pub struct PolygonBuilder {
    vertices: Vec<Point2D>,
    closed: bool,
}

impl PolygonBuilder {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            closed: false,
        }
    }

    pub fn closed(mut self) -> Self {
        self.closed = true;
        self
    }

    pub fn add_point(mut self, x: f32, y: f32) -> Self {
        self.vertices.push(Point2D::new(x, y));
        self
    }

    pub fn add_vertex(mut self, vertex: Point2D) -> Self {
        self.vertices.push(vertex);
        self
    }

    pub fn add_vertices(mut self, vertices: impl IntoIterator<Item = Point2D>) -> Self {
        self.vertices.extend(vertices);
        self
    }

    /// Fügt einen Kreis hinzu
    pub fn circle(mut self, center: Point2D, radius: f32, segments: usize) -> Self {
        let vertices = create_circle_vertices(center, radius, segments);
        self.vertices.extend(vertices);
        self
    }

    /// Fügt ein Rechteck hinzu
    pub fn rectangle(mut self, min: Point2D, max: Point2D) -> Self {
        self.vertices.extend([
            Point2D::new(min.x, min.y),
            Point2D::new(max.x, min.y),
            Point2D::new(max.x, max.y),
            Point2D::new(min.x, max.y),
        ]);
        self
    }

    /// Fügt einen Stern hinzu
    pub fn star(
        mut self,
        center: Point2D,
        inner_radius: f32,
        outer_radius: f32,
        points: usize,
    ) -> Self {
        let vertices = create_star_vertices(center, inner_radius, outer_radius, points);
        self.vertices.extend(vertices);
        self
    }

    /// Fügt eine Linie hinzu
    pub fn line_to(mut self, end: Point2D) -> Self {
        if let Some(&last) = self.vertices.last() {
            // Interpoliere zwischen letztem Punkt und end
            let direction = end - last;
            let distance = direction.length();
            let steps = (distance / 1.0).ceil() as usize; // Ein Punkt pro Einheit

            for i in 1..=steps {
                let t = i as f32 / steps as f32;
                self.vertices.push(last + direction * t);
            }
        } else {
            self.vertices.push(end);
        }
        self
    }

    /// Fügt eine Kurve hinzu (quadratische Bézier)
    pub fn curve_to(mut self, control: Point2D, end: Point2D, segments: usize) -> Self {
        if let Some(&start) = self.vertices.last() {
            for i in 1..=segments {
                let t = i as f32 / segments as f32;
                let point = quadratic_bezier(start, control, end, t);
                self.vertices.push(point);
            }
        } else {
            self.vertices.push(end);
        }
        self
    }

    /// Fügt eine kubische Bézier-Kurve hinzu
    pub fn cubic_curve_to(
        mut self,
        control1: Point2D,
        control2: Point2D,
        end: Point2D,
        segments: usize,
    ) -> Self {
        if let Some(&start) = self.vertices.last() {
            for i in 1..=segments {
                let t = i as f32 / segments as f32;
                let point = cubic_bezier(start, control1, control2, end, t);
                self.vertices.push(point);
            }
        } else {
            self.vertices.push(end);
        }
        self
    }

    /// Fügt einen Bogen hinzu
    pub fn arc_to(
        mut self,
        center: Point2D,
        radius: f32,
        start_angle: f32,
        end_angle: f32,
        segments: usize,
    ) -> Self {
        let angle_diff = end_angle - start_angle;

        for i in 0..=segments {
            let t = i as f32 / segments as f32;
            let angle = start_angle + angle_diff * t;
            let point = Point2D::new(
                center.x + radius * angle.cos(),
                center.y + radius * angle.sin(),
            );
            self.vertices.push(point);
        }
        self
    }

    /// Fügt eine Ellipse hinzu
    pub fn ellipse(mut self, center: Point2D, width: f32, height: f32, segments: usize) -> Self {
        let vertices = create_ellipse_vertices(center, width, height, segments);
        self.vertices.extend(vertices);
        self
    }

    /// Fügt ein regelmäßiges Polygon hinzu
    pub fn regular_polygon(mut self, center: Point2D, radius: f32, sides: usize) -> Self {
        let vertices = create_regular_polygon_vertices(center, radius, sides);
        self.vertices.extend(vertices);
        self
    }

    /// Skaliert alle bisherigen Vertices
    pub fn scale(mut self, factor: f32) -> Self {
        for vertex in &mut self.vertices {
            *vertex = *vertex * factor;
        }
        self
    }

    /// Skaliert mit verschiedenen Faktoren für X und Y
    pub fn scale_xy(mut self, factor_x: f32, factor_y: f32) -> Self {
        for vertex in &mut self.vertices {
            vertex.x *= factor_x;
            vertex.y *= factor_y;
        }
        self
    }

    /// Verschiebt alle Vertices um einen Offset
    pub fn translate(mut self, offset: Point2D) -> Self {
        for vertex in &mut self.vertices {
            *vertex = *vertex + offset;
        }
        self
    }

    /// Rotiert alle Vertices um einen Punkt
    pub fn rotate_around(mut self, center: Point2D, angle: f32) -> Self {
        let cos_a = angle.cos();
        let sin_a = angle.sin();

        for vertex in &mut self.vertices {
            let relative = *vertex - center;
            *vertex = center
                + Point2D::new(
                    relative.x * cos_a - relative.y * sin_a,
                    relative.x * sin_a + relative.y * cos_a,
                );
        }
        self
    }

    /// Spiegelt alle Vertices horizontal
    pub fn mirror_horizontal(mut self, axis_x: f32) -> Self {
        for vertex in &mut self.vertices {
            vertex.x = 2.0 * axis_x - vertex.x;
        }
        self
    }

    /// Spiegelt alle Vertices vertikal
    pub fn mirror_vertical(mut self, axis_y: f32) -> Self {
        for vertex in &mut self.vertices {
            vertex.y = 2.0 * axis_y - vertex.y;
        }
        self
    }

    pub fn build(self) -> MathResult<Polygon> {
        if self.closed {
            Polygon::closed(self.vertices)
        } else {
            Polygon::new(self.vertices)
        }
    }
}

impl Default for PolygonBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// Geometrie-Generatoren
pub fn create_circle_vertices(center: Point2D, radius: f32, segments: usize) -> Vec<Point2D> {
    let mut vertices = Vec::with_capacity(segments);

    for i in 0..segments {
        let angle = (i as f32 / segments as f32) * TAU;
        let x = center.x + radius * angle.cos();
        let y = center.y + radius * angle.sin();
        vertices.push(Point2D::new(x, y));
    }

    vertices
}

pub fn create_ellipse_vertices(
    center: Point2D,
    width: f32,
    height: f32,
    segments: usize,
) -> Vec<Point2D> {
    let mut vertices = Vec::with_capacity(segments);
    let half_width = width * 0.5;
    let half_height = height * 0.5;

    for i in 0..segments {
        let angle = (i as f32 / segments as f32) * TAU;
        let x = center.x + half_width * angle.cos();
        let y = center.y + half_height * angle.sin();
        vertices.push(Point2D::new(x, y));
    }

    vertices
}

pub fn create_star_vertices(
    center: Point2D,
    inner_radius: f32,
    outer_radius: f32,
    points: usize,
) -> Vec<Point2D> {
    let mut vertices = Vec::with_capacity(points * 2);

    for i in 0..(points * 2) {
        let angle = (i as f32 / (points * 2) as f32) * TAU;
        let radius = if i % 2 == 0 {
            outer_radius
        } else {
            inner_radius
        };

        let x = center.x + radius * angle.cos();
        let y = center.y + radius * angle.sin();
        vertices.push(Point2D::new(x, y));
    }

    vertices
}

pub fn create_regular_polygon_vertices(center: Point2D, radius: f32, sides: usize) -> Vec<Point2D> {
    create_circle_vertices(center, radius, sides)
}

/// Quadratische Bézier-Kurve
pub fn quadratic_bezier(p0: Point2D, p1: Point2D, p2: Point2D, t: f32) -> Point2D {
    let u = 1.0 - t;
    p0 * (u * u) + p1 * (2.0 * u * t) + p2 * (t * t)
}

/// Kubische Bézier-Kurve
pub fn cubic_bezier(p0: Point2D, p1: Point2D, p2: Point2D, p3: Point2D, t: f32) -> Point2D {
    let u = 1.0 - t;
    let tt = t * t;
    let uu = u * u;
    let uuu = uu * u;
    let ttt = tt * t;

    p0 * uuu + p1 * (3.0 * uu * t) + p2 * (3.0 * u * tt) + p3 * ttt
}

/// Erweiterte Geometrie-Generatoren
pub struct GeometryGenerators;

impl GeometryGenerators {
    /// Erstellt ein Rechteck mit abgerundeten Ecken
    pub fn rounded_rectangle(
        center: Point2D,
        width: f32,
        height: f32,
        corner_radius: f32,
        segments_per_corner: usize,
    ) -> Vec<Point2D> {
        let mut vertices = Vec::new();

        let half_w = width * 0.5;
        let half_h = height * 0.5;
        let r = corner_radius.min(half_w).min(half_h);

        // Ecken-Zentren
        let corners = [
            Point2D::new(center.x + half_w - r, center.y + half_h - r), // Oben rechts
            Point2D::new(center.x - half_w + r, center.y + half_h - r), // Oben links
            Point2D::new(center.x - half_w + r, center.y - half_h + r), // Unten links
            Point2D::new(center.x + half_w - r, center.y - half_h + r), // Unten rechts
        ];

        let start_angles = [0.0, TAU * 0.25, TAU * 0.5, TAU * 0.75];

        for (corner_center, &start_angle) in corners.iter().zip(start_angles.iter()) {
            for i in 0..segments_per_corner {
                let t = i as f32 / segments_per_corner as f32;
                let angle = start_angle + (TAU * 0.25) * t;
                vertices.push(Point2D::new(
                    corner_center.x + r * angle.cos(),
                    corner_center.y + r * angle.sin(),
                ));
            }
        }

        vertices
    }

    /// Erstellt eine Spirale
    pub fn spiral(
        center: Point2D,
        start_radius: f32,
        end_radius: f32,
        turns: f32,
        segments: usize,
    ) -> Vec<Point2D> {
        let mut vertices = Vec::with_capacity(segments);

        for i in 0..segments {
            let t = i as f32 / (segments - 1) as f32;
            let angle = turns * TAU * t;
            let radius = start_radius + (end_radius - start_radius) * t;

            vertices.push(Point2D::new(
                center.x + radius * angle.cos(),
                center.y + radius * angle.sin(),
            ));
        }

        vertices
    }

    /// Erstellt ein Polygon basierend auf einer Funktion
    pub fn function_polygon<F>(x_range: (f32, f32), func: F, segments: usize) -> Vec<Point2D>
    where
        F: Fn(f32) -> f32,
    {
        let mut vertices = Vec::with_capacity(segments);
        let (x_start, x_end) = x_range;

        for i in 0..segments {
            let t = i as f32 / (segments - 1) as f32;
            let x = x_start + (x_end - x_start) * t;
            let y = func(x);
            vertices.push(Point2D::new(x, y));
        }

        vertices
    }

    /// Erstellt ein Grid von Punkten
    pub fn grid(
        center: Point2D,
        width: f32,
        height: f32,
        cols: usize,
        rows: usize,
    ) -> Vec<Point2D> {
        let mut vertices = Vec::with_capacity(cols * rows);

        let cell_width = width / (cols - 1) as f32;
        let cell_height = height / (rows - 1) as f32;
        let start_x = center.x - width * 0.5;
        let start_y = center.y - height * 0.5;

        for row in 0..rows {
            for col in 0..cols {
                vertices.push(Point2D::new(
                    start_x + col as f32 * cell_width,
                    start_y + row as f32 * cell_height,
                ));
            }
        }

        vertices
    }

    /// Erstellt ein Herzform-Polygon
    pub fn heart(center: Point2D, size: f32, segments: usize) -> Vec<Point2D> {
        let mut vertices = Vec::with_capacity(segments);

        for i in 0..segments {
            let t = (i as f32 / segments as f32) * TAU;
            // Herz-Parametrierung
            let x = 16.0 * (t.sin().powi(3));
            let y =
                13.0 * t.cos() - 5.0 * (2.0 * t).cos() - 2.0 * (3.0 * t).cos() - (4.0 * t).cos();

            vertices.push(Point2D::new(
                center.x + x * size * 0.01,
                center.y + y * size * 0.01,
            ));
        }

        vertices
    }

    /// Erstellt ein Supershape (Gielis-Formel)
    pub fn supershape(
        center: Point2D,
        size: f32,
        m: f32,
        n1: f32,
        n2: f32,
        n3: f32,
        segments: usize,
    ) -> Vec<Point2D> {
        let mut vertices = Vec::with_capacity(segments);

        for i in 0..segments {
            let angle = (i as f32 / segments as f32) * TAU;

            // Gielis Supershape-Formel
            let part_a = ((m * angle * 0.25).cos().abs()).powf(n2);
            let part_b = ((m * angle * 0.25).sin().abs()).powf(n3);
            let r = (part_a + part_b).powf(-1.0 / n1);

            if r.is_finite() && r > 0.0 {
                let x = center.x + size * r * angle.cos();
                let y = center.y + size * r * angle.sin();
                vertices.push(Point2D::new(x, y));
            }
        }

        vertices
    }
}

/// Pattern-Generatoren
pub struct PatternGenerators;

impl PatternGenerators {
    /// Erstellt ein Hexagon-Grid Pattern
    pub fn hexagon_grid(center: Point2D, radius: f32, rings: usize) -> Vec<Point2D> {
        let mut vertices = Vec::new();

        // Zentrum
        vertices.push(center);

        // Für jeden Ring
        for ring in 1..=rings {
            let ring_radius = radius * ring as f32;
            let hexagons_in_ring = ring * 6;

            for i in 0..hexagons_in_ring {
                let angle = (i as f32 / hexagons_in_ring as f32) * TAU;
                vertices.push(Point2D::new(
                    center.x + ring_radius * angle.cos(),
                    center.y + ring_radius * angle.sin(),
                ));
            }
        }

        vertices
    }

    /// Erstellt ein Fibonacci-Spiral Pattern
    pub fn fibonacci_spiral(center: Point2D, radius: f32, points: usize) -> Vec<Point2D> {
        let mut vertices = Vec::with_capacity(points);
        let golden_ratio = (1.0 + 5.0_f32.sqrt()) * 0.5;

        for i in 0..points {
            let angle = i as f32 * golden_ratio * TAU;
            let r = radius * (i as f32 / points as f32).sqrt();

            vertices.push(Point2D::new(
                center.x + r * angle.cos(),
                center.y + r * angle.sin(),
            ));
        }

        vertices
    }

    /// Erstellt ein Voronoi-ähnliches zufälliges Pattern
    pub fn random_scattered(
        center: Point2D,
        radius: f32,
        points: usize,
        seed: Option<u64>,
    ) -> Vec<Point2D> {
        use rand::{Rng, SeedableRng};

        let mut rng = match seed {
            Some(s) => rand::rngs::StdRng::seed_from_u64(s),
            None => rand::rngs::StdRng::from_entropy(),
        };

        let mut vertices = Vec::with_capacity(points);

        for _ in 0..points {
            let angle = rng.gen_range(0.0..TAU);
            let r = radius * rng.r#gen::<f32>().sqrt(); // Gleichmäßige Verteilung im Kreis

            vertices.push(Point2D::new(
                center.x + r * angle.cos(),
                center.y + r * angle.sin(),
            ));
        }

        vertices
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polygon_builder_basic() {
        let polygon = PolygonBuilder::new()
            .add_point(0.0, 0.0)
            .add_point(1.0, 0.0)
            .add_point(1.0, 1.0)
            .add_point(0.0, 1.0)
            .closed()
            .build()
            .unwrap();

        assert_eq!(polygon.len(), 5); // 4 + geschlossener Punkt
        assert!(polygon.is_closed());
    }

    #[test]
    fn test_bezier_curves() {
        let start = Point2D::new(0.0, 0.0);
        let control = Point2D::new(1.0, 1.0);
        let end = Point2D::new(2.0, 0.0);

        let mid_point = quadratic_bezier(start, control, end, 0.5);
        assert!(mid_point.x > 0.0 && mid_point.x < 2.0);
        assert!(mid_point.y > 0.0);

        let cubic_point = cubic_bezier(start, control, control, end, 0.5);
        assert!(cubic_point.x > 0.0 && cubic_point.x < 2.0);
    }

    #[test]
    fn test_geometry_generators() {
        let circle = create_circle_vertices(Point2D::ZERO, 1.0, 8);
        assert_eq!(circle.len(), 8);

        let star = create_star_vertices(Point2D::ZERO, 0.5, 1.0, 5);
        assert_eq!(star.len(), 10); // 5 points * 2 (inner + outer)

        let heart = GeometryGenerators::heart(Point2D::ZERO, 1.0, 20);
        assert_eq!(heart.len(), 20);
    }

    #[test]
    fn test_transformations() {
        let original = PolygonBuilder::new()
            .rectangle(Point2D::new(-1.0, -1.0), Point2D::new(1.0, 1.0))
            .scale(2.0)
            .translate(Point2D::new(5.0, 5.0))
            .build()
            .unwrap();

        // Nach Skalierung und Translation sollten sich die Koordinaten entsprechend geändert haben
        let vertices = original.vertices();
        assert!(vertices.iter().all(|v| v.x >= 3.0 && v.x <= 7.0));
        assert!(vertices.iter().all(|v| v.y >= 3.0 && v.y <= 7.0));
    }
}
