// src/math/geometry/polygon/builder.rs

use crate::math::error::MathResult;
use crate::math::geometry::polygon::core::Polygon; // Korrekter Pfad zu Polygon
use crate::math::utils::constants::EPSILON;
use crate::math::utils::constants::TAU; // Korrekter Pfad zu TAU
use bevy::math::Vec2; // Vec2 direkt verwenden

/// Builder für komplexe Polygon-Erstellung durch Aneinanderreihung von geometrischen Formen.
/// Das Ergebnis der `build()` Methode ist ein `Polygon`.
pub struct PolygonBuilder {
    vertices: Vec<Vec2>,
    closed: bool,
}

impl PolygonBuilder {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            closed: false,
        }
    }

    /// Markiert das zu bauende Polygon als geschlossen.
    pub fn close_shape(mut self) -> Self {
        // Umbenannt von `closed` um Kollision mit Keyword zu vermeiden
        self.closed = true;
        self
    }

    pub fn add_point(mut self, x: f32, y: f32) -> Self {
        self.vertices.push(Vec2::new(x, y));
        self
    }

    pub fn add_vertex(mut self, vertex: Vec2) -> Self {
        self.vertices.push(vertex);
        self
    }

    pub fn add_vertices(mut self, vertices: impl IntoIterator<Item = Vec2>) -> Self {
        self.vertices.extend(vertices);
        self
    }

    /// Fügt einen Kreis hinzu
    pub fn circle(mut self, center: Vec2, radius: f32, segments: usize) -> Self {
        let vertices = ShapeGenerators::create_circle_vertices(center, radius, segments);
        self.vertices.extend(vertices);
        self
    }

    /// Fügt ein Rechteck hinzu
    pub fn rectangle(mut self, min: Vec2, max: Vec2) -> Self {
        self.vertices.extend([
            Vec2::new(min.x, min.y),
            Vec2::new(max.x, min.y),
            Vec2::new(max.x, max.y),
            Vec2::new(min.x, max.y),
        ]);
        self
    }

    /// Fügt einen Stern hinzu
    pub fn star(
        mut self,
        center: Vec2,
        inner_radius: f32,
        outer_radius: f32,
        points: usize,
    ) -> Self {
        let vertices =
            ShapeGenerators::create_star_vertices(center, inner_radius, outer_radius, points);
        self.vertices.extend(vertices);
        self
    }

    // ... (line_to, curve_to, cubic_curve_to, arc_to, ellipse, regular_polygon)
    // Diese Methoden bleiben ähnlich, müssen aber `ShapeGenerators` verwenden oder die Logik hierbehalten.
    // Beispiel für curve_to:
    pub fn curve_to(mut self, control: Vec2, end: Vec2, segments: usize) -> Self {
        if let Some(&start) = self.vertices.last() {
            for i in 1..=segments {
                let t = i as f32 / segments as f32;
                let point = ShapeGenerators::quadratic_bezier(start, control, end, t);
                self.vertices.push(point);
            }
        } else {
            // Wenn es keine Startpunkte gibt, füge Kontroll- und Endpunkte als Liniensegmente hinzu
            // oder entscheide, dass eine Kurve einen Startpunkt benötigt.
            // Hier fügen wir einfach den Endpunkt hinzu, um das Verhalten beizubehalten.
            self.vertices.push(control); // Eventuell nur 'end' hinzufügen oder Fehler werfen
            self.vertices.push(end);
        }
        self
    }

    // (scale, scale_xy, translate, rotate_around, mirror_horizontal, mirror_vertical)
    // bleiben logisch gleich, arbeiten auf self.vertices.

    pub fn scale(mut self, factor: f32) -> Self {
        for vertex in &mut self.vertices {
            *vertex *= factor;
        }
        self
    }

    pub fn scale_xy(mut self, factor_x: f32, factor_y: f32) -> Self {
        for vertex in &mut self.vertices {
            vertex.x *= factor_x;
            vertex.y *= factor_y;
        }
        self
    }

    pub fn translate(mut self, offset: Vec2) -> Self {
        for vertex in &mut self.vertices {
            *vertex += offset;
        }
        self
    }

    pub fn rotate_around(mut self, center: Vec2, angle_rad: f32) -> Self {
        let cos_a = angle_rad.cos();
        let sin_a = angle_rad.sin();

        for vertex in &mut self.vertices {
            let relative = *vertex - center;
            *vertex = center
                + Vec2::new(
                    relative.x * cos_a - relative.y * sin_a,
                    relative.x * sin_a + relative.y * cos_a,
                );
        }
        self
    }

    // ... weitere Transformationsmethoden ...

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

/// Geometrie-Generatoren, die Vec<Vec2> zurückgeben.
/// Diese könnten auch in `math/utils/shape_generators.rs` ausgelagert werden.
pub struct ShapeGenerators;

impl ShapeGenerators {
    pub fn create_circle_vertices(center: Vec2, radius: f32, segments: usize) -> Vec<Vec2> {
        let mut vertices = Vec::with_capacity(segments);
        for i in 0..segments {
            let angle = (i as f32 / segments as f32) * TAU;
            vertices.push(center + Vec2::new(radius * angle.cos(), radius * angle.sin()));
        }
        vertices
    }

    pub fn create_ellipse_vertices(
        center: Vec2,
        radius_x: f32, // Umbenannt von width für Klarheit
        radius_y: f32, // Umbenannt von height für Klarheit
        segments: usize,
    ) -> Vec<Vec2> {
        let mut vertices = Vec::with_capacity(segments);
        for i in 0..segments {
            let angle = (i as f32 / segments as f32) * TAU;
            vertices.push(center + Vec2::new(radius_x * angle.cos(), radius_y * angle.sin()));
        }
        vertices
    }

    pub fn create_star_vertices(
        center: Vec2,
        inner_radius: f32,
        outer_radius: f32,
        num_points: usize, // Umbenannt von points
    ) -> Vec<Vec2> {
        let mut vertices = Vec::with_capacity(num_points * 2);
        for i in 0..(num_points * 2) {
            let angle = (i as f32 / (num_points * 2) as f32) * TAU - TAU / 4.0; // Offset für Spitze nach oben
            let radius = if i % 2 == 0 {
                outer_radius
            } else {
                inner_radius
            };
            vertices.push(center + Vec2::new(radius * angle.cos(), radius * angle.sin()));
        }
        vertices
    }

    pub fn create_regular_polygon_vertices(center: Vec2, radius: f32, sides: usize) -> Vec<Vec2> {
        if sides < 3 {
            return Vec::new();
        } // Mindestens ein Dreieck
        Self::create_circle_vertices(center, radius, sides)
    }

    pub fn quadratic_bezier(p0: Vec2, p1: Vec2, p2: Vec2, t: f32) -> Vec2 {
        let u = 1.0 - t;
        p0 * (u * u) + p1 * (2.0 * u * t) + p2 * (t * t)
    }

    pub fn cubic_bezier(p0: Vec2, p1: Vec2, p2: Vec2, p3: Vec2, t: f32) -> Vec2 {
        let u = 1.0 - t;
        let tt = t * t;
        let uu = u * u;
        p0 * (uu * u) + p1 * (3.0 * uu * t) + p2 * (3.0 * u * tt) + p3 * (tt * t)
    }

    // Die erweiterten Generatoren (rounded_rectangle, spiral, etc.) können hier ebenfalls implementiert werden.
    // Beispiel für rounded_rectangle:
    pub fn rounded_rectangle(
        center: Vec2,
        width: f32,
        height: f32,
        corner_radius: f32,
        segments_per_corner: usize,
    ) -> Vec<Vec2> {
        let mut vertices = Vec::new();
        let half_w = width * 0.5;
        let half_h = height * 0.5;
        let r = corner_radius.min(half_w).min(half_h).max(0.0); // Sicherstellen, dass r nicht negativ ist

        if segments_per_corner == 0 || r < EPSILON {
            // Kein Radius oder keine Segmente -> normales Rechteck
            return vec![
                center + Vec2::new(-half_w, -half_h),
                center + Vec2::new(half_w, -half_h),
                center + Vec2::new(half_w, half_h),
                center + Vec2::new(-half_w, half_h),
            ];
        }

        // Eck-Zentren für die Bögen
        let corner_centers = [
            center + Vec2::new(half_w - r, half_h - r),  // Oben rechts
            center + Vec2::new(-half_w + r, half_h - r), // Oben links
            center + Vec2::new(-half_w + r, -half_h + r), // Unten links
            center + Vec2::new(half_w - r, -half_h + r), // Unten rechts
        ];

        // Startwinkel für die Bögen (mathematisch positiv, beginnend auf der positiven x-Achse)
        // Wir wollen bei der oberen rechten Ecke beginnen, die nach rechts und dann nach oben geht.
        // Also 1. Quadrant (0 bis PI/2), 2. Quadrant (PI/2 bis PI) etc.
        let start_angles = [
            0.0,        // Oben rechts (beginnt bei 0 rad, geht zu PI/2)
            TAU * 0.25, // Oben links (beginnt bei PI/2, geht zu PI)
            TAU * 0.5,  // Unten links (beginnt bei PI, geht zu 3PI/2)
            TAU * 0.75, // Unten rechts (beginnt bei 3PI/2, geht zu 2PI)
        ];

        for i in 0..4 {
            let corner_c = corner_centers[i];
            let start_angle_rad = start_angles[i];
            for j in 0..=segments_per_corner {
                // Inklusive Endpunkt des Bogens
                let t = j as f32 / segments_per_corner as f32;
                let angle = start_angle_rad + (TAU * 0.25) * t; // Jeder Bogen ist 90 Grad (PI/2 oder TAU/4)
                vertices.push(Vec2::new(
                    corner_c.x + r * angle.cos(),
                    corner_c.y + r * angle.sin(),
                ));
            }
        }
        // Duplikate entfernen, die durch ..=segments_per_corner entstehen könnten, wenn
        // der nächste Bogen am selben Punkt startet. Die Polygon::closed Logik wird das später behandeln.
        // Oder man macht `0..segments_per_corner` und fügt die geraden Linien dazwischen explizit ein.
        // Die aktuelle Logik fügt nur die Bogenpunkte hinzu, was zu einem "Blumen"-Polygon führt, wenn
        // die geraden Stücke nicht implizit durch die Polygonkonstruktion entstehen.
        // Um das korrekt zu machen, müssten wir die Reihenfolge ändern oder die geraden Stücke einfügen.

        // Korrigierte Logik für Rounded Rectangle:
        let mut correct_vertices = Vec::with_capacity(4 * (segments_per_corner + 1));
        // Unten rechts nach oben rechts
        for j in 0..=segments_per_corner {
            let t = j as f32 / segments_per_corner as f32;
            let angle = TAU * 0.75 + (TAU * 0.25) * t; // -PI/2 bis 0
            correct_vertices.push(corner_centers[3] + Vec2::new(r * angle.cos(), r * angle.sin()));
        }
        // Oben rechts nach oben links
        for j in 0..=segments_per_corner {
            let t = j as f32 / segments_per_corner as f32;
            let angle = 0.0 + (TAU * 0.25) * t; // 0 bis PI/2
            correct_vertices.push(corner_centers[0] + Vec2::new(r * angle.cos(), r * angle.sin()));
        }
        // Oben links nach unten links
        for j in 0..=segments_per_corner {
            let t = j as f32 / segments_per_corner as f32;
            let angle = TAU * 0.25 + (TAU * 0.25) * t; // PI/2 bis PI
            correct_vertices.push(corner_centers[1] + Vec2::new(r * angle.cos(), r * angle.sin()));
        }
        // Unten links nach unten rechts
        for j in 0..=segments_per_corner {
            let t = j as f32 / segments_per_corner as f32;
            let angle = TAU * 0.5 + (TAU * 0.25) * t; // PI bis 3PI/2
            correct_vertices.push(corner_centers[2] + Vec2::new(r * angle.cos(), r * angle.sin()));
        }
        // Dedup zur Sicherheit, falls Endpunkte der Bögen exakt aufeinander treffen
        correct_vertices.dedup_by(|a, b| a.distance_squared(*b) < 1e-9);
        correct_vertices
    }

    // Weitere Generatoren wie spiral, function_polygon, grid, heart, supershape
    // können hierhin verschoben oder neu implementiert werden,
    // wobei sie Vec<Vec2> zurückgeben.
}

// Die PatternGenerators könnten ähnlich hier platziert oder in ein
// eigenes Modul `math/utils/pattern_generators.rs` ausgelagert werden.
// Wichtig ist, dass sie `&mut impl Rng` oder `ResMut<SeedResource>` für Zufälligkeit nehmen.
