// src/math/algorithms/triangulation.rs

use crate::math::{
    error::{MathError, MathResult},
    utils::constants,
};
use bevy::math::Vec2;
use std::collections::HashMap; // Für Delaunay-Hilfsstrukturen

/// Verschiedene Algorithmen zur Triangulation eines einfachen Polygons.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TriangulationAlgorithm {
    /// Ear Clipping Algorithmus (O(n^2) Komplexität, robust für einfache Polygone).
    EarClipping,
    /// Delaunay Triangulation (O(n log n)), erzeugt Dreiecke, die die Delaunay-Bedingung erfüllen.
    /// Die Implementierung hier wird eine vereinfachte Version sein. Für eine volle
    /// Delaunay-Triangulation einer Punktmenge (nicht nur eines Polygons) sind
    /// spezialisierte Bibliotheken wie `spade` besser geeignet.
    DelaunayApprox, // Umbenannt, um klarzustellen, dass es eine Approximation ist
    /// Monotone Polygon Triangulation (O(n) nach initialer Sortierung/Partitionierung).
    /// Komplexer zu implementieren. (Noch nicht implementiert)
    Monotone,
    /// Fan Triangulation (O(n)), nur für konvexe Polygone geeignet.
    Fan,
}

impl Default for TriangulationAlgorithm {
    fn default() -> Self {
        TriangulationAlgorithm::EarClipping // Guter Allrounder für einfache Polygone
    }
}

/// Repräsentiert ein einzelnes Dreieck, definiert durch drei 2D-Punkte.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Triangle {
    pub a: Vec2,
    pub b: Vec2,
    pub c: Vec2,
}

impl Triangle {
    pub fn new(p_a: Vec2, p_b: Vec2, p_c: Vec2) -> Self {
        Self {
            a: p_a,
            b: p_b,
            c: p_c,
        }
    }

    /// Berechnet die (vorzeichenbehaftete) doppelte Fläche des Dreiecks.
    /// Das Vorzeichen gibt die Orientierung an (z.B. positiv für CCW).
    fn signed_area_doubled(&self) -> f32 {
        (self.b.x - self.a.x) * (self.c.y - self.a.y)
            - (self.c.x - self.a.x) * (self.b.y - self.a.y)
    }

    /// Berechnet die Fläche des Dreiecks.
    pub fn area(&self) -> f32 {
        0.5 * self.signed_area_doubled().abs()
    }

    /// Berechnet den Umkreismittelpunkt des Dreiecks.
    /// Gibt `None` zurück, wenn die Punkte kollinear sind.
    pub fn circumcenter(&self) -> Option<Vec2> {
        // Formel aus Wikipedia: Circumscribed_circle#Cartesian_coordinates
        let d = 2.0
            * (self.a.x * (self.b.y - self.c.y)
                + self.b.x * (self.c.y - self.a.y)
                + self.c.x * (self.a.y - self.b.y));
        if d.abs() < constants::EPSILON * constants::EPSILON {
            // Noch kleinere Toleranz für Division
            return None; // Kollinear
        }

        let a_sq = self.a.length_squared();
        let b_sq = self.b.length_squared();
        let c_sq = self.c.length_squared();

        let ux = (a_sq * (self.b.y - self.c.y)
            + b_sq * (self.c.y - self.a.y)
            + c_sq * (self.a.y - self.b.y))
            / d;
        let uy = (a_sq * (self.c.x - self.b.x)
            + b_sq * (self.a.x - self.c.x)
            + c_sq * (self.b.x - self.a.x))
            / d;

        Some(Vec2::new(ux, uy))
    }

    /// Prüft, ob ein Punkt innerhalb des Umkreises dieses Dreiecks liegt.
    /// `point` ist der zu prüfende Punkt.
    pub fn is_point_in_circumcircle(&self, point: Vec2) -> bool {
        if let Some(center) = self.circumcenter() {
            // Vergleiche quadratische Abstände, um sqrt zu vermeiden
            let radius_sq = center.distance_squared(self.a);
            let dist_to_point_sq = center.distance_squared(point);
            // Punkt ist drin, wenn dist < radius (oder dist_sq < radius_sq)
            // Eine kleine Toleranz für Punkte genau auf dem Kreis
            dist_to_point_sq < radius_sq - constants::EPSILON
        } else {
            // Kollineare Punkte haben keinen wohldefinierten Umkreis in diesem Sinne
            false
        }
    }

    /// Prüft, ob ein Punkt innerhalb dieses Dreiecks liegt (inklusive der Kanten).
    /// Verwendet baryzentrische Koordinaten oder Kreuzprodukt-Tests.
    pub fn contains_point(&self, point: Vec2) -> bool {
        // Kreuzprodukt-Methode: Alle Kreuzprodukte (p,a,b), (p,b,c), (p,c,a) müssen dasselbe Vorzeichen haben.
        let sign_ab = (point.x - self.b.x) * (self.a.y - self.b.y)
            - (self.a.x - self.b.x) * (point.y - self.b.y);
        let sign_bc = (point.x - self.c.x) * (self.b.y - self.c.y)
            - (self.b.x - self.c.x) * (point.y - self.c.y);
        let sign_ca = (point.x - self.a.x) * (self.c.y - self.a.y)
            - (self.c.x - self.a.x) * (point.y - self.a.y);

        // Toleranz für Punkte auf der Kante
        let on_edge_tolerance = -constants::EPSILON;
        let all_non_positive = sign_ab <= -on_edge_tolerance
            && sign_bc <= -on_edge_tolerance
            && sign_ca <= -on_edge_tolerance;
        let all_non_negative = sign_ab >= on_edge_tolerance
            && sign_bc >= on_edge_tolerance
            && sign_ca >= on_edge_tolerance;

        all_non_positive || all_non_negative
    }
}

/// Führt die Triangulation eines einfachen Polygons (definiert durch seine Eckpunkte) durch.
pub struct PolygonTriangulator {
    algorithm: TriangulationAlgorithm,
    tolerance: f32,
}

impl Default for PolygonTriangulator {
    fn default() -> Self {
        Self {
            algorithm: TriangulationAlgorithm::default(),
            tolerance: constants::EPSILON * 10.0,
        }
    }
}

impl PolygonTriangulator {
    pub fn new(algorithm: TriangulationAlgorithm) -> Self {
        Self {
            algorithm,
            ..Default::default()
        }
    }

    pub fn with_tolerance(mut self, tolerance: f32) -> Self {
        self.tolerance = tolerance.max(0.0);
        self
    }

    /// Trianguliert ein Polygon, das durch eine geordnete Liste von Eckpunkten `polygon_vertices` definiert ist.
    /// Die Punkte sollten ein einfaches Polygon (nicht selbst-überschneidend) in CCW- oder CW-Orientierung bilden.
    /// Der letzte Punkt sollte NICHT mit dem ersten identisch sein (d.h. keine explizite Schließung).
    pub fn triangulate_points(&self, polygon_vertices: &[Vec2]) -> MathResult<Vec<Triangle>> {
        if polygon_vertices.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: polygon_vertices.len(),
            });
        }
        // TODO: Optional: Validierung auf Selbstüberschneidungen oder Duplikate hier.

        match self.algorithm {
            TriangulationAlgorithm::EarClipping => self.ear_clipping(polygon_vertices),
            TriangulationAlgorithm::Fan => self.fan_triangulation(polygon_vertices),
            TriangulationAlgorithm::DelaunayApprox => self.delaunay_approximation(polygon_vertices),
            TriangulationAlgorithm::Monotone => Err(MathError::TriangulationFailed {
                reason: "Monotone polygon triangulation is not yet implemented.".to_string(),
            }),
        }
    }

    // --- Ear Clipping Algorithmus ---
    fn ear_clipping(&self, initial_vertices: &[Vec2]) -> MathResult<Vec<Triangle>> {
        let mut triangles = Vec::new();
        let mut polygon_points = initial_vertices.to_vec();
        let n_initial = polygon_points.len();

        if n_initial < 3 {
            return Ok(triangles);
        } // Nichts zu triangulieren

        // Bestimme Orientierung, um konsistente "Links"-Turns zu identifizieren
        let mut signed_area_sum = 0.0;
        for i in 0..n_initial {
            let p1 = polygon_points[i];
            let p2 = polygon_points[(i + 1) % n_initial];
            signed_area_sum += (p1.x * p2.y) - (p2.x * p1.y);
        }
        let is_ccw = signed_area_sum > 0.0; // True für CCW, False für CW

        let mut iteration_limit = n_initial * n_initial; // Sicherheitslimit gegen Endlosschleifen

        while polygon_points.len() > 3 && iteration_limit > 0 {
            iteration_limit -= 1;
            let mut ear_found_this_pass = false;
            let current_n = polygon_points.len();

            for i in 0..current_n {
                let idx_prev = (i + current_n - 1) % current_n;
                let idx_curr = i;
                let idx_next = (i + 1) % current_n;

                let p_prev = polygon_points[idx_prev];
                let p_curr = polygon_points[idx_curr];
                let p_next = polygon_points[idx_next];

                // 1. Ist der aktuelle Vertex konvex (bildet einen "Links"-Turn in CCW, oder "Rechts"-Turn in CW)?
                // Kreuzprodukt (p_curr - p_prev) x (p_next - p_prev)
                let cross_z = (p_curr.x - p_prev.x) * (p_next.y - p_prev.y)
                    - (p_curr.y - p_prev.y) * (p_next.x - p_prev.x);

                let is_convex_turn = if is_ccw {
                    cross_z > self.tolerance
                } else {
                    cross_z < -self.tolerance
                };

                if is_convex_turn {
                    // 2. Enthält das Dreieck (p_prev, p_curr, p_next) andere Polygonpunkte?
                    let candidate_ear = Triangle::new(p_prev, p_curr, p_next);
                    let mut no_other_points_inside = true;
                    for k in 0..current_n {
                        if k == idx_prev || k == idx_curr || k == idx_next {
                            continue;
                        }
                        if candidate_ear.contains_point(polygon_points[k]) {
                            no_other_points_inside = false;
                            break;
                        }
                    }

                    if no_other_points_inside {
                        triangles.push(candidate_ear);
                        polygon_points.remove(idx_curr); // Entferne das "Ohr"
                        ear_found_this_pass = true;
                        break; // Starte die Suche nach dem nächsten Ohr von vorn
                    }
                }
            }
            if !ear_found_this_pass {
                return Err(MathError::TriangulationFailed {
                    reason: "No ear found during Ear Clipping (polygon might be self-intersecting or have other issues).".to_string(),
                });
            }
        }
        if iteration_limit == 0 {
            return Err(MathError::TriangulationFailed {
                reason: "Ear Clipping iteration limit reached.".to_string(),
            });
        }

        // Das verbleibende Polygon ist ein Dreieck
        if polygon_points.len() == 3 {
            triangles.push(Triangle::new(
                polygon_points[0],
                polygon_points[1],
                polygon_points[2],
            ));
        }
        Ok(triangles)
    }

    // --- Fan Triangulation (nur für konvexe Polygone) ---
    fn fan_triangulation(&self, vertices: &[Vec2]) -> MathResult<Vec<Triangle>> {
        // TODO: Füge eine Konvexitätsprüfung hinzu, wenn dieser Algorithmus robust sein soll.
        //       Oder verlange, dass der Aufrufer die Konvexität sicherstellt.
        let mut triangles = Vec::new();
        let p0 = vertices[0];
        for i in 1..(vertices.len() - 1) {
            triangles.push(Triangle::new(p0, vertices[i], vertices[i + 1]));
        }
        Ok(triangles)
    }

    // --- Delaunay Approximation (sehr vereinfacht, nicht robust) ---
    fn delaunay_approximation(&self, vertices: &[Vec2]) -> MathResult<Vec<Triangle>> {
        // Dies ist KEINE korrekte Delaunay-Triangulation einer Punktmenge, sondern
        // versucht, die Kanten eines Polygons zu erhalten und dann zu triangulieren.
        // Eine echte Delaunay-Triangulation würde die Punkte als solche behandeln.
        // Für eine robuste Implementierung sollte `spade` verwendet werden.
        // Diese Version ist eher eine Fan-Triangulation oder EarClipping-Fallback.
        if vertices.len() < 3 {
            return Ok(vec![]);
        }

        // Fallback zu EarClipping für diese Approximation
        self.ear_clipping(vertices)
    }
}

/// Hilfsfunktionen und Analysen für Triangulationen.
pub struct TriangulationUtils;

impl TriangulationUtils {
    /// Berechnet die Gesamtfläche aller Dreiecke in einer Triangulation.
    pub fn total_triangulation_area(triangles: &[Triangle]) -> f32 {
        triangles.iter().map(|t| t.area()).sum()
    }

    /// Findet den Index des Dreiecks, das einen gegebenen Punkt enthält.
    /// Gibt `None` zurück, wenn kein Dreieck den Punkt enthält.
    pub fn find_triangle_containing_point(triangles: &[Triangle], point: Vec2) -> Option<usize> {
        triangles
            .iter()
            .position(|triangle| triangle.contains_point(point))
    }

    /// Konvertiert eine Liste von Dreiecken in eine indexierte Repräsentation (Vertices und Indices).
    /// Nützlich für das Rendern oder für Datenstrukturen, die indexierte Meshes erwarten.
    /// Entfernt Vertex-Duplikate.
    pub fn to_indexed_triangles(triangles: &[Triangle]) -> (Vec<Vec2>, Vec<[usize; 3]>) {
        let mut unique_vertices = Vec::new();
        let mut triangle_indices = Vec::new();
        // HashMap zur Verfolgung bereits hinzugefügter Vertices und ihrer Indizes.
        // Konvertiere Vec2 in eine hashbare Form (z.B. (i64, i64) durch Multiplikation und Cast).
        let mut vertex_to_index_map: HashMap<(i64, i64), usize> = HashMap::new();

        let to_key = |v: Vec2| {
            ((v.x * 1e6) as i64, (v.y * 1e6) as i64) // Skalierung zur Erhaltung von Präzision als Key
        };

        for triangle in triangles {
            let mut current_triangle_v_indices = [0; 3];
            for (i, &vertex_coord) in [triangle.a, triangle.b, triangle.c].iter().enumerate() {
                let key = to_key(vertex_coord);
                let index = *vertex_to_index_map.entry(key).or_insert_with(|| {
                    unique_vertices.push(vertex_coord);
                    unique_vertices.len() - 1
                });
                current_triangle_v_indices[i] = index;
            }
            triangle_indices.push(current_triangle_v_indices);
        }
        (unique_vertices, triangle_indices)
    }

    /// Extrahiert alle einzigartigen Kanten aus einer Triangulation.
    pub fn extract_unique_edges(triangles: &[Triangle]) -> Vec<(Vec2, Vec2)> {
        let mut edges_set = std::collections::HashSet::new();
        for triangle in triangles {
            let edges_in_triangle = [
                (triangle.a, triangle.b),
                (triangle.b, triangle.c),
                (triangle.c, triangle.a),
            ];
            for edge in edges_in_triangle {
                // Normalisiere Kante (kleinerer Punkt zuerst), um Duplikate wie (A,B) und (B,A) zu finden.
                let (p1, p2) = if (edge.0.x, edge.0.y) < (edge.1.x, edge.1.y) {
                    (edge.0, edge.1)
                } else {
                    (edge.1, edge.0)
                };
                let key1 = ((p1.x * 1e6) as i64, (p1.y * 1e6) as i64);
                let key2 = ((p2.x * 1e6) as i64, (p2.y * 1e6) as i64);
                edges_set.insert((key1, key2));
            }
        }
        edges_set
            .into_iter()
            .map(|(k1, k2)| {
                (
                    Vec2::new(k1.0 as f32 / 1e6, k1.1 as f32 / 1e6),
                    Vec2::new(k2.0 as f32 / 1e6, k2.1 as f32 / 1e6),
                )
            })
            .collect()
    }

    /// Validiert, ob die Gesamtfläche der Triangulation (annähernd) der Fläche des
    /// Originalpolygons entspricht. `original_polygon_area` muss vom Aufrufer bereitgestellt werden.
    pub fn validate_triangulation_area(
        triangles: &[Triangle],
        original_polygon_area: f32,
        tolerance_factor: f32, // z.B. 0.01 für 1% Toleranz
    ) -> bool {
        let triangulation_area = Self::total_triangulation_area(triangles);
        let area_diff = (triangulation_area - original_polygon_area).abs();
        // Relative Toleranz oder absolute Toleranz für sehr kleine Flächen
        area_diff <= (original_polygon_area * tolerance_factor).max(constants::EPSILON * 100.0)
    }
}
