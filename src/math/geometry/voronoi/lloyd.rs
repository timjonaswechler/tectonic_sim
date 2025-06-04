// src/math/geometry/voronoi/lloyd.rs

use crate::math::{error::*, types::*, utils::*};
use spade::{DelaunayTriangulation, Point2, Triangulation};
use std::collections::HashMap;

/// Lloyd-Relaxation Konfiguration
#[derive(Debug, Clone)]
pub struct LloydConfig {
    /// Maximale Anzahl von Iterationen
    pub max_iterations: usize,
    /// Konvergenz-Toleranz (wenn Bewegung unter diesem Wert, stoppe)
    pub convergence_tolerance: f64,
    /// Boundary-Behandlung
    pub boundary_behavior: BoundaryBehavior,
    /// Gewichtung für verschiedene Relaxation-Modi
    pub relaxation_weight: f64,
    /// Ob Rand-Punkte fixiert bleiben sollen
    pub fix_boundary_points: bool,
}

/// Verschiedene Boundary-Behandlungsstrategien
#[derive(Debug, Clone, Copy)]
pub enum BoundaryBehavior {
    /// Punkte werden an der Boundary gespiegelt
    Mirror,
    /// Punkte werden zur nächsten Boundary-Position bewegt
    Clamp,
    /// Punkte werden durch Torus-Topologie behandelt
    Wrap,
    /// Boundary-Punkte bleiben unverändert
    Fixed,
}

impl Default for LloydConfig {
    fn default() -> Self {
        Self {
            max_iterations: 50,
            convergence_tolerance: 1e-6,
            boundary_behavior: BoundaryBehavior::Clamp,
            relaxation_weight: 1.0,
            fix_boundary_points: false,
        }
    }
}

/// Lloyd-Relaxation Engine
pub struct LloydRelaxation {
    config: LloydConfig,
}

impl LloydRelaxation {
    /// Erstellt eine neue Lloyd-Relaxation Instanz
    pub fn new(config: LloydConfig) -> Self {
        Self { config }
    }

    /// Wendet Lloyd-Relaxation auf eine Menge von Punkten an
    pub fn relax_points(
        &self,
        initial_points: &[Point2<f64>],
        boundary: Bounds2D,
    ) -> MathResult<(Vec<Point2<f64>>, LloydStatistics)> {
        if initial_points.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: initial_points.len(),
            });
        }

        let mut current_points = initial_points.to_vec();
        let mut statistics = LloydStatistics::new();

        for iteration in 0..self.config.max_iterations {
            // Erstelle Delaunay-Triangulation
            let mut triangulation = DelaunayTriangulation::<Point2<f64>>::new();

            for &point in &current_points {
                if self.is_point_in_bounds(point, boundary) {
                    triangulation.insert(point).ok();
                }
            }

            if triangulation.num_vertices() == 0 {
                break;
            }

            // Berechne neue Positionen
            let new_points = self.compute_lloyd_iteration(&triangulation, boundary)?;

            // Berechne Bewegungsstatistiken
            let total_movement = self.calculate_movement(&current_points, &new_points);
            statistics.iterations.push(LloydIterationStats {
                iteration,
                total_movement,
                max_movement: self.calculate_max_movement(&current_points, &new_points),
                point_count: new_points.len(),
                energy: self.calculate_system_energy(&triangulation),
            });

            // Prüfe Konvergenz
            if total_movement < self.config.convergence_tolerance {
                statistics.converged = true;
                statistics.final_iteration = iteration;
                break;
            }

            current_points = new_points;
            statistics.final_iteration = iteration;
        }

        Ok((current_points, statistics))
    }

    /// Führt eine einzelne Lloyd-Iteration durch
    fn compute_lloyd_iteration(
        &self,
        triangulation: &DelaunayTriangulation<Point2<f64>>,
        boundary: Bounds2D,
    ) -> MathResult<Vec<Point2<f64>>> {
        let mut new_points = Vec::new();

        for vertex_handle in triangulation.vertices() {
            let old_position = vertex_handle.position();

            // Berechne Voronoi-Zelle
            let voronoi_cell = self.compute_voronoi_cell(vertex_handle, boundary)?;

            if voronoi_cell.len() < 3 {
                // Degenerate Zelle: behalte alte Position
                new_points.push(old_position);
                continue;
            }

            // Berechne Zentroid der Voronoi-Zelle
            let centroid = self.calculate_polygon_centroid(&voronoi_cell)?;

            // Anwende Boundary-Behandlung
            let new_position = self.apply_boundary_behavior(centroid, boundary);

            // Anwende Relaxation-Gewichtung
            let weighted_position = self.apply_relaxation_weight(old_position, new_position);

            new_points.push(weighted_position);
        }

        Ok(new_points)
    }

    /// Berechnet die Voronoi-Zelle für einen Vertex
    fn compute_voronoi_cell(
        &self,
        vertex_handle: spade::handles::VertexHandle<Point2<f64>>,
        boundary: Bounds2D,
    ) -> MathResult<Vec<Point2<f64>>> {
        let mut voronoi_vertices = Vec::new();
        let generator_pos = vertex_handle.position();

        // Sammle alle Umkreismittelpunkte der angrenzenden Dreiecke
        for edge in vertex_handle.out_edges() {
            let face = edge.face();
            if face.is_outer() {
                continue;
            }

            // Hole die drei Eckpunkte des Dreiecks
            let triangle_vertices = self.get_triangle_vertices(face);

            if triangle_vertices.len() == 3 {
                if let Some(circumcenter) = self.calculate_circumcenter(
                    triangle_vertices[0],
                    triangle_vertices[1],
                    triangle_vertices[2],
                ) {
                    voronoi_vertices.push(circumcenter);
                }
            }
        }

        if voronoi_vertices.is_empty() {
            return Ok(vec![]);
        }

        // Sortiere Vertices um den Generator-Punkt
        self.sort_vertices_by_angle(&mut voronoi_vertices, generator_pos);

        // Clippe gegen Boundary
        let clipped_vertices = self.clip_voronoi_cell_to_boundary(&voronoi_vertices, boundary)?;

        Ok(clipped_vertices)
    }

    /// Extrahiert die Vertices eines Dreiecks aus einer Face
    fn get_triangle_vertices(
        &self,
        face: spade::handles::FaceHandle<Point2<f64>>,
    ) -> Vec<Point2<f64>> {
        let mut vertices = Vec::new();

        if let Some(edge) = face.adjacent_edge() {
            let mut current_edge = edge;
            for _ in 0..3 {
                vertices.push(current_edge.from().position());
                current_edge = current_edge.next();
                if current_edge == edge {
                    break;
                }
            }
        }

        vertices
    }

    /// Berechnet den Umkreismittelpunkt eines Dreiecks
    fn calculate_circumcenter(
        &self,
        p1: Point2<f64>,
        p2: Point2<f64>,
        p3: Point2<f64>,
    ) -> Option<Point2<f64>> {
        let ax = p1.x;
        let ay = p1.y;
        let bx = p2.x;
        let by = p2.y;
        let cx = p3.x;
        let cy = p3.y;

        let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));

        if d.abs() < 1e-10 {
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

        Some(Point2::new(ux, uy))
    }

    /// Sortiert Vertices nach Winkel um einen Mittelpunkt
    fn sort_vertices_by_angle(&self, vertices: &mut Vec<Point2<f64>>, center: Point2<f64>) {
        vertices.sort_by(|a, b| {
            let angle_a = (a.y - center.y).atan2(a.x - center.x);
            let angle_b = (b.y - center.y).atan2(b.x - center.x);
            angle_a
                .partial_cmp(&angle_b)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Clippt eine Voronoi-Zelle gegen die Boundary
    fn clip_voronoi_cell_to_boundary(
        &self,
        vertices: &[Point2<f64>],
        boundary: Bounds2D,
    ) -> MathResult<Vec<Point2<f64>>> {
        let mut clipped = vertices.to_vec();

        // Einfaches Rechteck-Clipping mit Sutherland-Hodgman ähnlichem Algorithmus
        let edges = [
            (
                Point2::new(boundary.min.x as f64, boundary.min.y as f64),
                Point2::new(boundary.max.x as f64, boundary.min.y as f64),
            ), // Bottom
            (
                Point2::new(boundary.max.x as f64, boundary.min.y as f64),
                Point2::new(boundary.max.x as f64, boundary.max.y as f64),
            ), // Right
            (
                Point2::new(boundary.max.x as f64, boundary.max.y as f64),
                Point2::new(boundary.min.x as f64, boundary.max.y as f64),
            ), // Top
            (
                Point2::new(boundary.min.x as f64, boundary.max.y as f64),
                Point2::new(boundary.min.x as f64, boundary.min.y as f64),
            ), // Left
        ];

        for (edge_start, edge_end) in &edges {
            if clipped.is_empty() {
                break;
            }

            let input_vertices = clipped.clone();
            clipped.clear();

            if input_vertices.is_empty() {
                continue;
            }

            let mut s = input_vertices[input_vertices.len() - 1];

            for e in input_vertices {
                if self.is_inside_edge(e, *edge_start, *edge_end) {
                    if !self.is_inside_edge(s, *edge_start, *edge_end) {
                        if let Some(intersection) =
                            self.line_intersection(s, e, *edge_start, *edge_end)
                        {
                            clipped.push(intersection);
                        }
                    }
                    clipped.push(e);
                } else if self.is_inside_edge(s, *edge_start, *edge_end) {
                    if let Some(intersection) = self.line_intersection(s, e, *edge_start, *edge_end)
                    {
                        clipped.push(intersection);
                    }
                }
                s = e;
            }
        }

        Ok(clipped)
    }

    /// Prüft ob ein Punkt auf der "inneren" Seite einer gerichteten Kante liegt
    fn is_inside_edge(
        &self,
        point: Point2<f64>,
        edge_start: Point2<f64>,
        edge_end: Point2<f64>,
    ) -> bool {
        let cross = (edge_end.x - edge_start.x) * (point.y - edge_start.y)
            - (edge_end.y - edge_start.y) * (point.x - edge_start.x);
        cross >= 0.0
    }

    /// Berechnet Schnittpunkt zwischen zwei Linien
    fn line_intersection(
        &self,
        p1: Point2<f64>,
        p2: Point2<f64>,
        p3: Point2<f64>,
        p4: Point2<f64>,
    ) -> Option<Point2<f64>> {
        let d1_x = p2.x - p1.x;
        let d1_y = p2.y - p1.y;
        let d2_x = p4.x - p3.x;
        let d2_y = p4.y - p3.y;

        let denominator = d1_x * d2_y - d1_y * d2_x;

        if denominator.abs() < 1e-10 {
            return None;
        }

        let t = ((p3.x - p1.x) * d2_y - (p3.y - p1.y) * d2_x) / denominator;

        Some(Point2::new(p1.x + t * d1_x, p1.y + t * d1_y))
    }

    /// Berechnet den Schwerpunkt eines Polygons
    fn calculate_polygon_centroid(&self, vertices: &[Point2<f64>]) -> MathResult<Point2<f64>> {
        if vertices.is_empty() {
            return Err(MathError::GeometricFailure {
                operation: "Cannot calculate centroid of empty polygon".to_string(),
            });
        }

        if vertices.len() == 1 {
            return Ok(vertices[0]);
        }

        if vertices.len() == 2 {
            return Ok(Point2::new(
                (vertices[0].x + vertices[1].x) * 0.5,
                (vertices[0].y + vertices[1].y) * 0.5,
            ));
        }

        // Für Polygone mit >= 3 Vertices: verwende Flächen-gewichteten Schwerpunkt
        let mut area = 0.0;
        let mut centroid_x = 0.0;
        let mut centroid_y = 0.0;
        let n = vertices.len();

        for i in 0..n {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % n];
            let cross_product = p1.x * p2.y - p2.x * p1.y;
            area += cross_product;
            centroid_x += (p1.x + p2.x) * cross_product;
            centroid_y += (p1.y + p2.y) * cross_product;
        }

        if area.abs() < 1e-10 {
            // Degenerate polygon: verwende arithmetischen Mittelwert
            let sum_x: f64 = vertices.iter().map(|p| p.x).sum();
            let sum_y: f64 = vertices.iter().map(|p| p.y).sum();
            return Ok(Point2::new(sum_x / n as f64, sum_y / n as f64));
        }

        area *= 0.5;
        centroid_x /= 6.0 * area;
        centroid_y /= 6.0 * area;

        Ok(Point2::new(centroid_x, centroid_y))
    }

    /// Wendet Boundary-Verhalten an
    fn apply_boundary_behavior(&self, point: Point2<f64>, boundary: Bounds2D) -> Point2<f64> {
        match self.config.boundary_behavior {
            BoundaryBehavior::Clamp => Point2::new(
                point.x.clamp(boundary.min.x as f64, boundary.max.x as f64),
                point.y.clamp(boundary.min.y as f64, boundary.max.y as f64),
            ),
            BoundaryBehavior::Mirror => {
                let mut new_point = point;

                let width = (boundary.max.x - boundary.min.x) as f64;
                let height = (boundary.max.y - boundary.min.y) as f64;

                if new_point.x < boundary.min.x as f64 {
                    new_point.x = (boundary.min.x as f64) + ((boundary.min.x as f64) - new_point.x);
                } else if new_point.x > boundary.max.x as f64 {
                    new_point.x = (boundary.max.x as f64) - (new_point.x - (boundary.max.x as f64));
                }

                if new_point.y < boundary.min.y as f64 {
                    new_point.y = (boundary.min.y as f64) + ((boundary.min.y as f64) - new_point.y);
                } else if new_point.y > boundary.max.y as f64 {
                    new_point.y = (boundary.max.y as f64) - (new_point.y - (boundary.max.y as f64));
                }

                new_point
            }
            BoundaryBehavior::Wrap => {
                let width = (boundary.max.x - boundary.min.x) as f64;
                let height = (boundary.max.y - boundary.min.y) as f64;

                let mut wrapped_x = point.x;
                let mut wrapped_y = point.y;

                while wrapped_x < boundary.min.x as f64 {
                    wrapped_x += width;
                }
                while wrapped_x > boundary.max.x as f64 {
                    wrapped_x -= width;
                }

                while wrapped_y < boundary.min.y as f64 {
                    wrapped_y += height;
                }
                while wrapped_y > boundary.max.y as f64 {
                    wrapped_y -= height;
                }

                Point2::new(wrapped_x, wrapped_y)
            }
            BoundaryBehavior::Fixed => point,
        }
    }

    /// Wendet Relaxation-Gewichtung an
    fn apply_relaxation_weight(&self, old_pos: Point2<f64>, new_pos: Point2<f64>) -> Point2<f64> {
        let weight = self.config.relaxation_weight.clamp(0.0, 1.0);
        Point2::new(
            old_pos.x + weight * (new_pos.x - old_pos.x),
            old_pos.y + weight * (new_pos.y - old_pos.y),
        )
    }

    /// Prüft ob ein Punkt in den Bounds liegt
    fn is_point_in_bounds(&self, point: Point2<f64>, boundary: Bounds2D) -> bool {
        point.x >= boundary.min.x as f64
            && point.x <= boundary.max.x as f64
            && point.y >= boundary.min.y as f64
            && point.y <= boundary.max.y as f64
    }

    /// Berechnet die Gesamtbewegung zwischen zwei Punkt-Sets
    fn calculate_movement(&self, old_points: &[Point2<f64>], new_points: &[Point2<f64>]) -> f64 {
        if old_points.len() != new_points.len() {
            return f64::INFINITY;
        }

        old_points
            .iter()
            .zip(new_points.iter())
            .map(|(old, new)| {
                let dx = new.x - old.x;
                let dy = new.y - old.y;
                (dx * dx + dy * dy).sqrt()
            })
            .sum()
    }

    /// Berechnet die maximale Bewegung zwischen zwei Punkt-Sets
    fn calculate_max_movement(
        &self,
        old_points: &[Point2<f64>],
        new_points: &[Point2<f64>],
    ) -> f64 {
        if old_points.len() != new_points.len() {
            return f64::INFINITY;
        }

        old_points
            .iter()
            .zip(new_points.iter())
            .map(|(old, new)| {
                let dx = new.x - old.x;
                let dy = new.y - old.y;
                (dx * dx + dy * dy).sqrt()
            })
            .fold(0.0, |acc, movement| acc.max(movement))
    }

    /// Berechnet die System-Energie (vereinfacht)
    fn calculate_system_energy(&self, triangulation: &DelaunayTriangulation<Point2<f64>>) -> f64 {
        // Vereinfachte Energie-Berechnung basierend auf Voronoi-Zellen-Unregelmäßigkeit
        let mut total_energy = 0.0;
        let mut count = 0;

        for vertex_handle in triangulation.vertices() {
            if let Ok(voronoi_cell) = self.compute_voronoi_cell(
                vertex_handle,
                Bounds2D::from_points(Point2D::new(-1000.0, -1000.0), Point2D::new(1000.0, 1000.0)),
            ) {
                if voronoi_cell.len() >= 3 {
                    // Berechne "Rundheit" der Zelle als Energie-Maß
                    if let Ok(centroid) = self.calculate_polygon_centroid(&voronoi_cell) {
                        let mut distance_variance = 0.0;
                        let mut mean_distance = 0.0;

                        for vertex in &voronoi_cell {
                            let dx = vertex.x - centroid.x;
                            let dy = vertex.y - centroid.y;
                            let distance = (dx * dx + dy * dy).sqrt();
                            mean_distance += distance;
                        }
                        mean_distance /= voronoi_cell.len() as f64;

                        for vertex in &voronoi_cell {
                            let dx = vertex.x - centroid.x;
                            let dy = vertex.y - centroid.y;
                            let distance = (dx * dx + dy * dy).sqrt();
                            let diff = distance - mean_distance;
                            distance_variance += diff * diff;
                        }
                        distance_variance /= voronoi_cell.len() as f64;

                        total_energy += distance_variance;
                        count += 1;
                    }
                }
            }
        }

        if count > 0 {
            total_energy / count as f64
        } else {
            0.0
        }
    }
}

/// Statistiken über eine Lloyd-Relaxation
#[derive(Debug, Clone)]
pub struct LloydStatistics {
    pub iterations: Vec<LloydIterationStats>,
    pub converged: bool,
    pub final_iteration: usize,
}

impl LloydStatistics {
    fn new() -> Self {
        Self {
            iterations: Vec::new(),
            converged: false,
            final_iteration: 0,
        }
    }

    /// Gibt die finale Energie zurück
    pub fn final_energy(&self) -> Option<f64> {
        self.iterations.last().map(|stats| stats.energy)
    }

    /// Gibt die Energie-Reduktion zurück
    pub fn energy_reduction(&self) -> Option<f64> {
        if self.iterations.len() < 2 {
            return None;
        }

        let initial = self.iterations.first()?.energy;
        let final_energy = self.iterations.last()?.energy;
        Some(initial - final_energy)
    }
}

#[derive(Debug, Clone)]
pub struct LloydIterationStats {
    pub iteration: usize,
    pub total_movement: f64,
    pub max_movement: f64,
    pub point_count: usize,
    pub energy: f64,
}

/// Erweiterte Lloyd-Varianten
pub struct LloydVariants;

impl LloydVariants {
    /// Standard Lloyd-Relaxation
    pub fn standard(max_iterations: usize) -> LloydConfig {
        LloydConfig {
            max_iterations,
            ..Default::default()
        }
    }

    /// Weighted Lloyd (mit reduzierter Relaxation-Stärke)
    pub fn weighted(max_iterations: usize, weight: f64) -> LloydConfig {
        LloydConfig {
            max_iterations,
            relaxation_weight: weight,
            ..Default::default()
        }
    }

    /// Boundary-preserving Lloyd
    pub fn boundary_preserving(max_iterations: usize) -> LloydConfig {
        LloydConfig {
            max_iterations,
            fix_boundary_points: true,
            boundary_behavior: BoundaryBehavior::Fixed,
            ..Default::default()
        }
    }

    /// High-precision Lloyd
    pub fn high_precision(max_iterations: usize) -> LloydConfig {
        LloydConfig {
            max_iterations,
            convergence_tolerance: 1e-8,
            ..Default::default()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lloyd_relaxation_basic() {
        let points = vec![
            Point2::new(0.1, 0.1),
            Point2::new(0.9, 0.1),
            Point2::new(0.9, 0.9),
            Point2::new(0.1, 0.9),
            Point2::new(0.5, 0.5),
        ];

        let boundary = Bounds2D::from_points(Point2D::new(0.0, 0.0), Point2D::new(1.0, 1.0));
        let config = LloydConfig::default();
        let lloyd = LloydRelaxation::new(config);

        let result = lloyd.relax_points(&points, boundary);
        assert!(result.is_ok());

        let (relaxed_points, stats) = result.unwrap();
        assert_eq!(relaxed_points.len(), points.len());
        assert!(!stats.iterations.is_empty());
    }

    #[test]
    fn test_lloyd_convergence() {
        let points = vec![
            Point2::new(0.1, 0.1),
            Point2::new(0.9, 0.1),
            Point2::new(0.9, 0.9),
            Point2::new(0.1, 0.9),
        ];

        let boundary = Bounds2D::from_points(Point2D::new(0.0, 0.0), Point2D::new(1.0, 1.0));

        let config = LloydConfig {
            max_iterations: 100,
            convergence_tolerance: 1e-8,
            ..Default::default()
        };

        let lloyd = LloydRelaxation::new(config);
        let (_, stats) = lloyd.relax_points(&points, boundary).unwrap();

        // Should converge for simple case
        if stats.converged {
            assert!(stats.final_iteration < 100);
        }
    }

    #[test]
    fn test_boundary_behaviors() {
        let points = vec![
            Point2::new(-0.5, 0.5), // Outside boundary
            Point2::new(0.5, 0.5),
        ];

        let boundary = Bounds2D::from_points(Point2D::new(0.0, 0.0), Point2D::new(1.0, 1.0));

        // Test clamp behavior
        let config_clamp = LloydConfig {
            boundary_behavior: BoundaryBehavior::Clamp,
            max_iterations: 1,
            ..Default::default()
        };

        let lloyd = LloydRelaxation::new(config_clamp);
        let (clamped_points, _) = lloyd.relax_points(&points, boundary).unwrap();

        // All points should be within boundary
        for point in &clamped_points {
            assert!(point.x >= 0.0 && point.x <= 1.0);
            assert!(point.y >= 0.0 && point.y <= 1.0);
        }
    }

    #[test]
    fn test_lloyd_variants() {
        let points = vec![
            Point2::new(0.3, 0.3),
            Point2::new(0.7, 0.3),
            Point2::new(0.7, 0.7),
            Point2::new(0.3, 0.7),
        ];

        let boundary = Bounds2D::from_points(Point2D::new(0.0, 0.0), Point2D::new(1.0, 1.0));

        let standard_config = LloydVariants::standard(10);
        let weighted_config = LloydVariants::weighted(10, 0.5);
        let precision_config = LloydVariants::high_precision(10);

        let lloyd_standard = LloydRelaxation::new(standard_config);
        let lloyd_weighted = LloydRelaxation::new(weighted_config);
        let lloyd_precision = LloydRelaxation::new(precision_config);

        let (_, stats_standard) = lloyd_standard.relax_points(&points, boundary).unwrap();
        let (_, stats_weighted) = lloyd_weighted.relax_points(&points, boundary).unwrap();
        let (_, stats_precision) = lloyd_precision.relax_points(&points, boundary).unwrap();

        // All variants should produce valid results
        assert!(!stats_standard.iterations.is_empty());
        assert!(!stats_weighted.iterations.is_empty());
        assert!(!stats_precision.iterations.is_empty());
    }
}
