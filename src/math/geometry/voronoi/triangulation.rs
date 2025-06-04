// src/math/geometry/voronoi/triangulation.rs

use crate::math::{error::*, types::*, utils::*};
use spade::{DelaunayTriangulation, Point2, Triangulation};
use std::collections::{HashMap, HashSet};

/// Voronoi-Diagramm mit Delaunay-Triangulation
pub struct VoronoiTriangulation {
    triangulation: DelaunayTriangulation<Point2<f64>>,
    boundary: Option<Bounds2D>,
}

impl VoronoiTriangulation {
    /// Erstellt eine neue VoronoiTriangulation
    pub fn new() -> Self {
        Self {
            triangulation: DelaunayTriangulation::new(),
            boundary: None,
        }
    }

    /// Setzt die Boundary für das Voronoi-Diagramm
    pub fn with_boundary(mut self, boundary: Bounds2D) -> Self {
        self.boundary = Some(boundary);
        self
    }

    /// Fügt Punkte zur Triangulation hinzu
    pub fn add_points(&mut self, points: &[Point2<f64>]) -> MathResult<()> {
        for &point in points {
            if let Some(boundary) = &self.boundary {
                if self.is_point_in_boundary(point, *boundary) {
                    self.triangulation
                        .insert(point)
                        .map_err(|_| MathError::GeometricFailure {
                            operation: "Failed to insert point into triangulation".to_string(),
                        })?;
                }
            } else {
                self.triangulation
                    .insert(point)
                    .map_err(|_| MathError::GeometricFailure {
                        operation: "Failed to insert point into triangulation".to_string(),
                    })?;
            }
        }
        Ok(())
    }

    /// Fügt einen einzelnen Punkt hinzu
    pub fn add_point(&mut self, point: Point2<f64>) -> MathResult<()> {
        self.add_points(&[point])
    }

    /// Gibt alle Delaunay-Dreiecke zurück
    pub fn delaunay_triangles(&self) -> Vec<DelaunayTriangle> {
        let mut triangles = Vec::new();

        for face in self.triangulation.inner_faces() {
            if let Some(triangle) = self.face_to_triangle(face) {
                triangles.push(triangle);
            }
        }

        triangles
    }

    /// Gibt alle Voronoi-Zellen zurück
    pub fn voronoi_cells(&self) -> MathResult<Vec<VoronoiCell>> {
        let mut cells = Vec::new();

        for vertex_handle in self.triangulation.vertices() {
            let generator = vertex_handle.position();
            let cell_vertices = self.compute_voronoi_cell_vertices(vertex_handle)?;

            if cell_vertices.len() >= 3 {
                cells.push(VoronoiCell {
                    generator,
                    vertices: cell_vertices,
                    is_boundary_cell: self.is_boundary_cell(vertex_handle),
                });
            }
        }

        Ok(cells)
    }

    /// Berechnet die Voronoi-Kanten
    pub fn voronoi_edges(&self) -> MathResult<Vec<VoronoiEdge>> {
        let mut edges = Vec::new();
        let mut processed_edges = HashSet::new();

        for edge_handle in self.triangulation.undirected_edges() {
            let edge_id = self.edge_to_id(edge_handle);
            if processed_edges.contains(&edge_id) {
                continue;
            }
            processed_edges.insert(edge_id);

            // Hole die zwei angrenzenden Faces
            let face1 = edge_handle.face();
            let face2 = edge_handle.rev().face();

            if face1.is_outer() || face2.is_outer() {
                // Randkante - erstelle unendliche Voronoi-Kante
                if let Some(edge) = self.create_boundary_voronoi_edge(edge_handle) {
                    edges.push(edge);
                }
            } else {
                // Innere Kante - erstelle finite Voronoi-Kante
                if let Some(edge) = self.create_finite_voronoi_edge(edge_handle) {
                    edges.push(edge);
                }
            }
        }

        Ok(edges)
    }

    /// Findet die nächsten Nachbarn eines Punktes
    pub fn nearest_neighbors(&self, point: Point2<f64>) -> Vec<Point2<f64>> {
        let mut neighbors = Vec::new();

        if let Ok(vertex_handle) = self.triangulation.locate_vertex(point) {
            for edge in vertex_handle.out_edges() {
                neighbors.push(edge.to().position());
            }
        }

        neighbors
    }

    /// Findet das Delaunay-Dreieck das einen Punkt enthält
    pub fn locate_triangle(&self, point: Point2<f64>) -> Option<DelaunayTriangle> {
        match self.triangulation.locate(point) {
            Ok(spade::PositionInTriangulation::OnFace(face)) => self.face_to_triangle(face),
            _ => None,
        }
    }

    /// Berechnet Voronoi-Diagramm-Statistiken
    pub fn statistics(&self) -> MathResult<VoronoiStatistics> {
        let cells = self.voronoi_cells()?;
        let triangles = self.delaunay_triangles();
        let edges = self.voronoi_edges()?;

        let mut total_cell_area = 0.0;
        let mut cell_areas = Vec::new();
        let mut cell_perimeters = Vec::new();

        for cell in &cells {
            let area = self.calculate_polygon_area(&cell.vertices);
            let perimeter = self.calculate_polygon_perimeter(&cell.vertices);

            total_cell_area += area;
            cell_areas.push(area);
            cell_perimeters.push(perimeter);
        }

        let mean_area = if !cell_areas.is_empty() {
            total_cell_area / cell_areas.len() as f64
        } else {
            0.0
        };

        let area_variance = if cell_areas.len() > 1 {
            cell_areas
                .iter()
                .map(|&area| (area - mean_area).powi(2))
                .sum::<f64>()
                / (cell_areas.len() - 1) as f64
        } else {
            0.0
        };

        Ok(VoronoiStatistics {
            point_count: self.triangulation.num_vertices(),
            triangle_count: triangles.len(),
            voronoi_cell_count: cells.len(),
            voronoi_edge_count: edges.len(),
            total_area: total_cell_area,
            mean_cell_area: mean_area,
            cell_area_variance: area_variance,
            boundary_cell_count: cells.iter().filter(|c| c.is_boundary_cell).count(),
        })
    }

    // === Private Helper Methods ===

    fn is_point_in_boundary(&self, point: Point2<f64>, boundary: Bounds2D) -> bool {
        point.x >= boundary.min.x as f64
            && point.x <= boundary.max.x as f64
            && point.y >= boundary.min.y as f64
            && point.y <= boundary.max.y as f64
    }

    fn face_to_triangle(
        &self,
        face: spade::handles::FaceHandle<Point2<f64>>,
    ) -> Option<DelaunayTriangle> {
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

        if vertices.len() == 3 {
            Some(DelaunayTriangle {
                vertices: [vertices[0], vertices[1], vertices[2]],
                circumcenter: self.calculate_circumcenter(vertices[0], vertices[1], vertices[2]),
            })
        } else {
            None
        }
    }

    fn compute_voronoi_cell_vertices(
        &self,
        vertex_handle: spade::handles::VertexHandle<Point2<f64>>,
    ) -> MathResult<Vec<Point2<f64>>> {
        let mut circumcenters = Vec::new();

        // Sammle alle Umkreismittelpunkte der angrenzenden Dreiecke
        for edge in vertex_handle.out_edges() {
            let face = edge.face();
            if !face.is_outer() {
                if let Some(triangle) = self.face_to_triangle(face) {
                    if let Some(circumcenter) = triangle.circumcenter {
                        circumcenters.push(circumcenter);
                    }
                }
            }
        }

        if circumcenters.is_empty() {
            return Ok(vec![]);
        }

        // Sortiere Circumcenter um den Generator-Punkt
        let generator = vertex_handle.position();
        circumcenters.sort_by(|a, b| {
            let angle_a = (a.y - generator.y).atan2(a.x - generator.x);
            let angle_b = (b.y - generator.y).atan2(b.x - generator.x);
            angle_a
                .partial_cmp(&angle_b)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Clippe gegen Boundary falls vorhanden
        if let Some(boundary) = self.boundary {
            self.clip_voronoi_cell(&circumcenters, boundary)
        } else {
            Ok(circumcenters)
        }
    }

    fn clip_voronoi_cell(
        &self,
        vertices: &[Point2<f64>],
        boundary: Bounds2D,
    ) -> MathResult<Vec<Point2<f64>>> {
        // Vereinfachtes Sutherland-Hodgman Clipping
        let mut clipped = vertices.to_vec();

        let clip_edges = [
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

        for (edge_start, edge_end) in &clip_edges {
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

    fn is_boundary_cell(&self, vertex_handle: spade::handles::VertexHandle<Point2<f64>>) -> bool {
        // Eine Zelle ist eine Boundary-Zelle wenn sie mindestens eine unendliche Kante hat
        for edge in vertex_handle.out_edges() {
            if edge.face().is_outer() {
                return true;
            }
        }
        false
    }

    fn edge_to_id(&self, edge: spade::handles::UndirectedEdgeHandle<Point2<f64>>) -> (u64, u64) {
        let from_pos = edge.as_directed_edge().from().position();
        let to_pos = edge.as_directed_edge().to().position();

        let hash_from = self.point_to_hash(from_pos);
        let hash_to = self.point_to_hash(to_pos);

        if hash_from < hash_to {
            (hash_from, hash_to)
        } else {
            (hash_to, hash_from)
        }
    }

    fn point_to_hash(&self, point: Point2<f64>) -> u64 {
        let x_bits = (point.x * 1e9) as i64;
        let y_bits = (point.y * 1e9) as i64;
        ((x_bits as u64) << 32) | (y_bits as u64)
    }

    fn create_finite_voronoi_edge(
        &self,
        edge_handle: spade::handles::UndirectedEdgeHandle<Point2<f64>>,
    ) -> Option<VoronoiEdge> {
        let directed_edge = edge_handle.as_directed_edge();
        let face1 = directed_edge.face();
        let face2 = directed_edge.rev().face();

        if face1.is_outer() || face2.is_outer() {
            return None;
        }

        let triangle1 = self.face_to_triangle(face1)?;
        let triangle2 = self.face_to_triangle(face2)?;

        if let (Some(cc1), Some(cc2)) = (triangle1.circumcenter, triangle2.circumcenter) {
            Some(VoronoiEdge {
                start: cc1,
                end: Some(cc2),
                generator1: directed_edge.from().position(),
                generator2: directed_edge.to().position(),
                is_infinite: false,
            })
        } else {
            None
        }
    }

    fn create_boundary_voronoi_edge(
        &self,
        edge_handle: spade::handles::UndirectedEdgeHandle<Point2<f64>>,
    ) -> Option<VoronoiEdge> {
        let directed_edge = edge_handle.as_directed_edge();
        let inner_face = if directed_edge.face().is_outer() {
            directed_edge.rev().face()
        } else {
            directed_edge.face()
        };

        if inner_face.is_outer() {
            return None;
        }

        let triangle = self.face_to_triangle(inner_face)?;
        let circumcenter = triangle.circumcenter?;

        Some(VoronoiEdge {
            start: circumcenter,
            end: None, // Unendliche Kante
            generator1: directed_edge.from().position(),
            generator2: directed_edge.to().position(),
            is_infinite: true,
        })
    }

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
            return None;
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

    fn calculate_polygon_area(&self, vertices: &[Point2<f64>]) -> f64 {
        if vertices.len() < 3 {
            return 0.0;
        }

        let mut area = 0.0;
        let n = vertices.len();

        for i in 0..n {
            let j = (i + 1) % n;
            area += vertices[i].x * vertices[j].y;
            area -= vertices[j].x * vertices[i].y;
        }

        (area * 0.5).abs()
    }

    fn calculate_polygon_perimeter(&self, vertices: &[Point2<f64>]) -> f64 {
        if vertices.len() < 2 {
            return 0.0;
        }

        let mut perimeter = 0.0;
        let n = vertices.len();

        for i in 0..n {
            let j = (i + 1) % n;
            let dx = vertices[j].x - vertices[i].x;
            let dy = vertices[j].y - vertices[i].y;
            perimeter += (dx * dx + dy * dy).sqrt();
        }

        perimeter
    }
}

/// Delaunay-Dreieck
#[derive(Debug, Clone)]
pub struct DelaunayTriangle {
    pub vertices: [Point2<f64>; 3],
    pub circumcenter: Option<Point2<f64>>,
}

/// Voronoi-Zelle
#[derive(Debug, Clone)]
pub struct VoronoiCell {
    pub generator: Point2<f64>,
    pub vertices: Vec<Point2<f64>>,
    pub is_boundary_cell: bool,
}

impl VoronoiCell {
    /// Berechnet die Fläche der Zelle
    pub fn area(&self) -> f64 {
        if self.vertices.len() < 3 {
            return 0.0;
        }

        let mut area = 0.0;
        let n = self.vertices.len();

        for i in 0..n {
            let j = (i + 1) % n;
            area += self.vertices[i].x * self.vertices[j].y;
            area -= self.vertices[j].x * self.vertices[i].y;
        }

        (area * 0.5).abs()
    }

    /// Berechnet den Umfang der Zelle
    pub fn perimeter(&self) -> f64 {
        if self.vertices.len() < 2 {
            return 0.0;
        }

        let mut perimeter = 0.0;
        let n = self.vertices.len();

        for i in 0..n {
            let j = (i + 1) % n;
            let dx = self.vertices[j].x - self.vertices[i].x;
            let dy = self.vertices[j].y - self.vertices[i].y;
            perimeter += (dx * dx + dy * dy).sqrt();
        }

        perimeter
    }

    /// Prüft ob ein Punkt in der Zelle liegt
    pub fn contains_point(&self, point: Point2<f64>) -> bool {
        if self.vertices.len() < 3 {
            return false;
        }

        // Ray-casting algorithm
        let mut inside = false;
        let n = self.vertices.len();

        let mut j = n - 1;
        for i in 0..n {
            let vi = self.vertices[i];
            let vj = self.vertices[j];

            if ((vi.y > point.y) != (vj.y > point.y))
                && (point.x < (vj.x - vi.x) * (point.y - vi.y) / (vj.y - vi.y) + vi.x)
            {
                inside = !inside;
            }
            j = i;
        }

        inside
    }
}

/// Voronoi-Kante
#[derive(Debug, Clone)]
pub struct VoronoiEdge {
    pub start: Point2<f64>,
    pub end: Option<Point2<f64>>, // None für unendliche Kanten
    pub generator1: Point2<f64>,
    pub generator2: Point2<f64>,
    pub is_infinite: bool,
}

/// Voronoi-Diagramm Statistiken
#[derive(Debug, Clone)]
pub struct VoronoiStatistics {
    pub point_count: usize,
    pub triangle_count: usize,
    pub voronoi_cell_count: usize,
    pub voronoi_edge_count: usize,
    pub total_area: f64,
    pub mean_cell_area: f64,
    pub cell_area_variance: f64,
    pub boundary_cell_count: usize,
}

impl VoronoiStatistics {
    /// Berechnet die Uniformitäts-Metrik
    pub fn uniformity_coefficient(&self) -> f64 {
        if self.mean_cell_area <= 0.0 {
            return 0.0;
        }

        1.0 - (self.cell_area_variance.sqrt() / self.mean_cell_area)
    }

    /// Berechnet die Dichte (Punkte pro Flächeneinheit)
    pub fn point_density(&self) -> f64 {
        if self.total_area <= 0.0 {
            return 0.0;
        }

        self.point_count as f64 / self.total_area
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_voronoi_triangulation_basic() {
        let mut voronoi = VoronoiTriangulation::new();

        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.5, 1.0),
        ];

        voronoi.add_points(&points).unwrap();

        let triangles = voronoi.delaunay_triangles();
        assert_eq!(triangles.len(), 1);

        let cells = voronoi.voronoi_cells().unwrap();
        assert_eq!(cells.len(), 3);
    }

    #[test]
    fn test_voronoi_with_boundary() {
        let boundary = Bounds2D::from_points(Point2D::new(0.0, 0.0), Point2D::new(10.0, 10.0));
        let mut voronoi = VoronoiTriangulation::new().with_boundary(boundary);

        let points = vec![
            Point2::new(2.0, 2.0),
            Point2::new(8.0, 2.0),
            Point2::new(8.0, 8.0),
            Point2::new(2.0, 8.0),
            Point2::new(5.0, 5.0),
        ];

        voronoi.add_points(&points).unwrap();

        let cells = voronoi.voronoi_cells().unwrap();
        assert_eq!(cells.len(), 5);

        // Alle Zellen sollten positive Fläche haben
        for cell in &cells {
            assert!(cell.area() > 0.0);
        }
    }

    #[test]
    fn test_voronoi_statistics() {
        let mut voronoi = VoronoiTriangulation::new();

        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];

        voronoi.add_points(&points).unwrap();

        let stats = voronoi.statistics().unwrap();

        assert_eq!(stats.point_count, 4);
        assert!(stats.triangle_count >= 1);
        assert_eq!(stats.voronoi_cell_count, 4);
        assert!(stats.total_area > 0.0);
    }

    #[test]
    fn test_voronoi_cell_operations() {
        let cell = VoronoiCell {
            generator: Point2::new(0.5, 0.5),
            vertices: vec![
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 0.0),
                Point2::new(1.0, 1.0),
                Point2::new(0.0, 1.0),
            ],
            is_boundary_cell: false,
        };

        assert_eq!(cell.area(), 1.0);
        assert_eq!(cell.perimeter(), 4.0);
        assert!(cell.contains_point(Point2::new(0.5, 0.5)));
        assert!(!cell.contains_point(Point2::new(1.5, 0.5)));
    }

    #[test]
    fn test_nearest_neighbors() {
        let mut voronoi = VoronoiTriangulation::new();

        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        ];

        voronoi.add_points(&points).unwrap();

        let neighbors = voronoi.nearest_neighbors(Point2::new(0.0, 0.0));
        assert!(neighbors.len() >= 1);
    }
}
