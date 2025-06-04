// src/math/geometry/voronoi/merger.rs

use super::config::VoronoiConfig;
use super::triangulation::{VoronoiCell, VoronoiTriangulation};
use crate::math::{error::*, types::*, utils::*};
use spade::{DelaunayTriangulation, Point2, Triangulation};
use std::collections::{HashMap, HashSet};

/// Polygon-Merger für Voronoi-Zellen
pub struct PolygonMerger {
    config: VoronoiConfig,
    merge_strategy: MergeStrategy,
}

/// Verschiedene Strategien zum Zusammenfügen von Zellen
#[derive(Debug, Clone, Copy)]
pub enum MergeStrategy {
    /// Füge basierend auf Flächen-Kriterien zusammen
    AreaBased { min_area: f64, max_area: f64 },
    /// Füge basierend auf Nachbarschaft zusammen
    NeighborhoodBased { max_neighbors: usize },
    /// Füge die N größten Zellen zusammen
    LargestCells { count: usize },
    /// Füge Zellen mit ähnlichen Eigenschaften zusammen
    PropertyBased { similarity_threshold: f64 },
    /// Zufällige Zusammenfügung
    Random { target_count: usize },
}

impl Default for MergeStrategy {
    fn default() -> Self {
        Self::LargestCells { count: 5 }
    }
}

impl PolygonMerger {
    /// Erstellt einen neuen PolygonMerger
    pub fn new(config: &VoronoiConfig) -> Self {
        Self {
            config: config.clone(),
            merge_strategy: MergeStrategy::default(),
        }
    }

    /// Setzt die Merge-Strategie
    pub fn with_strategy(mut self, strategy: MergeStrategy) -> Self {
        self.merge_strategy = strategy;
        self
    }

    /// Führt Voronoi-Zellen zu Polygonen zusammen
    pub fn merge_cells(
        &self,
        triangulation: &DelaunayTriangulation<Point2<f64>>,
    ) -> MathResult<Vec<Point2<f64>>> {
        // Erstelle VoronoiTriangulation wrapper
        let mut voronoi = VoronoiTriangulation::new();

        // Extrahiere alle Punkte aus der Triangulation
        let points: Vec<Point2<f64>> = triangulation.vertices().map(|v| v.position()).collect();

        voronoi.add_points(&points)?;
        let cells = voronoi.voronoi_cells()?;

        if cells.is_empty() {
            return Err(MathError::InsufficientPoints {
                expected: 1,
                actual: 0,
            });
        }

        // Wähle Zellen basierend auf der Strategie aus
        let selected_cells = self.select_cells(&cells)?;

        // Führe die ausgewählten Zellen zusammen
        let merged_polygon = self.merge_selected_cells(&selected_cells)?;

        Ok(merged_polygon)
    }

    /// Wählt Zellen basierend auf der konfigurierten Strategie aus
    fn select_cells(&self, cells: &[VoronoiCell]) -> MathResult<Vec<VoronoiCell>> {
        match self.merge_strategy {
            MergeStrategy::AreaBased { min_area, max_area } => {
                self.select_by_area(cells, min_area, max_area)
            }
            MergeStrategy::NeighborhoodBased { max_neighbors } => {
                self.select_by_neighborhood(cells, max_neighbors)
            }
            MergeStrategy::LargestCells { count } => self.select_largest_cells(cells, count),
            MergeStrategy::PropertyBased {
                similarity_threshold,
            } => self.select_by_properties(cells, similarity_threshold),
            MergeStrategy::Random { target_count } => self.select_randomly(cells, target_count),
        }
    }

    /// Wählt Zellen basierend auf Flächenkriterien aus
    fn select_by_area(
        &self,
        cells: &[VoronoiCell],
        min_area: f64,
        max_area: f64,
    ) -> MathResult<Vec<VoronoiCell>> {
        let selected: Vec<VoronoiCell> = cells
            .iter()
            .filter(|cell| {
                let area = cell.area();
                area >= min_area && area <= max_area
            })
            .cloned()
            .collect();

        if selected.is_empty() {
            return Err(MathError::GeometricFailure {
                operation: "No cells found matching area criteria".to_string(),
            });
        }

        Ok(selected)
    }

    /// Wählt die größten Zellen aus
    fn select_largest_cells(
        &self,
        cells: &[VoronoiCell],
        count: usize,
    ) -> MathResult<Vec<VoronoiCell>> {
        let mut sorted_cells: Vec<VoronoiCell> = cells.to_vec();
        sorted_cells.sort_by(|a, b| {
            b.area()
                .partial_cmp(&a.area())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let target_count = count.min(sorted_cells.len()).max(1);
        Ok(sorted_cells.into_iter().take(target_count).collect())
    }

    /// Wählt Zellen basierend auf Nachbarschaft aus
    fn select_by_neighborhood(
        &self,
        cells: &[VoronoiCell],
        max_neighbors: usize,
    ) -> MathResult<Vec<VoronoiCell>> {
        // Vereinfachte Implementierung - wähle Zellen mit wenigen "Nachbarn" (basierend auf Vertex-Anzahl)
        let selected: Vec<VoronoiCell> = cells
            .iter()
            .filter(|cell| cell.vertices.len() <= max_neighbors)
            .cloned()
            .collect();

        if selected.is_empty() {
            // Fallback: nimm die ersten Zellen
            return Ok(cells.iter().take(5).cloned().collect());
        }

        Ok(selected)
    }

    /// Wählt Zellen basierend auf Eigenschaften aus
    fn select_by_properties(
        &self,
        cells: &[VoronoiCell],
        similarity_threshold: f64,
    ) -> MathResult<Vec<VoronoiCell>> {
        if cells.is_empty() {
            return Ok(vec![]);
        }

        // Berechne Durchschnitts-Eigenschaften
        let avg_area: f64 = cells.iter().map(|c| c.area()).sum::<f64>() / cells.len() as f64;
        let avg_perimeter: f64 =
            cells.iter().map(|c| c.perimeter()).sum::<f64>() / cells.len() as f64;

        let selected: Vec<VoronoiCell> = cells
            .iter()
            .filter(|cell| {
                let area_similarity = 1.0 - ((cell.area() - avg_area).abs() / avg_area.max(1e-6));
                let perimeter_similarity =
                    1.0 - ((cell.perimeter() - avg_perimeter).abs() / avg_perimeter.max(1e-6));
                let combined_similarity = (area_similarity + perimeter_similarity) * 0.5;

                combined_similarity >= similarity_threshold
            })
            .cloned()
            .collect();

        if selected.is_empty() {
            // Fallback: nimm die ersten Zellen
            return Ok(cells.iter().take(3).cloned().collect());
        }

        Ok(selected)
    }

    /// Wählt zufällige Zellen aus
    fn select_randomly(
        &self,
        cells: &[VoronoiCell],
        target_count: usize,
    ) -> MathResult<Vec<VoronoiCell>> {
        use rand::{SeedableRng, seq::SliceRandom};

        let mut rng = match self.config.seed {
            Some(seed) => rand::rngs::StdRng::seed_from_u64(seed),
            None => rand::rngs::StdRng::from_entropy(),
        };

        let count = target_count.min(cells.len()).max(1);
        let selected = cells.choose_multiple(&mut rng, count).cloned().collect();

        Ok(selected)
    }

    /// Führt ausgewählte Zellen zu einem zusammenhängenden Polygon zusammen
    fn merge_selected_cells(&self, cells: &[VoronoiCell]) -> MathResult<Vec<Point2<f64>>> {
        if cells.is_empty() {
            return Err(MathError::InsufficientPoints {
                expected: 1,
                actual: 0,
            });
        }

        if cells.len() == 1 {
            return Ok(cells[0].vertices.clone());
        }

        // Sammle alle Vertices
        let mut all_vertices = Vec::new();
        for cell in cells {
            all_vertices.extend(&cell.vertices);
        }

        if all_vertices.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: all_vertices.len(),
            });
        }

        // Berechne konvexe Hülle als Basis-Polygon
        let convex_hull = self.compute_convex_hull(&all_vertices)?;

        // Optional: Smoothing anwenden
        if self.config.smoothing_iterations > 0 {
            self.apply_smoothing(convex_hull)
        } else {
            Ok(convex_hull)
        }
    }

    /// Berechnet die konvexe Hülle von Punkten (Andrew's Algorithm)
    fn compute_convex_hull(&self, points: &[Point2<f64>]) -> MathResult<Vec<Point2<f64>>> {
        if points.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: points.len(),
            });
        }

        let mut sorted_points = points.to_vec();

        // Entferne Duplikate
        sorted_points.sort_by(|a, b| {
            a.x.partial_cmp(&b.x)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
        });
        sorted_points.dedup_by(|a, b| (a.x - b.x).abs() < 1e-10 && (a.y - b.y).abs() < 1e-10);

        if sorted_points.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: sorted_points.len(),
            });
        }

        // Andrew's Monotone Chain Algorithm
        let mut lower_hull = Vec::new();
        for point in &sorted_points {
            while lower_hull.len() >= 2 {
                let cross = self.cross_product(
                    lower_hull[lower_hull.len() - 2],
                    lower_hull[lower_hull.len() - 1],
                    *point,
                );
                if cross <= 0.0 {
                    lower_hull.pop();
                } else {
                    break;
                }
            }
            lower_hull.push(*point);
        }

        let mut upper_hull = Vec::new();
        for point in sorted_points.iter().rev() {
            while upper_hull.len() >= 2 {
                let cross = self.cross_product(
                    upper_hull[upper_hull.len() - 2],
                    upper_hull[upper_hull.len() - 1],
                    *point,
                );
                if cross <= 0.0 {
                    upper_hull.pop();
                } else {
                    break;
                }
            }
            upper_hull.push(*point);
        }

        // Entferne die letzten Punkte (Duplikate)
        lower_hull.pop();
        upper_hull.pop();

        // Kombiniere beide Hulls
        lower_hull.extend(upper_hull);

        if lower_hull.len() < 3 {
            return Err(MathError::GeometricFailure {
                operation: "Convex hull resulted in degenerate polygon".to_string(),
            });
        }

        Ok(lower_hull)
    }

    /// Berechnet das Kreuzprodukt für drei Punkte
    fn cross_product(&self, o: Point2<f64>, a: Point2<f64>, b: Point2<f64>) -> f64 {
        (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
    }

    /// Wendet Smoothing auf das Polygon an
    fn apply_smoothing(&self, polygon: Vec<Point2<f64>>) -> MathResult<Vec<Point2<f64>>> {
        let mut current = polygon;

        for _ in 0..self.config.smoothing_iterations {
            current = self.smooth_iteration(&current)?;
        }

        Ok(current)
    }

    /// Eine Iteration des Chaikin-Smoothing Algorithmus
    fn smooth_iteration(&self, polygon: &[Point2<f64>]) -> MathResult<Vec<Point2<f64>>> {
        if polygon.len() < 3 {
            return Ok(polygon.to_vec());
        }

        let mut smoothed = Vec::new();
        let n = polygon.len();

        for i in 0..n {
            let current = polygon[i];
            let next = polygon[(i + 1) % n];

            // Chaikin subdivision: ersetze jede Kante durch zwei neue Punkte
            let p1 = Point2::new(
                current.x + 0.25 * (next.x - current.x),
                current.y + 0.25 * (next.y - current.y),
            );
            let p2 = Point2::new(
                current.x + 0.75 * (next.x - current.x),
                current.y + 0.75 * (next.y - current.y),
            );

            smoothed.push(p1);
            smoothed.push(p2);
        }

        Ok(smoothed)
    }
}

/// Erweiterte Merge-Utilities
pub struct MergeUtils;

impl MergeUtils {
    /// Analysiert Voronoi-Zellen und schlägt optimale Merge-Parameter vor
    pub fn suggest_merge_strategy(cells: &[VoronoiCell]) -> MergeStrategy {
        if cells.is_empty() {
            return MergeStrategy::default();
        }

        let areas: Vec<f64> = cells.iter().map(|c| c.area()).collect();
        let mean_area = areas.iter().sum::<f64>() / areas.len() as f64;
        let area_variance = areas
            .iter()
            .map(|&area| (area - mean_area).powi(2))
            .sum::<f64>()
            / areas.len() as f64;

        // Wenn die Varianz gering ist, verwende area-based merging
        if area_variance < mean_area * 0.1 {
            MergeStrategy::AreaBased {
                min_area: mean_area * 0.8,
                max_area: mean_area * 1.2,
            }
        } else {
            // Sonst verwende largest cells
            let target_count = (cells.len() as f64 * 0.1).ceil() as usize + 1;
            MergeStrategy::LargestCells {
                count: target_count.max(3).min(cells.len()),
            }
        }
    }

    /// Berechnet Merge-Qualitäts-Metriken
    pub fn calculate_merge_quality(
        original_cells: &[VoronoiCell],
        merged_polygon: &[Point2<f64>],
    ) -> MergeQuality {
        let original_total_area: f64 = original_cells.iter().map(|c| c.area()).sum();
        let merged_area = Self::polygon_area(merged_polygon);

        let area_preservation = if original_total_area > 0.0 {
            merged_area / original_total_area
        } else {
            0.0
        };

        let complexity_reduction = if !original_cells.is_empty() {
            let original_vertices: usize = original_cells.iter().map(|c| c.vertices.len()).sum();
            let merged_vertices = merged_polygon.len();

            if original_vertices > 0 {
                1.0 - (merged_vertices as f64 / original_vertices as f64)
            } else {
                0.0
            }
        } else {
            0.0
        };

        MergeQuality {
            area_preservation,
            complexity_reduction,
            vertex_count_reduction: complexity_reduction,
            smoothness_score: Self::calculate_smoothness(merged_polygon),
        }
    }

    /// Berechnet die Fläche eines Polygons
    fn polygon_area(vertices: &[Point2<f64>]) -> f64 {
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

    /// Berechnet eine Smoothness-Metrik für ein Polygon
    fn calculate_smoothness(vertices: &[Point2<f64>]) -> f64 {
        if vertices.len() < 3 {
            return 0.0;
        }

        let mut total_curvature = 0.0;
        let n = vertices.len();

        for i in 0..n {
            let prev = vertices[(i + n - 1) % n];
            let curr = vertices[i];
            let next = vertices[(i + 1) % n];

            // Berechne Winkel-Änderung
            let v1 = Point2::new(curr.x - prev.x, curr.y - prev.y);
            let v2 = Point2::new(next.x - curr.x, next.y - curr.y);

            let v1_len = (v1.x * v1.x + v1.y * v1.y).sqrt();
            let v2_len = (v2.x * v2.x + v2.y * v2.y).sqrt();

            if v1_len > 1e-10 && v2_len > 1e-10 {
                let cos_angle = (v1.x * v2.x + v1.y * v2.y) / (v1_len * v2_len);
                let angle = cos_angle.clamp(-1.0, 1.0).acos();
                total_curvature += angle;
            }
        }

        // Normiere auf [0, 1], wobei 1 = sehr smooth
        let max_possible_curvature = n as f64 * std::f64::consts::PI;
        1.0 - (total_curvature / max_possible_curvature).min(1.0)
    }
}

/// Qualitäts-Metriken für Merge-Operationen
#[derive(Debug, Clone)]
pub struct MergeQuality {
    /// Wie gut die ursprüngliche Fläche erhalten wurde (0.0 bis 1.0)
    pub area_preservation: f64,
    /// Wie stark die Komplexität reduziert wurde (0.0 bis 1.0)
    pub complexity_reduction: f64,
    /// Vertex-Anzahl Reduktion
    pub vertex_count_reduction: f64,
    /// Smoothness des resultierenden Polygons (0.0 bis 1.0)
    pub smoothness_score: f64,
}

impl MergeQuality {
    /// Berechnet einen kombinierten Qualitäts-Score
    pub fn overall_score(&self) -> f64 {
        (self.area_preservation * 0.4
            + self.complexity_reduction * 0.3
            + self.smoothness_score * 0.3)
            .clamp(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::geometry::voronoi::config::VoronoiConfig;

    #[test]
    fn test_merge_strategy_suggestion() {
        // Erstelle Test-Zellen mit uniformer Größe
        let uniform_cells = vec![
            VoronoiCell {
                generator: Point2::new(0.0, 0.0),
                vertices: vec![
                    Point2::new(0.0, 0.0),
                    Point2::new(1.0, 0.0),
                    Point2::new(1.0, 1.0),
                    Point2::new(0.0, 1.0),
                ],
                is_boundary_cell: false,
            },
            VoronoiCell {
                generator: Point2::new(2.0, 0.0),
                vertices: vec![
                    Point2::new(1.0, 0.0),
                    Point2::new(2.0, 0.0),
                    Point2::new(2.0, 1.0),
                    Point2::new(1.0, 1.0),
                ],
                is_boundary_cell: false,
            },
        ];

        let strategy = MergeUtils::suggest_merge_strategy(&uniform_cells);

        match strategy {
            MergeStrategy::AreaBased { .. } => {
                // Erwartetes Verhalten für uniforme Zellen
            }
            _ => panic!("Expected AreaBased strategy for uniform cells"),
        }
    }

    #[test]
    fn test_convex_hull_computation() {
        let config = VoronoiConfig::default();
        let merger = PolygonMerger::new(&config);

        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.5, 0.5), // Interior point
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];

        let hull = merger.compute_convex_hull(&points).unwrap();

        // Should be 4 corners (excluding interior point)
        assert_eq!(hull.len(), 4);
    }

    #[test]
    fn test_merge_quality_calculation() {
        let original_cells = vec![VoronoiCell {
            generator: Point2::new(0.5, 0.5),
            vertices: vec![
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 0.0),
                Point2::new(1.0, 1.0),
                Point2::new(0.0, 1.0),
            ],
            is_boundary_cell: false,
        }];

        let merged_polygon = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];

        let quality = MergeUtils::calculate_merge_quality(&original_cells, &merged_polygon);

        assert!(quality.area_preservation > 0.9); // Should preserve most area
        assert!(quality.overall_score() > 0.5);
    }
}
