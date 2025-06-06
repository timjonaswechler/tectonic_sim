// src/math/point_distribution/voronoi/builder.rs

use crate::math::probability::SeedResource;
use crate::math::{
    algorithms::{
        point_relaxation::{BoundaryHandling, LloydConfig, LloydRelaxation},
        smoothing::{ChaikinSmoother, Smoothing}, // Beispiel-Smoother
                                                 // convex_hull::{ConvexHullComputer, ConvexHullAlgorithm}, // Falls benötigt für Boundary
    },
    error::{MathError, MathResult},
    point_distribution::voronoi::{
        config::VoronoiConfig,
        merger::{MergeStrategy, PolygonMerger},
        voronoi_diagram::VoronoiExtractor, // Unsere neuen Strukturen
    },
    types::{Bounds2D, bevy_to_spade_points},
};
use crate::next_random_range;
use bevy::math::Vec2;
use spade::{DelaunayTriangulation, Triangulation}; // Spade's Triangulation // Für Rng trait und thread_rng

/// Erzeugt Voronoi-basierte polygonale Formen.
/// Orchestriert die Schritte: initiale Punktgenerierung, Lloyd-Relaxation,
/// Delaunay-Triangulation, Voronoi-Zellextraktion, Zell-Merging und finale Glättung.
pub struct VoronoiBuilder {
    config: VoronoiConfig,
}

impl VoronoiBuilder {
    pub fn new(config: VoronoiConfig) -> MathResult<Self> {
        config.validate()?; // Validiert die übergebene Konfiguration
        Ok(Self { config })
    }

    /// Generiert ein oder mehrere Polygone (als `Vec<Vec2>`) basierend auf der Konfiguration.
    pub fn generate_polygons(
        &self,
        seed_resource: &mut SeedResource,
    ) -> MathResult<Vec<Vec<Vec2>>> {
        // 2. Bestimme Berechnungs-Boundary
        let bounds = self.config.calculation_bounds.unwrap_or_else(|| {
            // Standard-Boundary, wenn keine in der Config gesetzt ist
            // Diese sollte groß genug sein, um die initialen Punkte und den Padding-Faktor abzudecken.
            // Hier ein Beispiel, das angepasst werden muss:
            Bounds2D::from_points(Vec2::new(-100.0, -100.0), Vec2::new(100.0, 100.0))
            // Besser wäre es, die Boundary basierend auf der erwarteten Ausdehnung der Punkte zu wählen.
        });

        let padded_bounds = bounds.expand(
            (bounds.size().max_element() * self.config.boundary_padding_factor as f32).max(1.0),
        );

        // 3. Generiere initiale Punkte (Vec2)
        let initial_points_bevy = self.generate_initial_points(&padded_bounds, seed_resource)?;
        if initial_points_bevy.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: initial_points_bevy.len(),
            });
        }

        // 4. Wende Lloyd-Relaxation an (arbeitet mit Vec2)
        let relaxed_points_bevy = if self.config.lloyd_iterations > 0 {
            let lloyd_config = LloydConfig {
                max_iterations: self.config.lloyd_iterations,
                boundary_handling: BoundaryHandling::Clamp, // Oder aus VoronoiConfig ableiten
                ..Default::default()
            };
            let lloyd_relaxer = LloydRelaxation::new(lloyd_config);
            let (relaxed, _stats) =
                lloyd_relaxer.relax_points(&initial_points_bevy, &padded_bounds)?;
            relaxed
        } else {
            initial_points_bevy
        };
        if relaxed_points_bevy.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: relaxed_points_bevy.len(),
            });
        }

        // 5. Konvertiere zu SpadePoints und erstelle Delaunay-Triangulation
        let spade_points = bevy_to_spade_points(&relaxed_points_bevy);
        if spade_points.len() < 3 {
            // Sollte durch vorherige Checks nicht passieren
            return Err(MathError::TriangulationFailed {
                reason: "Not enough unique points for Delaunay Triangulation after conversion."
                    .to_string(),
            });
        }

        let delaunay_triangulation = DelaunayTriangulation::bulk_load_stable(spade_points)
            .map_err(|e| MathError::TriangulationFailed {
                reason: format!("Failed to create triangulation: {:?}", e),
            })?;
        if delaunay_triangulation.num_vertices() < 3 {
            return Err(MathError::TriangulationFailed {
                reason: "Delaunay Triangulation resulted in too few vertices.".to_string(),
            });
        }

        // 6. Extrahiere Voronoi-Zellen (Vec<VoronoiCell> mit Vec2 Vertices)
        // Die Zellen werden gegen die *originale* (nicht gepaddete) `bounds` geclippt,
        // um die Randeffekte des Paddings zu entfernen.
        let voronoi_cells =
            VoronoiExtractor::extract_cells(&delaunay_triangulation, Some(&bounds))?;
        if voronoi_cells.is_empty() {
            return Err(MathError::GeometricFailure {
                operation: "No Voronoi cells extracted.".to_string(),
            });
        }

        // 7. Wähle und merge Zellen (erzeugt Vec<Vec<Vec2>>)
        // Die MergeStrategy könnte aus der VoronoiConfig kommen oder hier hartkodiert sein.
        // Beispiel: Nimm die größte Zelle oder vereinige `target_cell_count` Zellen.
        let merge_strategy = MergeStrategy::LargestCells {
            count: self.config.target_cell_count.max(1),
        };
        // let merge_strategy = MergeStrategy::UnionAllSelected; // Um alle zu einer konvexen Hülle zu mergen

        let merger = PolygonMerger::new(merge_strategy);
        let merged_polygon_point_lists =
            merger.merge_selected_cells(&voronoi_cells, seed_resource)?;

        if merged_polygon_point_lists.is_empty() {
            return Err(MathError::GeometricFailure {
                operation: "Polygon merging resulted in no polygons.".to_string(),
            });
        }

        // 8. Wende Glättung auf jedes resultierende Polygon an
        let mut final_polygons = Vec::with_capacity(merged_polygon_point_lists.len());
        if self.config.smoothing_iterations > 0 {
            let smoother = ChaikinSmoother::new(self.config.smoothing_iterations); // Beispiel-Smoother
            for point_list in merged_polygon_point_lists {
                if point_list.len() >= 3 {
                    // Annahme: Die Punktliste ist für ein geschlossenes Polygon.
                    // ChaikinSmoother::smooth_points braucht `is_closed`.
                    // Wir nehmen an, dass der Merger geschlossene Polygone (ohne Duplikat am Ende) liefert.
                    let smoothed_list = smoother.smooth_points(&point_list, true)?; // true für geschlossen
                    if smoothed_list.len() >= 3 {
                        final_polygons.push(smoothed_list);
                    }
                } else if !point_list.is_empty() {
                    final_polygons.push(point_list); // Behalte Linien/Punkte, wenn sie entstehen
                }
            }
        } else {
            final_polygons = merged_polygon_point_lists
                .into_iter()
                .filter(|pl| pl.len() >= 3)
                .collect();
        }

        if final_polygons.is_empty() {
            return Err(MathError::GeometricFailure {
                operation: "No valid polygons remained after smoothing.".to_string(),
            });
        }

        Ok(final_polygons)
    }

    /// Generiert initiale Punkte innerhalb der gegebenen Grenzen.
    fn generate_initial_points(
        &self,
        generation_bounds: &Bounds2D,
        seed_resource: &mut SeedResource,
    ) -> MathResult<Vec<Vec2>> {
        if self.config.initial_point_count == 0 {
            return Ok(Vec::new());
        }
        if !generation_bounds.is_valid() {
            return Err(MathError::InvalidConfiguration {
                message: "Bounds for initial point generation are invalid.".to_string(),
            });
        }

        let mut points = Vec::with_capacity(self.config.initial_point_count);
        for _ in 0..self.config.initial_point_count {
            let point = Vec2::new(
                next_random_range!(
                    seed_resource,
                    generation_bounds.min.x,
                    generation_bounds.max.x
                ),
                next_random_range!(
                    seed_resource,
                    generation_bounds.min.y,
                    generation_bounds.max.y
                ),
            );
            points.push(point);
        }
        Ok(points)
    }
}
