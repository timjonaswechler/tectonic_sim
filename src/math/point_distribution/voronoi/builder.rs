// src/math/point_distribution/voronoi/builder.rs

use crate::debug::visualization::svg::create_merged_polygon_svg;
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
use bevy::log::{debug, error, info, warn};
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

    pub fn generate_polygons(
        &self,
        seed_resource: &mut SeedResource,
    ) -> MathResult<Vec<Vec<Vec2>>> {
        let step_description =
            |step_name: &str| format!("VoronoiBuilder::generate_polygons - Step: {}", step_name);

        let bounds = self.config.calculation_bounds.unwrap_or_else(|| {
            Bounds2D::from_points(Vec2::new(-100.0, -100.0), Vec2::new(100.0, 100.0))
        });
        debug!(
            "{:?}: Target Bounds: {:?}",
            step_description("BoundsDetermination"),
            bounds
        );

        let padded_bounds = bounds.expand(
            (bounds.size().max_element() * self.config.boundary_padding_factor as f32).max(1.0),
        );
        debug!(
            "{:?}: Padded Bounds for point generation/relaxation: {:?}",
            step_description("BoundsPadding"),
            padded_bounds
        );

        // 3. Generiere initiale Punkte (Vec2)
        let step_name = "InitialPointGeneration";
        let initial_points_bevy = self
            .generate_initial_points(&padded_bounds, seed_resource)
            .map_err(|e| MathError::GeometricFailure {
                // Beispiel: Mappe generischen Fehler auf spezifischeren Kontext
                operation: format!(
                    "{} - Failed in generate_initial_points. Bounds: {:?}, Error: {:?}",
                    step_description(step_name),
                    padded_bounds,
                    e
                ),
            })?;

        debug!(
            "{:?}: Initial points generated: {}",
            step_description(step_name),
            initial_points_bevy.len()
        );
        if initial_points_bevy.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: initial_points_bevy.len(),
            });
        }

        // 4. Wende Lloyd-Relaxation an (arbeitet mit Vec2)
        let step_name = "LloydRelaxation";
        let initial_points_bevy_for_debug = initial_points_bevy.clone();
        let relaxed_points_bevy = if self.config.lloyd_iterations > 0 {
            let lloyd_config = LloydConfig {
                max_iterations: self.config.lloyd_iterations,
                boundary_handling: BoundaryHandling::Clamp,
                convergence_tolerance_sq: 1e-8,
            };
            let lloyd_relaxer = LloydRelaxation::new(lloyd_config);
            debug!(
                "{:?}: Relaxing {} points within padded_bounds: {:?}",
                step_description(step_name),
                initial_points_bevy.len(),
                padded_bounds
            );

            let result = lloyd_relaxer
                .relax_points(&initial_points_bevy, &padded_bounds)
                .map_err(|e| MathError::GeometricFailure {
                    // Gib mehr Kontext!
                    operation: format!(
                        "{} - LloydRelaxation::relax_points failed. Error: {:?}",
                        step_description(step_name),
                        e
                    ),
                })?;
            let (relaxed, stats) = result;

            debug!(
                "{:?}: Lloyd relaxation stats: {:?}",
                step_description(step_name),
                stats
            );
            debug!(
                "{:?}: Points after relaxation: {}",
                step_description(step_name),
                relaxed.len()
            );
            relaxed
        } else {
            initial_points_bevy
        };

        if relaxed_points_bevy.len() < 3 {
            error!(
                "{:?}: Not enough points after relaxation ({} < 3). Config was: {:?}",
                step_description(step_name),
                relaxed_points_bevy.len(),
                self.config // Gibt die gesamte VoronoiConfig aus
            );
            if relaxed_points_bevy.len() < 10 {
                debug!("Relaxed points causing error: {:?}", relaxed_points_bevy);
            }
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: relaxed_points_bevy.len(),
            });
        }

        // 5. Konvertiere zu SpadePoints und erstelle Delaunay-Triangulation
        let step_name = "SpadeConversionAndTriangulation";
        let spade_points = bevy_to_spade_points(&relaxed_points_bevy);
        if spade_points.len() < 3 {
            error!(
                "{:?}: Not enough spade_points after conversion ({} < 3), from {} bevy_points. Relaxed points: {:?}",
                step_description(step_name),
                spade_points.len(),
                relaxed_points_bevy.len(),
                if relaxed_points_bevy.len() < 10 {
                    Some(&relaxed_points_bevy)
                } else {
                    None
                }
            );
            return Err(MathError::TriangulationFailed {
                // Besser, da es ein Triangulations-Setup-Fehler ist
                reason: format!(
                    "{}: Not enough unique points after conversion for Delaunay. Got {}, expected 3. Input bevy points: {}. Relaxed points (sample): {:?}",
                    step_description(step_name),
                    spade_points.len(),
                    relaxed_points_bevy.len(),
                    relaxed_points_bevy.iter().take(5).collect::<Vec<_>>()
                ),
            });
        }

        let delaunay_triangulation = DelaunayTriangulation::bulk_load_stable(spade_points.clone()) // Klonen, wenn spade_points noch für Debug gebraucht wird
        .map_err(|e| MathError::TriangulationFailed {
            reason: format!(
                "{}: Spade's bulk_load_stable failed: {:?}. Input point count: {}. Sample points: {:?}",
                step_description(step_name), e, spade_points.len(), spade_points.iter().take(5).collect::<Vec<_>>()
            ),
        })?;

        debug!(
            "{:?}: Delaunay triangulation created with {} vertices and {} faces.",
            step_description(step_name),
            delaunay_triangulation.num_vertices(),
            delaunay_triangulation.num_all_faces()
        );

        if delaunay_triangulation.num_vertices() < 3 {
            error!(
                "{:?}: Delaunay triangulation resulted in too few vertices ({} < 3). Initial spade points: {}. Sample points: {:?}",
                step_description(step_name),
                delaunay_triangulation.num_vertices(),
                spade_points.len(),
                spade_points.iter().take(5).collect::<Vec<_>>()
            );
            return Err(MathError::TriangulationFailed {
                reason: format!(
                    "{}: Delaunay triangulation resulted in {} vertices, expected at least 3. Input spade points: {}. Sample: {:?}",
                    step_description(step_name),
                    delaunay_triangulation.num_vertices(),
                    spade_points.len(),
                    spade_points.iter().take(5).collect::<Vec<_>>()
                ),
            });
        }
        //''''''''''################
        // Temporär für SVG-Debug: Roh-Voronoi-Zellen extrahieren *ohne* direktes Clipping
        // (Diese Implementierung für Roh-Zellen ist sehr rudimentär und zeigt nur Umkreismittelpunkte)
        let mut raw_cell_vertex_lists_for_svg: Vec<Vec<Vec2>> = Vec::new();
        for vertex_handle in delaunay_triangulation.vertices() {
            let generator_spade = vertex_handle.position();
            let mut cell_circumcenters_spade_raw: Vec<spade::Point2<f32>> = Vec::new(); // Verwende f32 für Spade
            for connected_edge in vertex_handle.out_edges() {
                if let Some(inner_face) = connected_edge.face().as_inner() {
                    cell_circumcenters_spade_raw.push(inner_face.circumcenter());
                }
            }
            if cell_circumcenters_spade_raw.len() >= 2 {
                cell_circumcenters_spade_raw.sort_unstable_by(|a, b| {
                    let angle_a = (a.y - generator_spade.y).atan2(a.x - generator_spade.x);
                    let angle_b = (b.y - generator_spade.y).atan2(b.x - generator_spade.x);
                    angle_a
                        .partial_cmp(&angle_b)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                raw_cell_vertex_lists_for_svg.push(
                    cell_circumcenters_spade_raw
                        .into_iter()
                        .map(|sp| Vec2::new(sp.x, sp.y))
                        .collect(),
                );
            } else if !cell_circumcenters_spade_raw.is_empty() {
                // Wenn nur ein Punkt, trotzdem als Linie mit sich selbst darstellen?
                raw_cell_vertex_lists_for_svg.push(
                    cell_circumcenters_spade_raw
                        .into_iter()
                        .map(|sp| Vec2::new(sp.x, sp.y))
                        .collect(),
                );
            }
        }
        // End Temporär SVG

        //''''''''''################
        // 6. Extrahiere Voronoi-Zellen
        let step_name = "VoronoiCellExtraction";
        let voronoi_cells = VoronoiExtractor::extract_cells(&delaunay_triangulation, Some(&bounds))
        .map_err(|e| MathError::GeometricFailure {
            operation: format!(
                "{}: Failed during VoronoiExtractor::extract_cells with clip_bounds {:?}. Triangulation had {} vertices. Error: {:?}",
                step_description(step_name), bounds, delaunay_triangulation.num_vertices(), e
            )
        })?;
        warn!(
            "{:?}: Extracted {} Voronoi cells.",
            step_description(step_name),
            voronoi_cells.len()
        );
        let voronoi_cells_for_merging = &voronoi_cells;
        if voronoi_cells.is_empty() {
            warn!(
                "{:?}: No Voronoi cells were extracted. Triangulation had {} vertices. Clip bounds: {:?}",
                step_description(step_name),
                delaunay_triangulation.num_vertices(),
                bounds
            );
            return Err(MathError::GeometricFailure {
                operation: format!(
                    "{}: No Voronoi cells extracted. Triangulation vertices: {}, Clip bounds: {:?}",
                    step_description(step_name),
                    delaunay_triangulation.num_vertices(),
                    bounds
                ),
            });
        }

        // 7. Wähle und merge Zellen
        let step_name = "CellMerging";
        let merge_strategy = MergeStrategy::RandomConnectedCells {
            min_vertices: 6,
            max_cells: 10,
        };
        debug!(
            "{:?}: Using merge strategy: {:?}",
            step_description(step_name),
            merge_strategy
        );

        let merger = PolygonMerger::new(merge_strategy);
        let merged_polygon_point_lists = merger.merge_selected_cells(&voronoi_cells_for_merging, seed_resource)
        .map_err(|e| MathError::GeometricFailure {
            operation: format!(
                "{}: PolygonMerger::merge_selected_cells failed. Strategy: {:?}, Cell count: {}. Error: {:?}",
                step_description(step_name), merge_strategy, voronoi_cells.len(), e
            )
        })?;
        debug!(
            "{:?}: Merged into {} polygon(s).",
            step_description(step_name),
            merged_polygon_point_lists.len()
        );

        if merged_polygon_point_lists.is_empty() {
            warn!(
                "{:?}: Polygon merging resulted in no polygons. Voronoi cells available: {}, Strategy: {:?}",
                step_description(step_name),
                voronoi_cells.len(),
                merge_strategy
            );
            return Err(MathError::GeometricFailure {
                operation: format!(
                    "{}: Polygon merging resulted in no polygons. Voronoi cells: {}, Strategy: {:?}",
                    step_description(step_name),
                    voronoi_cells.len(),
                    merge_strategy
                ),
            });
        }
        let initial_merged_count = merged_polygon_point_lists.len();
        // 8. Wende Glättung an
        let step_name = "Smoothing";
        let mut final_polygons = Vec::with_capacity(merged_polygon_point_lists.len());
        if self.config.smoothing_iterations > 0 {
            let smoother = ChaikinSmoother::new(self.config.smoothing_iterations); // Beispiel-Smoother
            debug!(
                "{:?}: Applying Chaikin smoothing with {} iterations.",
                step_description(step_name),
                self.config.smoothing_iterations
            );

            for (idx, point_list) in merged_polygon_point_lists.into_iter().enumerate() {
                if point_list.len() >= 3 {
                    let smoothed_list = smoother.smooth_points(&point_list, true).map_err(|e| {
                        MathError::GeometricFailure {
                            operation: format!(
                                "{}: Smoothing polygon #{} ({} points) failed. Error: {:?}",
                                step_description(step_name),
                                idx,
                                point_list.len(),
                                e
                            ),
                        }
                    })?;

                    if smoothed_list.len() >= 3 {
                        final_polygons.push(smoothed_list);
                    } else {
                        warn!(
                            "{}: Polygon #{} (originally {} points) degenerated to {} points after smoothing. Skipping.",
                            step_description(step_name),
                            idx,
                            point_list.len(),
                            smoothed_list.len()
                        );
                    }
                } else if !point_list.is_empty() {
                    warn!(
                        "{}: Polygon #{} ({} points) was too short to smooth, kept as is (if >=3 points).",
                        step_description(step_name),
                        idx,
                        point_list.len()
                    );
                    final_polygons.push(point_list); // Behalte Linien/Punkte, wenn sie entstehen und >=3 waren. Hier ist es aber schon zu spät für >=3.
                    // Besser wäre es, wenn <3 hier nicht mehr hinzugefügt wird.
                }
            }
        } else {
            debug!(
                "{:?}: Skipping smoothing (iterations set to 0).",
                step_description(step_name)
            );
            final_polygons = merged_polygon_point_lists
                .into_iter()
                .filter(|pl| pl.len() >= 3)
                .collect();
        }

        debug!(
            "{:?}: Produced {} final polygons.",
            step_description(step_name),
            final_polygons.len()
        );
        if final_polygons.is_empty() {
            warn!(
                "{:?}: No valid polygons remained after smoothing. Initial merged count: {}",
                step_description(step_name),
                initial_merged_count
            );
            return Err(MathError::GeometricFailure {
                operation: format!(
                    "{}: No valid polygons remained after smoothing. Input to smoothing: {} polygon(s).",
                    step_description(step_name),
                    if self.config.smoothing_iterations > 0 {
                        // merged_polygon_point_lists wurde konsumiert,
                        // wir bräuchten die Länge vorher oder müssen anders darauf zugreifen
                        "unknown (was > 0)".to_string()
                    } else {
                        initial_merged_count.to_string()
                    }
                ),
            });
        }
        #[cfg(debug_assertions)]
        {
            create_merged_polygon_svg(
                "output/polygons_scaled_to_fit.svg",
                &bounds,
                &final_polygons,
                1024.0,
                true,
            )
            .expect("Konnte SVG nicht erstellen");
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
