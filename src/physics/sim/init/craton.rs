use crate::debug::visualization::{normal_vector::NormalArrowVisual, polygon::PolygonVisual};
use crate::math::algorithms::point_relaxation::{
    BoundaryHandling as LloydBoundaryHandling, LloydConfig, LloydRelaxation,
};
use crate::math::geometry::metalballs::influence::{
    FieldInfluence,
    PolygonFalloff, // Influences
    PolygonInfluence,
};
use crate::math::geometry::metalballs::{
    Contour, // MarchingSquares and Contour
    GridConfig,
    MarchingSquares,
    Metaballs, // Metaballs struct itself
    MetaballsBuilder,
};
use crate::math::geometry::polygon::Polygon; // Custom Polygon struct
use crate::math::geometry::sphere::projection::{SphereProjectionType, SphereProjector};
use crate::math::geometry::sphere::sampling::{SphereSampler, SphereSamplingMethod};
use crate::math::point_distribution::voronoi::{VoronoiCell, VoronoiExtractor};
use crate::math::prelude::*;
use crate::math::probability::SeedResource; // For Vec2, MathError, MathResult, Bounds2D, constants etc.
use crate::next_random_range;
use crate::physics::sim::resources::SimulationParameters;
use bevy::prelude::*;
use spade::{DelaunayTriangulation, Triangulation as SpadeTriangulation}; // For bulk_load_stable and Triangulation trait // Assuming this path is correct // Assuming this macro now takes &mut SeedResource
use std::sync::Arc;

// Constants (seem fine)
const PROJECTED_PLANE_BOUNDS_SIZE: f32 = 2.0;
const METABALL_GRID_RESOLUTION: usize = 100;
const METABALL_FALLOFF_DISTANCE_FACTOR: f32 = 0.5;
const METABALL_STRENGTH: f32 = 1.0;
const METABALL_THRESHOLD: f32 = 0.5;

pub fn init_craton_voronoi_metaball(
    mut commands: Commands,
    sim_params: Res<SimulationParameters>,
    mut seed_resource: ResMut<SeedResource>, // Added SeedResource
) {
    info!("Initializing Craton generation using Voronoi and Metaballs (Refactored).");

    // 1. Bestimme ein Zentrum auf der Kugel für unsere Operationen
    let sphere_sampler =
        match SphereSampler::new(SphereSamplingMethod::Uniform, sim_params.planet_radius_km) {
            Ok(s) => s,
            Err(e) => {
                error!("Failed to create SphereSampler: {:?}", e);
                return;
            }
        };
    let area_center_on_sphere_3d = sphere_sampler.sample_single_point(&mut seed_resource);

    #[cfg(debug_assertions)]
    commands.spawn((
        NormalArrowVisual::new(
            area_center_on_sphere_3d, // Already scaled by sphere_sampler
            0.1 * sim_params.planet_radius_km,
            Color::RED,
        ),
        Name::new("Projection_Center_Debug"),
    ));

    // 2. Erstelle einen Kugelprojektor
    let sphere_projector = match SphereProjector::new(
        SphereProjectionType::Gnomonic {
            center_on_sphere: area_center_on_sphere_3d.normalize_or_zero(), // Gnomonic expects normalized center
        },
        sim_params.planet_radius_km,
    ) {
        Ok(p) => p.with_projection_scale(1.0), // Use with_projection_scale
        Err(e) => {
            error!("Failed to create SphereProjector: {:?}", e);
            return;
        }
    };

    // 3. Generiere initiale 2D-Punkte für Voronoi
    let bounds_2d =
        match Bounds2D::from_center_size(Vec2::ZERO, Vec2::splat(PROJECTED_PLANE_BOUNDS_SIZE)) {
            Ok(b) => b,
            Err(e) => {
                error!("Failed to create Bounds2D for Voronoi points: {:?}", e);
                return;
            }
        };

    let mut initial_bevy_points_for_voronoi: Vec<Vec2> = Vec::new();
    for _ in 0..100 {
        let x = next_random_range!(seed_resource, bounds_2d.min.x, bounds_2d.max.x);
        let y = next_random_range!(seed_resource, bounds_2d.min.y, bounds_2d.max.y);
        initial_bevy_points_for_voronoi.push(Vec2::new(x, y));
    }

    if initial_bevy_points_for_voronoi.len() < 3 {
        warn!(
            "Not enough initial points for Voronoi generation after filtering: {}. Expected at least 3.",
            initial_bevy_points_for_voronoi.len()
        );
        return;
    }

    // 4. Lloyd-Relaxation auf den 2D-Punkten (jetzt Vec2)
    let lloyd_config = LloydConfig {
        max_iterations: 10,
        convergence_tolerance_sq: 1e-4 * 1e-4, // Use squared tolerance from new API
        boundary_handling: LloydBoundaryHandling::Clamp,
    };
    let lloyd_relaxer = LloydRelaxation::new(lloyd_config);
    let (relaxed_bevy_points, _lloyd_stats) =
        match lloyd_relaxer.relax_points(&initial_bevy_points_for_voronoi, &bounds_2d) {
            Ok(result) => result,
            Err(e) => {
                warn!("Lloyd relaxation failed: {:?}", e);
                (initial_bevy_points_for_voronoi, Default::default())
            }
        };

    if relaxed_bevy_points.len() < 3 {
        warn!(
            "Not enough points for Voronoi after Lloyd relaxation: {}. Expected at least 3.",
            relaxed_bevy_points.len()
        );
        return;
    }

    // 5. Erzeuge das 2D-Voronoi-Diagramm
    // Konvertiere Vec2 zu SpadePoint für DelaunayTriangulation
    let spade_points = bevy_to_spade_points(&relaxed_bevy_points);
    if spade_points.len() < 3 {
        // Double check, though relaxed_bevy_points was checked
        warn!(
            "Not enough SpadePoints for Delaunay: {}",
            spade_points.len()
        );
        return;
    }

    let delaunay_triangulation = match DelaunayTriangulation::bulk_load_stable(spade_points) {
        Ok(dt) => dt,
        Err(e) => {
            error!("Failed to create Delaunay triangulation: {:?}", e);
            return;
        }
    };

    if delaunay_triangulation.num_vertices() < 3 {
        warn!(
            "Delaunay triangulation resulted in too few vertices: {}",
            delaunay_triangulation.num_vertices()
        );
        return;
    }

    let voronoi_cells: Vec<VoronoiCell> =
        match VoronoiExtractor::extract_cells(&delaunay_triangulation, Some(&bounds_2d)) {
            Ok(cells) => cells,
            Err(e) => {
                error!("Failed to compute Voronoi cells: {:?}", e);
                return;
            }
        };

    if voronoi_cells.is_empty() {
        warn!("No Voronoi cells generated.");
        return;
    }

    // 6. Erstelle PolygonInfluences für Metaballs aus den Voronoi-Zellen
    let mut metaball_sources: Vec<Arc<dyn FieldInfluence>> = Vec::new();
    for (idx, cell) in voronoi_cells.iter().enumerate() {
        if cell.vertices.len() < 3 {
            continue;
        }
        // cell.vertices ist bereits Vec<Vec2>

        let polygon_2d = match Polygon::closed(cell.vertices.clone()) {
            // Voronoi cells are closed
            Ok(p) => p,
            Err(e) => {
                warn!(
                    "Failed to create polygon from Voronoi cell {}: {:?}",
                    idx, e
                );
                continue;
            }
        };

        let cell_bounds = polygon_2d.bounds().unwrap_or(bounds_2d); // Fallback to main bounds
        let cell_avg_dim = (cell_bounds.width() + cell_bounds.height()) / 2.0;
        let falloff_distance = cell_avg_dim * METABALL_FALLOFF_DISTANCE_FACTOR;

        // Use PolygonInfluenceBuilder
        let influence_builder = PolygonInfluenceBuilder::new()
            .polygon(polygon_2d)
            .strength(METABALL_STRENGTH)
            .falloff_distance(falloff_distance.max(0.01))
            .falloff_type(PolygonFalloff::Gaussian {
                variance_factor: 0.3,
            }) // Example variance_factor
            .interior_only(true) // Typically for craton-like shapes
            .offset(0.0);

        if let Some(influence) = influence_builder.build() {
            metaball_sources.push(Arc::new(influence));
        } else {
            warn!("Failed to build PolygonInfluence for cell {}", idx);
        }
    }

    if metaball_sources.is_empty() {
        warn!("No metaball sources created from Voronoi cells.");
        return;
    }

    // 7. Erzeuge das Metaball-Feld und extrahiere Konturen
    let grid_config = GridConfig {
        width: METABALL_GRID_RESOLUTION,
        height: METABALL_GRID_RESOLUTION,
        cell_size: PROJECTED_PLANE_BOUNDS_SIZE / METABALL_GRID_RESOLUTION as f32,
        threshold: METABALL_THRESHOLD,
    };

    let metaballs_instance = MetaballsBuilder::new()
        .with_config(grid_config.clone()) // Use with_config for GridConfig
        .add_sources(metaball_sources) // add_sources takes an iterator of Arc<dyn FieldInfluence>
        .build();
    // .recompute_field() is called within build() if sources are present.

    // MarchingSquares::extract_contour_segments now uses an iterator
    // If you want all contours in a Vec:
    let contours_2d: Vec<Contour> = MarchingSquaresIterator::new(
        &metaballs_instance.field,
        metaballs_instance.config.threshold,
    )
    .collect();

    // 8. Projiziere die 2D-Konturen zurück auf die 3D-Kugel
    for (i, contour) in contours_2d.iter().enumerate() {
        if contour.vertices.len() < 3 {
            continue;
        }

        let mut world_vertices_3d: Vec<Vec3> = Vec::new();
        for point_2d_contour in &contour.vertices {
            // contour.vertices ist Vec<Vec2>
            match sphere_projector.unproject_point(*point_2d_contour) {
                // Pass by value if Vec2 is Copy
                Ok(mut point_3d_on_sphere) => {
                    point_3d_on_sphere =
                        point_3d_on_sphere.normalize_or_zero() * sim_params.planet_radius_km;
                    world_vertices_3d.push(point_3d_on_sphere);
                }
                Err(e) => {
                    warn!(
                        "Failed to unproject point {:?} from contour {}: {:?}",
                        point_2d_contour, i, e
                    );
                    world_vertices_3d.clear();
                    break;
                }
            }
        }

        if world_vertices_3d.len() >= 3 {
            let craton_name = format!("VoronoiCraton_{}", i);
            info!("Spawning visual for {}", craton_name);

            commands.spawn((
                PolygonVisual {
                    points: world_vertices_3d,
                    color: Color::rgba(0.0, 0.5 + (i % 5) as f32 * 0.1, 0.0, 0.7),
                },
                Name::new(craton_name),
            ));
        } else {
            warn!(
                "Skipping contour {} due to insufficient unprojected vertices ({}).",
                i,
                world_vertices_3d.len()
            );
        }
    }
    info!("Finished Craton generation using Voronoi and Metaballs (Refactored).");
}
