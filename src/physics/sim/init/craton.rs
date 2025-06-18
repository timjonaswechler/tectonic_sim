// src/physics/sim/init/craton.rs
use crate::math::{
    self,
    // Wichtig: Stelle sicher, dass PolygonProperties im Scope ist, wenn du .area() etc. aufrufst
    geometry::polygon::PolygonProperties,
    prelude::{
        Bounds2D,
        Polygon, // Polygon für 2D-Polygonobjekt
        SeedResource,
        SphereProjectionType,
        SphereProjector,
        SphereSampler,
        SphereSamplingMethod,
        Vector2DExt, // Für Polygon-Centroid, wenn nicht PolygonProperties verwendet wird
        VoronoiBuilder,
        VoronoiConfig,
    },
};
use crate::physics::sim::resources::SimulationParameters;
use bevy::prelude::*;
use std::f32::consts::PI;

// --- Konstanten für die Kraton-Generierung ---
const NUM_CRATONS_TO_GENERATE: usize = 5;
const VORONOI_INITIAL_POINT_COUNT: usize = 100; // Mehr Punkte für bessere Formen
const VORONOI_LLOYD_ITERATIONS: usize = 5;
const VORONOI_TARGET_CELL_COUNT: usize = 3;
const VORONOI_SMOOTHING_ITERATIONS: usize = 1;
const TARGET_CRATON_ANGULAR_EXTENT_DEGREES: f32 = 10.0;

// Größe der 2D "Leinwand" für die Voronoi-Generierung (z.B. -50 bis 50)
const INTERNAL_VORONOI_CALCULATION_EXTENT: f32 = 100.0;

#[derive(Debug, Component)]
pub struct ProjectedCraton {
    pub center_on_sphere: Vec3,
    pub vertices_3d: Vec<Vec3>,
    pub original_polygon_area_2d_scaled: f32, // Fläche des skalierten 2D-Polygons
}

pub fn init_craton_voronoi_metaball(
    mut commands: Commands,
    sim_params: Res<SimulationParameters>,
    mut seed_resource: ResMut<SeedResource>,
) {
    info!("System: Creating projected cratons...");
    let planet_radius = sim_params.planet_radius_km;
    let target_angular_extent_rad = TARGET_CRATON_ANGULAR_EXTENT_DEGREES.to_radians();

    // Dies sind die Bounds, in denen die 2D Voronoi-Berechnung intern stattfindet.
    // Das resultierende Polygon wird dann skaliert und unprojiziert.
    let internal_voronoi_bounds_2d = match Bounds2D::from_center_size(
        Vec2::ZERO,
        Vec2::splat(INTERNAL_VORONOI_CALCULATION_EXTENT),
    ) {
        Ok(bounds) => bounds,
        Err(e) => {
            error!("Failed to create internal 2D bounds for Voronoi: {:?}", e);
            return;
        }
    };

    for craton_idx in 0..NUM_CRATONS_TO_GENERATE {
        // 1. Wähle zufälligen Mittelpunkt auf der Kugel
        let sphere_sampler =
            SphereSampler::new(SphereSamplingMethod::Uniform, planet_radius).unwrap();
        let craton_center_on_sphere = sphere_sampler
            .sample_single_point(&mut seed_resource)
            .normalize_or_zero()
            * planet_radius;
        debug!(
            "Craton {}: Center on sphere: {:?}",
            craton_idx, &craton_center_on_sphere
        );

        // 2. Voronoi Konfiguration mit den *internen* Berechnungsbounds
        let voronoi_config = VoronoiConfig {
            initial_point_count: VORONOI_INITIAL_POINT_COUNT,
            lloyd_iterations: VORONOI_LLOYD_ITERATIONS,
            target_cell_count: VORONOI_TARGET_CELL_COUNT,
            smoothing_iterations: VORONOI_SMOOTHING_ITERATIONS,
            calculation_bounds: Some(internal_voronoi_bounds_2d.clone()), // Clip gegen diese internen Bounds
            boundary_padding_factor: 0.15, // Padding relativ zu internal_voronoi_bounds_2d
            seed: Some(seed_resource.seed as u64 + craton_idx as u64),
        };

        let voronoi_builder = match VoronoiBuilder::new(voronoi_config.clone()) {
            // Klonen, da config für Log gebraucht wird
            Ok(builder) => builder,
            Err(e) => {
                error!(
                    "Failed to create VoronoiBuilder for craton {}: {:?}",
                    craton_idx, e
                );
                continue;
            }
        };

        // 3. Generiere die 2D Voronoi-Polygone (Vertices sind im Bereich von internal_voronoi_bounds_2d)
        let generated_2d_polygon_groups =
            match voronoi_builder.generate_polygons(&mut seed_resource) {
                Ok(polys) => polys,
                Err(e) => {
                    // Hier sollten die Logs aus dem Builder mehr Infos geben
                    error!(
                        "Craton {}: Failed to generate Voronoi polygons. Config: {:?}, Error: {:?}",
                        craton_idx,
                        voronoi_config,
                        e // builder.config hinzugefügt
                    );
                    continue;
                }
            };

        if generated_2d_polygon_groups.is_empty() || generated_2d_polygon_groups[0].is_empty() {
            warn!(
                "Craton {}: Voronoi generation yielded no suitable 2D polygon.",
                craton_idx
            );
            continue;
        }
        let raw_2d_vertices = &generated_2d_polygon_groups[0];

        // 4. Skaliere das 2D Polygon, sodass seine "Ausdehnung" der Ziel-Bogenlänge auf der Kugel entspricht.
        //    Die Ziel-Bogenlänge für den *Radius* der Projektion (halbe Ausdehnung)
        let target_projection_radius_on_sphere = planet_radius * (target_angular_extent_rad / 2.0);

        // Bestimme die aktuelle "Ausdehnung" des rohen 2D Polygons.
        // Eine Möglichkeit: Nimm die Bounds und den max_element davon.
        let raw_2d_poly_obj = Polygon::new(raw_2d_vertices.clone()).unwrap(); // Für Bounds etc.
        let raw_2d_bounds = raw_2d_poly_obj.bounds().unwrap_or_else(|| {
            // Sollte nicht passieren, wenn raw_2d_vertices nicht leer ist
            Bounds2D::from_points(Vec2::ZERO, Vec2::ZERO)
        });
        // Aktueller "effektiver Radius" des Polygons von seinem Zentrum.
        // Wir nehmen hier die halbe Größe der Bounding Box als Annäherung.
        let current_max_extent_from_center_2d = raw_2d_bounds.size().max_element() / 2.0;

        let scale_factor_2d = if current_max_extent_from_center_2d > math::utils::constants::EPSILON
        {
            target_projection_radius_on_sphere / current_max_extent_from_center_2d
        } else {
            1.0 // Vermeide Division durch Null, falls Polygon degeneriert
        };

        // Polygon-Zentrum finden, um Skalierung darum durchzuführen (wenn nicht eh schon um 0,0)
        let raw_poly_centroid = raw_2d_poly_obj.geometric_centroid().unwrap_or(Vec2::ZERO); // geometric_centroid für Flächenschwerpunkt

        let scaled_2d_vertices: Vec<Vec2> = raw_2d_vertices
            .iter()
            .map(|&v| {
                (v - raw_poly_centroid) * scale_factor_2d // + raw_poly_centroid (wenn nicht auf Vec2::ZERO zentriert unprojeziert wird)
                // Für Lambert Azimuthal/Gnomonic, die vom 2D-Ursprung aus projizieren,
                // wollen wir die skalierten Punkte direkt um den 2D-Ursprung.
            })
            .collect();

        // 5. Projiziere (Unproject) die *skalierten* 2D-Vertices auf die Kugel
        let projector = SphereProjector::new(
            SphereProjectionType::LambertAzimuthal {
                center_on_sphere: craton_center_on_sphere,
            },
            planet_radius,
        )
        .unwrap()
        .with_projection_scale(1.0); // Wichtig: scale hier auf 1, da die Skalierung schon in scaled_2d_vertices steckt.

        let mut projected_3d_vertices = Vec::with_capacity(scaled_2d_vertices.len());
        for &vertex_2d_scaled in &scaled_2d_vertices {
            match projector.unproject_point(vertex_2d_scaled) {
                Ok(vertex_3d) => projected_3d_vertices.push(vertex_3d),
                Err(e) => {
                    warn!(
                        "Craton {}: Failed to unproject scaled point {:?} to sphere: {:?}. Skipping point.",
                        craton_idx, vertex_2d_scaled, e
                    );
                }
            }
        }

        if projected_3d_vertices.len() < 3 {
            warn!(
                "Craton {}: Projected 3D polygon has fewer than 3 vertices ({}). Skipping.",
                craton_idx,
                projected_3d_vertices.len()
            );
            continue;
        }

        // 6. Speichern
        let scaled_poly_2d_for_area = Polygon::new(scaled_2d_vertices)
            .unwrap_or_else(|_| Polygon::new(vec![Vec2::ZERO, Vec2::X, Vec2::Y]).unwrap());
        let craton_entity = commands
            .spawn(ProjectedCraton {
                center_on_sphere: craton_center_on_sphere,
                vertices_3d: projected_3d_vertices.clone(),
                original_polygon_area_2d_scaled: scaled_poly_2d_for_area.area(),
            })
            .id();

        debug!(
            "Craton {} (Entity: {:?}): Generated. 3D Vertices: {}, Scaled 2D Area: {:.2}, Target Proj Radius: {:.2}, ScaleFactor: {:.3}",
            craton_idx,
            craton_entity,
            projected_3d_vertices.len(),
            scaled_poly_2d_for_area.area(),
            target_projection_radius_on_sphere,
            scale_factor_2d
        );
    }
    info!("System: Finished creating projected cratons.");
}
