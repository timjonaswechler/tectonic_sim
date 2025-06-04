use crate::math::geometry::sphere::projection::project_normalized_2d_polygon_to_sphere;
use crate::math::geometry::voronoi::generate_normalized_voronoi_merged_shape;
use crate::{
    math::geometry::voronoi::VoronoiMergedShapeConfig,
    physics::sim::resources::SimulationParameters,
};
use bevy::prelude::Bundle;
use bevy::prelude::*;
use rand::{random_range, rng};

#[derive(Component, Debug, Clone)]
pub struct CratonVertices(pub Vec<Vec3>);

#[derive(Debug, Clone, Bundle)]
pub struct CratonBuilder {
    pub name: Name,
    pub vertices: CratonVertices,
}

impl CratonBuilder {
    pub fn new(
        name_str: String,
        center_of_craton: Vec3,
        sim_params: &Res<SimulationParameters>,
    ) -> Self {
        let mut craton_vertices_data = Vec::new();

        let config = VoronoiMergedShapeConfig {
            rng_seed: Some((sim_params.advanced_simulation_time_ma * 1000.0) as u64 + 1), // Anderer Seed
            ..Default::default()
        };

        let normalized_2d_polygon = generate_normalized_voronoi_merged_shape(&config);

        if !normalized_2d_polygon.is_empty() {
            let sphere_center_normalized = center_of_craton.normalize();
            let sphere_radius = sim_params.planet_radius_km;
            let max_angular_extent_deg =
                max_craton_angular_deg(config.target_connected_cell_group_size as f32); // Aus SimParams

            let projected_3d_vertices = project_normalized_2d_polygon_to_sphere(
                &normalized_2d_polygon,
                sphere_center_normalized,
                sphere_radius,
                max_angular_extent_deg,
            );
            if !projected_3d_vertices.is_empty() {
                craton_vertices_data = projected_3d_vertices;
            }
        }

        Self {
            name: Name::new(name_str),
            vertices: CratonVertices(craton_vertices_data),
        }
    }
}

fn max_craton_angular_deg(x: f32) -> f32 {
    let f = 0.002 * (x - 3.0).powi(2) + 0.05;
    let g = (0.002 * (x - 3.0).powi(2) + 0.05) * 1.5;
    random_range(f..g) * 200.0 // Skaliert auf auf werte vons
}
