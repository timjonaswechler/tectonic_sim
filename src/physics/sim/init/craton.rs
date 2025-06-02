use crate::debug::visualization::normal_vector::NormalArrowVisual;
use crate::debug::visualization::polygon::PolygonVisual;
use crate::math::sphere::algorithms::random_point_on_unit_sphere;
use crate::math::vector::functions::*;
use crate::physics::geology::tectonic::generation::craton_builder::CratonBuilder;
use crate::physics::sim::resources::SimulationParameters;
use bevy::prelude::*;

pub fn init_craton(mut commands: Commands, sim_params: Res<SimulationParameters>) {
    let craton_area_center = random_point_on_unit_sphere();
    #[cfg(debug_assertions)]
    {
        let normal_vector = NormalArrowVisual::new(
            craton_area_center * sim_params.planet_radius_km,
            0.05 * sim_params.planet_radius_km,
            Color::rgba(1.0, 0.0, 1.0, 0.5),
        );
        commands.spawn((normal_vector,));
    }
    for i in 0..sim_params.craton_num {
        let name_char = (b'A' + i as u8) as char;
        let craton_name = format!("Craton_{}", name_char);
        info!("Entity Name: {}", craton_name);

        let craton_center = random_point_within_distance(craton_area_center, Dist::Chord(1.0), 1.0);

        #[cfg(debug_assertions)]
        {
            let normal_vector = NormalArrowVisual::new(
                craton_center,
                0.05 * sim_params.planet_radius_km,
                Color::rgba(1.0, 1.0, 0.0, 0.5),
            );
            commands.spawn((
                normal_vector,
                Name::new(format!("NormalArrow_{}", name_char)),
            ));
        }

        let craton_builder_instance =
            CratonBuilder::new(craton_name.clone(), craton_center, &sim_params);
        let world_craton_vertices = craton_builder_instance.vertices.0.clone();

        let craton_entity_transform = Transform::from_translation(craton_center);

        commands.spawn((craton_builder_instance, craton_entity_transform));

        if world_craton_vertices.is_empty() {
            warn!("Craton {} has no vertices!", name_char);
            continue;
        } else {
        }

        commands.spawn(PolygonVisual {
            points: world_craton_vertices,
            color: Color::rgba(0.0, 1.0, 0.0, 0.5),
        });
    }
}
