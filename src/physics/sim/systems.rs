use super::resources::SimulationParameters;
use crate::debug::visualization::{normal_vector::NormalArrowVisual, polygon::PolygonVisual};
use crate::physics::geology::tectonic::feature::Craton;
use crate::physics::sim::time::resources::TickHistory;
use bevy::prelude::*; // Hinzufügen für clear_tick_history_system

// BEISPIEL: System für Plattentektonik-Dynamik
pub fn run_plate_dynamics_system(
    // Zugriff auf Weltzustand: Query für Platten-Entitäten, Ressourcen etc.
    // mut query_plates: Query<&mut PlateComponent>,
    sim_params: Res<SimulationParameters>, // Für time_step_ma etc.
                                           // mut geology_events: EventWriter<GeologyEvent>,
) {
    // Diese Funktion würde NUR laufen, wenn State::Running oder SingleStep im Paused-State.
    // Die Run-Condition wird in main.rs definiert.

    // Pseudo-Code:
    // for mut plate in query_plates.iter_mut() {
    //     plate.position += plate.velocity * sim_params.time_step_ma;
    //     // ... Kollisionserkennung, Interaktionen ...
    // }
    info!(
        "Running actual simulation dynamics for t_advanced = {:.2} (next step after this is {:.2})",
        sim_params.advanced_simulation_time_ma,
        sim_params.advanced_simulation_time_ma + sim_params.time_step_ma
    );
    // Dieses System modifiziert den "lebenden" Weltzustand.
}

// BEISPIEL: System, das das Ergebnis der Voronoi-Generierung verarbeitet
pub fn process_generated_shapes_system(
    mut commands: Commands,
    sim_params: Res<SimulationParameters>, // um den aktuellen Zeitstempel zu bekommen
                                           // query_existing_cratons: Query<Entity, With<CratonComponent>>, // um alte ggf. zu despawnen
                                           // mut gizmos: Gizmos, // Falls benötigt für Debug-Zeichnung
) {
    // Diese Funktion würde NUR laufen, wenn State::Running oder SingleStep im Paused-State.
    // Sie wird NACH run_plate_dynamics_system laufen.

    // let config = VoronoiMergedShapeConfig {
    //     rng_seed: Some((sim_params.advanced_simulation_time_ma * 1000.0) as u64 + 1), // Anderer Seed
    //     ..Default::default()
    // };
    // let normalized_2d_polygon = generate_normalized_voronoi_merged_shape(&config);

    // if !normalized_2d_polygon.is_empty() {
    //     let sphere_center_normalized = random_point_on_unit_sphere(); // Neuer zufälliger Punkt
    //     let sphere_radius = 1.0;
    //     let max_angular_extent_deg = sim_params.gen_max_craton_angular_extent_deg; // Aus SimParams

    //     let projected_3d_vertices = project_normalized_2d_polygon_to_sphere(
    //         &normalized_2d_polygon,
    //         sphere_center_normalized,
    //         sphere_radius,
    //         max_angular_extent_deg,
    //     );
    //     if !projected_3d_vertices.is_empty() {
    //         // ... CratonComponent spawnen ...
    //         info!("New Craton spawned at t_advanced = {:.2}", sim_params.advanced_simulation_time_ma);
    //     }
    // }
}

pub fn cleanup_simulation_entities_system(
    mut commands: Commands,
    craton_query: Query<Entity, With<Craton>>,
    polygon_query: Query<Entity, With<PolygonVisual>>,
    arrow_query: Query<Entity, With<NormalArrowVisual>>,
    // Füge hier weitere Queries für andere Simulations-Entitäten hinzu, falls nötig
) {
    info!("Cleaning up simulation entities...");
    for entity in craton_query.iter() {
        commands.entity(entity).despawn_recursive();
    }
    for entity in polygon_query.iter() {
        commands.entity(entity).despawn_recursive();
    }
    for entity in arrow_query.iter() {
        commands.entity(entity).despawn_recursive();
    }
    info!("Simulation entities cleaned up.");
}

pub fn clear_tick_history_system(mut history: ResMut<TickHistory>) {
    info!("Clearing TickHistory.");
    history.clear();
    // Der initiale Snapshot wird durch `create_initial_snapshot` hinzugefügt,
    // das nach dieser Systemausführung im Resetting-State aufgerufen wird.
}
