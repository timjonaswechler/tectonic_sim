use bevy::gizmos::GizmoPlugin;
// ./src/main.rs
use bevy::prelude::*;
use bevy_egui::EguiPlugin;
use bevy_panorbit_camera::PanOrbitCameraPlugin;

// Eigene Module deklarieren
pub mod debug;
pub mod math;
pub mod physics;
pub mod setup; // Wird für Kamera-Setup und Szene verwendet

// Importiere spezifische Elemente aus unseren Modulen
use bevy::gizmos::prelude::GizmoConfig;
use debug::{
    ui::simulation_control_ui_system,
    visualization::{
        normal_vector::draw_normal_arrows_system, polygon::draw_polygons_system,
        sphere_grid::draw_sphere_grid_gizmos,
    },
}; // Importiere SVG-Dumper
use physics::sim::init::craton::init_craton;
use physics::sim::resources::*; // Geänderter Pfad
use physics::sim::state::*;
use physics::sim::systems::*;
use physics::sim::time::resources::*;
use physics::sim::time::{
    resources::{SimulationSnapshot, TickHistory}, // Korrekte Pfade
    systems::*,
};
use setup::setup_scene;
fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(EguiPlugin)
        .add_plugins(PanOrbitCameraPlugin)
        .init_resource::<ExecuteSingleStepRequest>()
        .init_resource::<ExecuteStepBackwardRequest>()
        .init_resource::<SimulationParameters>()
        .insert_resource(TickHistory::new(1000)) // Max 1000 Snapshots
        .init_state::<SimulationState>()
        .add_systems(Startup, (setup_scene,).chain())
        // --- Initialisierungsphase ---
        .add_systems(
            OnEnter(SimulationState::Initializing),
            (
                // Systeme, die nur einmal zu Beginn der Initialisierung laufen
                setup_initial_world_conditions_system, // Generiert z.B. erste Platten
                init_craton,
            ),
        )
        .add_systems(
            Update,
            // Wechselt von Initializing zu Paused, nachdem Startup-Systeme gelaufen sind.
            (|mut next_state: ResMut<NextState<SimulationState>>,
              current_state: Res<State<SimulationState>>| {
                if *current_state == SimulationState::Initializing {
                    next_state.set(SimulationState::Paused);
                    info!("Initialization sequence complete, simulation is Paused.");
                }
            })
            .run_if(in_state(SimulationState::Initializing)),
        )
        .add_systems(
            Update,
            (
                // Block 1: UI und Anfrageverarbeitung
                simulation_control_ui_system,
                handle_simulation_requests_system, // Reagiert auf UI, setzt States/Flags
                handle_timeline_interaction_system, // Reagiert auf UI (Slider)
                // Block 2: Zeit vorrücken für Playback (wenn Running und nicht am Frontier)
                auto_advance_display_time_system.run_if(in_state(SimulationState::Running)),
                // Block 3: Zustand aus History laden (wenn nötig)
                load_state_from_history_system.run_if(should_load_from_history_condition),
                // Block 4: Simulationsschritt ausführen und Snapshot erstellen
                // Die Gruppe läuft, wenn:
                // - State ist Running ODER
                // - State ist Paused UND SingleStepRequest ist true
                // UND
                // - should_load_from_history_condition ist false (wir laden nicht gerade)
                (
                    run_plate_dynamics_system, // Deine Simulationslogik hier
                    record_snapshot_system // Erstellt Snapshot danach
                        .after(run_plate_dynamics_system), // explizite Ordnung innerhalb der Gruppe
                )
                    .run_if(
                        (in_state(SimulationState::Running).or_else(
                            in_state(SimulationState::Paused).and_then(
                                // Verwende resource_equals von Bevy oder deine eigene Implementierung
                                bevy::ecs::schedule::common_conditions::resource_equals(
                                    ExecuteSingleStepRequest(true),
                                ),
                            ),
                        ))
                        .and_then(not(should_load_from_history_condition)),
                    ),
            )
                .chain(), // .chain() auf das gesamte Tupel der Update-Systeme
        )
        .add_systems(
            Update,
            (
                draw_normal_arrows_system,
                draw_polygons_system,
                draw_sphere_grid_gizmos,
            ),
        )
        .run();
}

fn create_initial_snapshot(
    mut history: ResMut<TickHistory>,
    sim_params: Res<SimulationParameters>, // Nur lesend, um initiale Zeit zu bekommen
) {
    if history.is_empty() {
        let initial_snapshot = SimulationSnapshot {
            time_ma: sim_params.advanced_simulation_time_ma, // Sollte 0.0 sein
            dummy_value: 0.0,                                // Initialer Dummy-Wert
        };
        history.add_snapshot(initial_snapshot);
        info!("Initial snapshot at t=0 created for display.");
    }
}
