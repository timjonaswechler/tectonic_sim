// ./src/main.rs
use bevy::gizmos::GizmoPlugin;
use bevy::prelude::*;
use bevy_egui::EguiPlugin;
use bevy_inspector_egui::quick::WorldInspectorPlugin; // Import für den Inspector
use bevy_panorbit_camera::PanOrbitCameraPlugin;

// Eigene Module deklarieren
pub mod debug;
pub mod math;
pub mod physics;
pub mod setup; // Wird für Kamera-Setup und Szene verwendet

// Importiere spezifische Elemente aus unseren Modulen
use bevy::gizmos::prelude::GizmoConfig;
use debug::{
    ui::*,
    visualization::{
        normal_vector::draw_normal_arrows_system, polygon::draw_polygons_system,
        sphere_grid::draw_sphere_grid_gizmos,
    },
};
use egui_dock::{DockState, NodeIndex, Split, TabDestination, Tree};
use math::probability::seed::resource::SeedResource;
use physics::sim::{
    init::craton::init_craton,
    resources::*,
    state::*,
    systems::*,
    time::resources::*,
    time::{
        resources::{SimulationSnapshot, TickHistory},
        systems::*,
    },
};
use setup::setup_scene;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(EguiPlugin)
        .add_plugins(PanOrbitCameraPlugin)
        .add_plugins(WorldInspectorPlugin::new())
        .init_resource::<ExecuteSingleStepRequest>()
        .init_resource::<ExecuteStepBackwardRequest>()
        .init_resource::<SimulationParameters>()
        .insert_resource(TickHistory::new(1000)) // Max 1000 Snapshots
        .insert_resource(SeedResource::default())
        .init_state::<SimulationState>()
        .add_systems(Startup, (setup_scene, create_initial_snapshot).chain()) // create_initial_snapshot beim Start
        // --- Initialisierungsphase ---
        .add_systems(
            OnEnter(SimulationState::Initializing),
            (
                setup_initial_world_conditions_system,
                init_craton,
                // create_initial_snapshot wird jetzt in Startup aufgerufen, oder wir brauchen eine Variante für hier
            ),
        )
        .add_systems(
            Update,
            (|mut next_state: ResMut<NextState<SimulationState>>,
              current_state: Res<State<SimulationState>>| {
                if *current_state == SimulationState::Initializing {
                    next_state.set(SimulationState::Paused);
                    info!("Initialization sequence complete, simulation is Paused.");
                }
            })
            .run_if(in_state(SimulationState::Initializing)),
        )
        // --- Reset-Phase ---
        .add_systems(
            OnEnter(SimulationState::Resetting),
            (
                cleanup_simulation_entities_system,
                clear_tick_history_system,             // History leeren
                setup_initial_world_conditions_system, // SimParams (Zeiten) zurücksetzen
                init_craton,                           // Kratone neu generieren
                create_initial_snapshot,               // Initialen Snapshot für t=0 erstellen
            )
                .chain(), // Wichtig für die Reihenfolge!
        )
        .add_systems(
            Update,
            (|mut next_state: ResMut<NextState<SimulationState>>,
              current_state: Res<State<SimulationState>>| {
                if *current_state == SimulationState::Resetting {
                    next_state.set(SimulationState::Paused);
                    info!("Reset sequence complete, simulation is Paused.");
                }
            })
            .run_if(in_state(SimulationState::Resetting)),
        )
        // --- Haupt-Update-Logik ---
        .add_systems(
            Update,
            (
                // Block 1: UI und Anfrageverarbeitung
                simulation_control_ui_system,
                handle_simulation_requests_system,
                handle_timeline_interaction_system,
                // Block 2: Zeit vorrücken für Playback
                auto_advance_display_time_system.run_if(in_state(SimulationState::Running)),
                // Block 3: Zustand aus History laden
                load_state_from_history_system.run_if(should_load_from_history_condition),
                // Block 4: Simulationsschritt ausführen und Snapshot erstellen
                (
                    run_plate_dynamics_system,
                    record_snapshot_system.after(run_plate_dynamics_system),
                )
                    .run_if(
                        (in_state(SimulationState::Running).or_else(
                            in_state(SimulationState::Paused).and_then(
                                bevy::ecs::schedule::common_conditions::resource_equals(
                                    ExecuteSingleStepRequest(true),
                                ),
                            ),
                        ))
                        .and_then(not(should_load_from_history_condition)),
                    ),
            )
                .chain(),
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

// Dieses System wird beim Startup und beim Reset verwendet, um den initialen Snapshot zu erstellen.
fn create_initial_snapshot(
    mut history: ResMut<TickHistory>,
    sim_params: Res<SimulationParameters>,
) {
    // Stellen sicher, dass die History leer ist, bevor der erste Snapshot hinzugefügt wird
    // oder dass wir einen expliziten Reset-Kontext haben, in dem das okay ist.
    // Da `clear_tick_history_system` im Reset-Fall vorher läuft, ist `history.is_empty()` wieder wahr.
    if history.is_empty() {
        let initial_snapshot = SimulationSnapshot {
            time_ma: sim_params.advanced_simulation_time_ma, // Sollte 0.0 sein nach Reset/Init
            dummy_value: 0.0,
        };
        history.add_snapshot(initial_snapshot);
        info!(
            "Initial snapshot at t={} created for display.",
            sim_params.advanced_simulation_time_ma
        );
    } else {
        warn!(
            "Attempted to create initial snapshot, but history was not empty and not in reset context. Time: {}",
            sim_params.advanced_simulation_time_ma
        );
    }
}
