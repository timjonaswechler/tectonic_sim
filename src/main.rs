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

use debug::{
    ui::simulation_control_ui_system, visualization::normal_vector::draw_normal_arrows_system,
}; // Importiere SVG-Dumper
use physics::sim::resources::*; // Geänderter Pfad
use physics::sim::state::*;
use physics::sim::time::resources::*;

use physics::sim::time::{
    resources::{SimulationSnapshot, TickHistory}, // Korrekte Pfade
    systems::*,
};
use setup::setup_scene;
fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .init_resource::<ExecuteSingleStepRequest>()
        .init_resource::<ExecuteStepBackwardRequest>()
        .add_plugins(EguiPlugin)
        .add_plugins(PanOrbitCameraPlugin)
        .init_resource::<SimulationParameters>()
        .insert_resource(TickHistory::new(1000)) // Max 1000 Snapshots
        .init_state::<SimulationState>()
        .add_systems(
            Startup,
            (
                // Deine 3D-Szene etc.
                setup_scene,
                // ------------------------------------------
                create_initial_snapshot, // Läuft nach der 3D- und ggf. 2D-Initialisierung
            )
                .chain(),
        )
        // --- Initialisierungsphase ---
        .add_systems(
            OnEnter(SimulationState::Initializing),
            (
                // Systeme, die nur einmal zu Beginn der Initialisierung laufen
                setup_initial_world_conditions_system, // Generiert z.B. erste Platten
                create_initial_snapshot,               // Erstellt den allerersten Snapshot bei t=0
            ),
        )
        .add_systems(
            Update,
            // Systeme, die während der Initialisierung laufen, bis sie fertig sind
            // z.B. eine iterative Plattengenerierung, die mehrere Frames braucht
            // Dieses System muss dann den State zu Running/Paused wechseln
            finish_initialization_system.run_if(in_state(SimulationState::Initializing)),
        )
        // --- Hauptsimulationslogik (wenn Running oder Paused und Einzelschritt) ---
        // Diese Systeme ersetzen Teile des alten simulation_driver_system
        .add_systems(
            Update,
            (
                // Diese Reihenfolge ist wichtig!
                // 1. UI und Anfragen verarbeiten (setzt Parameter wie paused, step_request)
                simulation_control_ui_system, // Hier werden execute_..._request Flags gesetzt
                // 2. Anfragen auswerten und ggf. State-Änderungen oder Aktionen auslösen
                handle_simulation_requests_system, // Neues System, das auf Flags reagiert
                // 3. Timeline-Interaktion (kann auch display_time verändern)
                handle_timeline_interaction_system,
                // 4. Den eigentlichen Simulationsschritt ausführen, wenn nötig
                advance_simulation_step_system.run_if(
                    // Läuft, wenn Running ODER Paused UND ein Einzelschritt angefordert wurde
                    // Die Logik, ob wirklich ein neuer Schritt gemacht wird, ist *in* dem System
                    in_state(SimulationState::Running).or_else(
                        in_state(SimulationState::Paused).and_then(resource_exists_and_equals(
                            // Dieses Flag wird von handle_simulation_requests_system gesetzt
                            ExecuteSingleStepFlag(true),
                        )),
                    ),
                ),
                // 5. Zustand basierend auf History laden, wenn durch Timeline oder Step Backward getriggert
                load_state_from_history_system.run_if(
                    // Läuft, wenn die display_time sich geändert hat und nicht am Frontier ist
                    // ODER wenn ein Step Backward angefordert wurde
                    // Dieses System muss intelligent sein, um nicht ständig zu laden
                    should_load_from_history_condition,
                ),
            )
                .chain(), // Stelle sicher, dass sie in dieser Reihenfolge laufen
        )
        // Visualisierung
        .add_systems(Update, draw_normal_arrows_system)
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
