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
use physics::sim::resources::SimulationParameters; // Geänderter Pfad
use physics::sim::time::{
    resources::{SimulationSnapshot, TickHistory}, // Korrekte Pfade
    systems::{handle_timeline_interaction_system, simulation_driver_system},
};
use setup::setup_scene;
fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(EguiPlugin)
        .add_plugins(PanOrbitCameraPlugin)
        .init_resource::<SimulationParameters>()
        .insert_resource(TickHistory::new(1000)) // Max 1000 Snapshots
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
        // UI und Parameter-Anpassung zuerst
        .add_systems(Update, simulation_control_ui_system)
        // Dann die Logik für Timeline-Interaktion (Slider)
        .add_systems(
            Update,
            handle_timeline_interaction_system
                .after(simulation_control_ui_system) // Reagiert auf UI-Änderungen
                .before(simulation_driver_system), // Zustand setzen, bevor Sim-Driver läuft
        )
        // Dann der Haupt-Simulations-Treiber
        .add_systems(Update, simulation_driver_system)
        .add_systems(Update, draw_normal_arrows_system) // Füge das Zeichen-System hinzu
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
