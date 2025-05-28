use super::resources::{SimulationSnapshot, TickHistory};
use crate::debug::visualization::normal_vector::NormalArrowVisual;
use crate::debug::visualization::svg::create_debug_svg_for_normalized_polygon;
use crate::math::polygon::voronoi::*;
use crate::math::sphere::algorithms::random_point_on_unit_sphere;
use crate::math::sphere::projection::project_normalized_2d_polygon_to_sphere;
use crate::physics::geology::tectonic::components::CratonComponent; // Für Craton-Spawn
use crate::physics::sim::resources::SimulationParameters; // Importiere aus dem übergeordneten Sim-Modul
use bevy::prelude::*;
use geo::LineString as GeoLineString;
use geo::Point as GeoPoint; // Für geo::Polygon

pub fn simulation_driver_system(
    mut sim_params: ResMut<SimulationParameters>,
    mut history: ResMut<TickHistory>,
    mut commands: Commands,
    mut gizmos: Gizmos,
    // Query für Simulationsentitäten
) {
    if sim_params.execute_step_backward_request {
        if let Some(current_idx) = history.current_display_snapshot_index {
            if current_idx > 0 {
                let prev_idx = current_idx - 1;
                // Extract data from snapshot to drop the immutable borrow of `history`
                // before attempting a mutable borrow later.
                let snapshot_data = history
                    .get_snapshot_by_index(prev_idx)
                    .map(|s| (s.time_ma, s.dummy_value));

                if let Some((snapshot_time, snapshot_dummy)) = snapshot_data {
                    sim_params.display_time_ma = snapshot_time;
                    history.current_display_snapshot_index = Some(prev_idx);
                    // TODO: WELTZUSTAND AUF DATEN DES VORHERIGEN SNAPSHOTS LADEN (verwende snapshot_time, snapshot_dummy)
                    info!(
                        "Stepped backward. Displaying snapshot at: {:.3} Ma. Dummy: {:.2}",
                        snapshot_time, snapshot_dummy
                    );
                    sim_params.paused = true; // Nach einem Schritt zurück immer pausieren
                }
            }
        }
        sim_params.execute_step_backward_request = false; // Anfrage verarbeitet
        return; // Keine weitere Aktion in diesem Frame für Step Backward
    }

    // Bestimme, ob ein Vorwärts-Schritt ausgeführt werden soll
    let should_execute_forward_step = sim_params.execute_single_step_request || !sim_params.paused;

    if should_execute_forward_step {
        sim_params.execute_single_step_request = false; // Anfrage verarbeitet, falls es ein Single-Step war

        let at_simulation_frontier = history
            .current_display_snapshot_index
            .map_or(true, |idx| idx == history.len() - 1);

        if at_simulation_frontier {
            // --- NEUEN SIMULATIONSSCHRITT AUSFÜHREN ---
            //TODO: Hier die Logik für den neuen Simulationsschritt einfügen
            sim_params.advanced_simulation_time_ma += sim_params.time_step_ma;
            sim_params.display_time_ma = sim_params.advanced_simulation_time_ma;

            let new_dummy_value = (sim_params.advanced_simulation_time_ma * 10.0 % 100.0) as f32;
            info!(
                "Simulating new step. Time: {:.3} Ma. Dummy: {:.2}",
                sim_params.advanced_simulation_time_ma, new_dummy_value
            );

            let new_snapshot = SimulationSnapshot {
                time_ma: sim_params.advanced_simulation_time_ma,
                dummy_value: new_dummy_value,
            };
            history.add_snapshot(new_snapshot);
            // Wenn ein einzelner Schritt ausgeführt wurde (nicht im Play-Modus), pausieren wir danach
            if !(!sim_params.paused) {
                // d.h. es war ein single step
                sim_params.paused = true;
            }
        } else if !sim_params.paused {
            // Im "Play"-Modus und nicht am Frontier -> Historie abspielen
            // --- DURCH HISTORIE VORWÄRTS "SPIELEN" ---
            if let Some(current_idx) = history.current_display_snapshot_index {
                let next_idx = current_idx + 1;
                if next_idx < history.len() {
                    // Extract data from snapshot to drop the immutable borrow of `history`
                    // before attempting a mutable borrow later.
                    let snapshot_data = history
                        .get_snapshot_by_index(next_idx)
                        .map(|s| (s.time_ma, s.dummy_value));

                    if let Some((snapshot_time, snapshot_dummy)) = snapshot_data {
                        sim_params.display_time_ma = snapshot_time;
                        history.current_display_snapshot_index = Some(next_idx);
                        // TODO: WELTZUSTAND AUF `next_snapshot` DATEN LADEN (verwende snapshot_time, snapshot_dummy)
                        info!(
                            "Playing history forward. Displaying snapshot at: {:.3} Ma. Dummy: {:.2}",
                            snapshot_time, snapshot_dummy
                        );
                    }
                } else {
                    // Am Ende der Historie angekommen, aber Frontier-Check war falsch?
                    // Sollte nicht passieren, aber zur Sicherheit pausieren.
                    sim_params.paused = true;
                }
            }
        }
        // Wenn ein `execute_single_step_request` kam und wir NICHT am Frontier waren,
        // dann wollen wir auch nur einen Schritt in der Historie machen und pausieren.
        // Die obige Logik für `!sim_params.paused` deckt den "Play"-Modus ab.
        if let Some(current_idx) = history.current_display_snapshot_index {
            let next_idx = current_idx + 1;
            // Extract data from snapshot to drop the immutable borrow of `history`
            // before attempting a mutable borrow later.
            let snapshot_data = history
                .get_snapshot_by_index(next_idx)
                .map(|s| (s.time_ma, s.dummy_value));

            if let Some((snapshot_time, snapshot_dummy)) = snapshot_data {
                // Store values and update state after dropping immutable borrow
                sim_params.display_time_ma = snapshot_time;
                history.current_display_snapshot_index = Some(next_idx);
                // TODO: WELTZUSTAND AUF `next_snapshot` DATEN LADEN (verwende snapshot_time, snapshot_dummy)
                info!(
                    "Single Step through history. Displaying snapshot at: {:.3} Ma. Dummy: {:.2}",
                    snapshot_time, snapshot_dummy
                );
                sim_params.paused = true; // Nach einem Step in History pausieren
            }
        }
    }
}

pub fn handle_timeline_interaction_system(
    mut sim_params: ResMut<SimulationParameters>,
    mut history: ResMut<TickHistory>,
) {
    let current_snapshot_time = history
        .get_current_display_snapshot()
        .map_or(-1.0, |s| s.time_ma);

    if (sim_params.display_time_ma - current_snapshot_time).abs() > 0.00001 {
        if let Some(closest_idx) = history.find_closest_snapshot_index(sim_params.display_time_ma) {
            if history.current_display_snapshot_index != Some(closest_idx) {
                if let Some(snapshot_to_display) = history.get_snapshot_by_index(closest_idx) {
                    let snapshot_time = snapshot_to_display.time_ma;
                    let snapshot_dummy = snapshot_to_display.dummy_value;
                    // TODO: WELTZUSTAND AUF DIESEN SNAPSHOT ZURÜCKSETZEN/ANWENDEN
                    info!(
                        "Timeline scrub: Restoring state to snapshot at: {:.3} Ma. Dummy: {:.2}",
                        snapshot_time, snapshot_dummy
                    );
                    history.current_display_snapshot_index = Some(closest_idx);
                    sim_params.display_time_ma = snapshot_time;
                    sim_params.paused = true;
                }
            } else {
                if let Some(snapshot_to_display) = history.get_snapshot_by_index(closest_idx) {
                    sim_params.display_time_ma = snapshot_to_display.time_ma;
                }
            }
        } else if history.is_empty() {
            sim_params.display_time_ma = sim_params.advanced_simulation_time_ma;
        }
    }
}
