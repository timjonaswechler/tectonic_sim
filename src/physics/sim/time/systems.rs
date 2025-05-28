use super::super::state::SimulationState;
use super::resources::*;

use crate::physics::sim::resources::SimulationParameters; // Importiere aus dem übergeordneten Sim-Modul
use bevy::prelude::*;

// pub fn simulation_driver_system(
//     mut sim_params: ResMut<SimulationParameters>,
//     mut history: ResMut<TickHistory>,
// ) {
//     if sim_params.execute_step_backward_request {
//         if let Some(current_idx) = history.current_display_snapshot_index {
//             if current_idx > 0 {
//                 let prev_idx = current_idx - 1;
//                 // Extract data from snapshot to drop the immutable borrow of `history`
//                 // before attempting a mutable borrow later.
//                 let snapshot_data = history
//                     .get_snapshot_by_index(prev_idx)
//                     .map(|s| (s.time_ma, s.dummy_value));

//                 if let Some((snapshot_time, snapshot_dummy)) = snapshot_data {
//                     sim_params.display_time_ma = snapshot_time;
//                     history.current_display_snapshot_index = Some(prev_idx);
//                     // TODO: WELTZUSTAND AUF DATEN DES VORHERIGEN SNAPSHOTS LADEN (verwende snapshot_time, snapshot_dummy)
//                     info!(
//                         "Stepped backward. Displaying snapshot at: {:.3} Ma. Dummy: {:.2}",
//                         snapshot_time, snapshot_dummy
//                     );
//                     sim_params.paused = true; // Nach einem Schritt zurück immer pausieren
//                 }
//             }
//         }
//         sim_params.execute_step_backward_request = false; // Anfrage verarbeitet
//         return; // Keine weitere Aktion in diesem Frame für Step Backward
//     }

//     // Bestimme, ob ein Vorwärts-Schritt ausgeführt werden soll
//     let should_execute_forward_step = sim_params.execute_single_step_request || !sim_params.paused;

//     if should_execute_forward_step {
//         sim_params.execute_single_step_request = false; // Anfrage verarbeitet, falls es ein Single-Step war

//         let at_simulation_frontier = history
//             .current_display_snapshot_index
//             .map_or(true, |idx| idx == history.len() - 1);

//         if at_simulation_frontier {
//             // --- NEUEN SIMULATIONSSCHRITT AUSFÜHREN ---
//             //TODO: Hier die Logik für den neuen Simulationsschritt einfügen
//             sim_params.advanced_simulation_time_ma += sim_params.time_step_ma;
//             sim_params.display_time_ma = sim_params.advanced_simulation_time_ma;

//             let new_dummy_value = (sim_params.advanced_simulation_time_ma * 10.0 % 100.0) as f32;
//             info!(
//                 "Simulating new step. Time: {:.3} Ma. Dummy: {:.2}",
//                 sim_params.advanced_simulation_time_ma, new_dummy_value
//             );

//             let new_snapshot = SimulationSnapshot {
//                 time_ma: sim_params.advanced_simulation_time_ma,
//                 dummy_value: new_dummy_value,
//             };
//             history.add_snapshot(new_snapshot);
//             // Wenn ein einzelner Schritt ausgeführt wurde (nicht im Play-Modus), pausieren wir danach
//             if !(!sim_params.paused) {
//                 // d.h. es war ein single step
//                 sim_params.paused = true;
//             }
//         } else if !sim_params.paused {
//             // Im "Play"-Modus und nicht am Frontier -> Historie abspielen
//             // --- DURCH HISTORIE VORWÄRTS "SPIELEN" ---
//             if let Some(current_idx) = history.current_display_snapshot_index {
//                 let next_idx = current_idx + 1;
//                 if next_idx < history.len() {
//                     // Extract data from snapshot to drop the immutable borrow of `history`
//                     // before attempting a mutable borrow later.
//                     let snapshot_data = history
//                         .get_snapshot_by_index(next_idx)
//                         .map(|s| (s.time_ma, s.dummy_value));

//                     if let Some((snapshot_time, snapshot_dummy)) = snapshot_data {
//                         sim_params.display_time_ma = snapshot_time;
//                         history.current_display_snapshot_index = Some(next_idx);
//                         // TODO: WELTZUSTAND AUF `next_snapshot` DATEN LADEN (verwende snapshot_time, snapshot_dummy)
//                         info!(
//                             "Playing history forward. Displaying snapshot at: {:.3} Ma. Dummy: {:.2}",
//                             snapshot_time, snapshot_dummy
//                         );
//                     }
//                 } else {
//                     // Am Ende der Historie angekommen, aber Frontier-Check war falsch?
//                     // Sollte nicht passieren, aber zur Sicherheit pausieren.
//                     sim_params.paused = true;
//                 }
//             }
//         }
//         // Wenn ein `execute_single_step_request` kam und wir NICHT am Frontier waren,
//         // dann wollen wir auch nur einen Schritt in der Historie machen und pausieren.
//         // Die obige Logik für `!sim_params.paused` deckt den "Play"-Modus ab.
//         if let Some(current_idx) = history.current_display_snapshot_index {
//             let next_idx = current_idx + 1;
//             // Extract data from snapshot to drop the immutable borrow of `history`
//             // before attempting a mutable borrow later.
//             let snapshot_data = history
//                 .get_snapshot_by_index(next_idx)
//                 .map(|s| (s.time_ma, s.dummy_value));

//             if let Some((snapshot_time, snapshot_dummy)) = snapshot_data {
//                 // Store values and update state after dropping immutable borrow
//                 sim_params.display_time_ma = snapshot_time;
//                 history.current_display_snapshot_index = Some(next_idx);
//                 // TODO: WELTZUSTAND AUF `next_snapshot` DATEN LADEN (verwende snapshot_time, snapshot_dummy)
//                 info!(
//                     "Single Step through history. Displaying snapshot at: {:.3} Ma. Dummy: {:.2}",
//                     snapshot_time, snapshot_dummy
//                 );
//                 sim_params.paused = true; // Nach einem Step in History pausieren
//             }
//         }
//     }
// }

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

pub fn handle_simulation_requests_system(
    mut sim_params: ResMut<SimulationParameters>,
    mut next_state: ResMut<NextState<SimulationState>>,
    current_state: Res<State<SimulationState>>,
    mut exec_single_step: ResMut<ExecuteSingleStepRequest>,
    mut exec_step_backward: ResMut<ExecuteStepBackwardRequest>,
    // Optional: history Res<TickHistory> wenn für Entscheidungen gebraucht
) {
    // Play/Pause aus sim_params synchronisieren mit dem State
    if sim_params.paused && *current_state == SimulationState::Running {
        next_state.set(SimulationState::Paused);
        info!("Simulation Paused (State Change).");
    } else if !sim_params.paused && *current_state == SimulationState::Paused {
        next_state.set(SimulationState::Running);
        info!("Simulation Resumed (State Change).");
    }

    // Single Step Request
    if sim_params.execute_single_step_request {
        info!("Single step forward request received.");
        exec_single_step.0 = true; // Setze das globale Flag
        if *current_state == SimulationState::Running {
            // Wenn es schon läuft, wird der nächste Frame ein Schritt sein,
            // aber wir wollen danach pausieren. Setze sim_params.paused,
            // damit der nächste Aufruf dieses Systems in Paused wechselt.
            sim_params.paused = true;
        } else if *current_state == SimulationState::Paused {
            // Ist bereits pausiert, Flag reicht.
        }
        // Im Initializing State sollte ein single step ggf. was anderes bedeuten
        sim_params.execute_single_step_request = false; // Anfrage in sim_params verarbeitet
    } else {
        // Stelle sicher, dass das Flag zurückgesetzt wird, wenn keine Anfrage mehr da ist
        // und der Schritt (hoffentlich) verarbeitet wurde.
        // Das ist knifflig - advance_simulation_step_system sollte das Flag zurücksetzen.
    }

    // Step Backward Request
    if sim_params.execute_step_backward_request {
        info!("Step backward request received.");
        exec_step_backward.0 = true; // Setze das globale Flag
        if *current_state != SimulationState::Paused {
            next_state.set(SimulationState::Paused); // Nach einem Schritt zurück immer pausieren
        }
        sim_params.execute_step_backward_request = false; // Anfrage in sim_params verarbeitet
    } else {
        // exec_step_backward.0 wird von load_state_from_history_system zurückgesetzt
    }
}

pub fn advance_simulation_step_system(
    mut sim_params: ResMut<SimulationParameters>,
    mut history: ResMut<TickHistory>,
    current_state: Res<State<SimulationState>>,
    mut exec_single_step: ResMut<ExecuteSingleStepRequest>,
    // mut commands: Commands, mut gizmos: Gizmos, // Für Polygon-Generierung etc.
) {
    // Dieses System sollte nur laufen, wenn State::Running ODER (State::Paused UND exec_single_step.0)
    // Die run_if-Bedingung stellt das sicher.

    let at_simulation_frontier = history
        .current_display_snapshot_index
        .map_or(true, |idx| idx == history.len() - 1);

    if !at_simulation_frontier {
        // Wenn wir nicht am Frontier sind, soll dieses System keinen *neuen* Schritt berechnen.
        // Das Abspielen der Historie (wenn !paused) wird entweder hier oder in einem
        // separaten System gehandhabt. Für Klarheit: separates System.
        // Wenn es ein Single-Step war, der durch die History geht, macht das load_state_from_history.
        if exec_single_step.0 { // Es war ein single step, aber nicht am Frontier
            // Lasse load_state_from_history den Schritt machen und das Flag zurücksetzen.
        }
        return;
    }

    // Nur wenn am Frontier, einen neuen Schritt berechnen.
    info!("Advancing simulation step.");
    sim_params.advanced_simulation_time_ma += sim_params.time_step_ma;
    sim_params.display_time_ma = sim_params.advanced_simulation_time_ma;

    let new_dummy_value = (sim_params.advanced_simulation_time_ma * 10.0 % 100.0) as f32;
    // info!("New dummy value: {}", new_dummy_value);

    // TODO: Hier die eigentliche Simulationslogik für einen Schritt einfügen
    // (Polygon generieren, etc.)

    let new_snapshot = SimulationSnapshot {
        time_ma: sim_params.advanced_simulation_time_ma,
        dummy_value: new_dummy_value,
    };
    history.add_snapshot(new_snapshot);

    // Wenn es ein Einzelschritt war (angefordert durch Flag), oder wenn sim_params.paused
    // (was von handle_simulation_requests gesetzt wurde, wenn ein Play->SingleStep Übergang war),
    // dann setze sim_params.paused explizit (obwohl es schon sein sollte) und resette das Flag.
    if exec_single_step.0 {
        sim_params.paused = true; // Stellt sicher, dass es nach dem Single Step pausiert ist
        exec_single_step.0 = false; // Flag verarbeitet
        info!("Single step executed, simulation paused.");
    }
    // Wenn State::Running und kein exec_single_step, läuft es einfach weiter.
}

pub fn play_history_forward_system(
    mut sim_params: ResMut<SimulationParameters>,
    mut history: ResMut<TickHistory>,
    // Dieses System läuft nur, wenn State == Running
) {
    let at_simulation_frontier = history
        .current_display_snapshot_index
        .map_or(true, |idx| idx == history.len() - 1);

    if at_simulation_frontier {
        return; // Nichts zu tun, advance_simulation_step_system ist zuständig
    }

    // Im "Play"-Modus und nicht am Frontier -> Historie abspielen
    if let Some(current_idx) = history.current_display_snapshot_index {
        let next_idx = current_idx + 1;
        if next_idx < history.len() {
            if let Some(snapshot) = history.get_snapshot_by_index(next_idx) {
                sim_params.display_time_ma = snapshot.time_ma;
                // Der tatsächliche Ladevorgang des Zustands wird von load_state_from_history_system übernommen,
                // das auf Änderungen von display_time_ma oder Flags reagiert.
                // Wir setzen hier nur den Index, damit load_state_from_history weiß, was zu tun ist.
                history.current_display_snapshot_index = Some(next_idx);
                // info!("Playing history: Displaying snapshot at {:.3} Ma", snapshot.time_ma);
            }
        } else {
            // Sollte nicht passieren, da at_simulation_frontier das abfangen sollte.
            // sim_params.paused = true; // Sicherheitshalber
        }
    }
}

pub fn should_load_from_history_condition(
    sim_params: Res<SimulationParameters>,
    history: Res<TickHistory>,
    exec_single_step: Res<ExecuteSingleStepRequest>,
    exec_step_backward: Res<ExecuteStepBackwardRequest>,
) -> bool {
    if exec_step_backward.0 {
        return true;
    }

    let at_simulation_frontier = history
        .current_display_snapshot_index
        .map_or(true, |idx| idx == history.len() - 1);

    if exec_single_step.0 && !at_simulation_frontier {
        return true;
    }

    // Prüfe, ob sich display_time_ma signifikant vom aktuell geladenen Snapshot unterscheidet
    if let Some(current_display_idx) = history.current_display_snapshot_index {
        if let Some(current_snapshot) = history.get_snapshot_by_index(current_display_idx) {
            if (sim_params.display_time_ma - current_snapshot.time_ma).abs()
                > sim_params.time_step_ma * 0.1
            {
                // Toleranz
                return true;
            }
        } else {
            // Index ist Some, aber Snapshot nicht da -> inkonsistent, laden!
            return true;
        }
    } else if !history.is_empty() {
        // Kein Index, aber History hat was -> laden!
        return true;
    }
    false
}

pub fn load_state_from_history_system(
    mut sim_params: ResMut<SimulationParameters>,
    mut history: ResMut<TickHistory>,
    mut exec_single_step: ResMut<ExecuteSingleStepRequest>,
    mut exec_step_backward: ResMut<ExecuteStepBackwardRequest>,
    // mut commands: Commands, // Um Entitäten neu zu erstellen/modifizieren
    // Query für Simulationsentitäten, die aktualisiert werden müssen
) {
    let mut target_snapshot_index: Option<usize> = None;
    let mut pause_after_load = false;

    if exec_step_backward.0 {
        if let Some(current_idx) = history.current_display_snapshot_index {
            if current_idx > 0 {
                target_snapshot_index = Some(current_idx - 1);
                pause_after_load = true;
            } else {
                // Schon am ältesten Snapshot
                exec_step_backward.0 = false; // Kein Schritt möglich
                return;
            }
        } else if !history.is_empty() {
            // Kein Index, aber History da -> zum ersten
            target_snapshot_index = Some(0);
            pause_after_load = true;
        }
        exec_step_backward.0 = false; // Anfrage wird hier verarbeitet
    } else if exec_single_step.0 {
        let at_simulation_frontier = history
            .current_display_snapshot_index
            .map_or(true, |idx| idx == history.len() - 1);
        if !at_simulation_frontier {
            // Single Step durch die Historie
            if let Some(current_idx) = history.current_display_snapshot_index {
                if current_idx + 1 < history.len() {
                    target_snapshot_index = Some(current_idx + 1);
                    pause_after_load = true;
                } else {
                    // Sollte nicht passieren, wenn nicht am Frontier
                    exec_single_step.0 = false;
                    return;
                }
            } else if !history.is_empty() {
                // Kein Index, aber History da -> zum ersten
                target_snapshot_index = Some(0);
                pause_after_load = true;
            }
            exec_single_step.0 = false; // Anfrage wird hier verarbeitet
        } else {
            // Single Step am Frontier wird von advance_simulation_step_system gehandhabt
            // Dieses Flag wird dort zurückgesetzt. Hier nichts tun.
            return;
        }
    } else {
        // Getriggert durch Timeline-Scrubbing oder Play History
        target_snapshot_index = history.find_closest_snapshot_index(sim_params.display_time_ma);
        // Beim reinen Scrubbing nicht unbedingt pausieren, es sei denn, sim_params.paused ist schon true
        pause_after_load = sim_params.paused;
    }

    if let Some(load_idx) = target_snapshot_index {
        if history.current_display_snapshot_index == Some(load_idx) && !pause_after_load {
            // Der Ziel-Snapshot ist bereits geladen und wir müssen nicht extra pausieren.
            // Aber display_time_ma ggf. exakt angleichen.
            if let Some(snap) = history.get_snapshot_by_index(load_idx) {
                sim_params.display_time_ma = snap.time_ma;
            }
            return;
        }

        if let Some(snapshot_to_load) = history.get_snapshot_by_index(load_idx) {
            info!(
                "Loading state from history. Target Idx: {}, Time: {:.3} Ma, Dummy: {:.2}",
                load_idx, snapshot_to_load.time_ma, snapshot_to_load.dummy_value
            );
            // --- EIGENTLICHE LOGIK ZUM LADEN DES WELTZUSTANDS ---
            // Dies ist der komplexeste Teil:
            // 1. Bestehende Simulationsentitäten (Platten etc.) despawnen oder aktualisieren.
            // 2. Neue Entitäten basierend auf dem Snapshot-Daten erstellen.
            // Beispiel:
            // For entity_in_world in query_sim_entities.iter_mut() { commands.entity(entity_in_world.id()).despawn_recursive(); }
            // For saved_plate_state in snapshot_to_load.plate_states { spawn_plate_from_state(&mut commands, saved_plate_state); }
            // Hier nur mit dem Dummy-Wert:
            // commands.insert_resource(CurrentSimData { dummy: snapshot_to_load.dummy_value });

            sim_params.display_time_ma = snapshot_to_load.time_ma;
            history.current_display_snapshot_index = Some(load_idx);
            if pause_after_load {
                sim_params.paused = true;
            }
        } else {
            warn!("Konnte Snapshot mit Index {} nicht laden.", load_idx);
        }
    } else if history.is_empty()
        && (sim_params.display_time_ma - sim_params.advanced_simulation_time_ma).abs() > 1e-6
    {
        // History ist leer, aber display_time weicht ab -> auf advanced_simulation_time setzen (sollte 0 sein)
        sim_params.display_time_ma = sim_params.advanced_simulation_time_ma;
    }
}
pub fn finish_initialization_system(
    mut sim_params: ResMut<SimulationParameters>,
    mut history: ResMut<TickHistory>,
) {
    // Set initial display time to match advanced simulation time
    sim_params.display_time_ma = sim_params.advanced_simulation_time_ma;

    // Clear history
    history.clear();
}

pub fn setup_initial_world_conditions_system(
    mut commands: Commands,
    mut sim_params: ResMut<SimulationParameters>,
) {
    // Set initial simulation parameters
    sim_params.advanced_simulation_time_ma = 0.0;
    sim_params.display_time_ma = 0.0;
    sim_params.time_step_ma = 0.1; // Default time step
    sim_params.paused = true; // Start paused
    sim_params.execute_single_step_request = false;
    sim_params.execute_step_backward_request = false;

    info!("Initial world conditions set up.");
}
