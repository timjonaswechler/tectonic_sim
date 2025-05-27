use crate::physics::sim::resources::SimulationParameters; // Geänderter Pfad
use crate::physics::sim::time::resources::{SimulationSnapshot, TickHistory}; // Korrekte Pfade
use bevy::prelude::*;
use bevy_egui::{
    EguiContexts,
    egui::{self, ScrollArea, Slider, Window},
};

pub fn simulation_control_ui_system(
    mut contexts: EguiContexts,
    mut sim_params: ResMut<SimulationParameters>,
    mut history: ResMut<TickHistory>,
) {
    Window::new("Simulationssteuerung")
        .default_width(350.0)
        .show(contexts.ctx_mut(), |ui| {
            ui.heading("Globale Steuerung");
            ui.collapsing("Zeit & Geschwindigkeit", |ui| {
                ui.label(format!(
                    "Angezeigte Simulationszeit: {:.3} Ma",
                    sim_params.display_time_ma
                ));
                ui.label(format!(
                    "Fortgeschrittene Sim-Zeit: {:.3} Ma",
                    sim_params.advanced_simulation_time_ma
                ));

                if !history.is_empty() {
                    let num_snapshots = history.len();
                    let min_time = history.snapshots.front().map_or(0.0, |s| s.time_ma);
                    let max_time = history
                        .snapshots
                        .back()
                        .map_or(sim_params.advanced_simulation_time_ma, |s| s.time_ma);

                    let mut slider_time_value = sim_params.display_time_ma;
                    if ui
                        .add(
                            Slider::new(&mut slider_time_value, min_time..=max_time)
                                .text("Zeitleiste (Ma)")
                                .show_value(true)
                                .step_by(sim_params.time_step_ma.min(0.001)),
                        )
                        .changed()
                    {
                        sim_params.display_time_ma = slider_time_value;
                    }

                    ui.horizontal(|ui| {
                        if ui.button("⏮ Zum Ältesten").clicked() {
                            if let Some(snapshot) = history.snapshots.front() {
                                sim_params.display_time_ma = snapshot.time_ma;
                            }
                        }
                        if ui.button("Neuesten ⏭").clicked() {
                            if let Some(snapshot) = history.snapshots.back() {
                                sim_params.display_time_ma = snapshot.time_ma;
                            } else if history.is_empty() {
                                sim_params.display_time_ma = sim_params.advanced_simulation_time_ma;
                            }
                        }
                    });
                } else {
                    ui.label("Keine Historie aufgezeichnet.");
                }

                ui.separator();

                ui.add(
                    Slider::new(&mut sim_params.time_step_ma, 0.001..=1.0)
                        .logarithmic(true)
                        .text("Zeitschritt pro Tick (Ma)"),
                );
                // Entfernt: Slider für real_time_ms_per_tick

                ui.horizontal(|ui| {
                    if ui
                        .button(if sim_params.paused {
                            "▶ Play"
                        } else {
                            "⏸ Pause"
                        })
                        .clicked()
                    {
                        sim_params.paused = !sim_params.paused;
                    }
                    if ui.button("◀ Schritt Zurück").clicked() {
                        sim_params.execute_step_backward_request = true;
                    }
                    if ui.button("Schritt Vorwärts ▶").clicked() {
                        sim_params.execute_single_step_request = true;
                    }
                    if ui.button("↺ Reset").clicked() {
                        sim_params.advanced_simulation_time_ma = 0.0;
                        sim_params.display_time_ma = 0.0;
                        history.clear();
                        let initial_snapshot = SimulationSnapshot {
                            time_ma: 0.0,
                            dummy_value: 0.0,
                        };
                        history.add_snapshot(initial_snapshot);
                        sim_params.paused = true;
                        info!("Simulation Reset!");
                    }
                });
            });

            // Visualisierung & Debug Sektion bleibt gleich
            ui.collapsing("Visualisierung & Debug", |ui| {
                // ... (Code von vorher)
                ui.label(format!(
                    "Geologie Geschw.-Gizmo Skala: {:.2}",
                    sim_params.geology_vel_scale
                ));
                ui.add(Slider::new(&mut sim_params.geology_vel_scale, 0.1..=10.0));
                ui.label(format!(
                    "Vertex Geschw.-Gizmo Skala: {:.2}",
                    sim_params.vertex_vel_scale
                ));
                ui.add(Slider::new(&mut sim_params.vertex_vel_scale, 0.1..=10.0));
                ui.separator();
                ui.label("Grenztyp-Schwellenwerte (Faktoren der Rel.geschw.)");
                ui.add(
                    Slider::new(&mut sim_params.convergence_threshold_factor, 0.01..=0.99)
                        .text("Konvergenz Faktor"),
                );
                ui.add(
                    Slider::new(&mut sim_params.shear_threshold_factor, 0.01..=0.99)
                        .text("Scher Faktor"),
                );
                if let Some(snapshot) = history.get_current_display_snapshot() {
                    ui.label(format!(
                        "Angezeigter Snapshot Dummy Value: {:.2}",
                        snapshot.dummy_value
                    ));
                } else if history.is_empty() && sim_params.advanced_simulation_time_ma == 0.0 {
                    ui.label("Angezeigter Snapshot Dummy Value: 0.00 (Initial)");
                }
            });

            // Kamerasteuerung Info Sektion bleibt gleich
            ui.collapsing("Kamerasteuerung Info", |ui| {
                // ... (Code von vorher)
                ui.label("Rechte Maustaste + Ziehen: Orbit");
                ui.label("Mittlere Maustaste + Ziehen: Pan");
                ui.label("Mausrad: Zoom");
            });
        });
}
