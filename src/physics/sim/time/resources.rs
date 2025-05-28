use bevy::prelude::*;
use std::collections::VecDeque;

// ExecuteSingleStepFlag und ExecuteStepBackwardFlag müssen als Ressourcen definiert werden
#[derive(Resource, PartialEq, Debug, Clone, Copy)]
pub struct ExecuteSingleStepFlag(pub bool); // true, wenn ein Einzelschritt ausgeführt werden soll

impl Default for ExecuteSingleStepFlag {
    fn default() -> Self {
        Self(false)
    }
}

#[derive(Resource, PartialEq, Debug, Clone, Copy)]
pub struct ExecuteStepBackwardFlag(pub bool); // true, wenn ein Schritt zurück ausgeführt werden soll

impl Default for ExecuteStepBackwardFlag {
    fn default() -> Self {
        Self(false)
    }
}
// ExecuteSingleStepFlag und ExecuteStepBackwardFlag müssen als Ressourcen definiert werden
#[derive(Resource, Default, Debug, Clone, Copy, PartialEq)]
pub struct ExecuteSingleStepRequest(pub bool); // Ersetzt das Flag in SimParams
#[derive(Resource, Default, Debug, Clone, Copy, PartialEq)]
pub struct ExecuteStepBackwardRequest(pub bool); // Ersetzt das Flag in SimParams

// SimulationSnapshot definiert, was pro Zeitschritt gespeichert wird.
#[derive(Clone, Debug, Default)]
pub struct SimulationSnapshot {
    pub time_ma: f64,
    pub dummy_value: f32, // Wird später durch echte Zustandsdaten der Simulation ersetzt
                          // z.B. plate_states: Vec<PlateSaveState>,
                          //      world_temperature_grid: Vec<f32>, etc.
}

// TickHistory verwaltet die Sammlung von Snapshots.
#[derive(Resource, Debug)]
pub struct TickHistory {
    pub snapshots: VecDeque<SimulationSnapshot>, // Öffentlich für leichteren Zugriff, wenn nötig
    pub max_size: usize,
    /// Index des aktuell *angezeigten* Snapshots in der `snapshots` VecDeque.
    pub current_display_snapshot_index: Option<usize>,
}

impl TickHistory {
    pub fn new(max_size: usize) -> Self {
        Self {
            snapshots: VecDeque::with_capacity(max_size),
            max_size,
            current_display_snapshot_index: None,
        }
    }

    pub fn add_snapshot(&mut self, snapshot: SimulationSnapshot) {
        if self.snapshots.len() == self.max_size {
            self.snapshots.pop_front();
            // Index anpassen, wenn der älteste Snapshot entfernt wird und der Index davon betroffen war
            if let Some(idx) = self.current_display_snapshot_index.as_mut() {
                if *idx > 0 {
                    *idx -= 1;
                } else {
                    // Wenn der erste (index 0) entfernt wurde und angezeigt wurde
                    self.current_display_snapshot_index = if self.snapshots.is_empty() {
                        None
                    } else {
                        Some(0)
                    };
                }
            }
        }
        self.snapshots.push_back(snapshot);
        // Nach dem Hinzufügen eines neuen Snapshots wird dieser zum aktuell angezeigten
        self.current_display_snapshot_index = Some(self.snapshots.len() - 1);
    }

    pub fn get_snapshot_by_index(&self, index: usize) -> Option<&SimulationSnapshot> {
        self.snapshots.get(index)
    }

    pub fn get_current_display_snapshot(&self) -> Option<&SimulationSnapshot> {
        self.current_display_snapshot_index
            .and_then(|idx| self.snapshots.get(idx))
    }

    pub fn len(&self) -> usize {
        self.snapshots.len()
    }

    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }

    pub fn clear(&mut self) {
        self.snapshots.clear();
        self.current_display_snapshot_index = None;
    }

    /// Findet den Index des Snapshots, dessen `time_ma` am nächsten zu `target_time_ma` liegt.
    pub fn find_closest_snapshot_index(&self, target_time_ma: f64) -> Option<usize> {
        if self.snapshots.is_empty() {
            return None;
        }
        self.snapshots
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                (a.time_ma - target_time_ma)
                    .abs()
                    .partial_cmp(&(b.time_ma - target_time_ma).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(index, _)| index)
    }
}
