use bevy::prelude::*;

#[derive(States, Debug, Clone, Copy, Eq, PartialEq, Hash, Default)]
pub enum SimulationState {
    #[default] // Startzustand
    Initializing,
    Running,
    Paused,
    // Optional:
    // SteppingForward, // Für den Moment, in dem ein Einzelschritt ausgeführt wird
    // SteppingBackward,
    // Resetting,
}
