use bevy::prelude::*;

#[derive(Event, Debug, Clone)]
pub struct SeedChangedEvent {
    pub new_seed: u32,
}
