/// Represents an event that is triggered when the seed for random number generation changes.
///
/// This event carries the new seed value, allowing other parts of the application
/// to react to the change, for example, by re-initializing their own random number
/// generators or re-calculating procedurally generated content.
use bevy::prelude::*;

#[derive(Event, Debug, Clone)]
pub struct SeedChangedEvent {
    /// The new seed value.
    pub new_seed: u32,
}
impl SeedChangedEvent {
    /// Creates a new `SeedChangedEvent` with the specified seed value.
    ///
    /// # Arguments
    ///
    /// * `new_seed` - The new seed value to be set.
    pub fn new(new_seed: u32) -> Self {
        Self { new_seed }
    }
}
