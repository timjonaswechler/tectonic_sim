//! Provides a Bevy `Plugin` for integrating seed-based random number generation.
//!
//! This plugin initializes the `SeedResource` and sets up an event handler
//! for `SeedChangedEvent` to update the resource when the seed changes.

use super::{events::SeedChangedEvent, resource::SeedResource};
use bevy::prelude::*;

/// A Bevy `Plugin` responsible for managing the global `SeedResource`.
///
/// This plugin performs the following actions:
/// - Initializes the `SeedResource` if it doesn't already exist.
/// - Registers the `SeedChangedEvent`.
/// - Adds a system (`seed_changed_event_handler`) to listen for `SeedChangedEvent`
///   and update the `SeedResource` accordingly.
pub struct SeedPlugin;

impl Plugin for SeedPlugin {
    /// Builds the plugin by adding necessary resources, events, and systems to the Bevy `App`.
    fn build(&self, app: &mut App) {
        app.init_resource::<SeedResource>()
            .add_event::<SeedChangedEvent>()
            // The system to handle seed changes.
            // For Bevy 0.10+, .add_systems(Update, system) is used.
            // For older versions, .add_system(system) was used.
            .add_systems(Update, seed_changed_event_handler);
    }
}

/// Event handler system that listens for `SeedChangedEvent`.
///
/// When a `SeedChangedEvent` is detected, this system updates the `SeedResource`
/// with the new seed value provided in the event.
///
/// # Arguments
///
/// * `events`: An `EventReader` for `SeedChangedEvent` to consume incoming events.
/// * `seed_resource`: A mutable reference to the `SeedResource` to be updated.
fn seed_changed_event_handler(
    mut events: EventReader<SeedChangedEvent>,
    mut seed_resource: ResMut<SeedResource>,
) {
    for event in events.read() {
        info!(
            "SeedResource is being updated due to SeedChangedEvent: {}",
            event.new_seed
        );
        seed_resource.reset_with_new_seed(event.new_seed);
        // Additional actions could be performed here if needed after a seed change.
    }
}
