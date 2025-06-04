use super::{events::SeedChangedEvent, resource::SeedResource};
use bevy::prelude::*;

pub struct SeedPlugin;

impl Plugin for SeedPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<SeedResource>()
            .add_event::<SeedChangedEvent>()
            .add_system(seed_changed_event_handler);
    }
}

fn seed_changed_event_handler(mut events: EventReader<SeedChangedEvent>) {
    for event in events.iter() {
        info!("Seed geändert: {}", event.new_seed);
        // Hier ggf. weitere Aktionen
    }
}
