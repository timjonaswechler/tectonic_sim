use super::{events::SeedChangedEvent, resource::SeedResource};
use bevy::prelude::*;

pub struct SeedPlugin;

impl Plugin for SeedPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<SeedResource>()
            .add_event::<SeedChangedEvent>()
            // Für Bevy 0.10+ wäre es .add_systems(Update, seed_changed_event_handler);
            // Für ältere Versionen ist .add_system(seed_changed_event_handler) ok.
            .add_systems(Update, seed_changed_event_handler);
    }
}

fn seed_changed_event_handler(
    mut events: EventReader<SeedChangedEvent>,
    mut seed_resource: ResMut<SeedResource>, // Hinzugefügt
) {
    for event in events.read() {
        info!(
            "SeedResource wird aktualisiert wegen SeedChangedEvent: {}",
            event.new_seed
        );
        seed_resource.reset_with_new_seed(event.new_seed); // Aufruf der Methode
        // Hier ggf. weitere Aktionen
    }
}
