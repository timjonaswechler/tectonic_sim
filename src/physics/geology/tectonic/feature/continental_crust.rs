use bevy::prelude::*;

// Marker component
#[derive(Component, Debug, Clone)]
pub struct ContinentalCrust;

// Keine eigenen Felder mehr hier, diese werden durch
// common_components::PlateId, common_components::AppearsAt,
// bevy::prelude::Name, common_components::Vertices,
// common_components::Center, common_components::Area
// abgedeckt, wenn eine Entität mit ContinentalCrust erstellt wird.
