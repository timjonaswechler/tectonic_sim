use bevy::prelude::*;

// Hinweis: einen flowline geht vom continental rift bis zum continetal crust vertex
// alle Flowlines die zu einem Continent gehören kreuzen sich nicht
// generell dürfen flowlines nicht mit anderen Kontinttal platten in kontakt kommen.

#[derive(Component, Debug, Clone)]
pub struct Flowline; // Umbenannt von FlowlineMarker für Konsistenz

// Die eigentlichen Daten für die Flowline könnten in common_components::FlowlinePath
// oder direkt common_components::Vertices gespeichert werden, je nach Bedarf.
