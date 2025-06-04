use bevy::prelude::*;

// Gemeinsame Datenkomponenten

#[derive(Component, Debug, Clone, Copy)]
pub struct PlateId(pub u32); // Eindeutige ID der tektonischen Platte, zu der das Feature gehört

#[derive(Component, Debug, Clone, Copy)]
pub struct AppearsAt(pub f64); // Zeitpunkt des Erscheinens in Millionen Jahren

#[derive(Component, Debug, Clone)]
pub struct Vertices(pub Vec<Vec3>); // Liste von Eckpunkten, die das Feature definieren (Polygon, Polylinie)

#[derive(Component, Debug, Clone, Copy)]
pub struct Center(pub Vec3); // Zentraler Punkt des Features (kann berechnet oder explizit gesetzt werden)

#[derive(Component, Debug, Clone, Copy)]
pub struct Area(pub f64); // Fläche des Features (relevant für polygonale Features)

// Beispiel für eine spezifischere Geometriekomponente, falls benötigt
#[derive(Component, Debug, Clone)]
pub struct FlowlinePath {
    pub start_rift_entity: Option<Entity>, // Referenz zur Rift-Entität (optional)
    pub end_crust_vertex_index: Option<usize>, // Index des Vertex im ContinentalCrust (optional)
    pub path_vertices: Vec<Vec3>,
}
