// src/physics/geology/tectonic/components.rs
use bevy::prelude::*;
use geo::Polygon; // Wir verwenden f32 für Konsistenz mit Bevy

#[derive(Component, Debug, Clone)] // Clone ist nützlich für Snapshots
pub struct CratonComponent {
    pub id: u32,
    pub name: String, // Z.B. "Craton A"

    /// Die ursprüngliche 2D-Form des Cratons, wie von `spade` generiert (planar).
    /// Die Koordinaten sind relativ zum 2D-Centroid dieser Form.
    pub shape_2d: Polygon<f32>,

    /// Der 3D-Punkt auf der Kugeloberfläche, an dem das Zentrum des Cratons platziert wird.
    pub center_on_sphere: Vec3,

    /// Die auf die Kugel projizierten 3D-Eckpunkte der Craton-Form.
    /// Die Reihenfolge der Vertices definiert das Polygon auf der Kugel.
    pub vertices_on_sphere: Vec<Vec3>,
    // Später könnten hier noch Eigenschaften wie "Starrheit" oder "Alter" hinzukommen.
}
