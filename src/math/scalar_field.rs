// src/math/scalar_field.rs

use bevy::math::Vec2;

/// Trait für ein zweidimensionales Skalarfeld.
/// Ermöglicht es Algorithmen wie Marching Squares, auf verschiedenen
/// Feldimplementierungen zu operieren.
pub trait ScalarField2D {
    /// Gibt die Breite des Feldes in Zellen zurück.
    fn width(&self) -> usize;

    /// Gibt die Höhe des Feldes in Zellen zurück.
    fn height(&self) -> usize;

    /// Gibt die Größe einer einzelnen Zelle in Weltkoordinaten zurück.
    fn cell_size(&self) -> f32;

    /// Gibt den Skalarwert an der Zelle (x_idx, y_idx) zurück.
    /// Indizes sind zellenbasiert (0-basiert).
    /// Sollte 0.0 oder einen passenden Default-Wert zurückgeben, wenn die
    /// Koordinaten außerhalb der Grenzen liegen oder der Index ungültig ist.
    fn get_value(&self, x_idx: usize, y_idx: usize) -> f32;

    /// Konvertiert eine Gitterzell-Koordinate (Index) in Weltkoordinaten (linke obere Ecke der Zelle).
    /// Standardimplementierung.
    fn cell_to_world(&self, x_idx: usize, y_idx: usize) -> Vec2 {
        Vec2::new(
            x_idx as f32 * self.cell_size(),
            y_idx as f32 * self.cell_size(),
        )
    }

    // Optional: Methoden zum Setzen von Werten oder Zugriff auf Rohdaten, falls nötig,
    // aber für Marching Squares `get_value` und Dimensionen sind primär.
}

// Implementierung für unser existierendes `MetaballField`
// Dies wird in der MetaballField-Datei selbst geschehen, aber hier als Idee.
/*
use crate::math::geometry::metaballs::field::MetaballField; // Pfad anpassen

impl ScalarField2D for MetaballField {
    fn width(&self) -> usize {
        self.width
    }
    fn height(&self) -> usize {
        self.height
    }
    fn cell_size(&self) -> f32 {
        self.cell_size
    }
    fn get_value(&self, x_idx: usize, y_idx: usize) -> f32 {
        self.get(x_idx, y_idx) // Verwendet die bestehende `get`-Methode
    }
}
*/
