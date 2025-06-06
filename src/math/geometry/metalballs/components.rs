// src/math/geometry/metalballs/components.rs

/// Gemeinsame Konfiguration für Grid-basierte Berechnungen wie Metaballs oder Marching Squares.
#[derive(Debug, Clone)]
pub struct GridConfig {
    /// Breite des Grids in Zellen.
    pub width: usize,
    /// Höhe des Grids in Zellen.
    pub height: usize,
    /// Größe einer einzelnen Zelle in Weltkoordinaten.
    pub cell_size: f32,
    /// Schwellenwert, der für die Iso-Kontur-Extraktion (z.B. Marching Squares) verwendet wird.
    pub threshold: f32,
}

impl Default for GridConfig {
    fn default() -> Self {
        Self {
            width: 100,
            height: 100,
            cell_size: 1.0, // Geändert von 8.0 für eine eventuell intuitivere Skalierung
            threshold: 0.5, // Üblicher Wert für normalisierte Felder [0,1] oder symmetrische Felder
        }
    }
}

impl GridConfig {
    /// Erstellt eine neue Grid-Konfiguration.
    pub fn new(width: usize, height: usize, cell_size: f32, threshold: f32) -> Self {
        Self {
            width,
            height,
            cell_size,
            threshold,
        }
    }
}
