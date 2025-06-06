// src/math/geometry/metalballs/field.rs

/// Repräsentiert ein zweidimensionales Skalarfeld, das die summierten Einflüsse
/// von Metaball-Quellen an jedem Gitterpunkt speichert.
/// Die Daten werden zeilenweise (row-major) gespeichert.
#[derive(Debug, Clone, PartialEq)]
pub struct MetaballField {
    data: Vec<f32>,     // Die Skalarwerte des Feldes.
    pub width: usize,   // Breite des Feldes in Zellen.
    pub height: usize,  // Höhe des Feldes in Zellen.
    pub cell_size: f32, // Größe einer Zelle in Weltkoordinaten.
}

impl MetaballField {
    /// Erstellt ein neues, mit Nullen initialisiertes MetaballField.
    pub fn new(width: usize, height: usize, cell_size: f32) -> Self {
        if width == 0 || height == 0 {
            // Vermeide ein leeres `data` Array, wenn width oder height 0 ist,
            // was zu Panics bei `idx` führen könnte, wenn nicht sorgfältig geprüft.
            // Besser ist es, ein Feld der Größe 0x0 oder 1x0 etc. zu erlauben,
            // aber die Operationen müssen damit umgehen können.
            // Fürs Erste: leeres Feld, wenn eine Dimension 0 ist.
            Self {
                data: Vec::new(),
                width: 0,
                height: 0,
                cell_size,
            }
        } else {
            Self {
                data: vec![0.0; width * height],
                width,
                height,
                cell_size,
            }
        }
    }

    /// Berechnet den linearen Index für gegebene 2D-Koordinaten (x, y).
    #[inline]
    fn idx(&self, x: usize, y: usize) -> Option<usize> {
        if x < self.width && y < self.height {
            Some(y * self.width + x)
        } else {
            None // Außerhalb der Grenzen
        }
    }

    /// Gibt den Skalarwert an der Zelle (x, y) zurück.
    /// Gibt `0.0` zurück, wenn die Koordinaten außerhalb der Grenzen liegen.
    pub fn get(&self, x: usize, y: usize) -> f32 {
        match self.idx(x, y) {
            Some(index) => self.data.get(index).copied().unwrap_or(0.0), // unwrap_or für den Fall, dass data leer ist
            None => 0.0,                                                 // Außerhalb der Grenzen
        }
    }

    /// Setzt den Skalarwert an der Zelle (x, y).
    /// Macht nichts, wenn die Koordinaten außerhalb der Grenzen liegen.
    pub fn set(&mut self, x: usize, y: usize, value: f32) {
        if let Some(index) = self.idx(x, y) {
            if let Some(val_ref) = self.data.get_mut(index) {
                *val_ref = value;
            }
        }
    }

    /// Gibt die Daten des Feldes als Slice zurück.
    pub fn data(&self) -> &[f32] {
        &self.data
    }
}
