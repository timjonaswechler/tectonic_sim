// src/math/algorithms/metaballs/field.rs
use crate::math::scalar_field::ScalarField2D; // Annahme: dieser Trait existiert noch

/// Repräsentiert ein zweidimensionales Skalarfeld, das die summierten Einflüsse
/// von Einflussquellen an jedem Gitterpunkt speichert.
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

    #[inline]
    fn idx(&self, x: usize, y: usize) -> Option<usize> {
        if x < self.width && y < self.height {
            Some(y * self.width + x)
        } else {
            None
        }
    }

    pub fn get(&self, x: usize, y: usize) -> f32 {
        match self.idx(x, y) {
            Some(index) => self.data.get(index).copied().unwrap_or(0.0),
            None => 0.0,
        }
    }

    pub fn set(&mut self, x: usize, y: usize, value: f32) {
        if let Some(index) = self.idx(x, y) {
            if let Some(val_ref) = self.data.get_mut(index) {
                *val_ref = value;
            }
        }
    }

    pub fn data(&self) -> &[f32] {
        &self.data
    }
}

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
        self.get(x_idx, y_idx)
    }
}
