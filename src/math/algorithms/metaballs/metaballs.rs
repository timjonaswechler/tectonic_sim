// src/math/algorithms/metaballs/metaballs.rs

// components::GridConfig wird nicht mehr direkt importiert
use crate::math::algorithms::metaballs::{field::MetaballField, influence::FieldInfluence};
use bevy::math::Vec2;
use std::sync::Arc;

/// Hauptstruktur für das Feld-Summierungs-System.
#[derive(Debug, Clone)]
pub struct Metaballs {
    // Felder von GridConfig (außer ggf. Threshold, wenn der nur zur Weitergabe an MarchingSquares dient)
    pub width: usize,   // Grid-Breite in Zellen
    pub height: usize,  // Grid-Höhe in Zellen
    pub cell_size: f32, // Größe einer Zelle in Weltkoordinaten
    pub threshold: f32, // Schwellenwert für die Konturextraktion
    pub sources: Vec<Arc<dyn FieldInfluence>>,
    pub field: MetaballField, // Das eigentliche Feld
}

impl Metaballs {
    // Angepasste `from_builder_data` oder `new` Methode
    pub fn from_builder_data(
        width: usize,
        height: usize,
        cell_size: f32,
        threshold: f32,
        sources: Vec<Arc<dyn FieldInfluence>>,
    ) -> Self {
        let field = MetaballField::new(width, height, cell_size);

        let mut metaballs_instance = Self {
            width,
            height,
            cell_size,
            threshold,
            sources,
            field,
        };

        if !metaballs_instance.sources.is_empty() {
            metaballs_instance.recompute_field();
        }

        metaballs_instance
    }

    pub fn recompute_field(&mut self) {
        // Überprüfen, ob sich die Feld-Dimensionen seit der letzten Berechnung geändert haben
        // Dies ist notwendig, wenn der Benutzer die `width`, `height`, `cell_size` der `Metaballs`-Instanz
        // nach der Erstellung direkt modifizieren könnte (was bei `pub` Feldern möglich ist).
        // Wenn diese Felder `pub(super)` oder privat mit Settern wären, könnte diese Prüfung
        // entfallen oder anders gehandhabt werden.
        if self.field.width != self.width
            || self.field.height != self.height
            || (self.field.cell_size - self.cell_size).abs()
                > crate::math::utils::constants::EPSILON
        {
            self.field = MetaballField::new(self.width, self.height, self.cell_size);
        }

        if self.sources.is_empty() {
            // Feld ggf. leeren (oder sicherstellen, dass es mit Nullen initialisiert ist)
            self.field = MetaballField::new(self.width, self.height, self.cell_size);
            return;
        }

        for y_idx in 0..self.height {
            for x_idx in 0..self.width {
                let world_x = (x_idx as f32 + 0.5) * self.cell_size;
                let world_y = (y_idx as f32 + 0.5) * self.cell_size;
                let current_point = Vec2::new(world_x, world_y);

                let total_influence: f32 = self // Typ-Annotation für Klarheit
                    .sources
                    .iter()
                    .map(|source| source.influence_at(current_point))
                    .sum();

                self.field.set(x_idx, y_idx, total_influence);
            }
        }
    }

    pub fn add_source<I: FieldInfluence + 'static>(&mut self, source: I, recompute: bool) {
        self.sources.push(Arc::new(source));
        if recompute {
            self.recompute_field();
        }
    }

    pub fn clear_sources(&mut self, recompute: bool) {
        self.sources.clear();
        if recompute {
            self.recompute_field();
        }
    }

    pub fn get_field_value(&self, x_idx: usize, y_idx: usize) -> f32 {
        self.field.get(x_idx, y_idx)
    }

    // Methode, um den Threshold an Marching Squares weiterzugeben
    pub fn get_threshold(&self) -> f32 {
        self.threshold
    }
}
