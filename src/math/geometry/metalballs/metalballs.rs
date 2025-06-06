// src/math/geometry/metalballs/metalballs.rs

use crate::math::geometry::metalballs::{
    components::GridConfig, field::MetaballField, influence::FieldInfluence,
};
use bevy::math::Vec2; // Für die Konvertierung von Zellkoordinaten zu Weltkoordinaten
use std::sync::Arc;

/// Hauptstruktur für das Metaballs-System.
/// Enthält die Konfiguration, die Einflussquellen und das resultierende Skalarfeld.
#[derive(Debug, Clone)] // Clone ist möglich, wenn FieldInfluence-Arcs cloneable sind
pub struct Metaballs {
    pub config: GridConfig,
    pub sources: Vec<Arc<dyn FieldInfluence>>,
    pub field: MetaballField,
}

impl Metaballs {
    /// Erstellt eine neue `Metaballs`-Instanz aus Builder-Daten.
    /// Diese Methode wird typischerweise vom `MetaballsBuilder` aufgerufen.
    pub fn from_builder_data(config: GridConfig, sources: Vec<Arc<dyn FieldInfluence>>) -> Self {
        // Erstelle ein leeres Feld passend zur Konfiguration.
        // Die `new` Methode von MetaballField sollte mit width/height=0 umgehen können,
        // aber der Builder stellt jetzt sicher, dass sie >= 1 sind.
        let field = MetaballField::new(config.width, config.height, config.cell_size);

        let mut metaballs_instance = Self {
            config,
            sources,
            field,
        };

        // Berechne das Feld initial, wenn Quellen vorhanden sind.
        if !metaballs_instance.sources.is_empty() {
            metaballs_instance.recompute_field();
        }

        metaballs_instance
    }

    /// Berechnet die Werte des Skalarfeldes basierend auf den aktuellen Einflussquellen neu.
    pub fn recompute_field(&mut self) {
        // Stelle sicher, dass das Feld die korrekten Dimensionen hat, falls sich die Config geändert hat
        // (obwohl der Builder das Feld initialisiert).
        if self.field.width != self.config.width
            || self.field.height != self.config.height
            || (self.field.cell_size - self.config.cell_size).abs()
                > crate::math::utils::constants::EPSILON
        {
            self.field =
                MetaballField::new(self.config.width, self.config.height, self.config.cell_size);
        }

        if self.sources.is_empty() {
            // Wenn keine Quellen, Feld mit Nullen füllen (oder was auch immer der Default ist)
            self.field =
                MetaballField::new(self.config.width, self.config.height, self.config.cell_size);
            return;
        }

        for y_idx in 0..self.config.height {
            for x_idx in 0..self.config.width {
                // Berechne die Weltkoordinaten des aktuellen Zellzentrums
                // (oder der oberen linken Ecke, je nach Konvention für `influence_at`)
                // Nehmen wir Zellzentrum für genauere Ergebnisse mit Falloffs.
                let world_x = (x_idx as f32 + 0.5) * self.config.cell_size;
                let world_y = (y_idx as f32 + 0.5) * self.config.cell_size;
                let current_point = Vec2::new(world_x, world_y);

                let total_influence = self
                    .sources
                    .iter()
                    .map(|source| source.influence_at(current_point))
                    .sum(); // Einfache Addition aller Einflüsse

                self.field.set(x_idx, y_idx, total_influence);
            }
        }
    }

    /// Fügt eine neue Einflussquelle hinzu und berechnet das Feld optional neu.
    pub fn add_source<I: FieldInfluence + 'static>(&mut self, source: I, recompute: bool) {
        self.sources.push(Arc::new(source));
        if recompute {
            self.recompute_field();
        }
    }

    /// Entfernt alle Einflussquellen.
    pub fn clear_sources(&mut self, recompute: bool) {
        self.sources.clear();
        if recompute {
            self.recompute_field(); // Wird das Feld mit Nullen füllen
        }
    }

    /// Gibt den Wert des Feldes an einer bestimmten Zellkoordinate zurück.
    pub fn get_field_value(&self, x_idx: usize, y_idx: usize) -> f32 {
        self.field.get(x_idx, y_idx)
    }
}
