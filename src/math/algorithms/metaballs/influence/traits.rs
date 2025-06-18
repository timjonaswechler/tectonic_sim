// src/math/geometry/metalballs/influence.rs

use super::composite::{CombinationMode, CompositeInfluence};
use crate::math::types::*;
use bevy::prelude::*;
use std::sync::Arc;

/// Trait für alle Objekte, die ein Skalarfeld beeinflussen können
pub trait FieldInfluence: Send + Sync + std::fmt::Debug {
    /// Berechnet den Einfluss an einem gegebenen Punkt
    fn influence_at(&self, point: Vec2) -> f32;

    /// Optionale Bounding Box für Optimierungen
    fn bounding_box(&self) -> Option<Bounds2D> {
        None
    }

    /// Prüft ob ein Punkt möglicherweise beeinflusst wird (grobe Schätzung)
    fn might_influence(&self, point: Vec2) -> bool {
        if let Some(bounds) = self.bounding_box() {
            bounds.contains_point(point)
        } else {
            true // Konservative Schätzung
        }
    }

    /// Gibt die maximale Reichweite des Einflusses zurück
    fn max_influence_distance(&self) -> Option<f32> {
        None
    }

    /// Gibt den Namen/Typ der Influence zurück für Debugging
    fn influence_type(&self) -> &'static str {
        "Unknown"
    }
}

/// Hilfsfunktionen für Field Influence
pub struct InfluenceUtils;

impl InfluenceUtils {
    /// Erstellt einen optimierten Composite-Einfluss basierend auf räumlicher Nähe
    pub fn create_spatial_composite(
        influences: Vec<Arc<dyn FieldInfluence>>,
        mode: CombinationMode,
    ) -> CompositeInfluence {
        // TODO: Hier könnte eine räumliche Datenstruktur (Quadtree) implementiert werden
        // für bessere Performance bei vielen Einflüssen

        let mut composite = CompositeInfluence::new(mode);
        for influence in influences {
            composite.influences.push(influence);
        }
        composite
    }

    /// Optimiert eine Liste von Einflüssen durch Removal von überlappenden oder schwachen Einflüssen
    pub fn optimize_influences(
        influences: Vec<Arc<dyn FieldInfluence>>,
        min_strength_threshold: f32,
    ) -> Vec<Arc<dyn FieldInfluence>> {
        influences
            .into_iter()
            .filter(|influence| {
                // Einfache Heuristik: teste Einfluss am eigenen Zentrum
                if let Some(bounds) = influence.bounding_box() {
                    let center = bounds.center();
                    let center_vec2 = Vec2::new(center.x, center.y);
                    influence.influence_at(center_vec2) >= min_strength_threshold
                } else {
                    true // Behalte Einflüsse ohne Bounds
                }
            })
            .collect()
    }
}
