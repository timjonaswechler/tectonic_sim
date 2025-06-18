use super::FieldInfluence;
use crate::math::types::Bounds2D; // Unser Bounds2D
use bevy::math::Vec2;
use std::sync::Arc; // Für Arc<dyn FieldInfluence>

/// Definiert, wie die Einflüsse mehrerer Quellen kombiniert werden.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CombinationMode {
    Add,      // Addiert alle Einflüsse.
    Multiply, // Multipliziert alle Einflüsse.
    Maximum,  // Nimmt den maximalen Einflusswert.
    Minimum,  // Nimmt den minimalen Einflusswert.
    Average,  // Berechnet den Durchschnitt der Einflüsse.
}

/// Kombiniert die Einflüsse mehrerer `FieldInfluence`-Quellen.
#[derive(Debug)] // Manuelles Debug, da Vec<Arc<dyn ...>> nicht direkt ableitbar ist
pub struct CompositeInfluence {
    pub influences: Vec<Arc<dyn FieldInfluence>>,
    pub combination_mode: CombinationMode,
    cached_bounds: Option<Bounds2D>,
    cached_max_influence_distance: Option<f32>,
}

impl CompositeInfluence {
    /// Erstellt eine neue `CompositeInfluence`-Instanz.
    pub fn new(mode: CombinationMode) -> Self {
        Self {
            influences: Vec::new(),
            combination_mode: mode,
            cached_bounds: None,
            cached_max_influence_distance: None,
        }
    }

    /// Fügt eine neue Einflussquelle hinzu (als `Arc`).
    pub fn add_influence_arc(&mut self, influence_arc: Arc<dyn FieldInfluence>) {
        self.influences.push(influence_arc);
        self.invalidate_caches(); // Caches müssen neu berechnet werden
    }

    /// Fügt eine neue Einflussquelle hinzu (nimmt Ownership und erstellt intern `Arc`).
    pub fn add_influence<I: FieldInfluence + 'static>(&mut self, influence: I) {
        self.add_influence_arc(Arc::new(influence));
    }

    /// Fluent-Interface zum Hinzufügen einer Einflussquelle.
    pub fn with_influence<I: FieldInfluence + 'static>(mut self, influence: I) -> Self {
        self.add_influence(influence);
        self
    }

    /// Fluent-Interface zum Hinzufügen eines `Arc`-basierten Einflusses.
    pub fn with_influence_arc(mut self, influence_arc: Arc<dyn FieldInfluence>) -> Self {
        self.add_influence_arc(influence_arc);
        self
    }

    fn invalidate_caches(&mut self) {
        self.cached_bounds = None;
        self.cached_max_influence_distance = None;
    }

    fn calculate_and_cache_bounds(&mut self) {
        if self.cached_bounds.is_some() {
            return;
        }

        let mut combined_bounds: Option<Bounds2D> = None;
        for influence in &self.influences {
            if let Some(bounds) = influence.bounding_box() {
                combined_bounds = match combined_bounds {
                    Some(existing) => Some(existing.union(&bounds)),
                    None => Some(bounds),
                };
            } else {
                // Wenn ein Einfluss keine Bounds hat, kann der Composite auch keine sinnvollen haben
                combined_bounds = None;
                break;
            }
        }
        self.cached_bounds = combined_bounds;
    }

    fn calculate_and_cache_max_distance(&mut self) {
        if self.cached_max_influence_distance.is_some() {
            return;
        }

        self.cached_max_influence_distance = self
            .influences
            .iter()
            .filter_map(|i| i.max_influence_distance())
            .fold(None, |acc, dist| {
                Some(match acc {
                    None => dist,
                    Some(current) => current.max(dist), // Oder eine andere Logik je nach CombinationMode
                })
            });
    }
}

impl FieldInfluence for CompositeInfluence {
    fn influence_at(&self, point: Vec2) -> f32 {
        if self.influences.is_empty() {
            return 0.0;
        }

        // Optional: Früher Ausstieg, wenn Punkt außerhalb der kombinierten Bounding Box.
        // Benötigt aber &mut self für den Cache oder eine &self-Version von get_bounds.
        // if let Some(bounds) = self.bounding_box() { // bounding_box() nimmt &self
        //     if !bounds.contains_point(point) {
        //         return 0.0;
        //     }
        // }

        let values: Vec<f32> = self
            .influences
            .iter()
            .map(|influence| influence.influence_at(point))
            .collect();

        match self.combination_mode {
            CombinationMode::Add => values.iter().sum(),
            CombinationMode::Multiply => values.iter().product(), // Vorsicht bei vielen Quellen, kann schnell 0 oder sehr groß werden
            CombinationMode::Maximum => values.iter().fold(f32::NEG_INFINITY, |acc, &x| acc.max(x)),
            CombinationMode::Minimum => values.iter().fold(f32::INFINITY, |acc, &x| acc.min(x)),
            CombinationMode::Average => {
                if values.is_empty() {
                    0.0
                } else {
                    values.iter().sum::<f32>() / values.len() as f32
                }
            }
        }
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        // Da bounding_box() &self nimmt, können wir den Cache nicht hier füllen.
        // Der Aufrufer müsste `calculate_and_cache_bounds(&mut self)` vorher aufrufen,
        // oder wir akzeptieren, dass es hier immer neu berechnet wird, wenn nicht schon gecacht.
        // Für eine einfache &self-Implementierung ohne direkten Cache-Zugriff:
        if let Some(cached) = self.cached_bounds {
            return Some(cached);
        }
        // Neuberechnung (ohne Caching hier, da &self)
        let mut combined_bounds: Option<Bounds2D> = None;
        for influence in &self.influences {
            if let Some(bounds) = influence.bounding_box() {
                combined_bounds = match combined_bounds {
                    Some(existing) => Some(existing.union(&bounds)),
                    None => Some(bounds),
                };
            } else {
                return None; // Wenn ein Teil keine Bounds hat, hat der Composite auch keine.
            }
        }
        combined_bounds
    }

    fn max_influence_distance(&self) -> Option<f32> {
        if let Some(cached) = self.cached_max_influence_distance {
            return Some(cached);
        }
        // Neuberechnung
        self.influences
            .iter()
            .filter_map(|i| i.max_influence_distance())
            .fold(None, |acc, dist| {
                Some(match acc {
                    None => dist,
                    Some(current) => current.max(dist), // Typischerweise das Maximum relevant
                })
            })
    }

    fn influence_type(&self) -> &'static str {
        "CompositeInfluence"
    }
}
