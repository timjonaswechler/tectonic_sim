// src/math/point_distribution/voronoi/config.rs

use crate::math::{
    error::{MathError, MathResult},
    types::Bounds2D,
}; // Bounds2D für optionale Boundary in LloydConfig

/// Konfiguration für die Erzeugung eines Voronoi-basierten Polygons.
#[derive(Debug, Clone)]
pub struct VoronoiConfig {
    /// Anzahl der initialen Punkte, die für die Voronoi-Generierung verwendet werden.
    pub initial_point_count: usize,
    /// Anzahl der Lloyd-Relaxationsiterationen zur Verbesserung der Punktverteilung.
    pub lloyd_iterations: usize,
    /// Zielanzahl der Zellen/Regionen, die im finalen Polygon angestrebt werden (für den Merger).
    pub target_cell_count: usize,
    /// Anzahl der Glättungsiterationen für das finale Polygon.
    pub smoothing_iterations: usize,
    /// Faktor, wie weit die Voronoi-Berechnung über die eigentliche Ziel-Boundary hinausgehen soll,
    /// um Randeffekte zu minimieren (z.B. 0.1 für 10% größere Berechnungs-Boundary).
    pub boundary_padding_factor: f64, // f64, um mit Spade konsistent zu sein, falls nötig, sonst f32
    /// Optionaler Seed für die Zufallszahlengenerierung.
    pub seed: Option<u64>,
    // Optionale explizite Boundary für die Voronoi-Berechnung und Lloyd-Relaxation
    pub calculation_bounds: Option<Bounds2D>,
}

impl Default for VoronoiConfig {
    fn default() -> Self {
        Self {
            initial_point_count: 200, // Weniger als vorher, da Relaxation und Merging teuer sein können
            lloyd_iterations: 5,      // Weniger Iterationen für schnellere Generierung
            target_cell_count: 1,     // Ziel ist oft *ein* zentrales Polygon
            smoothing_iterations: 2,
            boundary_padding_factor: 0.15,
            seed: None,
            calculation_bounds: None, // Standardmäßig keine feste Boundary
        }
    }
}

impl VoronoiConfig {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_points(mut self, count: usize) -> Self {
        self.initial_point_count = count.max(3); // Mindestens 3 Punkte
        self
    }

    pub fn with_lloyd_iterations(mut self, iterations: usize) -> Self {
        self.lloyd_iterations = iterations;
        self
    }

    pub fn with_target_cells(mut self, count: usize) -> Self {
        self.target_cell_count = count.max(1);
        self
    }

    pub fn with_smoothing_iterations(mut self, iterations: usize) -> Self {
        self.smoothing_iterations = iterations;
        self
    }

    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    pub fn with_bounds(mut self, bounds: Bounds2D) -> Self {
        self.calculation_bounds = Some(bounds);
        self
    }

    pub fn validate(&self) -> MathResult<()> {
        if self.initial_point_count < 3 {
            return Err(MathError::InvalidConfiguration {
                message: "Voronoi generation requires at least 3 initial points.".to_string(),
            });
        }
        if self.target_cell_count == 0 {
            return Err(MathError::InvalidConfiguration {
                message: "Target cell count for Voronoi merging must be greater than 0."
                    .to_string(),
            });
        }
        if !(0.0..=1.0).contains(&self.boundary_padding_factor) { // Kann auch > 1 sein, wenn es Sinn macht
            // Keine strikte obere Grenze, aber typischerweise klein
        }
        Ok(())
    }
}

/// Konfiguration für den Lloyd-Relaxationsalgorithmus.
#[derive(Debug, Clone)]
pub struct LloydConfig {
    /// Maximale Anzahl von Iterationen für die Relaxation.
    pub max_iterations: usize,
    /// Konvergenz-Toleranz: Wenn die durchschnittliche Punktbewegung unter diesen Wert fällt, stoppt die Iteration.
    pub convergence_tolerance: f32,
    /// Definiert, wie Punkte an den Grenzen des Relaxationsbereichs behandelt werden.
    pub boundary_handling: BoundaryHandling,
    // Optional: Eine explizite Boundary für die Relaxation.
    // Wenn None, wird sie ggf. aus den Punkten abgeleitet oder ist unendlich.
    // Diese wird typischerweise vom VoronoiBuilder gesetzt.
    // pub fixed_boundary: Option<Bounds2D>, // Entfernt, da dies eher Laufzeitinfo ist
}

/// Strategien für die Behandlung von Punkten an den Grenzen während der Lloyd-Relaxation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryHandling {
    /// Punkte werden an der Grenze gespiegelt.
    Mirror,
    /// Punkte werden auf die nächstgelegene Grenzposition geklemmt.
    Clamp,
    /// Punkte, die eine Grenze überschreiten, erscheinen auf der gegenüberliegenden Seite (Torus-Topologie).
    Wrap,
    /// Grenzpunkte bleiben fixiert (nicht implementiert im generischen Lloyd, wäre spezifisch für Polygon-Relaxation).
    Fixed,
}

impl Default for LloydConfig {
    fn default() -> Self {
        Self {
            max_iterations: 15,          // Standardmäßig weniger Iterationen
            convergence_tolerance: 1e-4, // Etwas lockere Toleranz
            boundary_handling: BoundaryHandling::Clamp,
        }
    }
}
