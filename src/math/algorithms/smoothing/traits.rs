// src/math/algorithms/smoothing/traits.rs

use crate::math::{error::MathResult, geometry::polygon::core::Polygon}; // Polygon wird noch für estimate_iterations etc. benötigt
use bevy::math::Vec2;
// Die anderen Smoothing-Traits (Adaptive, QualityControlled etc.) bleiben erstmal hier,
// ihre Signaturen müssen ggf. auch angepasst werden, wenn sie `Polygon` verwenden.

/// Haupt-Trait für Algorithmen, die eine Sequenz von 2D-Punkten glätten.
pub trait Smoothing {
    /// Glättet eine Sequenz von Punkten und gibt die geglätteten Punkte zurück.
    /// `is_closed` gibt an, ob die Punktsequenz einen geschlossenen Pfad darstellt.
    fn smooth_points(&self, points: &[Vec2], is_closed: bool) -> MathResult<Vec<Vec2>>;

    /// Hilfsmethode, um ein Polygon direkt zu glätten.
    /// Dies ist eine Convenience-Funktion, die den Trait nutzt.
    fn smooth_polygon(&self, polygon: &Polygon) -> MathResult<Polygon> {
        if polygon.is_empty() {
            // Für leere Polygone das Original zurückgeben oder einen Fehler, je nach Definition
            return Polygon::new(vec![]); // Erstellt ein leeres, offenes Polygon
        }
        let smoothed_vertices = self.smooth_points(polygon.vertices(), polygon.is_closed())?;

        if polygon.is_closed() {
            Polygon::closed(smoothed_vertices)
        } else {
            Polygon::new(smoothed_vertices)
        }
    }

    /// Gibt an, ob der Smoother mehrere Iterationen unterstützt.
    /// Wenn true, kann `smooth_points` mehrmals aufgerufen werden, um den Effekt zu verstärken.
    fn supports_iterations(&self) -> bool {
        false // Standardmäßig nicht unterstützt
    }

    /// Schätzt eine sinnvolle Anzahl von Iterationen für ein gegebenes Polygon.
    /// Diese Methode könnte komplexer sein und Polygon-Eigenschaften analysieren.
    fn estimate_iterations_for_polygon(&self, polygon: &Polygon) -> usize {
        let _ = polygon; // polygon wird momentan nicht für die Schätzung verwendet
        1 // Standardmäßig eine Iteration
    }
}

// --- Die anderen Traits (AdaptiveSmoothing, QualityControlledSmoothing etc.) ---
// müssen ebenfalls überprüft werden. Wenn sie `Polygon` in ihren Signaturen verwenden,
// müssen wir entscheiden, ob sie auch auf `&[Vec2]` umgestellt werden oder
// spezifisch für Polygone bleiben.

// Beispiel für QualityControlledSmoothing Anpassung:
/*
use crate::math::utils::constants;

pub trait QualityControlledSmoothing: Smoothing {
    fn smooth_with_quality_control(
        &self,
        points: &[Vec2],
        is_closed: bool,
        target_quality: f32, // Ein Metrikwert [0,1]
        max_iterations: usize,
    ) -> MathResult<SmoothingOutput>; // SmoothingOutput statt SmoothingResult für Klarheit

    fn evaluate_smoothness_quality(&self, original_points: &[Vec2], smoothed_points: &[Vec2]) -> SmoothingQualityMetrics;
}

#[derive(Debug, Clone)]
pub struct SmoothingOutput {
    pub points: Vec<Vec2>,
    pub quality: SmoothingQualityMetrics,
    pub iterations_used: usize,
    pub converged: bool,
}

#[derive(Debug, Clone, Default)]
pub struct SmoothingQualityMetrics {
    pub smoothness_score: f32,      // z.B. basierend auf Krümmungsänderung
    pub length_preservation: f32,   // Wie gut die ursprüngliche Länge erhalten blieb
    pub vertex_reduction_ratio: f32,// Verhältnis der Vertex-Anzahl
}

impl SmoothingQualityMetrics {
    pub fn overall_score(&self) -> f32 {
        // Beispielhafte Gewichtung
        (self.smoothness_score * 0.6 + self.length_preservation * 0.4).clamp(0.0, 1.0)
    }
}
*/
// Die alten SmoothingParameters, AlgorithmParameters etc. können hier vorerst bleiben,
// werden aber ggf. angepasst, wenn die Smoother-Implementierungen überarbeitet werden.
// Fürs Erste konzentrieren wir uns auf den Haupt-Smoothing-Trait.

// Bestehende Parameter-Strukturen (gekürzt, müssen ggf. angepasst werden):
#[derive(Debug, Clone)]
pub struct SmoothingParameters {
    pub strength: f32,
    pub iterations: usize,
    pub algorithm_params: AlgorithmParameters,
}
impl Default for SmoothingParameters {
    fn default() -> Self {
        Self {
            strength: 1.0,
            iterations: 1,
            algorithm_params: AlgorithmParameters::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub enum AlgorithmParameters {
    Chaikin { cut_ratio: f32 },
    CatmullRom { segments: usize, tension: f32 },
    Bezier { ratio: f32, samples: usize },
    Generic,
}
impl Default for AlgorithmParameters {
    fn default() -> Self {
        Self::Generic
    }
}
