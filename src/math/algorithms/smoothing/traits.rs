// src/math/algorithms/smoothing/traits.rs

use crate::math::{error::MathResult, geometry::polygon::core::Polygon};
use bevy::math::Vec2;

/// Haupt-Trait für Algorithmen, die eine Sequenz von 2D-Punkten glätten.
///
/// Dieser Trait definiert die grundlegende Schnittstelle, die von allen
/// Glättungsalgorithmen implementiert werden muss. Er ermöglicht es,
/// verschiedene Glättungsstrategien austauschbar zu verwenden.
pub trait Smoothing {
    /// Glättet eine Sequenz von Punkten und gibt die geglätteten Punkte zurück.
    ///
    /// # Arguments
    ///
    /// * `points` - Ein Slice von `Vec2`-Punkten, die den zu glättenden Linienzug darstellen.
    /// * `is_closed` - Ein boolescher Wert, der angibt, ob der Linienzug geschlossen ist
    ///   (d.h. der letzte Punkt ist mit dem ersten verbunden).
    ///
    /// # Returns
    ///
    /// Ein `MathResult`, das entweder einen `Vec<Vec2>` mit den geglätteten Punkten
    /// oder einen `MathError` im Fehlerfall enthält (z.B. wenn nicht genügend
    /// Punkte für den Algorithmus vorhanden sind).
    fn smooth_points(&self, points: &[Vec2], is_closed: bool) -> MathResult<Vec<Vec2>>;

    /// Hilfsmethode, um ein `Polygon` direkt zu glätten.
    ///
    /// Dies ist eine Convenience-Funktion, die `smooth_points` intern verwendet,
    /// indem sie die Vertices und den Geschlossenheitsstatus des Polygons extrahiert.
    ///
    /// # Arguments
    ///
    /// * `polygon` - Eine Referenz auf das zu glättende `Polygon`.
    ///
    /// # Returns
    ///
    /// Ein `MathResult`, das entweder ein neues, geglättetes `Polygon`
    /// oder einen `MathError` enthält.
    fn smooth_polygon(&self, polygon: &Polygon) -> MathResult<Polygon> {
        if polygon.is_empty() {
            // Für leere Polygone das Original zurückgeben oder einen Fehler, je nach Definition.
            // Hier wird ein leeres, offenes Polygon zurückgegeben.
            return Polygon::new(vec![]);
        }
        let smoothed_vertices = self.smooth_points(polygon.vertices(), polygon.is_closed())?;

        if polygon.is_closed() {
            Polygon::closed(smoothed_vertices)
        } else {
            Polygon::new(smoothed_vertices)
        }
    }

    /// Gibt an, ob der spezifische Glättungsalgorithmus mehrere Iterationen unterstützt.
    ///
    /// Einige Algorithmen (z.B. Chaikin) können iterativ angewendet werden, um den
    /// Glättungseffekt zu verstärken. Andere (z.B. einmalige Bézier-Interpolation)
    /// sind nicht für iterative Anwendung gedacht.
    ///
    /// # Returns
    ///
    /// `true`, wenn mehrere Iterationen sinnvoll sind, andernfalls `false`.
    /// Standardmäßig wird `false` zurückgegeben.
    fn supports_iterations(&self) -> bool {
        false
    }

    /// Schätzt eine sinnvolle Anzahl von Iterationen für ein gegebenes Polygon,
    /// falls der Algorithmus Iterationen unterstützt.
    ///
    /// Diese Methode kann von Implementierungen überschrieben werden, um eine
    /// intelligentere Schätzung basierend auf den Eigenschaften des Polygons
    /// (z.B. Komplexität, Vertex-Anzahl) zu ermöglichen.
    ///
    /// # Arguments
    ///
    /// * `polygon` - Das Polygon, für das die Iterationsanzahl geschätzt werden soll.
    ///
    /// # Returns
    ///
    /// Eine geschätzte Anzahl von Iterationen. Standardmäßig wird 1 zurückgegeben.
    fn estimate_iterations_for_polygon(&self, polygon: &Polygon) -> usize {
        let _ = polygon; // polygon wird in der Standardimplementierung nicht verwendet
        1 // Standardmäßig eine Iteration
    }
}

// Zukünftige Erweiterungen für adaptive oder qualitätsgesteuerte Glättung könnten hier folgen.
// Beispiel für QualityControlledSmoothing Anpassung:
/*
use crate::math::utils::constants;

/// Ein Trait für Glättungsalgorithmen, die Qualitätskontrollen ermöglichen.
pub trait QualityControlledSmoothing: Smoothing {
    /// Glättet Punkte unter Berücksichtigung eines Qualitätsziels.
    ///
    /// # Arguments
    ///
    /// * `points` - Der zu glättende Linienzug.
    /// * `is_closed` - Ob der Linienzug geschlossen ist.
    /// * `target_quality` - Ein Metrikwert [0,1], der die gewünschte Glättungsqualität angibt.
    /// * `max_iterations` - Die maximale Anzahl von Iterationen, um das Qualitätsziel zu erreichen.
    ///
    /// # Returns
    ///
    /// Ein `MathResult` mit einem `SmoothingOutput`, das die geglätteten Punkte und Qualitätsmetriken enthält.
    fn smooth_with_quality_control(
        &self,
        points: &[Vec2],
        is_closed: bool,
        target_quality: f32,
        max_iterations: usize,
    ) -> MathResult<SmoothingOutput>;

    /// Evaluiert die Qualität der Glättung.
    ///
    /// # Arguments
    ///
    /// * `original_points` - Der ursprüngliche Linienzug.
    /// * `smoothed_points` - Der geglättete Linienzug.
    ///
    /// # Returns
    ///
    /// `SmoothingQualityMetrics`, die verschiedene Aspekte der Glättung bewerten.
    fn evaluate_smoothness_quality(&self, original_points: &[Vec2], smoothed_points: &[Vec2]) -> SmoothingQualityMetrics;
}

/// Enthält das Ergebnis einer Glättungsoperation zusammen mit Qualitätsmetriken.
#[derive(Debug, Clone)]
pub struct SmoothingOutput {
    pub points: Vec<Vec2>,
    pub quality: SmoothingQualityMetrics,
    pub iterations_used: usize,
    pub converged: bool,
}

/// Metriken zur Bewertung der Qualität einer Glättungsoperation.
#[derive(Debug, Clone, Default)]
pub struct SmoothingQualityMetrics {
    /// Ein Wert, der die Glattheit repräsentiert (z.B. basierend auf Krümmungsänderung).
    pub smoothness_score: f32,
    /// Ein Wert, der angibt, wie gut die ursprüngliche Länge des Linienzugs erhalten blieb.
    pub length_preservation: f32,
    /// Das Verhältnis der Vertex-Anzahl vor und nach der Glättung.
    pub vertex_reduction_ratio: f32,
}

impl SmoothingQualityMetrics {
    /// Berechnet einen Gesamtqualitätswert basierend auf den einzelnen Metriken.
    pub fn overall_score(&self) -> f32 {
        // Beispielhafte Gewichtung der Metriken
        (self.smoothness_score * 0.6 + self.length_preservation * 0.4).clamp(0.0, 1.0)
    }
}
*/

/// Allgemeine Parameter für Glättungsoperationen, die von verschiedenen Algorithmen genutzt werden können.
///
/// Diese Struktur kann verwendet werden, um Glättungsparameter auf einer höheren Ebene zu definieren,
/// bevor sie an spezifische Algorithmen weitergegeben werden.
#[derive(Debug, Clone)]
pub struct SmoothingParameters {
    /// Ein allgemeiner Wert für die "Stärke" der Glättung, dessen Interpretation
    /// vom spezifischen Algorithmus abhängen kann.
    pub strength: f32,
    /// Die Anzahl der Iterationen, die der Glättungsalgorithmus durchführen soll,
    /// falls der Algorithmus iterativ ist.
    pub iterations: usize,
    /// Algorithmusspezifische Parameter.
    pub algorithm_params: AlgorithmParameters,
}

impl Default for SmoothingParameters {
    /// Erstellt Standard-Glättungsparameter.
    ///
    /// Standardwerte:
    /// - `strength`: 1.0
    /// - `iterations`: 1
    /// - `algorithm_params`: `AlgorithmParameters::Generic`
    fn default() -> Self {
        Self {
            strength: 1.0,
            iterations: 1,
            algorithm_params: AlgorithmParameters::default(),
        }
    }
}

/// Eine Enumeration für algorithmusspezifische Parameter.
///
/// Dies ermöglicht es, Parameter für verschiedene Glättungsalgorithmen
/// in einer einzigen Struktur (`SmoothingParameters`) zu kapseln.
#[derive(Debug, Clone)]
pub enum AlgorithmParameters {
    /// Parameter für den Chaikin-Algorithmus.
    Chaikin {
        /// Das Schnittverhältnis, siehe [`ChaikinSmoother::cut_ratio`].
        cut_ratio: f32,
    },
    /// Parameter für den Catmull-Rom-Algorithmus.
    CatmullRom {
        /// Anzahl der Segmente pro Kante, siehe [`CatmullRomSmoother::segments_per_edge`].
        segments: usize,
        /// Der Tension-Parameter (Alpha), siehe [`CatmullRomSmoother::tension`].
        tension: f32,
    },
    /// Parameter für den Bézier-Algorithmus.
    Bezier {
        /// Verhältnis für Kontrollpunkte, siehe [`BezierSmoother::control_point_ratio`].
        ratio: f32,
        /// Samples pro Kurve, siehe [`BezierSmoother::samples_per_curve`].
        samples: usize,
    },
    /// Ein generischer Fall, wenn keine spezifischen Parameter benötigt werden oder
    /// der Algorithmus seine Standardwerte verwendet.
    Generic,
}

impl Default for AlgorithmParameters {
    /// Standardmäßig werden generische Parameter zurückgegeben.
    fn default() -> Self {
        Self::Generic
    }
}
