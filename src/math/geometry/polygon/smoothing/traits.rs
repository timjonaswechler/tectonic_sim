// src/math/geometry/polygon/smoothing/traits.rs

use super::super::Polygon;
use crate::math::{error::*, types::*, utils::*};

/// Haupt-Trait für Polygon-Glättung
pub trait Smoothing {
    /// Glättet ein Polygon und gibt das Ergebnis zurück
    fn smooth(&self, polygon: &Polygon) -> MathResult<Polygon>;

    /// Glättet ein Polygon in-place
    fn smooth_mut(&self, polygon: &mut Polygon) -> MathResult<()> {
        let smoothed = self.smooth(polygon)?;
        *polygon = smoothed;
        Ok(())
    }

    /// Gibt an, ob der Smoother mehrere Iterationen unterstützt
    fn supports_iterations(&self) -> bool {
        false
    }

    /// Schätzt die optimale Anzahl von Iterationen für ein Polygon
    fn estimate_iterations(&self, polygon: &Polygon) -> usize {
        let _ = polygon;
        1
    }
}

/// Trait für adaptive Glättung basierend auf Polygon-Eigenschaften
pub trait AdaptiveSmoothing: Smoothing {
    /// Analysiert ein Polygon und passt Glättungs-Parameter an
    fn analyze_and_smooth(&mut self, polygon: &Polygon) -> MathResult<Polygon>;

    /// Berechnet optimale Parameter für ein gegebenes Polygon
    fn calculate_optimal_parameters(&self, polygon: &Polygon) -> SmoothingParameters;

    /// Wendet berechnete Parameter an
    fn apply_parameters(&mut self, params: SmoothingParameters);
}

/// Trait für Glättung mit Qualitätskontrolle
pub trait QualityControlledSmoothing: Smoothing {
    /// Glättet mit Qualitätsprüfung und Rückgabe
    fn smooth_with_quality_control(
        &self,
        polygon: &Polygon,
        target_quality: f32,
    ) -> MathResult<SmoothingResult>;

    /// Bewertet die Qualität einer Glättung
    fn evaluate_quality(&self, original: &Polygon, smoothed: &Polygon) -> SmoothingQuality;

    /// Führt iterative Glättung mit Qualitätskontrolle durch
    fn iterative_smooth_with_control(
        &self,
        polygon: &Polygon,
        max_iterations: usize,
        target_quality: f32,
    ) -> MathResult<SmoothingResult> {
        let mut current = polygon.clone();
        let mut iteration = 0;
        let mut best_result = SmoothingResult {
            polygon: current.clone(),
            quality: self.evaluate_quality(polygon, &current),
            iterations_used: 0,
            converged: false,
        };

        while iteration < max_iterations {
            let smoothed = self.smooth(&current)?;
            let quality = self.evaluate_quality(polygon, &smoothed);

            iteration += 1;

            // Prüfe ob Qualität erreicht wurde
            if quality.overall_score() >= target_quality {
                return Ok(SmoothingResult {
                    polygon: smoothed,
                    quality,
                    iterations_used: iteration,
                    converged: true,
                });
            }

            // Prüfe ob Verbesserung stagniert
            if quality.overall_score() <= best_result.quality.overall_score() {
                break;
            }

            best_result = SmoothingResult {
                polygon: smoothed.clone(),
                quality,
                iterations_used: iteration,
                converged: false,
            };

            current = smoothed;
        }

        Ok(best_result)
    }
}

/// Trait für progressive Glättung
pub trait ProgressiveSmoothing: Smoothing {
    /// Führt schrittweise Glättung durch mit Zwischenergebnissen
    fn progressive_smooth(&self, polygon: &Polygon, steps: usize) -> MathResult<Vec<Polygon>>;

    /// Interpoliert zwischen Original und vollständig geglättetem Polygon
    fn interpolated_smooth(&self, polygon: &Polygon, factor: f32) -> MathResult<Polygon>;
}

/// Trait für Glättung mit Constraints
pub trait ConstrainedSmoothing: Smoothing {
    /// Glättet unter Beibehaltung bestimmter Vertices
    fn smooth_with_fixed_vertices(
        &self,
        polygon: &Polygon,
        fixed_indices: &[usize],
    ) -> MathResult<Polygon>;

    /// Glättet unter Beibehaltung der ungefähren Fläche
    fn smooth_preserving_area(&self, polygon: &Polygon, tolerance: f32) -> MathResult<Polygon>;

    /// Glättet unter Beibehaltung der Bounding Box
    fn smooth_preserving_bounds(&self, polygon: &Polygon) -> MathResult<Polygon>;
}

/// Parameter für Glättungs-Algorithmen
#[derive(Debug, Clone)]
pub struct SmoothingParameters {
    /// Glättungs-Stärke (0.0 bis 1.0)
    pub strength: f32,
    /// Anzahl der Iterationen
    pub iterations: usize,
    /// Konvergenz-Toleranz
    pub tolerance: f32,
    /// Algorithmus-spezifische Parameter
    pub algorithm_params: AlgorithmParameters,
}

impl Default for SmoothingParameters {
    fn default() -> Self {
        Self {
            strength: 0.5,
            iterations: 3,
            tolerance: 1e-4,
            algorithm_params: AlgorithmParameters::default(),
        }
    }
}

/// Algorithmus-spezifische Parameter
#[derive(Debug, Clone)]
pub enum AlgorithmParameters {
    /// Chaikin-spezifische Parameter
    Chaikin {
        cut_ratio: f32,
        preserve_endpoints: bool,
    },
    /// Catmull-Rom spezifische Parameter
    CatmullRom {
        segments_per_edge: usize,
        tension: f32,
        alpha: f32,
    },
    /// Bézier-spezifische Parameter
    Bezier {
        control_point_ratio: f32,
        samples_per_curve: usize,
    },
    /// Laplacian-spezifische Parameter
    Laplacian {
        lambda: f32,
        preserve_boundary: bool,
    },
    /// Keine spezifischen Parameter
    Generic,
}

impl Default for AlgorithmParameters {
    fn default() -> Self {
        Self::Generic
    }
}

/// Qualitätsbewertung für Glättungs-Ergebnisse
#[derive(Debug, Clone)]
pub struct SmoothingQuality {
    /// Smoothness-Score (0.0 bis 1.0, höher = glatter)
    pub smoothness: f32,
    /// Erhaltung der ursprünglichen Form (0.0 bis 1.0)
    pub shape_preservation: f32,
    /// Erhaltung der Fläche (0.0 bis 1.0)
    pub area_preservation: f32,
    /// Vertex-Anzahl Veränderung (ratio)
    pub vertex_count_ratio: f32,
    /// Feature-Erhaltung (Ecken, wichtige Punkte)
    pub feature_preservation: f32,
}

impl SmoothingQuality {
    /// Berechnet einen kombinierten Qualitäts-Score
    pub fn overall_score(&self) -> f32 {
        (self.smoothness * 0.3
            + self.shape_preservation * 0.25
            + self.area_preservation * 0.2
            + self.feature_preservation * 0.25)
            .clamp(0.0, 1.0)
    }

    /// Prüft ob die Qualität über einem Schwellwert liegt
    pub fn meets_threshold(&self, threshold: f32) -> bool {
        self.overall_score() >= threshold
    }
}

/// Ergebnis einer Glättungs-Operation
#[derive(Debug, Clone)]
pub struct SmoothingResult {
    /// Das geglättete Polygon
    pub polygon: Polygon,
    /// Qualitätsbewertung
    pub quality: SmoothingQuality,
    /// Anzahl verwendeter Iterationen
    pub iterations_used: usize,
    /// Ob der Algorithmus konvergiert ist
    pub converged: bool,
}

/// Utility-Funktionen für Smoothing-Implementierungen
pub struct SmoothingUtils;

impl SmoothingUtils {
    /// Berechnet die Krümmung an einem Vertex
    pub fn calculate_curvature(prev: Point2D, current: Point2D, next: Point2D) -> f32 {
        let v1 = current - prev;
        let v2 = next - current;

        let v1_len = v1.length();
        let v2_len = v2.length();

        if v1_len < constants::EPSILON || v2_len < constants::EPSILON {
            return 0.0;
        }

        let cos_angle = v1.dot(v2) / (v1_len * v2_len);
        let angle = cos_angle.clamp(-1.0, 1.0).acos();

        // Normalisiere auf [0, 1]
        angle / std::f32::consts::PI
    }

    /// Berechnet die Smoothness eines Polygons
    pub fn calculate_smoothness(polygon: &Polygon) -> f32 {
        let vertices = polygon.vertices();
        if vertices.len() < 3 {
            return 0.0;
        }

        let mut total_curvature = 0.0;
        let n = vertices.len();

        for i in 0..n {
            let prev = vertices[(i + n - 1) % n];
            let curr = vertices[i];
            let next = vertices[(i + 1) % n];

            total_curvature += Self::calculate_curvature(prev, curr, next);
        }

        // Invertiere, da niedrige Krümmung = höhere Smoothness
        1.0 - (total_curvature / n as f32).min(1.0)
    }

    /// Berechnet die Shape-Preservation zwischen zwei Polygonen
    pub fn calculate_shape_preservation(original: &Polygon, smoothed: &Polygon) -> f32 {
        use crate::math::geometry::polygon::PolygonProperties;

        let original_area = original.area();
        let smoothed_area = smoothed.area();

        if original_area < constants::EPSILON {
            return if smoothed_area < constants::EPSILON {
                1.0
            } else {
                0.0
            };
        }

        let area_ratio = (smoothed_area / original_area).min(original_area / smoothed_area);

        // Zusätzlich: Vergleiche Perimeter
        let original_perimeter = original.perimeter();
        let smoothed_perimeter = smoothed.perimeter();

        let perimeter_ratio = if original_perimeter > constants::EPSILON {
            (smoothed_perimeter / original_perimeter).min(original_perimeter / smoothed_perimeter)
        } else {
            1.0
        };

        (area_ratio + perimeter_ratio) * 0.5
    }

    /// Detektiert wichtige Features (Ecken, Wendepunkte)
    pub fn detect_features(polygon: &Polygon, curvature_threshold: f32) -> Vec<usize> {
        let vertices = polygon.vertices();
        let mut features = Vec::new();

        for i in 0..vertices.len() {
            let prev = vertices[(i + vertices.len() - 1) % vertices.len()];
            let curr = vertices[i];
            let next = vertices[(i + 1) % vertices.len()];

            let curvature = Self::calculate_curvature(prev, curr, next);
            if curvature > curvature_threshold {
                features.push(i);
            }
        }

        features
    }

    /// Berechnet optimale Glättungs-Parameter basierend auf Polygon-Eigenschaften
    pub fn suggest_parameters(polygon: &Polygon) -> SmoothingParameters {
        use crate::math::geometry::polygon::PolygonProperties;

        let vertex_count = polygon.len();
        let area = polygon.area();
        let perimeter = polygon.perimeter();
        let smoothness = Self::calculate_smoothness(polygon);

        // Basis-Parameter basierend auf Polygon-Größe
        let base_iterations = match vertex_count {
            0..=10 => 1,
            11..=50 => 2,
            51..=200 => 3,
            _ => 4,
        };

        // Anpassung basierend auf aktueller Smoothness
        let iterations = if smoothness < 0.3 {
            base_iterations + 2 // Mehr Glättung für eckige Polygone
        } else if smoothness > 0.8 {
            1 // Wenig Glättung für bereits glatte Polygone
        } else {
            base_iterations
        };

        // Stärke basierend auf Komplexität
        let complexity_ratio = if perimeter > 0.0 {
            area / (perimeter * perimeter)
        } else {
            0.0
        };

        let strength = if complexity_ratio < 0.01 {
            0.3 // Schwache Glättung für elongierte Polygone
        } else {
            0.6 // Stärkere Glättung für kompakte Polygone
        };

        SmoothingParameters {
            strength,
            iterations,
            tolerance: 1e-4,
            algorithm_params: AlgorithmParameters::Chaikin {
                cut_ratio: 0.25,
                preserve_endpoints: false,
            },
        }
    }
}
