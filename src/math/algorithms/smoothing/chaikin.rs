// src/math/algorithms/smoothing/chaikin.rs

use crate::math::{
    algorithms::smoothing::traits::Smoothing,
    error::{MathError, MathResult},
    utils::constants,
};
use bevy::math::Vec2;

/// Implementiert den Chaikin-Algorithmus zur Subdivision und Glättung von Linienzügen.
///
/// Der Chaikin-Algorithmus ist ein iteratives Verfahren, das einen Linienzug glättet,
/// indem es die Ecken "abschneidet". In jeder Iteration werden für jede Kante des
/// Linienzugs zwei neue Punkte erzeugt, die die ursprünglichen Punkte ersetzen.
/// Dies führt zu einer Verdopplung der Punktanzahl (ungefähr) und einer
/// glatteren Kurve.
#[derive(Debug, Clone)]
pub struct ChaikinSmoother {
    /// Anzahl der anzuwendenden Glättungsiterationen.
    ///
    /// Mehr Iterationen führen zu einer glatteren Kurve, aber auch zu einem
    /// stärkeren "Schrumpfen" des ursprünglichen Linienzugs und einer höheren
    /// Anzahl von Punkten.
    pub iterations: usize,
    /// Schnittverhältnis für neue Punkte auf Kanten (typischerweise 0.25 für Chaikin).
    ///
    /// Dieses Verhältnis bestimmt, wo die neuen Punkte auf einer Kante platziert werden.
    /// Ein Wert von `r` bedeutet, dass neue Punkte bei `r` und `1-r` der Kantenlänge
    /// von einem Endpunkt aus gesehen platziert werden. Für den Standard-Chaikin-Algorithmus
    /// ist dieser Wert 0.25, was bedeutet, dass die neuen Punkte bei 25% und 75%
    /// der Kantenlänge liegen. Der Wert muss im Bereich (0, 0.5) liegen.
    pub cut_ratio: f32,
    /// Ob die exakten Endpunkte eines offenen Linienzugs beibehalten werden sollen.
    ///
    /// Wenn `true`, bleiben der erste und letzte Punkt eines offenen Linienzugs
    /// während des Glättungsprozesses unverändert. Bei geschlossenen Linienzügen
    /// hat diese Einstellung keine Auswirkung.
    pub preserve_endpoints_for_open_paths: bool,
}

impl Default for ChaikinSmoother {
    /// Erstellt eine Standardinstanz von `ChaikinSmoother`.
    ///
    /// Standardwerte:
    /// - `iterations`: 1
    /// - `cut_ratio`: 0.25
    /// - `preserve_endpoints_for_open_paths`: true
    fn default() -> Self {
        Self {
            iterations: 1,
            cut_ratio: 0.25,
            preserve_endpoints_for_open_paths: true,
        }
    }
}

impl ChaikinSmoother {
    /// Erstellt eine neue Instanz von `ChaikinSmoother` mit einer gegebenen Anzahl von Iterationen.
    ///
    /// # Arguments
    ///
    /// * `iterations` - Die Anzahl der Glättungsiterationen. Muss mindestens 1 sein.
    pub fn new(iterations: usize) -> Self {
        Self {
            iterations: iterations.max(1), // Mindestens eine Iteration
            ..Default::default()
        }
    }

    /// Setzt das Schnittverhältnis für den Algorithmus.
    ///
    /// # Arguments
    ///
    /// * `ratio` - Das neue Schnittverhältnis. Wird auf den Bereich (epsilon, 0.5 - epsilon) geklemmt.
    pub fn with_cut_ratio(mut self, ratio: f32) -> Self {
        self.cut_ratio = ratio.clamp(0.0 + constants::EPSILON, 0.5 - constants::EPSILON); // Muss zwischen (0, 0.5) liegen
        self
    }

    /// Legt fest, ob die Endpunkte offener Pfade erhalten bleiben sollen.
    ///
    /// # Arguments
    ///
    /// * `preserve` - `true`, um Endpunkte zu erhalten, `false` andernfalls.
    pub fn preserve_endpoints(mut self, preserve: bool) -> Self {
        self.preserve_endpoints_for_open_paths = preserve;
        self
    }

    /// Führt eine einzelne Iteration des Chaikin-Algorithmus auf einem gegebenen Satz von Punkten durch.
    ///
    /// # Arguments
    ///
    /// * `points` - Der Linienzug, der geglättet werden soll.
    /// * `is_closed` - `true`, wenn der Linienzug geschlossen ist, `false` andernfalls.
    /// * `preserve_open_endpoints` - `true`, wenn die Endpunkte eines offenen Linienzugs erhalten bleiben sollen.
    ///
    /// # Returns
    ///
    /// Ein neuer Linienzug, der das Ergebnis einer Chaikin-Iteration darstellt.
    fn chaikin_iteration_on_points(
        &self,
        points: &[Vec2],
        is_closed: bool,
        preserve_open_endpoints: bool,
    ) -> Vec<Vec2> {
        let n = points.len();
        if n < 2 {
            // Nicht genug Punkte für Kanten
            return points.to_vec();
        }
        if is_closed && n < 3 {
            // Geschlossener Pfad braucht mind. 3 Punkte
            return points.to_vec();
        }

        let mut smoothed_points = Vec::new();
        // Kapazität abschätzen: jede Kante erzeugt 2 Punkte, plus ggf. Endpunkte
        smoothed_points.reserve(n * 2);

        if !is_closed && preserve_open_endpoints {
            smoothed_points.push(points[0]); // Ersten Punkt bei offenen Pfaden beibehalten
        }

        let num_segments = if is_closed { n } else { n - 1 };

        for i in 0..num_segments {
            let p1 = points[i];
            let p2 = points[(i + 1) % n]; // %n für Umlauf bei geschlossenen Pfaden

            let point_q = p1.lerp(p2, self.cut_ratio);
            let point_r = p1.lerp(p2, 1.0 - self.cut_ratio);

            if is_closed {
                smoothed_points.push(point_q);
                smoothed_points.push(point_r);
            } else {
                // Bei offenen Pfaden:
                if i == 0 && !preserve_open_endpoints {
                    // Wenn Endpunkte nicht erhalten werden, den originalen Startpunkt hinzufügen,
                    // bevor die neuen Punkte der ersten Kante kommen.
                    smoothed_points.push(points[0]);
                }

                smoothed_points.push(point_q);

                // Den zweiten Punkt (point_r) nur hinzufügen, wenn es nicht das letzte Segment ist UND Endpunkte erhalten bleiben,
                // ODER wenn Endpunkte nicht erhalten bleiben (dann immer beide hinzufügen).
                if i < num_segments - 1 || !preserve_open_endpoints {
                    smoothed_points.push(point_r);
                }

                if i == num_segments - 1 && !preserve_open_endpoints {
                    // Wenn Endpunkte nicht erhalten werden, den originalen Endpunkt hinzufügen,
                    // nachdem die neuen Punkte der letzten Kante hinzugefügt wurden.
                    smoothed_points.push(points[n - 1]);
                }
            }
        }

        if !is_closed && preserve_open_endpoints && n > 1 {
            // Sicherstellen, dass der letzte Originalpunkt enthalten ist, wenn er erhalten bleiben soll.
            // Dies ist notwendig, da die Schleife point_r für das letzte Segment möglicherweise nicht hinzufügt.
            if smoothed_points.last().map_or(true, |last_v| {
                last_v.distance_squared(points[n - 1]) > constants::EPSILON_SQUARED
            }) {
                // Nur hinzufügen, wenn der letzte Punkt nicht bereits der Endpunkt ist.
                smoothed_points.push(points[n - 1]);
            } else if let Some(last_v) = smoothed_points.last_mut() {
                // Wenn der letzte Punkt sehr nah am Endpunkt ist, ihn exakt setzen.
                *last_v = points[n - 1];
            }
        }
        smoothed_points
    }
}

impl Smoothing for ChaikinSmoother {
    /// Glättet eine Sequenz von Punkten unter Verwendung des Chaikin-Algorithmus.
    ///
    /// Führt die konfigurierte Anzahl von Iterationen des Algorithmus durch.
    ///
    /// # Arguments
    ///
    /// * `points` - Der Linienzug (Slice von `Vec2`), der geglättet werden soll.
    /// * `is_closed` - `true`, wenn der Linienzug geschlossen ist, `false` andernfalls.
    ///
    /// # Returns
    ///
    /// Ein `MathResult` mit dem geglätteten Linienzug (`Vec<Vec2>`) oder einem `MathError`,
    /// falls nicht genügend Punkte vorhanden sind.
    fn smooth_points(&self, points: &[Vec2], is_closed: bool) -> MathResult<Vec<Vec2>> {
        let n = points.len();
        if (is_closed && n < 3) || (!is_closed && n < 2) {
            return Err(MathError::InsufficientPoints {
                expected: if is_closed { 3 } else { 2 },
                actual: n,
            });
        }

        let mut current_vertices = points.to_vec();

        for _ in 0..self.iterations {
            current_vertices = self.chaikin_iteration_on_points(
                &current_vertices,
                is_closed,
                self.preserve_endpoints_for_open_paths,
            );
            // Wenn nach einer Iteration nicht mehr genug Punkte für eine weitere sinnvolle Iteration da sind
            if (is_closed && current_vertices.len() < 3)
                || (!is_closed && current_vertices.len() < 2)
            {
                break;
            }
        }
        Ok(current_vertices)
    }

    /// Gibt an, ob der Chaikin-Smoother mehrere Iterationen unterstützt.
    ///
    /// # Returns
    ///
    /// `true`, da Chaikin ein iterativer Algorithmus ist.
    fn supports_iterations(&self) -> bool {
        true
    }

    // estimate_iterations_for_polygon kann die Standardimplementierung aus dem Trait verwenden.
    // Eine spezifischere Implementierung könnte hier basierend auf der gewünschten Glattheit
    // oder dem Detailerhalt erfolgen, eventuell unter Berücksichtigung der Polygon-Komplexität.
}
