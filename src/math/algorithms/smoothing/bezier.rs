// src/math/algorithms/smoothing/bezier.rs

use crate::math::{
    algorithms::smoothing::traits::Smoothing, // Der angepasste Trait
    error::{MathError, MathResult},
    utils::constants,
};
use bevy::math::Vec2;

/// Glättet Linienzüge unter Verwendung von kubischen Bézier-Kurven.
/// Für jedes Segment des Originalpfades wird eine Bézier-Kurve generiert,
/// wobei die Kontrollpunkte basierend auf den benachbarten Punkten berechnet werden.
#[derive(Debug, Clone)]
pub struct BezierSmoother {
    /// Verhältnis für den Abstand der Kontrollpunkte von den Ankerpunkten
    /// (typischerweise zwischen 0.1 und 0.5). Ein kleinerer Wert führt zu
    /// engeren Kurven, ein größerer zu weicheren, weiter ausladenden Kurven.
    pub control_point_ratio: f32,
    /// Anzahl der Punkte, die pro Bézier-Kurvensegment generiert werden.
    pub samples_per_curve: usize,
    /// Spannungs-Parameter (oft als 'tension' bezeichnet, typisch 0.0 bis 1.0).
    /// Beeinflusst, wie stark die Kurve "durch" die Kontrollpunkte gezogen wird.
    /// 0.0 = Standard Catmull-Rom-ähnliche Kontrollpunkte,
    /// 1.0 = Straffere Kurven, die näher an den Originalsegmenten bleiben.
    /// (Hier anders interpretiert als im Original: tension von 0 bis 1, wo 0 die meiste "Rundung" gibt)
    pub tension: f32,
    /// Ob die exakten Endpunkte eines offenen Linienzugs beibehalten werden sollen.
    pub preserve_endpoints_for_open_paths: bool,
}

impl Default for BezierSmoother {
    fn default() -> Self {
        Self {
            control_point_ratio: 0.33, // Ein Drittel ist ein gängiger Wert
            samples_per_curve: 10,
            tension: 0.0, // Keine zusätzliche Spannung, Standardverhalten
            preserve_endpoints_for_open_paths: true,
        }
    }
}

impl BezierSmoother {
    pub fn new(control_point_ratio: f32, samples_per_curve: usize) -> Self {
        Self {
            control_point_ratio: control_point_ratio.clamp(0.05, 0.95),
            samples_per_curve: samples_per_curve.max(2), // Mindestens Start- und Endpunkt pro Kurve
            ..Default::default()
        }
    }

    pub fn with_tension(mut self, tension: f32) -> Self {
        self.tension = tension.clamp(0.0, 1.0); // tension [0,1]
        self
    }

    pub fn preserve_endpoints(mut self, preserve: bool) -> Self {
        self.preserve_endpoints_for_open_paths = preserve;
        self
    }

    /// Berechnet die beiden Kontrollpunkte für eine kubische Bézier-Kurve,
    /// die das Segment zwischen `p1` (aktueller Punkt) und `p2` (nächster Punkt) glättet.
    /// `p0` ist der Punkt vor `p1`, `p3` ist der Punkt nach `p2`.
    /// Diese werden verwendet, um die Tangenten an `p1` und `p2` zu schätzen.
    fn calculate_cubic_bezier_control_points(
        &self,
        p0: Vec2, // Punkt vor dem Start der Kurve
        p1: Vec2, // Startpunkt der Kurve (Anker)
        p2: Vec2, // Endpunkt der Kurve (Anker)
        p3: Vec2, // Punkt nach dem Ende der Kurve
    ) -> (Vec2, Vec2) {
        // (Kontrollpunkt bei p1, Kontrollpunkt bei p2)
        // Tangente bei p1 (zeigt von p0 zu p2)
        let tangent1 = (p2 - p0).normalize_or_zero();
        // Tangente bei p2 (zeigt von p1 zu p3)
        let tangent2 = (p3 - p1).normalize_or_zero();

        // Länge der Segmente
        let dist1 = p1.distance(p2);

        // Skalierungsfaktor für Kontrollpunkte
        let scale = dist1 * self.control_point_ratio * (1.0 - self.tension);

        let cp1 = p1 + tangent1 * scale; // Kontrollpunkt nach p1
        let cp2 = p2 - tangent2 * scale; // Kontrollpunkt vor p2

        (cp1, cp2)
    }
}

impl Smoothing for BezierSmoother {
    fn smooth_points(&self, points: &[Vec2], is_closed: bool) -> MathResult<Vec<Vec2>> {
        let n = points.len();
        if (is_closed && n < 3) || (!is_closed && n < 2) {
            return Err(MathError::InsufficientPoints {
                expected: if is_closed { 3 } else { 2 },
                actual: n,
            });
        }
        if n <= 2 && !is_closed {
            // Linie mit 2 Punkten oder Einzelpunkt
            return Ok(points.to_vec());
        }
        if n <= 3 && is_closed {
            // Dreieck
            if self.preserve_endpoints_for_open_paths && !is_closed {
                // Sollte nicht passieren bei n<=3
                return Ok(points.to_vec());
            }
            // Für ein Dreieck keine Bézier-Glättung, oder es wird zu einer Linie/Punkt
            // Es sei denn, samples_per_curve ist sehr hoch.
            // Hier geben wir das Original zurück oder eine angenäherte Kurve.
            // Vorerst das Original, da Bézier-Kontrollpunkte für ein Dreieck schwierig sind.
            // return Ok(points.to_vec()); // Oder eine spezifische Logik für kleine Polygone
        }

        let mut smoothed_path = Vec::new();
        // Bei offenen Pfaden den ersten Punkt direkt hinzufügen, wenn er erhalten bleiben soll
        if !is_closed && self.preserve_endpoints_for_open_paths {
            smoothed_path.push(points[0]);
        }

        let num_segments = if is_closed { n } else { n - 1 };

        for i in 0..num_segments {
            let p1_idx = i;
            let p2_idx = (i + 1) % n; // Nächster Punkt, bei geschlossenem Pfad Umlauf

            let p1 = points[p1_idx];
            let p2 = points[p2_idx];

            // Bestimme p0 und p3 für die Tangentenberechnung
            let p0 = if is_closed {
                points[(p1_idx + n - 1) % n] // Vorheriger Punkt mit Umlauf
            } else {
                // Für offenen Pfad: p0 für erstes Segment extrapolieren oder p1 verwenden
                if p1_idx == 0 {
                    points[p1_idx] - (points[p1_idx.saturating_add(1).min(n - 1)] - points[p1_idx])
                } else {
                    points[p1_idx - 1]
                }
            };

            let p3 = if is_closed {
                points[(p2_idx + 1) % n] // Punkt nach p2 mit Umlauf
            } else {
                // Für offenen Pfad: p3 für letztes Segment extrapolieren oder p2 verwenden
                if p2_idx == n - 1 {
                    points[p2_idx] + (points[p2_idx] - points[p2_idx.saturating_sub(1)])
                } else {
                    points[p2_idx + 1]
                }
            };

            let (cp1, cp2) = self.calculate_cubic_bezier_control_points(p0, p1, p2, p3);

            // Sample die kubische Bézier-Kurve zwischen p1 und p2
            // Der erste Punkt (t=0) der Kurve ist p1.
            // Der letzte Punkt (t=1) der Kurve ist p2.
            let num_samples_for_this_curve = if !is_closed
                && self.preserve_endpoints_for_open_paths
                && (i == 0 || i == num_segments - 1)
            {
                // Für erstes/letztes Segment bei Erhaltung der Endpunkte ggf. weniger Samples oder andere Logik
                self.samples_per_curve
            } else {
                self.samples_per_curve
            };

            for j in 0..num_samples_for_this_curve {
                let t = j as f32 / (num_samples_for_this_curve - 1).max(1) as f32; // t von 0 bis 1

                // Bei geschlossenen Pfaden: den Startpunkt jeder Kurve (außer der ersten) überspringen,
                // da er mit dem Endpunkt der vorherigen Kurve identisch ist.
                // Bei offenen Pfaden:
                //  - wenn preserve_endpoints: erster Punkt [0] ist schon drin. Für i=0 nur j>0 hinzufügen.
                //  -                      letzter Punkt [n-1] wird am Ende hinzugefügt. Für i=num_segments-1, bis j < num_samples-1 gehen.
                //  - wenn NICHT preserve_endpoints: Alle Punkte der Kurven hinzufügen.

                let add_this_point = if is_closed {
                    i > 0 || j > 0 // Für die erste Kurve alle Punkte, für spätere nur j>0
                } else {
                    // Offener Pfad
                    if self.preserve_endpoints_for_open_paths {
                        (i == 0 && j > 0)
                            || (i > 0 && i < num_segments - 1)
                            || (i > 0
                                && i == num_segments - 1
                                && j < num_samples_for_this_curve - 1)
                    } else {
                        true // Alle Punkte hinzufügen
                    }
                };

                if add_this_point
                    || (i == 0 && j == 0 && (!is_closed || !self.preserve_endpoints_for_open_paths))
                {
                    // Ersten Punkt immer hinzufügen (außer bei preserve_endpoints und i>0)
                    let point_on_curve = cubic_bezier_evaluate(p1, cp1, cp2, p2, t);
                    smoothed_path.push(point_on_curve);
                }
            }
        }

        // Bei offenen Pfaden, wenn Endpunkte erhalten werden sollen und der letzte Punkt
        // der letzten Kurve noch nicht der tatsächliche Endpunkt ist.
        if !is_closed && self.preserve_endpoints_for_open_paths && n > 0 {
            if smoothed_path.is_empty() || smoothed_path.last() != Some(&points[n - 1]) {
                // Falls der letzte Punkt der letzten Kurve nicht hinzugefügt wurde, jetzt den Endpunkt hinzufügen.
                // Dies ist der Fall, wenn die Schleife für j bis < num_samples-1 ging.
                // Oder wenn die Kurve selbst nicht exakt am Endpunkt endet (numerische Fehler).
                // Besser: Den letzten Punkt explizit setzen, wenn er erhalten werden soll.
                let path_length = smoothed_path.len();
                if let Some(last_smoothed) = smoothed_path.last_mut() {
                    if last_smoothed.distance_squared(points[n - 1]) > constants::EPSILON_SQUARED {
                        // Wenn der letzte berechnete Punkt nicht der Endpunkt ist, den Endpunkt hinzufügen.
                        // Dies kann passieren, wenn die samples_per_curve nicht perfekt passen.
                        // Entferne ggf. den letzten berechneten Punkt, wenn er zu nah ist, um Duplikate zu vermeiden.
                        let path_long_enough = path_length > 1;
                        if path_long_enough
                            && last_smoothed.distance_squared(points[n - 1])
                                < (self.control_point_ratio
                                    * points[n - 2].distance(points[n - 1])
                                    * 0.1)
                                    .powi(2)
                        {
                            // smoothed_path.pop(); // Optional: den "fast" letzten Punkt entfernen
                        }
                        smoothed_path.push(points[n - 1]);
                    } else {
                        *last_smoothed = points[n - 1]; // Exakt setzen
                    }
                } else {
                    // Falls smoothed_path leer war (nur bei n=1 und preserve_endpoints der Fall, was oben abgefangen wird)
                    smoothed_path.push(points[n - 1]);
                }
            }
        } else if is_closed && n > 0 {
            // Stelle sicher, dass bei geschlossenen Pfaden der erste und letzte Punkt übereinstimmen,
            // wenn das Polygon tatsächlich als geschlossene Schleife interpretiert werden soll.
            // Polygon::closed() wird dies später handhaben.
        }

        Ok(smoothed_path)
    }

    fn supports_iterations(&self) -> bool {
        false // Mehrere Iterationen von Bézier-Smoothing sind unüblich und führen oft zu starker Schrumpfung.
    }
}

/// Evaluiert eine kubische Bézier-Kurve an Parameter t.
fn cubic_bezier_evaluate(p0: Vec2, p1: Vec2, p2: Vec2, p3: Vec2, t: f32) -> Vec2 {
    let t_clamped = t.clamp(0.0, 1.0);
    let u = 1.0 - t_clamped;
    let uu = u * u;
    let uuu = uu * u;
    let tt = t_clamped * t_clamped;
    let ttt = tt * t_clamped;

    p0 * uuu + p1 * (3.0 * uu * t_clamped) + p2 * (3.0 * u * tt) + p3 * ttt
}
