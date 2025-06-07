// src/math/algorithms/smoothing/bezier.rs

use crate::math::{
    algorithms::smoothing::traits::Smoothing,
    error::{MathError, MathResult},
    utils::constants,
};
use bevy::math::Vec2;

/// Glättet Linienzüge unter Verwendung von kubischen Bézier-Kurven.
///
/// Für jedes Segment des Originalpfades wird eine kubische Bézier-Kurve generiert.
/// Die Kontrollpunkte dieser Kurven werden basierend auf den benachbarten Punkten
/// des Originalpfades berechnet, um einen glatten Übergang zwischen den Segmenten
/// zu gewährleisten.
#[derive(Debug, Clone)]
pub struct BezierSmoother {
    /// Verhältnis für den Abstand der Kontrollpunkte von den Ankerpunkten.
    ///
    /// Dieser Wert liegt typischerweise zwischen 0.1 und 0.5. Ein kleinerer Wert
    /// führt zu engeren Kurven, die näher an den Originalsegmenten bleiben,
    /// während ein größerer Wert zu weicheren, weiter ausladenden Kurven führt.
    pub control_point_ratio: f32,
    /// Anzahl der Punkte, die pro Bézier-Kurvensegment generiert werden.
    ///
    /// Eine höhere Anzahl führt zu einer detaillierteren und glatteren Darstellung
    /// der Kurve, erhöht aber auch die Gesamtanzahl der Punkte im geglätteten Pfad.
    pub samples_per_curve: usize,
    /// Spannungs-Parameter (oft als 'tension' bezeichnet, typisch 0.0 bis 1.0).
    ///
    /// Dieser Parameter beeinflusst, wie stark die Kurve "durch" die impliziten
    /// Kontrollpunkte gezogen wird, die aus den Nachbarpunkten abgeleitet werden.
    /// - `0.0`: Erzeugt Kontrollpunkte, die zu Standard-Catmull-Rom-ähnlichen Kurven führen können,
    ///   was oft weichere Übergänge bedeutet.
    /// - `1.0`: Führt zu strafferen Kurven, die näher an den ursprünglichen Liniensegmenten bleiben,
    ///   indem die Kontrollpunkte näher an den Ankerpunkten liegen oder die Tangenten kürzer sind.
    /// Die Interpretation hier ist, dass eine höhere `tension` die Länge der Tangentenvektoren
    /// (und damit den Einfluss der Nachbarpunkte auf die Krümmung) reduziert.
    pub tension: f32,
    /// Ob die exakten Endpunkte eines offenen Linienzugs beibehalten werden sollen.
    ///
    /// Wenn `true`, bleiben der erste und letzte Punkt eines offenen Linienzugs
    /// während des Glättungsprozesses unverändert. Bei geschlossenen Linienzügen
    /// hat diese Einstellung in der Regel keine direkte Auswirkung auf die expliziten
    /// Start-/Endpunkte, da die Kurve ohnehin durch die Originalpunkte verläuft
    /// (die Bézier-Kurven starten und enden an den Originalpunkten).
    pub preserve_endpoints_for_open_paths: bool,
}

impl Default for BezierSmoother {
    /// Erstellt eine Standardinstanz von `BezierSmoother`.
    ///
    /// Standardwerte:
    /// - `control_point_ratio`: 0.33 (ein Drittel ist ein gängiger Wert)
    /// - `samples_per_curve`: 10
    /// - `tension`: 0.0 (keine zusätzliche Spannung, Standardverhalten)
    /// - `preserve_endpoints_for_open_paths`: true
    fn default() -> Self {
        Self {
            control_point_ratio: 0.33,
            samples_per_curve: 10,
            tension: 0.0,
            preserve_endpoints_for_open_paths: true,
        }
    }
}

impl BezierSmoother {
    /// Erstellt eine neue Instanz von `BezierSmoother`.
    ///
    /// # Arguments
    ///
    /// * `control_point_ratio` - Das Verhältnis für den Abstand der Kontrollpunkte.
    ///   Wird auf den Bereich [0.05, 0.95] geklemmt.
    /// * `samples_per_curve` - Die Anzahl der Samples pro Kurvensegment. Muss mindestens 2 sein.
    pub fn new(control_point_ratio: f32, samples_per_curve: usize) -> Self {
        Self {
            control_point_ratio: control_point_ratio.clamp(0.05, 0.95),
            samples_per_curve: samples_per_curve.max(2), // Mindestens Start- und Endpunkt pro Kurve
            ..Default::default()
        }
    }

    /// Setzt den Tension-Parameter für den Smoother.
    ///
    /// # Arguments
    ///
    /// * `tension` - Der neue Tension-Wert. Wird auf den Bereich [0.0, 1.0] geklemmt.
    pub fn with_tension(mut self, tension: f32) -> Self {
        self.tension = tension.clamp(0.0, 1.0);
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

    /// Berechnet die beiden Kontrollpunkte für eine kubische Bézier-Kurve,
    /// die das Segment zwischen `p1` (aktueller Punkt) und `p2` (nächster Punkt) glättet.
    /// `p0` ist der Punkt vor `p1`, `p3` ist der Punkt nach `p2`.
    /// Diese werden verwendet, um die Tangenten an `p1` und `p2` zu schätzen.
    ///
    /// # Arguments
    /// * `p0` - Der Punkt vor dem Startpunkt der Kurve (`p1`).
    /// * `p1` - Der Startpunkt (Anker) der zu definierenden Bézier-Kurve.
    /// * `p2` - Der Endpunkt (Anker) der zu definierenden Bézier-Kurve.
    /// * `p3` - Der Punkt nach dem Endpunkt der Kurve (`p2`).
    ///
    /// # Returns
    /// Ein Tupel `(cp1, cp2)`:
    /// * `cp1`: Der erste Kontrollpunkt (assoziiert mit `p1`).
    /// * `cp2`: Der zweite Kontrollpunkt (assoziiert mit `p2`).
    fn calculate_cubic_bezier_control_points(
        &self,
        p0: Vec2,
        p1: Vec2,
        p2: Vec2,
        p3: Vec2,
    ) -> (Vec2, Vec2) {
        // Tangente bei p1 (zeigt von p0 zu p2, beeinflusst die Ausrichtung der Kurve nach p1)
        let tangent1 = (p2 - p0).normalize_or_zero();
        // Tangente bei p2 (zeigt von p1 zu p3, beeinflusst die Ausrichtung der Kurve vor p2)
        let tangent2 = (p3 - p1).normalize_or_zero();

        // Distanz zwischen den Ankerpunkten des aktuellen Segments
        let dist_p1_p2 = p1.distance(p2);

        // Skalierungsfaktor für die Länge der Tangentenvektoren, die die Kontrollpunkte definieren.
        // Beeinflusst durch control_point_ratio und tension.
        // Eine höhere tension reduziert den `scale`, was die Kurve "straffer" macht.
        let scale = dist_p1_p2 * self.control_point_ratio * (1.0 - self.tension);

        // Erster Kontrollpunkt: Startet bei p1 und bewegt sich entlang der Tangente tangent1.
        let cp1 = p1 + tangent1 * scale;
        // Zweiter Kontrollpunkt: Startet bei p2 und bewegt sich entgegen der Tangente tangent2.
        let cp2 = p2 - tangent2 * scale;

        (cp1, cp2)
    }
}

impl Smoothing for BezierSmoother {
    /// Glättet eine Sequenz von Punkten unter Verwendung von kubischen Bézier-Kurven.
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
        // Für sehr kurze Linienzüge (z.B. 2 Punkte offen, 3 Punkte geschlossen als Dreieck)
        // ist eine Bézier-Glättung oft nicht notwendig oder führt zu unerwünschten Ergebnissen,
        // wenn nicht genügend Kontext für Tangenten vorhanden ist.
        if n <= 2 && !is_closed {
            // Linie mit 2 Punkten oder Einzelpunkt
            return Ok(points.to_vec());
        }
        // Ein Dreieck (n=3, is_closed) könnte man speziell behandeln oder wie unten fortfahren.
        // Momentan wird es wie andere Polygone behandelt.

        let mut smoothed_path = Vec::new();
        // Kapazität abschätzen: (Anzahl Segmente) * samples_per_curve
        let num_segments_approx = if is_closed { n } else { n - 1 };
        smoothed_path.reserve(num_segments_approx * self.samples_per_curve);

        // Bei offenen Pfaden den ersten Punkt direkt hinzufügen, wenn er erhalten bleiben soll.
        // Die erste Bézier-Kurve startet dann ebenfalls an diesem Punkt.
        if !is_closed && self.preserve_endpoints_for_open_paths {
            smoothed_path.push(points[0]);
        }

        let num_segments_to_process = if is_closed { n } else { n - 1 };

        for i in 0..num_segments_to_process {
            let p1_idx = i;
            let p2_idx = (i + 1) % n; // Nächster Punkt, bei geschlossenem Pfad Umlauf

            let p1 = points[p1_idx]; // Startpunkt des aktuellen Bézier-Segments
            let p2 = points[p2_idx]; // Endpunkt des aktuellen Bézier-Segments

            // Bestimme p0 (Punkt vor p1) und p3 (Punkt nach p2) für die Tangentenberechnung.
            // Für Endsegmente offener Pfade werden diese Punkte extrapoliert,
            // um plausible Tangenten zu erhalten.
            let p0 = if is_closed {
                points[(p1_idx + n - 1) % n] // Vorheriger Punkt mit Umlauf
            } else {
                // Offener Pfad:
                if p1_idx == 0 {
                    // Erster Punkt: Extrapoliere p0 zurück von p1 und p2
                    // p0 = p1 - (p2 - p1) = 2*p1 - p2
                    points[p1_idx] - (points[p1_idx.saturating_add(1).min(n - 1)] - points[p1_idx])
                } else {
                    points[p1_idx - 1] // Normaler vorheriger Punkt
                }
            };

            let p3 = if is_closed {
                points[(p2_idx + 1) % n] // Punkt nach p2 mit Umlauf
            } else {
                // Offener Pfad:
                if p2_idx == n - 1 {
                    // Letzter Punkt (p2 ist der letzte Originalpunkt): Extrapoliere p3 vorwärts von p1 und p2
                    // p3 = p2 + (p2 - p1)
                    points[p2_idx] + (points[p2_idx] - points[p2_idx.saturating_sub(1).min(n - 1)])
                } else {
                    points[p2_idx + 1] // Normaler nächster Punkt
                }
            };

            let (cp1, cp2) = self.calculate_cubic_bezier_control_points(p0, p1, p2, p3);

            // Sample die kubische Bézier-Kurve zwischen p1 und p2.
            // Der erste Punkt (t=0) der Kurve ist p1.
            // Der letzte Punkt (t=1) der Kurve ist p2.
            let num_samples = self.samples_per_curve.max(2); // Mindestens Start- und Endpunkt

            for j in 0..num_samples {
                let t = j as f32 / (num_samples - 1) as f32; // t von 0 bis 1

                // Logik zum Hinzufügen von Punkten, um Duplikate an Segmentgrenzen zu vermeiden:
                // - Bei geschlossenen Pfaden: Füge den Startpunkt (t=0) jeder Kurve nur hinzu,
                //   wenn es die allererste Kurve ist (i=0). Für spätere Kurven wird t=0 übersprungen,
                //   da er mit dem Endpunkt (t=1) der vorherigen Kurve identisch ist.
                // - Bei offenen Pfaden:
                //   - Wenn `preserve_endpoints_for_open_paths` true ist:
                //     - Der allererste Originalpunkt (points[0]) wurde bereits hinzugefügt.
                //     - Für die erste Kurve (i=0): Füge Punkte für t > 0 hinzu.
                //     - Für mittlere Kurven (i > 0 und nicht die letzte): Füge Punkte für t > 0 hinzu.
                //     - Für die letzte Kurve: Füge Punkte für t > 0 hinzu, aber der allerletzte Punkt
                //       (points[n-1]) wird separat am Ende gehandhabt, um Genauigkeit zu gewährleisten.
                //   - Wenn `preserve_endpoints_for_open_paths` false ist:
                //     - Füge für die erste Kurve (i=0) alle Punkte (t >= 0) hinzu.
                //     - Für spätere Kurven (i > 0): Füge Punkte für t > 0 hinzu.

                let mut add_this_point = true;

                if is_closed {
                    if i > 0 && j == 0 {
                        // Überspringe Startpunkt späterer Kurven
                        add_this_point = false;
                    }
                } else {
                    // Offener Pfad
                    if self.preserve_endpoints_for_open_paths {
                        if j == 0 {
                            // Überspringe Startpunkt jeder Kurve (da p1=points[0] oder p1=Ende der Vor-Kurve)
                            add_this_point = false;
                        }
                        if i == num_segments_to_process - 1 && j == num_samples - 1 {
                            // Überspringe den letzten Punkt der letzten Kurve, da points[n-1] explizit gesetzt wird.
                            add_this_point = false;
                        }
                    } else {
                        // Nicht preserve_endpoints
                        if i > 0 && j == 0 {
                            // Überspringe Startpunkt späterer Kurven
                            add_this_point = false;
                        }
                    }
                }
                // Sonderfall: Erster Punkt des allerersten Segments, wenn nicht schon durch preserve_endpoints hinzugefügt
                if i == 0 && j == 0 && (is_closed || !self.preserve_endpoints_for_open_paths) {
                    add_this_point = true;
                }

                if add_this_point {
                    let point_on_curve = cubic_bezier_evaluate(p1, cp1, cp2, p2, t);
                    smoothed_path.push(point_on_curve);
                }
            }
        }

        // Bei offenen Pfaden und Erhaltung der Endpunkte: Stelle sicher, dass der letzte
        // Originalpunkt exakt der letzte Punkt im geglätteten Pfad ist.
        if !is_closed && self.preserve_endpoints_for_open_paths && n > 0 {
            // Überprüfe, ob der letzte hinzugefügte Punkt bereits der Endpunkt ist.
            // Wenn nicht, oder wenn der Pfad leer ist (sollte bei n>0 nicht passieren), füge ihn hinzu.
            if smoothed_path.last().map_or(true, |last_p| {
                last_p.distance_squared(points[n - 1]) > constants::EPSILON_SQUARED * 0.01
            }) {
                // Wenn der letzte Punkt nicht exakt der Endpunkt ist, kann es sein, dass er durch die
                // `add_this_point`-Logik übersprungen wurde oder numerische Ungenauigkeiten bestehen.
                // Wir entfernen den potenziell "fast" letzten Punkt und setzen den exakten Endpunkt.
                if smoothed_path.last().is_some() { // Nur pop, wenn etwas da ist
                    // Optional: Wenn der letzte Punkt sehr nah ist, könnte man ihn entfernen, um Duplikate zu vermeiden.
                    // Hier wird er aber eher durch den exakten Endpunkt ersetzt, falls er nicht schon exakt ist.
                }
                smoothed_path.push(points[n - 1]);
            } else if let Some(last_p) = smoothed_path.last_mut() {
                // Wenn der letzte Punkt sehr nah ist, setze ihn exakt auf den Endpunkt.
                *last_p = points[n - 1];
            }
        }
        Ok(smoothed_path)
    }

    /// Gibt an, ob der Bézier-Smoother mehrere Iterationen unterstützt.
    ///
    /// # Returns
    ///
    /// `false`, da mehrere Iterationen von Bézier-Smoothing typischerweise
    /// nicht angewendet werden und zu starker Schrumpfung oder unerwünschten
    /// Artefakten führen können. Die Glättung erfolgt in einem Durchgang.
    fn supports_iterations(&self) -> bool {
        false
    }
}

/// Evaluiert eine kubische Bézier-Kurve am Parameter `t`.
///
/// Die Kurve wird definiert durch vier Punkte:
/// `p0` (Startpunkt), `p1` (erster Kontrollpunkt), `p2` (zweiter Kontrollpunkt) und `p3` (Endpunkt).
///
/// # Arguments
///
/// * `p0` - Startpunkt der Kurve.
/// * `p1` - Erster Kontrollpunkt.
/// * `p2` - Zweiter Kontrollpunkt.
/// * `p3` - Endpunkt der Kurve.
/// * `t` - Interpolationsparameter, typischerweise im Bereich [0, 1].
///         Wird intern auf [0, 1] geklemmt.
///
/// # Returns
///
/// Der Punkt auf der Bézier-Kurve für den gegebenen Parameter `t`.
fn cubic_bezier_evaluate(p0: Vec2, p1: Vec2, p2: Vec2, p3: Vec2, t: f32) -> Vec2 {
    let t_clamped = t.clamp(0.0, 1.0);
    let u = 1.0 - t_clamped;
    let uu = u * u;
    let uuu = uu * u;
    let tt = t_clamped * t_clamped;
    let ttt = tt * t_clamped;

    // Formel: (1-t)^3 * P0 + 3*(1-t)^2*t * P1 + 3*(1-t)*t^2 * P2 + t^3 * P3
    p0 * uuu + p1 * (3.0 * uu * t_clamped) + p2 * (3.0 * u * tt) + p3 * ttt
}
