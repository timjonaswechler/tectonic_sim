// src/math/algorithms/smoothing/catmull_rom.rs

use crate::math::{
    algorithms::smoothing::traits::Smoothing,
    error::{MathError, MathResult},
    utils::constants,
};
use bevy::math::Vec2;

/// Glättet einen Linienzug unter Verwendung von Catmull-Rom Splines.
///
/// Catmull-Rom Splines sind eine Familie von kubischen interpolierenden Splines,
/// die so formuliert sind, dass die Tangente an jedem Punkt `Pi` durch die Punkte
/// `P(i-1)` und `P(i+1)` bestimmt wird. Die Kurve verläuft durch jeden gegebenen
/// Kontrollpunkt.
/// Diese Implementierung unterstützt verschiedene Parametrisierungen (Uniform, Centripetal, Chordal)
/// über den `tension` (alpha) Parameter.
#[derive(Debug, Clone)]
pub struct CatmullRomSmoother {
    /// Anzahl der zu generierenden Punkte zwischen jedem Paar von Originalpunkten.
    ///
    /// Mehr Segmente führen zu einer glatteren Kurve, aber auch zu mehr Punkten.
    /// Muss mindestens 1 sein.
    pub segments_per_edge: usize,
    /// Tension-Parameter (entspricht Alpha in der Catmull-Rom-Parametrisierung), typischerweise [0, 1].
    ///
    /// Dieser Parameter steuert die Form der Spline-Kurve durch die Parametrisierung der Knotenpunkte:
    /// - `0.0`: **Uniform Catmull-Rom**. Die Knotenpunkte sind gleichmäßig verteilt. Kann zu Schleifen
    ///   oder Überschwingern bei ungleichmäßigen Punktabständen führen.
    /// - `0.5`: **Centripetal Catmull-Rom**. Die Knotenpunktabstände basieren auf der Quadratwurzel
    ///   der Distanz zwischen den Punkten. Dies ist oft die beste Wahl, da es Schleifen und
    ///   Überschwinger verhindert und eine natürlichere Kurve erzeugt.
    /// - `1.0`: **Chordal Catmull-Rom**. Die Knotenpunktabstände basieren direkt auf der Distanz
    ///   zwischen den Punkten.
    ///
    /// Der hier als `tension` bezeichnete Parameter ist das Alpha in der Formel
    /// `ti+1 = ti + |Pi+1 - Pi|^alpha`.
    pub tension: f32,
    /// Ob die exakten Endpunkte eines offenen Linienzugs beibehalten werden sollen.
    ///
    /// Bei Catmull-Rom Splines verläuft die Kurve definitionsgemäß durch die Kontrollpunkte.
    /// Diese Einstellung beeinflusst hauptsächlich, wie die "Phantom"-Punkte an den Enden
    /// eines offenen Pfades behandelt werden, was die Form der Kurve an den Enden bestimmt.
    /// Wenn `true`, wird versucht, die Tangenten an den Enden so zu wählen, dass die Kurve
    /// "natürlich" an den Endpunkten endet, oft durch Duplizieren der Endpunkte als
    /// imaginäre Nachbarpunkte.
    pub preserve_endpoints_for_open_paths: bool,
}

impl Default for CatmullRomSmoother {
    /// Erstellt eine Standardinstanz von `CatmullRomSmoother`.
    ///
    /// Standardwerte:
    /// - `segments_per_edge`: 10
    /// - `tension`: 0.5 (Centripetal Catmull-Rom, eine gute allgemeine Wahl)
    /// - `preserve_endpoints_for_open_paths`: true
    fn default() -> Self {
        Self {
            segments_per_edge: 10,
            tension: 0.5, // Centripetal ist ein guter Standard
            preserve_endpoints_for_open_paths: true,
        }
    }
}

impl CatmullRomSmoother {
    /// Erstellt eine neue Instanz von `CatmullRomSmoother`.
    ///
    /// # Arguments
    ///
    /// * `segments_per_edge` - Die Anzahl der zu generierenden Punkte zwischen jedem
    ///   Paar von Originalpunkten. Muss mindestens 1 sein.
    pub fn new(segments_per_edge: usize) -> Self {
        Self {
            segments_per_edge: segments_per_edge.max(1), // Mindestens 1 Segment
            ..Default::default()
        }
    }

    /// Setzt den Tension-Parameter (Alpha für die Knotenparametrisierung).
    ///
    /// # Arguments
    ///
    /// * `tension_alpha` - Der neue Tension-Wert (Alpha).
    ///   - 0.0 für Uniform,
    ///   - 0.5 für Centripetal,
    ///   - 1.0 für Chordal.
    ///   Wird auf den Bereich [0.0, 1.0] geklemmt.
    pub fn with_tension(mut self, tension_alpha: f32) -> Self {
        self.tension = tension_alpha.clamp(0.0, 1.0);
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

    /// Interpoliert einen Punkt auf einem Catmull-Rom Spline-Segment.
    ///
    /// Das Segment wird durch die vier Kontrollpunkte `p0`, `p1`, `p2`, `p3` definiert,
    /// wobei die eigentliche Kurve zwischen `p1` und `p2` verläuft.
    /// Der Parameter `t_segment` interpoliert von `p1` (bei `t_segment`=0) nach `p2` (bei `t_segment`=1).
    ///
    /// # Arguments
    /// * `p0` - Kontrollpunkt vor `p1`.
    /// * `p1` - Erster Kontrollpunkt des Segments (Kurve startet hier).
    /// * `p2` - Zweiter Kontrollpunkt des Segments (Kurve endet hier).
    /// * `p3` - Kontrollpunkt nach `p2`.
    /// * `t_segment` - Interpolationsparameter innerhalb des Segments [0, 1] zwischen `p1` und `p2`.
    ///
    /// # Returns
    /// Der interpolierte Punkt auf der Kurve.
    fn catmull_rom_interpolate_segment(
        &self,
        p0: Vec2,
        p1: Vec2,
        p2: Vec2,
        p3: Vec2,
        t_segment: f32, // Parameter innerhalb des aktuellen Segments p1-p2 [0,1]
    ) -> Vec2 {
        // Knotenparameter ti basierend auf der gewählten Parametrisierung (Uniform, Centripetal, Chordal)
        // t_i+1 = t_i + |P_i+1 - P_i|^alpha, wobei alpha = self.tension
        // Wir brauchen t0, t1, t2, t3 für die Interpolation zwischen P1 und P2.
        let knot_t0 = 0.0f32;
        let knot_t1 = knot_t0
            + (p1 - p0)
                .length()
                .powf(self.tension)
                .max(constants::EPSILON);
        let knot_t2 = knot_t1
            + (p2 - p1)
                .length()
                .powf(self.tension)
                .max(constants::EPSILON);
        let knot_t3 = knot_t2
            + (p3 - p2)
                .length()
                .powf(self.tension)
                .max(constants::EPSILON);
        // Hinweis: Originalcode verwendete length_squared().powf(tension * 0.5), was äquivalent zu length().powf(tension) ist.
        // max(EPSILON) um Division durch Null bei kollinearen Punkten zu vermeiden.

        // Der globale Parameter t für die Interpolation liegt im Intervall [knot_t1, knot_t2].
        // Wir mappen t_segment (von 0..1) auf dieses Intervall.
        let t_global = knot_t1 + (knot_t2 - knot_t1) * t_segment;

        // Catmull-Rom Spline Interpolationsformel (basierend auf Neville's Algorithmus oder direkter Formel)
        // A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
        // A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
        // A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3
        // B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
        // B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3
        // C  = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2

        // Hilfsfunktion für sichere Division, um NaN bei gleichen Knotenwerten zu vermeiden.
        let safe_div = |numerator: f32, denominator: f32| -> f32 {
            if denominator.abs() < constants::EPSILON {
                // fast null
                0.0 // oder eine andere passende Behandlung, z.B. Grenzwert
            } else {
                numerator / denominator
            }
        };

        // Wenn Knoten zusammenfallen (z.B. t0=t1), müssen die Formeln sorgfältig behandelt werden.
        // Oft bedeutet dies, dass der Beitrag dieses Terms 0 ist oder der Punkt direkt verwendet wird.
        // Die safe_div behandelt dies teilweise.

        let a1_factor0 = safe_div(knot_t1 - t_global, knot_t1 - knot_t0);
        let a1_factor1 = safe_div(t_global - knot_t0, knot_t1 - knot_t0);
        let a1 = p0 * a1_factor0 + p1 * a1_factor1;

        let a2_factor1 = safe_div(knot_t2 - t_global, knot_t2 - knot_t1);
        let a2_factor2 = safe_div(t_global - knot_t1, knot_t2 - knot_t1);
        let a2 = p1 * a2_factor1 + p2 * a2_factor2;

        let a3_factor2 = safe_div(knot_t3 - t_global, knot_t3 - knot_t2);
        let a3_factor3 = safe_div(t_global - knot_t2, knot_t3 - knot_t2);
        let a3 = p2 * a3_factor2 + p3 * a3_factor3;

        // Wenn t0=t1, dann ist a1 = p1 (wenn t_global >= t1). Wenn t1=t2, a2 = p2.
        // Dies ist wichtig, wenn Punkte dupliziert werden (z.B. an Enden).
        // Die Formel sollte robust sein, wenn t_global exakt auf einem Knoten liegt.
        // Wenn t_global = knot_t1, dann A1=P1, A2=P1.
        // Wenn t_global = knot_t2, dann A2=P2, A3=P2.

        let b1_factor1 = safe_div(knot_t2 - t_global, knot_t2 - knot_t0);
        let b1_factor2 = safe_div(t_global - knot_t0, knot_t2 - knot_t0);
        let b1 = a1 * b1_factor1 + a2 * b1_factor2;

        let b2_factor2 = safe_div(knot_t3 - t_global, knot_t3 - knot_t1);
        let b2_factor3 = safe_div(t_global - knot_t1, knot_t3 - knot_t1);
        let b2 = a2 * b2_factor2 + a3 * b2_factor3;

        let c_factor1 = safe_div(knot_t2 - t_global, knot_t2 - knot_t1);
        let c_factor2 = safe_div(t_global - knot_t1, knot_t2 - knot_t1);
        let c = b1 * c_factor1 + b2 * c_factor2;

        // Wenn t_global = knot_t1, dann c = p1.
        // Wenn t_global = knot_t2, dann c = p2.
        // Dies stellt sicher, dass die Kurve durch p1 und p2 interpoliert.
        if (knot_t2 - knot_t1).abs() < constants::EPSILON {
            // p1 und p2 sind quasi derselbe Punkt im Parameterraum
            return p1; // oder p2, sie sind sehr nah
        }

        c
    }
}

impl Smoothing for CatmullRomSmoother {
    /// Glättet eine Sequenz von Punkten unter Verwendung von Catmull-Rom Splines.
    ///
    /// # Arguments
    ///
    /// * `points` - Der Linienzug (Slice von `Vec2`), der geglättet werden soll.
    /// * `is_closed` - `true`, wenn der Linienzug geschlossen ist, `false` andernfalls.
    ///
    /// # Returns
    ///
    /// Ein `MathResult` mit dem geglätteten Linienzug (`Vec<Vec2>`) oder einem `MathError`,
    /// falls nicht genügend Punkte vorhanden sind (mindestens 2 für offen, 3 für geschlossen).
    fn smooth_points(&self, points: &[Vec2], is_closed: bool) -> MathResult<Vec<Vec2>> {
        let n = points.len();

        // Benötigt mindestens 2 Punkte für einen offenen Pfad (um 1 Segment zu bilden)
        // und 3 Punkte für einen geschlossenen Pfad (um 3 Segmente zu bilden).
        // Die Catmull-Rom Interpolation für ein Segment P1-P2 benötigt P0, P1, P2, P3.
        if (!is_closed && n < 2) || (is_closed && n < 3) {
            return Err(MathError::InsufficientPoints {
                expected: if is_closed { 3 } else { 2 },
                actual: n,
            });
        }
        // Wenn nur 2 Punkte für offenen Pfad oder 3 für geschlossenen (Dreieck) vorhanden sind,
        // ist eine "Glättung" im Sinne von Kurvenbildung zwischen Punkten oft nicht das Ziel,
        // da die Kurve ohnehin durch die Punkte geht. Man könnte die Originalpunkte zurückgeben.
        // Die aktuelle Implementierung wird jedoch versuchen, Kurven zu bilden.

        let mut smoothed_path = Vec::new();
        let num_output_points_hint = n * self.segments_per_edge;
        smoothed_path.reserve(num_output_points_hint);

        if is_closed {
            // Für geschlossene Pfade iterieren wir über jeden Punkt als Start eines Segments P1-P2.
            for i in 0..n {
                let p0 = points[(i + n - 1) % n]; // Vorheriger Punkt (zyklisch)
                let p1 = points[i]; // Startpunkt des aktuellen Segments
                let p2 = points[(i + 1) % n]; // Endpunkt des aktuellen Segments
                let p3 = points[(i + 2) % n]; // Nächster Punkt nach P2 (zyklisch)

                // Generiere `segments_per_edge` Punkte *zwischen* p1 und p2 (exklusive p2).
                // Der Punkt p1 (t_segment=0) wird hier hinzugefügt.
                // Der Punkt p2 (t_segment=1) wird als Startpunkt des nächsten Segments hinzugefügt.
                for j in 0..self.segments_per_edge {
                    let t_segment = j as f32 / self.segments_per_edge as f32; // t von 0 bis knapp vor 1
                    smoothed_path
                        .push(self.catmull_rom_interpolate_segment(p0, p1, p2, p3, t_segment));
                }
            }
            // Bei geschlossenen Pfaden stellt Polygon::closed() sicher, dass der Pfad korrekt geschlossen wird.
            // Die obige Logik erzeugt N * segments_per_edge Punkte.
        } else {
            // Offener Pfad
            if n < 2 {
                // Sollte durch die Prüfung oben abgedeckt sein
                return Ok(points.to_vec());
            }

            // Für offene Pfade müssen "Phantom"-Punkte an den Enden für P0 und P3 bereitgestellt werden.
            // Iteriere über jedes Segment von points[i] zu points[i+1].
            for i in 0..(n - 1) {
                let p1 = points[i]; // Startpunkt des aktuellen Segments
                let p2 = points[i + 1]; // Endpunkt des aktuellen Segments

                // Bestimme P0
                let p0 = if i == 0 {
                    // Für das erste Segment: P0 wird oft als P1 (oder 2*P1-P2) gewählt, um die Tangente zu definieren.
                    // Wenn preserve_endpoints_for_open_paths, ist eine flache Tangente (P0=P1) üblich.
                    // Andernfalls Extrapolation.
                    if self.preserve_endpoints_for_open_paths {
                        p1 // P0 = P1 -> Tangente bei P1 ist (P2-P1)
                    } else {
                        p1 - (p2 - p1) // P0 = 2*P1 - P2 -> Tangente bei P1 ist (P2-P0)/2 = (P2 - (2P1-P2))/2 = P2-P1
                    }
                } else {
                    points[i - 1] // Normaler vorheriger Punkt
                };

                // Bestimme P3
                let p3 = if i == n - 2 {
                    // Für das letzte Segment (von P(n-2) zu P(n-1)): P3 wird oft als P(n-1) (oder 2*P(n-1)-P(n-2)) gewählt.
                    if self.preserve_endpoints_for_open_paths {
                        p2 // P3 = P2 -> Tangente bei P2 ist (P2-P1)
                    } else {
                        p2 + (p2 - p1) // P3 = 2*P2 - P1 -> Tangente bei P2 ist (P3-P1)/2 = P2-P1
                    }
                } else {
                    points[i + 2] // Normaler nächster Punkt
                };

                // Füge den Startpunkt P1 des Segments hinzu, außer wenn es sich um P1 eines
                // vorherigen Segments handelt (was hier nicht der Fall ist, da wir explizit P1 hinzufügen).
                // Bei Catmull-Rom geht die Kurve durch P1.
                if i == 0 {
                    // Nur für das allererste Segment den Startpunkt P1 explizit hinzufügen.
                    smoothed_path.push(p1);
                }

                // Generiere Punkte für das aktuelle Segment (P1 nach P2).
                // Wir wollen `segments_per_edge` neue Punkte *nach* P1 bis *einschließlich* P2.
                // Also `segments_per_edge` Schritte von t > 0 bis t = 1.0.
                // Die Schleife geht von j=1 bis segments_per_edge.
                let num_samples_to_add = self.segments_per_edge;

                for j in 1..=num_samples_to_add {
                    // Von t > 0 bis t = 1.0
                    let t_segment = j as f32 / num_samples_to_add as f32;

                    // Wenn es das letzte Segment ist und der letzte Punkt (t_segment=1.0)
                    // und Endpunkte erhalten bleiben sollen, dann direkt den Original-Endpunkt P2 nehmen.
                    if i == n - 2
                        && j == num_samples_to_add
                        && self.preserve_endpoints_for_open_paths
                    {
                        // Dieser Punkt ist points[n-1]. Er wird am Ende explizit hinzugefügt oder überschrieben.
                        // Um Duplikate zu vermeiden, fügen wir ihn hier nicht hinzu, wenn er der letzte ist.
                        // Stattdessen stellen wir sicher, dass der letzte Punkt des Pfades points[n-1] ist.
                        // Die aktuelle Logik fügt ihn hinzu und korrigiert ihn dann ggf.
                        smoothed_path
                            .push(self.catmull_rom_interpolate_segment(p0, p1, p2, p3, t_segment));
                    } else {
                        smoothed_path
                            .push(self.catmull_rom_interpolate_segment(p0, p1, p2, p3, t_segment));
                    }
                }
            }

            // Stelle sicher, dass der letzte Originalpunkt exakt der letzte Punkt im Pfad ist,
            // wenn Endpunkte erhalten bleiben sollen.
            if self.preserve_endpoints_for_open_paths && n > 0 {
                if let Some(last_v) = smoothed_path.last_mut() {
                    if last_v.distance_squared(points[n - 1]) > constants::EPSILON_SQUARED * 0.01 {
                        // Wenn der interpolierte letzte Punkt nicht exakt der Original-Endpunkt ist,
                        // könnte man ihn entweder ersetzen oder den Originalpunkt zusätzlich hinzufügen.
                        // Hier wird er ersetzt, um die Punktzahl nicht unnötig zu erhöhen.
                        *last_v = points[n - 1];
                    } else {
                        // Wenn er schon sehr nah dran ist, exakt setzen.
                        *last_v = points[n - 1];
                    }
                } else if n == 1 && smoothed_path.is_empty() {
                    // Spezialfall: Nur ein Punkt im Original
                    smoothed_path.push(points[0]);
                }
                // Falls der Pfad aus irgendeinem Grund kürzer ist und der letzte Punkt fehlt:
                if smoothed_path.is_empty() && n > 0 {
                    // Sollte nicht passieren bei n>=2
                    smoothed_path.push(points[n - 1]);
                }
            }
        }
        Ok(smoothed_path)
    }

    /// Gibt an, ob der Catmull-Rom-Smoother mehrere Iterationen unterstützt.
    ///
    /// # Returns
    ///
    /// `false`, da Catmull-Rom-Splines interpolierend sind und eine einmalige Anwendung
    /// die Punkte bereits auf die gewünschte Kurve legt. Mehrere Iterationen sind
    /// typischerweise nicht sinnvoll oder verändern die Kurve nicht weiter, es sei denn,
    /// die Kontrollpunkte selbst würden modifiziert.
    fn supports_iterations(&self) -> bool {
        false
    }
}
