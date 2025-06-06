// src/math/algorithms/smoothing/catmull_rom.rs

use crate::math::{
    algorithms::smoothing::traits::Smoothing,
    error::{MathError, MathResult},
    utils::constants,
};
use bevy::math::Vec2;

/// Glättet einen Linienzug unter Verwendung von Catmull-Rom Splines.
/// Die Kurve interpoliert durch jeden gegebenen Punkt.
#[derive(Debug, Clone)]
pub struct CatmullRomSmoother {
    /// Anzahl der zu generierenden Punkte zwischen jedem Paar von Originalpunkten.
    pub segments_per_edge: usize,
    /// Tension-Parameter (typischerweise [0, 1]).
    /// 0.0: Standard Catmull-Rom (relativ rund).
    /// 0.5: Centripetal Catmull-Rom (verhindert oft Schleifen und Überschwinger).
    /// 1.0: Chordal Catmull-Rom.
    /// Ein anderer Tension-Begriff (oft auch Alpha genannt) wird hier nicht direkt verwendet,
    /// sondern die alpha-Parametrisierung ist implizit durch die Wahl der `tension` (0.0, 0.5, 1.0).
    /// Wir nennen es hier `tension` wie im Originalcode, aber beachten die Bedeutung.
    pub tension: f32, // Entspricht dem 'alpha' im Centripetal/Chordal Kontext
    /// Ob die exakten Endpunkte eines offenen Linienzugs beibehalten werden sollen.
    /// Bei Catmull-Rom ist dies oft Standard, da die Kurve durch die Punkte geht.
    pub preserve_endpoints_for_open_paths: bool,
}

impl Default for CatmullRomSmoother {
    fn default() -> Self {
        Self {
            segments_per_edge: 10,
            tension: 0.5, // Centripetal ist ein guter Standard
            preserve_endpoints_for_open_paths: true,
        }
    }
}

impl CatmullRomSmoother {
    pub fn new(segments_per_edge: usize) -> Self {
        Self {
            segments_per_edge: segments_per_edge.max(1), // Mindestens 1 Segment
            ..Default::default()
        }
    }

    /// Setzt den Tension-Parameter (alpha für Parametrisierung).
    /// 0.0 = Uniform, 0.5 = Centripetal, 1.0 = Chordal.
    pub fn with_tension(mut self, tension_alpha: f32) -> Self {
        self.tension = tension_alpha.clamp(0.0, 1.0);
        self
    }

    pub fn preserve_endpoints(mut self, preserve: bool) -> Self {
        self.preserve_endpoints_for_open_paths = preserve;
        self
    }

    /// Interpoliert einen Punkt auf einem Catmull-Rom Spline-Segment.
    /// `p0, p1, p2, p3` sind die vier Kontrollpunkte, die das Segment zwischen `p1` und `p2` definieren.
    /// `t_segment` ist der Interpolationsparameter innerhalb des Segments [0, 1].
    fn catmull_rom_interpolate_segment(
        &self,
        p0: Vec2,
        p1: Vec2,
        p2: Vec2,
        p3: Vec2,
        t_segment: f32, // Parameter innerhalb des aktuellen Segments p1-p2
    ) -> Vec2 {
        // Parameterisierung der Knotenpunkte t0, t1, t2, t3
        // t0 = 0
        // t1 = t0 + |p1 - p0|^alpha
        // t2 = t1 + |p2 - p1|^alpha
        // t3 = t2 + |p3 - p2|^alpha
        // Wobei alpha = self.tension ist.

        let t0 = 0.0;
        let t1 = t0
            + (p1 - p0)
                .length_squared()
                .powf(self.tension * 0.5)
                .max(constants::EPSILON); // Längen mit alpha potenziert
        let t2 = t1
            + (p2 - p1)
                .length_squared()
                .powf(self.tension * 0.5)
                .max(constants::EPSILON);
        let t3 = t2
            + (p3 - p2)
                .length_squared()
                .powf(self.tension * 0.5)
                .max(constants::EPSILON);

        // Der interpolierte Punkt liegt auf dem Segment zwischen t1 und t2.
        // Skaliere t_segment (0..1) auf den Bereich [t1, t2]
        let t = t1 + (t2 - t1) * t_segment;

        // Formel von https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline
        // A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
        // A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
        // A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3
        // B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
        // B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3
        // C  = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2  <- Fehler in Wiki-Formel, sollte t2-t1 im Nenner sein für B1/B2?
        // Die Wiki Formel für C ist: C = (t2 - t) / (t2 - t1) * B1 + (t - t1) / (t2 - t1) * B2
        // Das ist korrekt und zielt auf das Segment zwischen P1 und P2 ab.

        // Hilfsfunktion für sichere Division
        let safe_div = |num: f32, den: f32| -> f32 {
            if den.abs() < constants::EPSILON {
                0.0
            } else {
                num / den
            }
        };

        let a1 = p0 * safe_div(t1 - t, t1 - t0) + p1 * safe_div(t - t0, t1 - t0);
        let a2 = p1 * safe_div(t2 - t, t2 - t1) + p2 * safe_div(t - t1, t2 - t1);
        let a3 = p2 * safe_div(t3 - t, t3 - t2) + p3 * safe_div(t - t2, t3 - t2);

        let b1 = a1 * safe_div(t2 - t, t2 - t0) + a2 * safe_div(t - t0, t2 - t0);
        let b2 = a2 * safe_div(t3 - t, t3 - t1) + a3 * safe_div(t - t1, t3 - t1);

        let c = b1 * safe_div(t2 - t, t2 - t1) + b2 * safe_div(t - t1, t2 - t1);
        c
    }
}

impl Smoothing for CatmullRomSmoother {
    fn smooth_points(&self, points: &[Vec2], is_closed: bool) -> MathResult<Vec<Vec2>> {
        let n = points.len();

        if (is_closed && n < 3) || (!is_closed && n < 2) {
            return Err(MathError::InsufficientPoints {
                expected: if is_closed { 3 } else { 2 },
                actual: n,
            });
        }
        // Für Catmull-Rom brauchen wir mindestens 2 Segmente (3 Punkte für offen, 3 für geschlossen)
        // um sinnvolle Kontrollpunkte p0,p1,p2,p3 zu haben.
        // Bei 2 Punkten (offen) oder 3 Punkten (geschlossen, die ein Dreieck bilden)
        // ist eine Catmull-Rom-Kurve oft nicht gut definiert oder führt zu keiner Glättung.
        if (!is_closed && n < 3) || (is_closed && n < 3) {
            // Hier könnte man einfach die Originalpunkte zurückgeben oder einen Fehler.
            // Da Catmull-Rom durch die Punkte geht, ist das Original schon "gesmoothed".
            return Ok(points.to_vec());
        }

        let mut smoothed_path = Vec::new();
        let num_output_points_hint = n * self.segments_per_edge;
        smoothed_path.reserve(num_output_points_hint);

        if is_closed {
            for i in 0..n {
                // Iteriere über jeden Originalpunkt als Start eines Segments
                let p0 = points[(i + n - 1) % n]; // Vorheriger Punkt
                let p1 = points[i]; // Aktueller Punkt (Start des Segments)
                let p2 = points[(i + 1) % n]; // Nächster Punkt (Ende des Segments)
                let p3 = points[(i + 2) % n]; // Punkt nach dem nächsten

                // Füge p1 hinzu, wenn es das allererste Segment ist (oder wenn preserve_endpoints false wäre)
                // Bei Catmull-Rom geht die Kurve durch p1 und p2. Wir sampeln zwischen p1 und p2.
                // Der erste Punkt p1 des Segments wird hinzugefügt.
                // Der letzte Punkt p2 des Segments wird NICHT hier hinzugefügt,
                // sondern als p1 des nächsten Segments.
                // Ausnahme: Das allererste p1 muss hinzugefügt werden.
                // Und für das letzte Segment muss der letzte Punkt (p2) hinzugefügt werden.
                // Diese Logik ist für geschlossene Pfade einfacher:
                // Wir fügen alle interpolierten Punkte hinzu. Polygon::closed() kümmert sich ums Schließen.

                for j in 0..self.segments_per_edge {
                    // `segments_per_edge` Punkte *zwischen* p1 und p2
                    let t_segment = j as f32 / self.segments_per_edge as f32;
                    smoothed_path
                        .push(self.catmull_rom_interpolate_segment(p0, p1, p2, p3, t_segment));
                }
            }
            // Der letzte Punkt der letzten Kurve (der p2 des letzten Segments) fehlt noch,
            // wenn die obige Schleife nur bis `segments_per_edge -1` ginge.
            // Da sie bis `segments_per_edge` geht, ist der Endpunkt des Segments (bei t_segment=1.0)
            // auch enthalten, wenn j = segments_per_edge.
            // Wenn wir `0..self.segments_per_edge` machen, dann ist t_segment nie 1.0.
            // Die `interpolate_segment` interpoliert zwischen p1 und p2.
            // Wenn wir also `segments_per_edge` Samples haben wollen, die p1 enthalten, aber p2 ausschließen
            // (da p2 der Start des nächsten Segments ist), dann `j in 0..self.segments_per_edge`.
            // Der allerletzte Punkt des Polygons muss dann manuell hinzugefügt werden.
            // Bei geschlossenen ist es einfacher: Polygon::closed() kümmert sich drum.
            // Die aktuelle Logik generiert N * segments_per_edge Punkte.
            // Polygon::closed wird den ersten und letzten verbinden.
        } else {
            // Offener Pfad
            if n < 2 {
                return Ok(points.to_vec());
            } // Linie oder Punkt

            // Für offene Pfade brauchen wir "Phantom"-Punkte an den Enden.
            // p0 für das erste Segment und p3 für das letzte Segment.
            for i in 0..(n - 1) {
                // Iteriere über jedes Segment
                let p1 = points[i];
                let p2 = points[i + 1];

                let p0 = if i == 0 {
                    if self.preserve_endpoints_for_open_paths {
                        points[i] // Erster Kontrollpunkt ist der Startpunkt selbst (flache Tangente)
                    // oder extrapolieren: points[i] - (points[i+1] - points[i])
                    } else {
                        points[i] - (points[i + 1] - points[i]) // Extrapoliert
                    }
                } else {
                    points[i - 1]
                };

                let p3 = if i == n - 2 {
                    // Letztes Segment (von n-2 zu n-1)
                    if self.preserve_endpoints_for_open_paths {
                        points[i + 1] // Letzter Kontrollpunkt ist der Endpunkt selbst
                    // oder extrapolieren: points[i+1] + (points[i+1] - points[i])
                    } else {
                        points[i + 1] + (points[i + 1] - points[i]) // Extrapoliert
                    }
                } else {
                    points[i + 2]
                };

                // Für das erste Segment (i=0), füge den Startpunkt p1 hinzu, wenn Endpunkte erhalten bleiben
                // oder wenn es der allererste Punkt ist.
                if i == 0 && self.preserve_endpoints_for_open_paths {
                    smoothed_path.push(p1);
                }

                // Generiere Punkte für das aktuelle Segment (p1 nach p2)
                // Wenn preserve_endpoints, dann ist der erste Punkt (p1) schon drin für i=0.
                // Wir wollen `segments_per_edge` neue Punkte zwischen p1 und p2 generieren.
                // Inklusive p2, wenn es nicht der allerletzte Punkt des gesamten Pfades ist und preserve_endpoints gilt.
                let num_samples_in_segment = self.segments_per_edge + 1; // Inkl. p1 und p2

                for j in 0..num_samples_in_segment {
                    let t_segment = j as f32 / (num_samples_in_segment - 1) as f32; // t von 0 bis 1

                    // Überspringe p1, wenn es nicht das erste Segment ist (oder wenn Endpunkte nicht erhalten werden)
                    // weil p1 der p2 des vorherigen Segments war.
                    // Wenn Endpunkte erhalten werden:
                    //  - Für i=0: p1 (points[0]) ist schon drin. Füge Punkte für t > 0 hinzu.
                    //  - Für i > 0: p1 (points[i]) war der p2 des vorherigen Segments. Füge Punkte für t > 0 hinzu.
                    //  - Für den letzten Punkt p2 (points[n-1]): Dieser wird speziell behandelt.

                    if j == 0 && i > 0 {
                        // p1 dieses Segments ist p2 des vorherigen, schon hinzugefügt
                        continue;
                    }
                    if j == 0 && i == 0 && !self.preserve_endpoints_for_open_paths {
                        // Wenn Endpunkte nicht erhalten, füge ersten interpolierten Punkt hinzu
                        smoothed_path
                            .push(self.catmull_rom_interpolate_segment(p0, p1, p2, p3, t_segment));
                        continue;
                    }
                    if j > 0 {
                        // Alle Punkte außer dem ersten (t=0) des Segments
                        if i == n - 2
                            && j == num_samples_in_segment - 1
                            && self.preserve_endpoints_for_open_paths
                        {
                            // Dies ist der letzte Punkt des gesamten Pfades (points[n-1]), nicht interpolieren, sondern Original nehmen
                            // Wird unten explizit hinzugefügt
                            continue;
                        }
                        smoothed_path
                            .push(self.catmull_rom_interpolate_segment(p0, p1, p2, p3, t_segment));
                    }
                }
            }
            // Füge den letzten Originalpunkt hinzu, wenn Endpunkte erhalten bleiben sollen
            if self.preserve_endpoints_for_open_paths && n > 0 {
                if smoothed_path.last().map_or(true, |v| {
                    v.distance_squared(points[n - 1]) > constants::EPSILON_SQUARED
                }) {
                    smoothed_path.push(points[n - 1]);
                } else if let Some(last_v) = smoothed_path.last_mut() {
                    // Stelle sicher, dass er exakt ist
                    *last_v = points[n - 1];
                }
            }
        }
        Ok(smoothed_path)
    }

    fn supports_iterations(&self) -> bool {
        false // Mehrere Iterationen von Catmull-Rom sind unüblich.
    }
}
