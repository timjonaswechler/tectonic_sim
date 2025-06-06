// src/math/algorithms/smoothing/chaikin.rs

use crate::math::{
    algorithms::smoothing::traits::Smoothing, // Der angepasste Trait
    error::{MathError, MathResult},
    utils::constants, // Für EPSILON
};
use bevy::math::Vec2;

/// Implementiert den Chaikin-Algorithmus zur Subdivision und Glättung von Linienzügen.
#[derive(Debug, Clone)]
pub struct ChaikinSmoother {
    /// Anzahl der anzuwendenden Glättungsiterationen.
    pub iterations: usize,
    /// Schnittverhältnis für neue Punkte auf Kanten (typischerweise 0.25 für Chaikin).
    /// Ein Wert von 0.25 bedeutet, neue Punkte werden bei 25% und 75% der Kantenlänge gesetzt.
    pub cut_ratio: f32,
    /// Ob die Endpunkte eines offenen Linienzugs beibehalten werden sollen.
    pub preserve_endpoints_for_open_paths: bool,
}

impl Default for ChaikinSmoother {
    fn default() -> Self {
        Self {
            iterations: 1,
            cut_ratio: 0.25,
            preserve_endpoints_for_open_paths: true,
        }
    }
}

impl ChaikinSmoother {
    pub fn new(iterations: usize) -> Self {
        Self {
            iterations: iterations.max(1), // Mindestens eine Iteration
            ..Default::default()
        }
    }

    pub fn with_cut_ratio(mut self, ratio: f32) -> Self {
        self.cut_ratio = ratio.clamp(0.0 + constants::EPSILON, 0.5 - constants::EPSILON); // Muss zwischen (0, 0.5) liegen
        self
    }

    pub fn preserve_endpoints(mut self, preserve: bool) -> Self {
        self.preserve_endpoints_for_open_paths = preserve;
        self
    }

    /// Führt eine einzelne Iteration des Chaikin-Algorithmus durch.
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

            // Zwei neue Punkte auf der Kante (p1, p2)
            // q_i = (1-ratio)*p_i + ratio*p_{i+1}
            // r_i = ratio*p_i + (1-ratio)*p_{i+1}
            // Vertauscht im Originalcode: (1-ratio) ist der größere Anteil.
            // Standard Chaikin: 1/4 und 3/4 => ratio = 0.25
            // q_i: näher an p1, r_i: näher an p2 (wenn ratio < 0.5)
            let point_q = p1.lerp(p2, self.cut_ratio);
            let point_r = p1.lerp(p2, 1.0 - self.cut_ratio);

            // Für Chaikin (Standard cut_ratio=0.25):
            // Erster neuer Punkt ist bei 0.75*P_i + 0.25*P_{i+1} (wenn man von P_i ausgeht)
            // Zweiter neuer Punkt ist bei 0.25*P_i + 0.75*P_{i+1}
            // Deine alte Logik: p1 + (p2 - p1) * ratio  => (1-ratio)p1 + ratio*p2
            //                  p1 + (p2 - p1) * (1-ratio) => ratio*p1 + (1-ratio)p2
            // Das ist korrekt.

            if is_closed {
                smoothed_points.push(point_q);
                smoothed_points.push(point_r);
            } else {
                // Bei offenen Pfaden:
                // Für das erste Segment (i=0), wenn Endpunkte nicht erhalten werden, beide Punkte hinzufügen.
                // Für mittlere Segmente, beide Punkte.
                // Für das letzte Segment, wenn Endpunkte nicht erhalten, beide Punkte.
                // Mit preserve_endpoints:
                //  - erster Punkt schon drin
                //  - letztes Segment (i = n-2): p1 ist points[n-2], p2 ist points[n-1]
                //    point_q wird hinzugefügt. point_r wird zum neuen Endpunkt (points[n-1])
                if i == 0 && !preserve_open_endpoints {
                    smoothed_points.push(points[0]);
                } // Startpunkt falls nicht preserved

                smoothed_points.push(point_q);

                // Nur point_r hinzufügen, wenn es nicht das vorletzte Segment ist,
                // bei dem der Endpunkt (p2) erhalten bleibt.
                if i < num_segments - 1 || !preserve_open_endpoints {
                    smoothed_points.push(point_r);
                }

                if i == num_segments - 1 && !preserve_open_endpoints {
                    smoothed_points.push(points[n - 1]);
                } // Endpunkt falls nicht preserved
            }
        }

        if !is_closed && preserve_open_endpoints && n > 1 {
            // Letzten Punkt bei offenen Pfaden beibehalten, falls noch nicht durch point_r abgedeckt
            if smoothed_points.last() != Some(&points[n - 1]) {
                // Wenn das letzte `point_r` (aus Segment n-2 zu n-1) nicht der Endpunkt war, füge Endpunkt hinzu.
                // Das passiert, wenn preserve_open_endpoints true ist und i < num_segments -1 nicht zutraf.
                // Korrekterweise:
                // Wenn preserve_open_endpoints, dann ist der erste Punkt schon drin.
                // Die Schleife geht bis num_segments.
                // Das letzte p2 ist points[n-1]. Der letzte point_r ist vom Segment (n-2, n-1).
                // Dieser point_r sollte eigentlich nicht der Endpunkt sein.
                // Die Logik ist etwas knifflig.
                // Alternative für offene Pfade mit Erhaltung:
                // 1. smoothed_points.push(points[0])
                // 2. for i in 0..n-1:
                //    p1 = points[i], p2 = points[i+1]
                //    q = p1.lerp(p2, ratio)
                //    r = p1.lerp(p2, 1-ratio)
                //    if i > 0: smoothed_points.push(q) // q von vorheriger Kante ist schon drin
                //    else: smoothed_points.push(q) // Für die erste Kante
                //    smoothed_points.push(r)
                // 3. smoothed_points.push(points[n-1])
                // Diese Logik muss angepasst werden, um Duplikate zu vermeiden.

                // Einfachere Logik für preserve_endpoints_for_open_paths:
                // Die Iteration erzeugt eine neue Punktmenge. Wenn `preserve_endpoints_for_open_paths` true ist,
                // werden der erste und letzte Punkt der *Originalmenge* als erster und letzter Punkt
                // der *neuen Menge* gesetzt, nachdem die Iteration auf den *inneren* Punkten lief.
                // Die aktuelle Implementierung mit `smoothed_points.push(points[0]);` am Anfang ist ein Teil davon.
                // Am Ende: `smoothed_points.push(points[n-1]);`
                // Dazwischen die Logik für die inneren Kanten.

                // Korrigierte Logik für offenen Pfad mit Endpunkterhaltung:
                if n > 1
                    && smoothed_points.last().map_or(true, |last_v| {
                        last_v.distance_squared(points[n - 1]) > constants::EPSILON_SQUARED
                    })
                {
                    smoothed_points.push(points[n - 1]);
                }
            }
        }

        // Bei geschlossenen Polygonen: Der erste Punkt der Iteration (aus Kante 0->1)
        // und der letzte Punkt (aus Kante n-1 -> 0) sollten zusammenfallen oder nahe sein.
        // Die Polygon::closed() Methode behandelt das Schließen später.
        smoothed_points
    }
}

impl Smoothing for ChaikinSmoother {
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

    fn supports_iterations(&self) -> bool {
        true
    }

    // estimate_iterations_for_polygon kann die Standardimplementierung aus dem Trait verwenden
    // oder hier spezifischer werden, z.B. basierend auf der gewünschten Glattheit vs. Detailerhalt.
}
