// src/math/algorithms/clipping.rs

use crate::math::{
    error::{MathError, MathResult},
    types::Bounds2D, // Unser Bounds2D mit Vec2
    utils::constants,
};
use bevy::math::Vec2;

/// Algorithmen zur Durchführung von Polygon-Clipping-Operationen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClippingAlgorithm {
    /// Sutherland-Hodgman Algorithmus (für konvexe Clipper-Polygone).
    SutherlandHodgman,
    // WeilerAtherton, // Komplexer, für beliebige Polygone (hier nicht implementiert)
    // LiangBarsky,    // Für Linien-Clipping gegen Rechtecke (hier nicht für Polygone)
    // CohenSutherland // Für Linien-Clipping gegen Rechtecke (hier nicht für Polygone)
}

impl Default for ClippingAlgorithm {
    fn default() -> Self {
        ClippingAlgorithm::SutherlandHodgman
    }
}

/// Führt Polygon-Clipping-Operationen durch.
pub struct PolygonClipper {
    algorithm: ClippingAlgorithm,
    tolerance: f32,
}

impl Default for PolygonClipper {
    fn default() -> Self {
        Self {
            algorithm: ClippingAlgorithm::default(),
            tolerance: constants::EPSILON * 10.0,
        }
    }
}

impl PolygonClipper {
    pub fn new(algorithm: ClippingAlgorithm) -> Self {
        Self {
            algorithm,
            ..Default::default()
        }
    }

    pub fn with_tolerance(mut self, tolerance: f32) -> Self {
        self.tolerance = tolerance.max(0.0);
        self
    }

    /// Clippt ein Subjekt-Polygon (Liste von Punkten) gegen ein Clipper-Polygon (Liste von Punkten).
    /// Beide Polygone werden als geschlossen und im Gegen-Uhrzeigersinn (CCW) orientiert angenommen
    /// für den Sutherland-Hodgman Algorithmus.
    /// Gibt eine Liste von Punktlisten zurück, die die resultierenden geclippten Polygone darstellen.
    /// Derzeit wird nur Sutherland-Hodgman unterstützt, der ein einzelnes Polygon (oder kein Polygon) zurückgibt.
    pub fn clip_polygon_points(
        &self,
        subject_polygon_points: &[Vec2],
        clipper_polygon_points: &[Vec2],
    ) -> MathResult<Vec<Vec<Vec2>>> {
        if subject_polygon_points.len() < 3 || clipper_polygon_points.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: subject_polygon_points
                    .len()
                    .min(clipper_polygon_points.len()),
            });
        }

        match self.algorithm {
            ClippingAlgorithm::SutherlandHodgman => {
                let result_points =
                    self.sutherland_hodgman_clip(subject_polygon_points, clipper_polygon_points)?;
                if result_points.len() >= 3 {
                    Ok(vec![result_points])
                } else {
                    Ok(Vec::new()) // Kein gültiges Polygon Ergebnis
                }
            } // Andere Algorithmen könnten hier implementiert werden
        }
    }

    /// Clippt ein Subjekt-Polygon (Liste von Punkten) gegen ein achsenparalleles Rechteck (Bounds2D).
    pub fn clip_polygon_against_rectangle(
        &self,
        subject_polygon_points: &[Vec2],
        clip_bounds: &Bounds2D,
    ) -> MathResult<Vec<Vec<Vec2>>> {
        if subject_polygon_points.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: subject_polygon_points.len(),
            });
        }
        if !clip_bounds.is_valid() {
            return Err(MathError::InvalidConfiguration {
                message: "Clip bounds are invalid.".to_string(),
            });
        }

        // Konvertiere Bounds2D zu einer Polygon-Punkteliste für Sutherland-Hodgman
        let clipper_rect_points = [
            clip_bounds.min,
            Vec2::new(clip_bounds.max.x, clip_bounds.min.y),
            clip_bounds.max,
            Vec2::new(clip_bounds.min.x, clip_bounds.max.y),
        ];

        // Sutherland-Hodgman kann auch für Rechteck-Clipping verwendet werden.
        // Spezifische Linien-Clipper wie Cohen-Sutherland oder Liang-Barsky sind hier nicht
        // direkt für Polygon-Clipping anwendbar, sondern für einzelne Liniensegmente.
        match self.algorithm {
            ClippingAlgorithm::SutherlandHodgman => {
                let result_points =
                    self.sutherland_hodgman_clip(subject_polygon_points, &clipper_rect_points)?;
                if result_points.len() >= 3 {
                    Ok(vec![result_points])
                } else {
                    Ok(Vec::new())
                }
            }
        }
    }

    // --- Sutherland-Hodgman Algorithmus Implementierung ---
    // Nimmt an, dass das clipper_polygon konvex ist und die Punkte in CCW-Reihenfolge sind.
    fn sutherland_hodgman_clip(
        &self,
        subject_points: &[Vec2],
        clipper_points: &[Vec2], // Muss konvex und CCW sein
    ) -> MathResult<Vec<Vec2>> {
        let mut current_clipped_points = subject_points.to_vec();

        for i in 0..clipper_points.len() {
            if current_clipped_points.len() < 3 && !current_clipped_points.is_empty() {
                // Weniger als 3 Punkte können kein Polygon mehr bilden
                current_clipped_points.clear(); // Kein Ergebnis mehr möglich
                break;
            }
            if current_clipped_points.is_empty() {
                break;
            }

            let clip_edge_p1 = clipper_points[i];
            let clip_edge_p2 = clipper_points[(i + 1) % clipper_points.len()]; // Nächster Punkt, mit Umlauf

            let input_for_this_edge = current_clipped_points.clone();
            current_clipped_points.clear();

            if input_for_this_edge.is_empty() {
                continue;
            }

            let mut s = input_for_this_edge.last().unwrap().clone(); // Startpunkt des aktuellen Segments des Subjektpolygons

            for &e in &input_for_this_edge {
                // Endpunkt des aktuellen Segments
                let s_is_inside = self.is_point_inside_clip_edge(s, clip_edge_p1, clip_edge_p2);
                let e_is_inside = self.is_point_inside_clip_edge(e, clip_edge_p1, clip_edge_p2);

                if e_is_inside {
                    // Fall 1 (Ende innen) oder Fall 4 (beide innen)
                    if !s_is_inside {
                        // Fall 1: S außen, E innen -> Kante betritt Clipping-Bereich
                        if let Some(intersection) =
                            self.line_segment_intersection(s, e, clip_edge_p1, clip_edge_p2)
                        {
                            current_clipped_points.push(intersection);
                        } else {
                            // Sollte nicht passieren, wenn s außen und e innen ist, und Kanten nicht kollinear sind
                            // oder ein Fehler in is_point_inside_clip_edge oder line_segment_intersection
                        }
                    }
                    current_clipped_points.push(e); // E ist innen, also hinzufügen
                } else if s_is_inside {
                    // Fall 2: S innen, E außen -> Kante verlässt Clipping-Bereich
                    if let Some(intersection) =
                        self.line_segment_intersection(s, e, clip_edge_p1, clip_edge_p2)
                    {
                        current_clipped_points.push(intersection);
                    } else {
                        // Wie oben
                    }
                }
                // Fall 3: S außen, E außen -> nichts hinzufügen
                s = e; // Nächster Startpunkt ist aktueller Endpunkt
            }
        }
        Ok(current_clipped_points)
    }

    /// Prüft, ob ein Punkt `p_test` auf der "inneren" Seite der Clipping-Kante (clip_p1 -> clip_p2) liegt.
    /// Annahme: Clipper-Polygon ist CCW orientiert. "Innen" ist dann links von der Kante.
    #[inline]
    fn is_point_inside_clip_edge(&self, p_test: Vec2, clip_p1: Vec2, clip_p2: Vec2) -> bool {
        // Kreuzprodukt: (clip_p2.x - clip_p1.x) * (p_test.y - clip_p1.y) - (clip_p2.y - clip_p1.y) * (p_test.x - clip_p1.x)
        // Positiv -> p_test ist links von der Kante (oder auf der Kante, wenn CCW)
        // Negativ -> p_test ist rechts von der Kante
        // Null    -> p_test ist kollinear
        let cross_product = (clip_p2.x - clip_p1.x) * (p_test.y - clip_p1.y)
            - (clip_p2.y - clip_p1.y) * (p_test.x - clip_p1.x);
        cross_product >= -self.tolerance // Inklusive Punkte auf der Kante (oder sehr nahe dran)
    }

    /// Berechnet den Schnittpunkt zweier Liniensegmente (p1-p2) und (p3-p4).
    /// Gibt `Some(intersection_point)` zurück, wenn sie sich schneiden, sonst `None`.
    fn line_segment_intersection(&self, p1: Vec2, p2: Vec2, p3: Vec2, p4: Vec2) -> Option<Vec2> {
        let d1 = p2 - p1; // Vektor von p1 zu p2
        let d2 = p4 - p3; // Vektor von p3 zu p4

        let denominator = d1.x * d2.y - d1.y * d2.x; // Kreuzprodukt d1 x d2

        if denominator.abs() < self.tolerance * constants::EPSILON {
            // Linien sind parallel oder kollinear
            // TODO: Kollineare Überlappung behandeln, falls erforderlich.
            // Für einfaches Clipping wird dies oft als kein Schnittpunkt gewertet.
            return None;
        }

        let t_num = (p3.x - p1.x) * d2.y - (p3.y - p1.y) * d2.x;
        let u_num = (p3.x - p1.x) * d1.y - (p3.y - p1.y) * d1.x;
        // Für Schnittpunkt:
        // t_num = (p1-p3).cross_product(d2)
        // u_num = (p1-p3).cross_product(d1)
        // Beachte Vorzeichen: d1.cross(d2) = d1.x*d2.y - d1.y*d2.x
        // (p1-p3).cross(d2) = (p1.x-p3.x)*d2.y - (p1.y-p3.y)*d2.x
        // (p1-p3).cross(d1) = (p1.x-p3.x)*d1.y - (p1.y-p3.y)*d1.x
        // t = ( (p1.x-p3.x)*d2.y - (p1.y-p3.y)*d2.x ) / denominator
        // u = ( (p1.x-p3.x)*d1.y - (p1.y-p3.y)*d1.x ) / denominator <- hier ist das Vorzeichen in der Formel wichtig für u
        // Standardformel für u ist oft: u = - (d1.x*(p1.y-p3.y) - d1.y*(p1.x-p3.x)) / denominator;

        let t = t_num / denominator;
        let u = u_num / denominator; // Wenn u_num wie oben, dann u = u_num / denominator

        // Prüfe, ob der Schnittpunkt auf beiden Segmenten liegt (0 <= t <= 1 und 0 <= u <= 1)
        // Mit Toleranz für Endpunkte
        let t_on_segment = t >= -self.tolerance && t <= 1.0 + self.tolerance;
        let u_on_segment = u >= -self.tolerance && u <= 1.0 + self.tolerance;

        if t_on_segment && u_on_segment {
            Some(p1 + d1 * t.clamp(0.0, 1.0)) // Clamp t, um innerhalb des Segments zu bleiben bei Ungenauigkeiten
        } else {
            None
        }
    }
}

/// Spezifische Clipping-Operationen als Utility-Funktionen.
pub struct ClippingOperations;

impl ClippingOperations {
    /// Clippt eine Punkteliste (Polygon) gegen ein Rechteck.
    pub fn clip_points_to_rectangle(
        subject_points: &[Vec2],
        clip_bounds: &Bounds2D,
        algorithm: ClippingAlgorithm, // Erlaube Algorithmuswahl
        tolerance: Option<f32>,
    ) -> MathResult<Vec<Vec<Vec2>>> {
        let mut clipper = PolygonClipper::new(algorithm);
        if let Some(tol) = tolerance {
            clipper = clipper.with_tolerance(tol);
        }
        clipper.clip_polygon_against_rectangle(subject_points, clip_bounds)
    }

    // clip_to_circle und andere spezifische Formen würden hier eigene Logik benötigen
    // oder eine Approximation des Clippers durch ein Polygon.
}
