// src/math/geometry/polygon/validation.rs

use crate::math::geometry::polygon::{
    core::Polygon,
    properties::{Orientation, PolygonProperties},
};
use crate::math::types::Vector2DExt; // Für Vec2 Methoden wie length(), dot(), cross_product()
use crate::math::utils::constants; // Für EPSILON, nearly_equal, etc.
use bevy::math::Vec2;

/// Verschiedene Validierungsebenen für Polygone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationLevel {
    Basic,         // Nur grundlegende Checks (z.B. Punkt-Anzahl).
    Standard,      // Mittlere Checks (z.B. Selbstüberschneidungen, Degenerierung).
    Comprehensive, // Umfassende Checks (z.B. Orientierung, Winding, Konvexität).
}

/// Ergebnis einer Polygon-Validierung.
#[derive(Debug, Clone)]
pub struct ValidationReport {
    // Umbenannt von ValidationResult, um Konflikt mit Result<T,E> zu vermeiden
    pub is_valid: bool,
    pub errors: Vec<ValidationError>,
    pub warnings: Vec<ValidationWarning>,
    pub suggestions: Vec<ValidationSuggestion>,
}

impl ValidationReport {
    fn new() -> Self {
        Self {
            is_valid: true, // Annahme: gültig, bis Fehler gefunden werden
            errors: Vec::new(),
            warnings: Vec::new(),
            suggestions: Vec::new(),
        }
    }

    fn add_error(&mut self, error: ValidationError) {
        self.is_valid = false;
        self.errors.push(error);
    }

    fn add_warning(&mut self, warning: ValidationWarning) {
        self.warnings.push(warning);
    }

    fn add_suggestion(&mut self, suggestion: ValidationSuggestion) {
        self.suggestions.push(suggestion);
    }
}

/// Mögliche Fehler, die ein Polygon ungültig machen.
#[derive(Debug, Clone, PartialEq)]
pub enum ValidationError {
    InsufficientVertices {
        count: usize,
        minimum: usize,
        is_closed: bool,
    },
    SelfIntersection {
        edge1_indices: (usize, usize),
        edge2_indices: (usize, usize),
        intersection_point: Vec2,
    },
    DegenerateEdge {
        vertex_index1: usize,
        vertex_index2: usize,
        length: f32,
    },
    InvalidVertex {
        vertex_index: usize,
        reason: String,
    },
    NonFiniteVertex {
        vertex_index: usize,
        vertex: Vec2,
    },
    DuplicateConsecutiveVertices {
        vertex_index: usize,
    },
}

/// Warnungen bezüglich potenzieller Probleme, die das Polygon nicht unbedingt ungültig machen.
#[derive(Debug, Clone, PartialEq)]
pub enum ValidationWarning {
    NearlyDegenerateEdge {
        vertex_index1: usize,
        vertex_index2: usize,
        length: f32,
    },
    NearlyCollinearVertices {
        indices: Vec<usize>,
    },
    SmallArea {
        area: f32,
    },
    LargeCoordinateValues {
        vertex_index: usize,
        vertex: Vec2,
    },
    InconsistentWinding {
        total_turn_rad: f32,
        expected_turn_rad: f32,
    },
    ClockwiseOrientation, // Wenn CCW erwartet wird
}

/// Vorschläge zur Verbesserung oder Korrektur des Polygons.
#[derive(Debug, Clone, PartialEq)]
pub enum ValidationSuggestion {
    RemoveVertex { index: usize, reason: String },
    MergeVertices { indices: Vec<usize>, reason: String },
    SimplifyPolygon { reason: String },
    ChangeOrientation,
    ApplySmoothing { iterations: usize },
    CheckSelfIntersections,
}

/// Validator für Polygone.
pub struct PolygonValidator {
    level: ValidationLevel,
    tolerance: f32,
    area_tolerance_factor: f32, // Faktor für Flächen-bezogene Toleranzen
}

impl PolygonValidator {
    pub fn new(level: ValidationLevel) -> Self {
        Self {
            level,
            tolerance: constants::EPSILON * 10.0, // Etwas größere Toleranz für Geometrie
            area_tolerance_factor: 100.0, // Für Flächenvergleiche (area < tolerance * area_tolerance_factor)
        }
    }

    pub fn with_tolerance(mut self, tolerance: f32) -> Self {
        self.tolerance = tolerance;
        self
    }

    pub fn validate(&self, polygon: &Polygon) -> ValidationReport {
        let mut report = ValidationReport::new();

        // --- Basic Validation ---
        self.validate_vertex_count(polygon, &mut report);
        self.validate_vertex_values(polygon, &mut report);
        self.validate_consecutive_duplicates(polygon, &mut report);

        if !report.is_valid || self.level == ValidationLevel::Basic {
            return report;
        }

        // --- Standard Validation (benötigt gültige Basis) ---
        self.validate_edge_lengths(polygon, &mut report);

        if polygon.vertices().len() >= 3 {
            // Self-intersection braucht mind. 3 Kanten (oder 4 Punkte für nicht-adjazente)
            self.validate_self_intersections(polygon, &mut report);
        }

        if !report.is_valid || self.level == ValidationLevel::Standard {
            return report;
        }

        // --- Comprehensive Validation (benötigt gültige Standard-Eigenschaften) ---
        if polygon.vertices().len() >= 3 {
            // Eigenschaften wie Area, Orientation brauchen mind. 3 Punkte
            self.validate_area(polygon, &mut report);
            self.validate_orientation(polygon, &mut report);
            self.validate_winding_consistency(polygon, &mut report);
            // Konvexitätsprüfung könnte hier auch hin, ist aber eher eine Eigenschaft als ein Fehler
        }

        report
    }

    fn validate_vertex_count(&self, polygon: &Polygon, report: &mut ValidationReport) {
        let count = polygon.len();
        // Ein "Polygon" als Fläche braucht mind. 3 *unterschiedliche* Punkte.
        // Eine offene Linienkette braucht mind. 2 Punkte.
        let min_points = if polygon.is_closed() { 3 } else { 2 };

        if count < min_points {
            report.add_error(ValidationError::InsufficientVertices {
                count,
                minimum: min_points,
                is_closed: polygon.is_closed(),
            });
        }
        // Spezifischer Check für geschlossene Polygone mit Duplikat am Ende
        if polygon.is_closed() && count >= 1 {
            let unique_vertices_count =
                if polygon.vertices().first() == polygon.vertices().last() && count > 1 {
                    count - 1
                } else {
                    count
                };
            if unique_vertices_count < 3 {
                report.add_error(ValidationError::InsufficientVertices {
                    count: unique_vertices_count, // Zeige die Anzahl einzigartiger Punkte
                    minimum: 3,
                    is_closed: true,
                });
            }
        }
    }

    fn validate_vertex_values(&self, polygon: &Polygon, report: &mut ValidationReport) {
        for (i, &vertex) in polygon.vertices().iter().enumerate() {
            if !vertex.x.is_finite() || !vertex.y.is_finite() {
                report.add_error(ValidationError::NonFiniteVertex {
                    vertex_index: i,
                    vertex,
                });
            }
            if vertex.x.abs() > 1e7 || vertex.y.abs() > 1e7 {
                // Noch größer als vorher
                report.add_warning(ValidationWarning::LargeCoordinateValues {
                    vertex_index: i,
                    vertex,
                });
            }
        }
    }

    fn validate_consecutive_duplicates(&self, polygon: &Polygon, report: &mut ValidationReport) {
        let vertices = polygon.vertices();
        if vertices.len() < 2 {
            return;
        }

        for i in 0..vertices.len() - 1 {
            let p1 = vertices[i];
            let p2 = vertices[i + 1];
            if p1.distance_squared(p2) < self.tolerance * self.tolerance {
                report.add_error(ValidationError::DuplicateConsecutiveVertices { vertex_index: i });
                report.add_suggestion(ValidationSuggestion::RemoveVertex {
                    index: i + 1,
                    reason: "Consecutive duplicate".to_string(),
                });
            }
        }
        // Bei geschlossenen Polygonen, die nicht explizit den Endpunkt duplizieren,
        // könnte man auch ersten und letzten prüfen (aber das macht Polygon::closed schon).
    }

    fn validate_edge_lengths(&self, polygon: &Polygon, report: &mut ValidationReport) {
        let vertices = polygon.vertices();
        let n = vertices.len();
        if n < 2 {
            return;
        }

        let num_edges = if polygon.is_closed() {
            if n > 0 && vertices.first() == vertices.last() {
                n - 1
            } else {
                n
            }
        } else {
            n - 1
        };
        if num_edges == 0 && n > 0 {
            return;
        } // Einzelner Punkt oder leeres Polygon

        for i in 0..num_edges {
            let p1_idx = i;
            let p2_idx = (i + 1) % n; // %n für Umlauf bei geschlossenen
            let p1 = vertices[p1_idx];
            let p2 = vertices[p2_idx];

            // Überspringe das letzte Segment eines offenen Polygons, wenn es nur aus dem letzten Punkt besteht.
            if !polygon.is_closed() && i == n - 1 {
                continue;
            }

            let length = p1.distance(p2);
            if length < self.tolerance {
                report.add_error(ValidationError::DegenerateEdge {
                    vertex_index1: p1_idx,
                    vertex_index2: p2_idx,
                    length,
                });
                report.add_suggestion(ValidationSuggestion::MergeVertices {
                    indices: vec![p1_idx, p2_idx],
                    reason: "Degenerate edge".to_string(),
                });
            } else if length < self.tolerance * 10.0 {
                report.add_warning(ValidationWarning::NearlyDegenerateEdge {
                    vertex_index1: p1_idx,
                    vertex_index2: p2_idx,
                    length,
                });
            }
        }
    }

    fn validate_area(&self, polygon: &Polygon, report: &mut ValidationReport) {
        if !polygon.is_closed() {
            return;
        } // Fläche ist primär für geschlossene Polygone relevant
        let area = polygon.area(); // PolygonProperties Trait muss im Scope sein
        if area < self.tolerance * self.tolerance * self.area_tolerance_factor {
            report.add_warning(ValidationWarning::SmallArea { area });
            if polygon.len() > 3 {
                report.add_suggestion(ValidationSuggestion::SimplifyPolygon {
                    reason: "Very small area".to_string(),
                });
            }
        }
    }

    fn validate_orientation(&self, polygon: &Polygon, report: &mut ValidationReport) {
        if !polygon.is_closed() {
            return;
        }
        match polygon.orientation() {
            // PolygonProperties Trait muss im Scope sein
            Orientation::Clockwise => {
                report.add_warning(ValidationWarning::ClockwiseOrientation);
                report.add_suggestion(ValidationSuggestion::ChangeOrientation);
            }
            Orientation::Collinear => {
                // Bereits durch validate_area oder edge_lengths abgedeckt, aber ggf. spezifische Warnung
            }
            Orientation::CounterClockwise => {} // Erwünschte Orientierung
        }
    }

    fn line_segment_intersection(
        p1: Vec2,
        p2: Vec2,
        p3: Vec2,
        p4: Vec2,
        tolerance: f32,
    ) -> Option<Vec2> {
        let d1 = p2 - p1; // s1_x, s1_y
        let d2 = p4 - p3; // s2_x, s2_y

        let denominator = -d2.x * d1.y + d1.x * d2.y;

        if denominator.abs() < tolerance {
            return None; // Linien sind parallel oder kollinear
        }

        let s = (-d1.y * (p1.x - p3.x) + d1.x * (p1.y - p3.y)) / denominator;
        let t = (d2.x * (p1.y - p3.y) - d2.y * (p1.x - p3.x)) / denominator;

        if s >= -tolerance && s <= 1.0 + tolerance && t >= -tolerance && t <= 1.0 + tolerance {
            // Kollisionspunkt
            Some(Vec2::new(p1.x + (t * d1.x), p1.y + (t * d1.y)))
        } else {
            None // Kein Schnittpunkt auf den Segmenten
        }
    }

    fn validate_self_intersections(&self, polygon: &Polygon, report: &mut ValidationReport) {
        let vertices = polygon.vertices();
        let n = vertices.len();

        // Mindestens 4 Vertices für nicht-adjazente Kanten-Schnittpunkte bei einfachen Polygonen
        let min_verts_for_check = if polygon.is_closed() {
            if n > 0 && vertices.first() == vertices.last() {
                4
            } else {
                3
            } // Wenn dupliziert, brauchen wir 4 phys. Punkte
        } else {
            3 // Für eine offene Kette mit 3 Punkten gibt es nur 2 Kanten, kein Selbstschnitt möglich
        };
        if n < min_verts_for_check {
            return;
        }

        let num_edges = if polygon.is_closed() {
            if n > 0 && vertices.first() == vertices.last() {
                n - 1
            } else {
                n
            }
        } else {
            n - 1
        };
        if num_edges < 2 {
            return;
        }

        for i in 0..num_edges {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % n]; // %n für Umlauf

            // Vergleiche Kante (i, i+1) mit allen nachfolgenden Kanten (j, j+1),
            // die nicht adjazent sind.
            // Für geschlossene Polygone: starte j bei i+2, bis num_edges-1 (wenn i=0, dann j bis num_edges-2)
            // Für offene Polygone: starte j bei i+2
            let j_start = i + 2;
            let j_end = if polygon.is_closed() && i == 0 {
                num_edges - 1
            } else {
                num_edges
            };

            for j in j_start..j_end {
                if (j + 1) % n == i {
                    continue;
                } // Überspringe adjazente Kante bei Umlauf

                let p3 = vertices[j];
                let p4 = vertices[(j + 1) % n];

                if let Some(intersect_pt) =
                    Self::line_segment_intersection(p1, p2, p3, p4, self.tolerance)
                {
                    // Stelle sicher, dass der Schnittpunkt nicht einer der Endpunkte ist (außer bei "echten" Überlappungen)
                    let is_endpoint_intersection = p1.distance_squared(intersect_pt)
                        < self.tolerance * self.tolerance
                        || p2.distance_squared(intersect_pt) < self.tolerance * self.tolerance
                        || p3.distance_squared(intersect_pt) < self.tolerance * self.tolerance
                        || p4.distance_squared(intersect_pt) < self.tolerance * self.tolerance;

                    // Ein "echter" Schnittpunkt ist, wenn er nicht nur ein gemeinsamer Vertex ist.
                    // Diese Bedingung ist komplex, für einfache Polygone ohne duplizierte Vertices
                    // reicht oft die Prüfung, dass die Kanten nicht adjazent sind.
                    // Hier vereinfacht: Wenn ein Schnittpunkt gefunden wird, der nicht nur ein Endpunkt ist,
                    // dann ist es ein Problem. Die Adjazenzprüfung oben sollte die meisten trivialen Fälle abfangen.
                    if !is_endpoint_intersection || (p1 == p4 || p2 == p3) {
                        // Erlaube "Touching corners" nicht als Fehler
                        // Prüfe ob der Schnittpunkt wirklich "zwischen" den Endpunkten liegt,
                        // und nicht nur auf einem Endpunkt.
                        let t_for_edge1 =
                            (intersect_pt - p1).dot(p2 - p1) / (p2 - p1).length_squared();
                        let t_for_edge2 =
                            (intersect_pt - p3).dot(p4 - p3) / (p4 - p3).length_squared();

                        let is_true_intersection = t_for_edge1 > self.tolerance
                            && t_for_edge1 < 1.0 - self.tolerance
                            && t_for_edge2 > self.tolerance
                            && t_for_edge2 < 1.0 - self.tolerance;

                        if is_true_intersection {
                            report.add_error(ValidationError::SelfIntersection {
                                edge1_indices: (i, (i + 1) % n),
                                edge2_indices: (j, (j + 1) % n),
                                intersection_point: intersect_pt,
                            });
                            // Ein Fehler reicht oft schon aus
                            return;
                        }
                    }
                }
            }
        }
    }

    fn validate_winding_consistency(&self, polygon: &Polygon, report: &mut ValidationReport) {
        let vertices = polygon.vertices();
        let n = vertices.len();
        if n < 3 {
            return;
        }

        let num_unique_points =
            if polygon.is_closed() && n > 0 && vertices.first() == vertices.last() {
                n - 1
            } else {
                n
            };
        if num_unique_points < 3 {
            return;
        }

        let mut total_turn = 0.0;
        for i in 0..num_unique_points {
            let p_prev = vertices[(i + num_unique_points - 1) % num_unique_points];
            let p_curr = vertices[i];
            let p_next = vertices[(i + 1) % num_unique_points];

            let v1 = p_curr - p_prev;
            let v2 = p_next - p_curr;

            // Sicherstellen, dass Vektoren nicht null sind
            if v1.length_squared() > self.tolerance * self.tolerance
                && v2.length_squared() > self.tolerance * self.tolerance
            {
                // angle_between gibt den Winkel zwischen -PI und PI.
                // Wir brauchen das Kreuzprodukt für die Richtung des Turns.
                let cross_z = v1.cross_product(v2); // Vector2DExt
                let dot = v1.dot(v2);
                total_turn += cross_z.atan2(dot);
            }
        }

        let _expected_turn = if polygon.is_closed() {
            constants::TAU
        } else {
            0.0
        }; // Für offene Polygone ist die Summe der Turns nicht zwingend 0

        // Für geschlossene Polygone sollte die Summe der Außenwinkel 2PI sein (oder -2PI je nach Orientierung).
        // Die Summe der Innenwinkel ist (n-2)*PI. Ein Turn nach links ist positiv.
        // Wenn CCW, Summe der Turns = 2PI. Wenn CW, Summe der Turns = -2PI.
        if polygon.is_closed() {
            if (total_turn.abs() - constants::TAU).abs()
                > self.tolerance * num_unique_points as f32 * 0.1
            {
                // Größere Toleranz
                report.add_warning(ValidationWarning::InconsistentWinding {
                    total_turn_rad: total_turn,
                    expected_turn_rad: if total_turn > 0.0 {
                        constants::TAU
                    } else {
                        -constants::TAU
                    },
                });
            }
        }
    }
}

/// Schnellvalidierungsfunktionen.
pub struct QuickValidation;

impl QuickValidation {
    pub fn is_valid_basic(polygon: &Polygon) -> bool {
        PolygonValidator::new(ValidationLevel::Basic)
            .validate(polygon)
            .is_valid
    }

    pub fn has_self_intersections(polygon: &Polygon) -> bool {
        let report = PolygonValidator::new(ValidationLevel::Standard).validate(polygon);
        report
            .errors
            .iter()
            .any(|e| matches!(e, ValidationError::SelfIntersection { .. }))
    }

    pub fn get_validation_report(polygon: &Polygon, level: ValidationLevel) -> ValidationReport {
        PolygonValidator::new(level).validate(polygon)
    }
}
