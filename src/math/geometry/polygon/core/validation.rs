// src/math/geometry/polygon/core/validation.rs

use super::super::Polygon;
use crate::math::{types::*, utils::*};

/// Verschiedene Validierungsebenen
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationLevel {
    /// Nur grundlegende Checks (Punkt-Anzahl)
    Basic,
    /// Mittlere Checks (Selbstüberschneidungen, Degenerierung)
    Standard,
    /// Umfassende Checks (Orientierung, Winding, Konvexität)
    Comprehensive,
}

/// Validation-Ergebnisse
#[derive(Debug, Clone)]
pub struct ValidationResult {
    pub is_valid: bool,
    pub errors: Vec<ValidationError>,
    pub warnings: Vec<ValidationWarning>,
    pub suggestions: Vec<ValidationSuggestion>,
}

#[derive(Debug, Clone)]
pub enum ValidationError {
    InsufficientVertices { count: usize, minimum: usize },
    SelfIntersection { point1: Point2D, point2: Point2D },
    DegenerateEdge { vertex_index: usize },
    InvalidVertex { vertex_index: usize, reason: String },
    TopologicalError { description: String },
}

#[derive(Debug, Clone)]
pub enum ValidationWarning {
    NearlyDegenerateEdge { vertex_index: usize, length: f32 },
    NearlyCollinearVertices { indices: Vec<usize> },
    IrregularSpacing { max_ratio: f32 },
    PotentialNumericalInstability { reason: String },
}

#[derive(Debug, Clone)]
pub enum ValidationSuggestion {
    RemoveVertex { index: usize },
    MergeVertices { indices: Vec<usize> },
    SimplifyPolygon { target_vertex_count: usize },
    ChangeOrientation,
    ApplySmoothing { iterations: usize },
}

/// Polygon-Validator
pub struct PolygonValidator {
    level: ValidationLevel,
    tolerance: f32,
    check_self_intersections: bool,
    check_orientation: bool,
    check_winding_number: bool,
}

impl PolygonValidator {
    /// Erstellt einen Standard-Validator
    pub fn new(level: ValidationLevel) -> Self {
        Self {
            level,
            tolerance: constants::EPSILON * 1000.0, // Etwas größere Toleranz für praktische Anwendung
            check_self_intersections: matches!(
                level,
                ValidationLevel::Standard | ValidationLevel::Comprehensive
            ),
            check_orientation: matches!(level, ValidationLevel::Comprehensive),
            check_winding_number: matches!(level, ValidationLevel::Comprehensive),
        }
    }

    /// Setzt die Toleranz für numerische Vergleiche
    pub fn with_tolerance(mut self, tolerance: f32) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Aktiviert/deaktiviert Selbstüberschneidungs-Check
    pub fn check_self_intersections(mut self, check: bool) -> Self {
        self.check_self_intersections = check;
        self
    }

    /// Validiert ein Polygon
    pub fn validate(&self, polygon: &Polygon) -> ValidationResult {
        let mut errors = Vec::new();
        let mut warnings = Vec::new();
        let mut suggestions = Vec::new();

        // Basic validation
        self.validate_vertex_count(polygon, &mut errors);
        self.validate_vertices(polygon, &mut errors, &mut warnings);

        if self.level == ValidationLevel::Basic {
            return ValidationResult {
                is_valid: errors.is_empty(),
                errors,
                warnings,
                suggestions,
            };
        }

        // Standard validation
        self.validate_edges(polygon, &mut errors, &mut warnings, &mut suggestions);

        if self.check_self_intersections {
            self.validate_self_intersections(polygon, &mut errors);
        }

        if self.level == ValidationLevel::Standard {
            return ValidationResult {
                is_valid: errors.is_empty(),
                errors,
                warnings,
                suggestions,
            };
        }

        // Comprehensive validation
        if self.check_orientation {
            self.validate_orientation(polygon, &mut warnings, &mut suggestions);
        }

        if self.check_winding_number {
            self.validate_winding_consistency(polygon, &mut warnings);
        }

        self.validate_geometric_properties(polygon, &mut warnings, &mut suggestions);

        ValidationResult {
            is_valid: errors.is_empty(),
            errors,
            warnings,
            suggestions,
        }
    }

    // === Private Validation Methods ===

    fn validate_vertex_count(&self, polygon: &Polygon, errors: &mut Vec<ValidationError>) {
        let count = polygon.len();
        let minimum = if polygon.is_closed() { 4 } else { 3 }; // Geschlossen: 3 unique + 1 duplicate

        if count < minimum {
            errors.push(ValidationError::InsufficientVertices { count, minimum });
        }
    }

    fn validate_vertices(
        &self,
        polygon: &Polygon,
        errors: &mut Vec<ValidationError>,
        warnings: &mut Vec<ValidationWarning>,
    ) {
        for (i, &vertex) in polygon.vertices().iter().enumerate() {
            // Check for NaN or infinite coordinates
            if !vertex.x.is_finite() || !vertex.y.is_finite() {
                errors.push(ValidationError::InvalidVertex {
                    vertex_index: i,
                    reason: "Non-finite coordinates".to_string(),
                });
            }

            // Check for extremely large coordinates (potential overflow)
            if vertex.x.abs() > 1e6 || vertex.y.abs() > 1e6 {
                warnings.push(ValidationWarning::PotentialNumericalInstability {
                    reason: format!("Very large coordinates at vertex {}: {:?}", i, vertex),
                });
            }
        }
    }

    fn validate_edges(
        &self,
        polygon: &Polygon,
        errors: &mut Vec<ValidationError>,
        warnings: &mut Vec<ValidationWarning>,
        suggestions: &mut Vec<ValidationSuggestion>,
    ) {
        let vertices = polygon.vertices();
        let n = vertices.len();

        if n < 2 {
            return;
        }

        let mut edge_lengths = Vec::new();

        for i in 0..n {
            let j = if polygon.is_closed() {
                (i + 1) % (n - 1) // Skip duplicate last vertex
            } else if i == n - 1 {
                break;
            } else {
                i + 1
            };

            let edge_length = vertices[i].distance(vertices[j]);
            edge_lengths.push(edge_length);

            // Degenerate edge check
            if edge_length < self.tolerance {
                if edge_length == 0.0 {
                    errors.push(ValidationError::DegenerateEdge { vertex_index: i });
                    suggestions.push(ValidationSuggestion::RemoveVertex { index: j });
                } else {
                    warnings.push(ValidationWarning::NearlyDegenerateEdge {
                        vertex_index: i,
                        length: edge_length,
                    });
                    suggestions.push(ValidationSuggestion::MergeVertices {
                        indices: vec![i, j],
                    });
                }
            }
        }

        // Check for irregular edge spacing
        if !edge_lengths.is_empty() {
            let min_length = edge_lengths.iter().fold(f32::INFINITY, |a, &b| a.min(b));
            let max_length = edge_lengths.iter().fold(0.0, |a, &b| a.max(b));

            if min_length > 0.0 {
                let ratio = max_length / min_length;
                if ratio > 100.0 {
                    warnings.push(ValidationWarning::IrregularSpacing { max_ratio: ratio });
                    suggestions.push(ValidationSuggestion::ApplySmoothing { iterations: 1 });
                }
            }
        }

        // Check for nearly collinear vertices
        self.check_collinear_vertices(polygon, warnings, suggestions);
    }

    fn check_collinear_vertices(
        &self,
        polygon: &Polygon,
        warnings: &mut Vec<ValidationWarning>,
        suggestions: &mut Vec<ValidationSuggestion>,
    ) {
        let vertices = polygon.vertices();
        let n = vertices.len();

        if n < 3 {
            return;
        }

        let mut collinear_groups = Vec::new();
        let mut current_group = Vec::new();

        for i in 0..n {
            let prev = vertices[(i + n - 1) % n];
            let curr = vertices[i];
            let next = vertices[(i + 1) % n];

            let v1 = curr - prev;
            let v2 = next - curr;

            let cross = v1.x * v2.y - v1.y * v2.x;

            if cross.abs() < self.tolerance * 10.0 {
                if current_group.is_empty() {
                    current_group.push(i.saturating_sub(1));
                }
                current_group.push(i);
            } else {
                if current_group.len() >= 3 {
                    collinear_groups.push(current_group.clone());
                }
                current_group.clear();
            }
        }

        if current_group.len() >= 3 {
            collinear_groups.push(current_group);
        }

        for group in collinear_groups {
            warnings.push(ValidationWarning::NearlyCollinearVertices {
                indices: group.clone(),
            });
            if group.len() > 3 {
                // Suggest removing middle vertices
                let middle_indices = group[1..group.len() - 1].to_vec();
                for &index in &middle_indices {
                    suggestions.push(ValidationSuggestion::RemoveVertex { index });
                }
            }
        }
    }

    fn validate_self_intersections(&self, polygon: &Polygon, errors: &mut Vec<ValidationError>) {
        let vertices = polygon.vertices();
        let n = vertices.len();

        if n < 4 {
            return;
        }

        // Check all edge pairs for intersections
        for i in 0..n {
            let i_next = (i + 1) % n;

            for j in (i + 2)..n {
                let j_next = (j + 1) % n;

                // Don't check adjacent edges
                if j_next == i {
                    continue;
                }

                if let Some(intersection) = self.line_intersection(
                    vertices[i],
                    vertices[i_next],
                    vertices[j],
                    vertices[j_next],
                ) {
                    errors.push(ValidationError::SelfIntersection {
                        point1: vertices[i],
                        point2: intersection,
                    });
                }
            }
        }
    }

    fn line_intersection(
        &self,
        p1: Point2D,
        p2: Point2D,
        p3: Point2D,
        p4: Point2D,
    ) -> Option<Point2D> {
        let d1 = p2 - p1;
        let d2 = p4 - p3;

        let denominator = d1.x * d2.y - d1.y * d2.x;

        if denominator.abs() < self.tolerance {
            return None; // Lines are parallel
        }

        let t = ((p3.x - p1.x) * d2.y - (p3.y - p1.y) * d2.x) / denominator;
        let u = ((p3.x - p1.x) * d1.y - (p3.y - p1.y) * d1.x) / denominator;

        // Check if intersection is within both line segments
        if t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0 {
            Some(p1 + d1 * t)
        } else {
            None
        }
    }

    fn validate_orientation(
        &self,
        polygon: &Polygon,
        warnings: &mut Vec<ValidationWarning>,
        suggestions: &mut Vec<ValidationSuggestion>,
    ) {
        use crate::math::geometry::polygon::PolygonProperties;

        let orientation = polygon.orientation();

        // Most applications expect counter-clockwise orientation
        if matches!(
            orientation,
            crate::math::geometry::polygon::core::properties::Orientation::Clockwise
        ) {
            warnings.push(ValidationWarning::PotentialNumericalInstability {
                reason: "Polygon has clockwise orientation, might be intended as counter-clockwise"
                    .to_string(),
            });
            suggestions.push(ValidationSuggestion::ChangeOrientation);
        }
    }

    fn validate_winding_consistency(
        &self,
        polygon: &Polygon,
        warnings: &mut Vec<ValidationWarning>,
    ) {
        // Check if the polygon has consistent winding by examining turn angles
        let vertices = polygon.vertices();
        let n = vertices.len();

        if n < 3 {
            return;
        }

        let mut total_turn = 0.0;

        for i in 0..n {
            let prev = vertices[(i + n - 1) % n];
            let curr = vertices[i];
            let next = vertices[(i + 1) % n];

            let v1 = curr - prev;
            let v2 = next - curr;

            if v1.length() > self.tolerance && v2.length() > self.tolerance {
                let turn_angle = v1.cross(v2).atan2(v1.dot(v2));
                total_turn += turn_angle;
            }
        }

        // For a simple polygon, total turn should be ±2π
        let expected_turn = if polygon.is_closed() { TAU } else { 0.0 };
        let turn_error = (total_turn.abs() - expected_turn).abs();

        if turn_error > self.tolerance * 100.0 {
            warnings.push(ValidationWarning::PotentialNumericalInstability {
                reason: format!(
                    "Inconsistent winding: total turn = {:.3}, expected = {:.3}",
                    total_turn, expected_turn
                ),
            });
        }
    }

    fn validate_geometric_properties(
        &self,
        polygon: &Polygon,
        warnings: &mut Vec<ValidationWarning>,
        suggestions: &mut Vec<ValidationSuggestion>,
    ) {
        use crate::math::geometry::polygon::PolygonProperties;

        // Check for very small area (might indicate degenerate polygon)
        let area = polygon.area();
        if area < self.tolerance * 100.0 {
            warnings.push(ValidationWarning::PotentialNumericalInstability {
                reason: format!("Very small polygon area: {:.2e}", area),
            });
            suggestions.push(ValidationSuggestion::SimplifyPolygon {
                target_vertex_count: 3,
            });
        }

        // Check perimeter to area ratio (elongated polygons might cause issues)
        let perimeter = polygon.perimeter();
        if perimeter > 0.0 && area > 0.0 {
            let compactness = 4.0 * PI * area / (perimeter * perimeter);
            if compactness < 0.01 {
                warnings.push(ValidationWarning::PotentialNumericalInstability {
                    reason: format!("Very elongated polygon (compactness = {:.3})", compactness),
                });
                suggestions.push(ValidationSuggestion::ApplySmoothing { iterations: 2 });
            }
        }
    }
}

/// Quick validation functions for common use cases
pub struct QuickValidation;

impl QuickValidation {
    /// Schnelle Prüfung ob Polygon grundlegend gültig ist
    pub fn is_valid_basic(polygon: &Polygon) -> bool {
        let validator = PolygonValidator::new(ValidationLevel::Basic);
        validator.validate(polygon).is_valid
    }

    /// Prüft auf Selbstüberschneidungen
    pub fn has_self_intersections(polygon: &Polygon) -> bool {
        let validator =
            PolygonValidator::new(ValidationLevel::Standard).check_self_intersections(true);
        let result = validator.validate(polygon);

        result
            .errors
            .iter()
            .any(|e| matches!(e, ValidationError::SelfIntersection { .. }))
    }

    /// Prüft auf degenerierte Kanten
    pub fn has_degenerate_edges(polygon: &Polygon, tolerance: f32) -> Vec<usize> {
        let validator = PolygonValidator::new(ValidationLevel::Standard).with_tolerance(tolerance);
        let result = validator.validate(polygon);

        result
            .errors
            .iter()
            .filter_map(|e| {
                if let ValidationError::DegenerateEdge { vertex_index } = e {
                    Some(*vertex_index)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Gibt Verbesserungsvorschläge für ein Polygon
    pub fn get_improvement_suggestions(polygon: &Polygon) -> Vec<ValidationSuggestion> {
        let validator = PolygonValidator::new(ValidationLevel::Comprehensive);
        validator.validate(polygon).suggestions
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::geometry::polygon::Polygon;

    #[test]
    fn test_basic_validation() {
        // Valid polygon
        let valid_polygon = Polygon::closed(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(0.0, 1.0),
        ])
        .unwrap();

        assert!(QuickValidation::is_valid_basic(&valid_polygon));

        // Invalid polygon (too few vertices)
        let invalid_polygon =
            Polygon::new(vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)]).unwrap();

        assert!(!QuickValidation::is_valid_basic(&invalid_polygon));
    }

    #[test]
    fn test_degenerate_edge_detection() {
        let polygon = Polygon::closed(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(1.0, 0.0), // Duplicate point
            Point2D::new(1.0, 1.0),
            Point2D::new(0.0, 1.0),
        ])
        .unwrap();

        let degenerate_edges = QuickValidation::has_degenerate_edges(&polygon, 1e-6);
        assert!(!degenerate_edges.is_empty());
    }

    #[test]
    fn test_self_intersection_detection() {
        // Figure-8 polygon (self-intersecting)
        let self_intersecting = Polygon::closed(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(0.0, 1.0),
        ])
        .unwrap();

        assert!(QuickValidation::has_self_intersections(&self_intersecting));

        // Simple square (no self-intersection)
        let simple_square = Polygon::closed(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(0.0, 1.0),
        ])
        .unwrap();

        assert!(!QuickValidation::has_self_intersections(&simple_square));
    }

    #[test]
    fn test_comprehensive_validation() {
        let polygon = Polygon::closed(vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(0.0, 1.0),
        ])
        .unwrap();

        let validator = PolygonValidator::new(ValidationLevel::Comprehensive);
        let result = validator.validate(&polygon);

        assert!(result.is_valid);
        assert!(result.errors.is_empty());
    }
}
