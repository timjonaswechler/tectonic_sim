// src/math/algorithms/metaballs/influence/polygon_source.rs
// Alter Pfad: .../geometry/metaballs/influence/polygon.rs

use crate::math::algorithms::metaballs::influence::FieldInfluence; // Pfad angepasst
use crate::math::{
    geometry::polygon::{core::Polygon, properties::PolygonProperties}, // Pfad für Polygon evtl. prüfen/anpassen
    types::Bounds2D,
    utils::{constants, numerical},
};
use bevy::math::Vec2;

// ... (Rest des Codes von PolygonInfluence und PolygonInfluenceBuilder wie oben) ...
/// Definiert, wie der Einfluss eines Polygons mit dem Abstand abnimmt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PolygonFalloff {
    Linear,
    /// Normalisierter quadratischer Abfall: `strength * (1 - normalized_distance)^2`. Erreicht 0 am `falloff_distance`.
    QuadraticSmooth,
    /// Algebraischer Abfall, ähnlich Metaballs: `strength / ( (scaled_distance)^2 + 1)`.
    AlgebraicInverseSq {
        /// Faktor, durch den `distance_abs` geteilt wird, bevor es quadriert wird.
        distance_scale: f32,
    },
    Exponential {
        /// Steuert, wie schnell der Einfluss abfällt. Größere Werte = schnellerer Abfall.
        steepness: f32,
    },
    Gaussian {
        /// Beeinflusst die "Breite" der Gausskurve relativ zur `falloff_distance`.
        variance_factor: f32,
    },
    /// Harter Übergang. Wert ist `strength` bis zur Hälfte der `falloff_distance`, danach 0.
    Step, // Feld midpoint_factor entfernt
    Sigmoid {
        /// Steuert die Steilheit des Übergangs. Größere Werte = steilerer Übergang.
        steepness: f32,
        /// Relativer Punkt (0.0 bis 1.0) innerhalb der `falloff_distance`, um den der Übergang zentriert ist.
        midpoint_factor: f32,
    },
}

impl Default for PolygonFalloff {
    fn default() -> Self {
        PolygonFalloff::Linear
    }
}

/// Definiert einen Einflussbereich, der von der Form eines Polygons ausgeht.
#[derive(Clone)]
pub struct PolygonInfluence {
    pub polygon: Polygon, // Annahme: Pfad zu Polygon ist crate::math::geometry::polygon::core::Polygon
    pub strength: f32,
    pub falloff_distance: f32,
    pub falloff_type: PolygonFalloff,
    pub interior_only: bool,
    pub offset: f32,
    cached_bounds_expanded: Option<Bounds2D>,
}

impl PolygonInfluence {
    #[allow(clippy::too_many_arguments)]
    fn new(
        polygon: Polygon,
        strength: f32,
        falloff_distance: f32,
        falloff_type: PolygonFalloff,
        interior_only: bool,
        offset: f32,
    ) -> Self {
        let mut new_self = Self {
            polygon,
            strength,
            falloff_distance: falloff_distance.max(constants::EPSILON),
            falloff_type,
            interior_only,
            offset,
            cached_bounds_expanded: None,
        };
        new_self.recalculate_cached_bounds();
        new_self
    }

    pub fn recalculate_cached_bounds(&mut self) {
        if let Some(poly_bounds) = self.polygon.bounds() {
            self.cached_bounds_expanded = Some(poly_bounds.expand(self.falloff_distance));
        } else {
            self.cached_bounds_expanded = None;
        }
    }

    fn signed_distance_to_polygon(&self, point: Vec2) -> f32 {
        if self.polygon.vertices().is_empty() {
            return f32::INFINITY;
        }
        let distance_to_boundary = self.distance_to_boundary(point);

        if self.polygon.contains_point(point) {
            -distance_to_boundary
        } else {
            distance_to_boundary
        }
    }

    fn distance_to_boundary(&self, point: Vec2) -> f32 {
        let vertices = self.polygon.vertices();
        let n = vertices.len();

        if n == 0 {
            return f32::INFINITY;
        }
        if n == 1 {
            return point.distance(vertices[0]);
        }

        let mut min_dist_sq = f32::INFINITY;
        let num_segments = if n > 1 && vertices.first() == vertices.last() {
            n - 1
        } else {
            n
        };

        if num_segments == 0 {
            return if n == 1 {
                point.distance(vertices[0])
            } else {
                f32::INFINITY
            };
        }

        for i in 0..num_segments {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % n];
            min_dist_sq = min_dist_sq.min(point_to_line_segment_distance_sq(point, p1, p2));
        }
        min_dist_sq.sqrt()
    }

    fn apply_falloff(&self, distance_abs: f32) -> f32 {
        if distance_abs >= self.falloff_distance - constants::EPSILON {
            return 0.0;
        }
        if self.falloff_distance < constants::EPSILON {
            return if distance_abs < constants::EPSILON {
                self.strength
            } else {
                0.0
            };
        }

        let normalized_distance = (distance_abs / self.falloff_distance).clamp(0.0, 1.0);

        match self.falloff_type {
            PolygonFalloff::Linear => self.strength * (1.0 - normalized_distance),
            PolygonFalloff::QuadraticSmooth => self.strength * (1.0 - normalized_distance).powi(2),
            PolygonFalloff::AlgebraicInverseSq { distance_scale } => {
                let effective_dist = distance_abs / distance_scale.max(constants::EPSILON);
                self.strength / (effective_dist * effective_dist + 1.0)
            }
            PolygonFalloff::Exponential { steepness } => {
                self.strength * (-normalized_distance * steepness.max(0.01)).exp()
            }
            PolygonFalloff::Gaussian { variance_factor } => {
                let sigma = self.falloff_distance * variance_factor.max(0.01);
                self.strength * (-distance_abs.powi(2) / (2.0 * sigma.powi(2))).exp()
            }
            PolygonFalloff::Step => {
                if normalized_distance <= 0.5 {
                    self.strength
                } else {
                    0.0
                }
            }
            PolygonFalloff::Sigmoid {
                steepness,
                midpoint_factor,
            } => {
                let midpoint_abs = self.falloff_distance * midpoint_factor.clamp(0.0, 1.0);
                let x_val = (distance_abs - midpoint_abs)
                    * (steepness.max(0.01) / self.falloff_distance.max(constants::EPSILON));
                self.strength * (1.0 - numerical::sigmoid(x_val))
            }
        }
    }
}

fn point_to_line_segment_distance_sq(point: Vec2, seg_a: Vec2, seg_b: Vec2) -> f32 {
    let l2 = seg_a.distance_squared(seg_b);
    if l2 < constants::EPSILON_SQUARED {
        return point.distance_squared(seg_a);
    }
    let t = ((point - seg_a).dot(seg_b - seg_a) / l2).clamp(0.0, 1.0);
    let projection = seg_a + t * (seg_b - seg_a);
    point.distance_squared(projection)
}

impl std::fmt::Debug for PolygonInfluence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PolygonInfluence")
            .field("polygon_vertices", &self.polygon.vertices().len())
            .field("strength", &self.strength)
            .field("falloff_distance", &self.falloff_distance)
            .field("falloff_type", &self.falloff_type)
            .field("interior_only", &self.interior_only)
            .field("offset", &self.offset)
            .field(
                "cached_bounds_expanded",
                &self.cached_bounds_expanded.is_some(),
            )
            .finish()
    }
}

impl FieldInfluence for PolygonInfluence {
    fn influence_at(&self, point: Vec2) -> f32 {
        if let Some(bounds) = self.cached_bounds_expanded {
            if !bounds.contains_point(point) {
                return 0.0;
            }
        }

        let signed_dist = self.signed_distance_to_polygon(point);

        if signed_dist == f32::INFINITY {
            return 0.0;
        }

        if self.interior_only && signed_dist > constants::EPSILON {
            return 0.0;
        }

        let mut influence_val = if self.interior_only {
            let dist_from_boundary_abs = signed_dist.abs();
            if dist_from_boundary_abs <= self.falloff_distance {
                self.apply_falloff(dist_from_boundary_abs)
            } else {
                self.strength
            }
        } else {
            self.apply_falloff(signed_dist.abs())
        };

        if influence_val > constants::EPSILON
            || (self.interior_only && signed_dist <= constants::EPSILON)
        {
            influence_val += self.offset;
        }

        influence_val.max(0.0)
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        self.cached_bounds_expanded
    }

    fn max_influence_distance(&self) -> Option<f32> {
        if let Some(bounds) = self.cached_bounds_expanded {
            let half_size = bounds.size() * 0.5;
            Some(half_size.length() + bounds.center().length())
        } else {
            Some(self.falloff_distance)
        }
    }

    fn influence_type(&self) -> &'static str {
        "PolygonInfluence"
    }
}

pub struct PolygonInfluenceBuilder {
    polygon: Option<Polygon>,
    strength: f32,
    falloff_distance: f32,
    falloff_type: PolygonFalloff,
    interior_only: bool,
    offset: f32,
}

impl Default for PolygonInfluenceBuilder {
    fn default() -> Self {
        Self {
            polygon: None,
            strength: 1.0,
            falloff_distance: 10.0,
            falloff_type: PolygonFalloff::default(),
            interior_only: false,
            offset: 0.0,
        }
    }
}

impl PolygonInfluenceBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn polygon(mut self, polygon: Polygon) -> Self {
        self.polygon = Some(polygon);
        self
    }

    pub fn vertices(mut self, vertices: Vec<Vec2>) -> Self {
        match Polygon::closed(vertices) {
            Ok(polygon) => self.polygon = Some(polygon),
            Err(e) => {
                eprintln!(
                    "Failed to create polygon from vertices for PolygonInfluenceBuilder: {:?}. Polygon will be None.",
                    e
                );
                self.polygon = None;
            }
        }
        self
    }

    pub fn strength(mut self, strength: f32) -> Self {
        self.strength = strength;
        self
    }

    pub fn falloff_distance(mut self, distance: f32) -> Self {
        self.falloff_distance = distance.max(constants::EPSILON);
        self
    }

    pub fn falloff_type(mut self, falloff_type: PolygonFalloff) -> Self {
        self.falloff_type = falloff_type;
        self
    }

    pub fn interior_only(mut self, interior_only: bool) -> Self {
        self.interior_only = interior_only;
        self
    }

    pub fn offset(mut self, offset: f32) -> Self {
        self.offset = offset;
        self
    }

    pub fn build(self) -> Option<PolygonInfluence> {
        self.polygon.map(|p| {
            PolygonInfluence::new(
                p,
                self.strength,
                self.falloff_distance,
                self.falloff_type,
                self.interior_only,
                self.offset,
            )
        })
    }
}
