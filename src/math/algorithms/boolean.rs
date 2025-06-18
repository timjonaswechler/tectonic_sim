// src/math/algorithms/boolean.rs

//! # Polygon Boolean Operations Module
//!
//! This module provides functionalities to perform boolean (set-theoretic) operations
//! on polygons, such as union, intersection, difference, and XOR.
//! Polygons are represented as slices of 2D points (`&[Vec2]`).
//!
//! The current implementation relies on clipping algorithms for some operations and
//! may have simplifications or placeholders for more complex robust algorithms.

use crate::math::{
    algorithms::clipping::{ClippingAlgorithm, PolygonClipper},
    error::{MathError, MathResult},
    types::Bounds2D, // Our Bounds2D using Vec2
    utils::constants,
};
use bevy::math::Vec2;
// For this module to operate on polygons, it ideally works with point lists.
// Conversion to/from a dedicated `Polygon` struct should happen externally.
// Helper functions operating on `&[Vec2]` are used here.

/// Defines the type of boolean operation to be performed on polygons.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOpType {
    /// Union of two polygons (A ∪ B). Combines both polygons into one.
    Union,
    /// Intersection of two polygons (A ∩ B). Results in the common area.
    Intersection,
    /// Difference between two polygons (A - B). Subtracts polygon B from A.
    Difference,
    /// Symmetrical Difference (XOR) of two polygons (A ⊕ B).
    /// Results in the areas that are in A or B, but not in both.
    Xor,
}

/// Performs boolean operations (union, intersection, difference, XOR) on polygons
/// represented as lists of 2D points.
///
/// Note: The robustness of these operations, especially Union and Difference,
/// can be highly dependent on the complexity of the input polygons and the
/// underlying algorithms. The current version may use simplified approaches.
pub struct PolygonBooleanOps {
    clipping_algorithm: ClippingAlgorithm, // Algorithm used for intersection and parts of other ops
    tolerance: f32,                        // Tolerance for floating-point comparisons
}

impl Default for PolygonBooleanOps {
    /// Creates a default `PolygonBooleanOps` instance.
    ///
    /// Default values:
    /// - `clipping_algorithm`: `ClippingAlgorithm::SutherlandHodgman`
    /// - `tolerance`: `constants::EPSILON * 10.0`
    fn default() -> Self {
        Self {
            clipping_algorithm: ClippingAlgorithm::SutherlandHodgman,
            tolerance: constants::EPSILON * 10.0,
        }
    }
}

impl PolygonBooleanOps {
    /// Creates a new `PolygonBooleanOps` instance with default settings.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the clipping algorithm to be used internally.
    ///
    /// # Arguments
    /// * `algo` - The `ClippingAlgorithm` to use.
    pub fn with_clipping_algorithm(mut self, algo: ClippingAlgorithm) -> Self {
        self.clipping_algorithm = algo;
        self
    }

    /// Sets the tolerance for floating-point comparisons.
    ///
    /// # Arguments
    /// * `tolerance` - The desired tolerance value. Must be non-negative.
    pub fn with_tolerance(mut self, tolerance: f32) -> Self {
        self.tolerance = tolerance.max(0.0);
        self
    }

    /// Executes the specified boolean operation between two polygons.
    ///
    /// `poly_a_points` and `poly_b_points` are slices of `Vec2` representing the vertices
    /// of the polygons, assumed to be ordered (e.g., counter-clockwise).
    ///
    /// # Arguments
    /// * `poly_a_points` - Vertices of the first polygon (subject polygon).
    /// * `poly_b_points` - Vertices of the second polygon (clip polygon or operand).
    /// * `operation` - The `BooleanOpType` to perform.
    ///
    /// # Returns
    /// A `MathResult` containing a `Vec<Vec<Vec2>>`, where each inner `Vec<Vec2>`
    /// represents a resulting polygon. Boolean operations can result in multiple
    /// disjoint polygons. Returns an error if input polygons have fewer than 3 points.
    pub fn execute(
        &self,
        poly_a_points: &[Vec2],
        poly_b_points: &[Vec2],
        operation: BooleanOpType,
    ) -> MathResult<Vec<Vec<Vec2>>> {
        if poly_a_points.len() < 3 || poly_b_points.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: poly_a_points.len().min(poly_b_points.len()),
            });
        }
        // TODO: Validation for simple polygons (no self-intersections) would be beneficial here.

        match operation {
            BooleanOpType::Intersection => self.intersection(poly_a_points, poly_b_points),
            BooleanOpType::Union => self.union(poly_a_points, poly_b_points),
            BooleanOpType::Difference => self.difference(poly_a_points, poly_b_points),
            BooleanOpType::Xor => self.xor(poly_a_points, poly_b_points),
        }
    }

    /// Computes the intersection of Polygon A and Polygon B (A ∩ B).
    ///
    /// This implementation primarily uses a clipping approach. If the `clipping_algorithm`
    /// is Sutherland-Hodgman, `poly_b_points` (the clipper) is assumed to be convex.
    /// For non-convex clippers with Sutherland-Hodgman, the result may be incorrect.
    ///
    /// # Returns
    /// A `MathResult` containing a list of polygons (usually one, or none if no intersection).
    fn intersection(
        &self,
        poly_a_points: &[Vec2],
        poly_b_points: &[Vec2],
    ) -> MathResult<Vec<Vec<Vec2>>> {
        let clipper = PolygonClipper::new(self.clipping_algorithm).with_tolerance(self.tolerance);
        let result_a_clipped_by_b = clipper.clip_polygon_points(poly_a_points, poly_b_points)?;
        Ok(result_a_clipped_by_b)
    }

    /// Computes the union of Polygon A and Polygon B (A ∪ B).
    ///
    /// **Note:** This is a highly simplified and often incorrect placeholder implementation.
    /// Robust polygon union is a complex algorithm (e.g., Greiner-Hormann, Vatti).
    /// The current version may return the convex hull of all points or handle only
    /// trivial cases like containment or disjointness.
    ///
    /// # Returns
    /// A `MathResult` containing a list of polygons.
    fn union(&self, poly_a_points: &[Vec2], poly_b_points: &[Vec2]) -> MathResult<Vec<Vec<Vec2>>> {
        let bounds_a = Bounds2D::from_points_iter(poly_a_points.iter().copied());
        let bounds_b = Bounds2D::from_points_iter(poly_b_points.iter().copied());

        if let (Some(ba), Some(bb)) = (bounds_a, bounds_b) {
            if !ba.intersects(&bb) {
                // Quick AABB test
                // Further check for actual segment intersections if AABBs don't intersect
                if find_line_segment_intersections(poly_a_points, poly_b_points, self.tolerance)
                    .is_empty()
                {
                    return Ok(vec![poly_a_points.to_vec(), poly_b_points.to_vec()]); // Disjoint
                }
            }
        }

        if polygon_contains_all_points(poly_a_points, poly_b_points, self.tolerance) {
            return Ok(vec![poly_a_points.to_vec()]); // A contains B
        }
        if polygon_contains_all_points(poly_b_points, poly_a_points, self.tolerance) {
            return Ok(vec![poly_b_points.to_vec()]); // B contains A
        }

        // Placeholder: Compute convex hull of all points. This is generally NOT a correct union.
        // TODO: Replace with a robust union algorithm.
        let mut all_points = poly_a_points.to_vec();
        all_points.extend_from_slice(poly_b_points);

        let hull_computer = super::convex_hull::ConvexHullComputer::new(
            super::convex_hull::ConvexHullAlgorithm::AndrewMonotone,
        )
        .with_tolerance(self.tolerance);
        let hull = hull_computer.compute_hull_points(&all_points)?;
        Ok(vec![hull])
    }

    /// Computes the difference of Polygon A and Polygon B (A - B).
    /// This results in the parts of Polygon A that are not in Polygon B.
    ///
    /// **Note:** This implementation is a simplified approach, potentially using
    /// clipping with an inverted clipper. It is generally correct only if `poly_b_points`
    /// (the subtrahend) is convex. Robust difference for arbitrary polygons is complex.
    ///
    /// # Returns
    /// A `MathResult` containing a list of polygons. The result can be multiple disjoint polygons.
    fn difference(
        &self,
        poly_a_points: &[Vec2], // The polygon from which to subtract
        poly_b_points: &[Vec2], // The polygon to subtract
    ) -> MathResult<Vec<Vec<Vec2>>> {
        if polygon_contains_all_points(poly_b_points, poly_a_points, self.tolerance) {
            // B contains A: A - B = {} (empty set)
            return Ok(Vec::new());
        }

        // Simplified approach: Clip A against the "reversed" B.
        // This approximates A INTERSECT (COMPLEMENT B) if B is convex.
        let clipper = PolygonClipper::new(self.clipping_algorithm).with_tolerance(self.tolerance);
        let mut reversed_b_points = poly_b_points.to_vec();
        reversed_b_points.reverse(); // Reversing order changes "inside" for Sutherland-Hodgman

        let result = clipper.clip_polygon_points(poly_a_points, &reversed_b_points)?;
        Ok(result)
        // IMPORTANT: This is a simplified difference, relying on specific behavior of S-H with reversed clipper.
        // Not robust for general concave polygons.
        // TODO: Replace with a robust difference algorithm (e.g., Weiler-Atherton, Greiner-Hormann).
    }

    /// Computes the symmetrical difference (XOR) of Polygon A and Polygon B.
    /// A ⊕ B = (A - B) ∪ (B - A).
    ///
    /// This operation relies on the `difference` and `union` operations.
    /// Its correctness is therefore dependent on their correctness.
    ///
    /// # Returns
    /// A `MathResult` containing a list of polygons.
    fn xor(&self, poly_a_points: &[Vec2], poly_b_points: &[Vec2]) -> MathResult<Vec<Vec<Vec2>>> {
        // A XOR B = (A \ B) U (B \ A)
        let diff_a_minus_b = self.difference(poly_a_points, poly_b_points)?;
        let diff_b_minus_a = self.difference(poly_b_points, poly_a_points)?;

        // The union of these two differences.
        // This simple concatenation is only correct if (A-B) and (B-A) are disjoint
        // and their union doesn't require merging. A proper union of the results might be needed.
        let mut result_polygons = diff_a_minus_b;
        result_polygons.extend(diff_b_minus_a);
        // TODO: Potentially merge/clean up result_polygons if they overlap/touch.
        // For now, we assume they form the correct set of disjoint polygons for XOR.
        Ok(result_polygons)
    }
}

// --- Helper functions for polygon point lists ---

/// Checks if a point is inside a polygon using the ray casting algorithm.
///
/// # Arguments
/// * `point` - The `Vec2` point to test.
/// * `poly_points` - A slice of `Vec2` representing the polygon's vertices in order.
///   The polygon should be simple (not self-intersecting).
/// * `tolerance` - Floating point tolerance for comparisons.
///
/// # Returns
/// `true` if the point is inside the polygon, `false` otherwise.
fn point_in_polygon(point: Vec2, poly_points: &[Vec2], tolerance: f32) -> bool {
    let n = poly_points.len();
    if n < 3 {
        return false; // Not a valid polygon
    }

    let mut inside = false;
    let mut p1 = poly_points[n - 1]; // Start with the last vertex to form the first edge with poly_points[0]

    for i in 0..n {
        let p2 = poly_points[i];
        // Check if the horizontal ray from `point` intersects the edge (p1, p2)
        if (p1.y > point.y) != (p2.y > point.y) {
            // Edge crosses the ray's y-level
            // Compute the x-coordinate of the intersection of the ray and the line extending the edge
            if (p2.y - p1.y).abs() > tolerance {
                // Avoid division by zero for horizontal edges
                let x_intersection = (point.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;
                if x_intersection > point.x {
                    // Intersection is to the right of the point
                    inside = !inside;
                }
            } else if p1.y == point.y { // Horizontal edge at the same y-level as the point
                // This case requires careful handling for robustness, often ignored or handled by specific rules.
                // If the point is on a horizontal segment, standard ray casting might be ambiguous.
                // For simplicity, we only count clear crossings of non-horizontal segments.
            }
        }
        p1 = p2; // Move to the next edge
    }
    inside
}

/// Checks if all points of `inner_poly_points` are contained within `outer_poly_points`.
///
/// # Arguments
/// * `outer_poly_points` - Vertices of the containing polygon.
/// * `inner_poly_points` - Vertices of the polygon to check for containment.
/// * `tolerance` - Floating point tolerance for `point_in_polygon` checks.
///
/// # Returns
/// `true` if all points of the inner polygon are inside the outer polygon, `false` otherwise.
fn polygon_contains_all_points(
    outer_poly_points: &[Vec2],
    inner_poly_points: &[Vec2],
    tolerance: f32,
) -> bool {
    if outer_poly_points.len() < 3 {
        return false; // Outer polygon is not valid
    }
    if inner_poly_points.is_empty() {
        return true; // An empty inner polygon is considered contained
    }

    for &inner_point in inner_poly_points {
        if !point_in_polygon(inner_point, outer_poly_points, tolerance) {
            return false; // Found a point of inner_poly not in outer_poly
        }
    }
    true // All points of inner_poly are in outer_poly
}

/// Finds all intersection points between the edges of two polygons.
///
/// # Arguments
/// * `poly_a` - Vertices of the first polygon.
/// * `poly_b` - Vertices of the second polygon.
/// * `tolerance` - Floating point tolerance for intersection and deduplication.
///
/// # Returns
/// A `Vec<Vec2>` of unique intersection points.
fn find_line_segment_intersections(poly_a: &[Vec2], poly_b: &[Vec2], tolerance: f32) -> Vec<Vec2> {
    let mut intersections = Vec::new();
    let n_a = poly_a.len();
    let n_b = poly_b.len();

    if n_a < 2 || n_b < 2 {
        // Need at least one segment
        return intersections;
    }

    for i in 0..n_a {
        let p1a = poly_a[i];
        let p2a = poly_a[(i + 1) % n_a]; // Loop back for the last segment
        for j in 0..n_b {
            let p1b = poly_b[j];
            let p2b = poly_b[(j + 1) % n_b]; // Loop back for the last segment

            // Line segment intersection logic (from PolygonClipper or similar)
            let d1 = p2a - p1a;
            let d2 = p2b - p1b;
            let denominator = d1.x * d2.y - d1.y * d2.x; // d1.perp_dot(d2)

            if denominator.abs() > tolerance * constants::EPSILON {
                // Lines are not parallel
                let t_num = (p1b - p1a).perp_dot(d2); // (p1b.x - p1a.x) * d2.y - (p1b.y - p1a.y) * d2.x
                let u_num = (p1b - p1a).perp_dot(d1); // (p1b.x - p1a.x) * d1.y - (p1b.y - p1a.y) * d1.x

                let t = t_num / denominator;
                let u = u_num / denominator;

                // Check if intersection point lies on both segments
                if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
                    intersections.push(p1a + d1 * t);
                }
            }
            // TODO: Handle collinear overlapping segments if necessary.
        }
    }

    // Remove duplicate intersection points
    if !intersections.is_empty() {
        intersections.sort_by(|a, b| {
            a.x.partial_cmp(&b.x)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
        });
        intersections.dedup_by(|a, b| a.distance_squared(*b) < tolerance * tolerance);
    }
    intersections
}
