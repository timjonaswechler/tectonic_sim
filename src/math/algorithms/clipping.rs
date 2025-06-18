// src/math/algorithms/clipping.rs

//! # Polygon Clipping Module
//!
//! This module provides structures and algorithms for clipping polygons.
//! Clipping is the process of determining the portion of a polygon (the "subject")
//! that lies within another polygon (the "clipper") or a rectangular region.
//!
//! Currently, the primary supported algorithm is Sutherland-Hodgman, which is
//! suitable for convex clipper polygons.

use crate::math::{
    error::{MathError, MathResult},
    types::Bounds2D, // Custom Bounds2D struct using Vec2
    utils::constants,
};
use bevy::math::Vec2;

/// Specifies the algorithm to be used for polygon clipping operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClippingAlgorithm {
    /// The Sutherland-Hodgman algorithm.
    /// This algorithm clips a subject polygon against a convex clipper polygon.
    /// It processes the subject polygon against each edge of the clipper polygon sequentially.
    SutherlandHodgman,
    // WeilerAtherton, // More complex, handles arbitrary polygons (not implemented here)
    // LiangBarsky,    // Primarily for line clipping against rectangles (not for polygon-polygon here)
    // CohenSutherland // Primarily for line clipping against rectangles (not for polygon-polygon here)
}

impl Default for ClippingAlgorithm {
    /// Returns the default clipping algorithm.
    ///
    /// Default: `ClippingAlgorithm::SutherlandHodgman`
    fn default() -> Self {
        ClippingAlgorithm::SutherlandHodgman
    }
}

/// Performs polygon clipping operations using a specified algorithm.
pub struct PolygonClipper {
    algorithm: ClippingAlgorithm,
    tolerance: f32, // Tolerance for floating-point comparisons
}

impl Default for PolygonClipper {
    /// Creates a default `PolygonClipper` instance.
    ///
    /// Default values:
    /// - `algorithm`: `ClippingAlgorithm::SutherlandHodgman`
    /// - `tolerance`: `constants::EPSILON * 10.0`
    fn default() -> Self {
        Self {
            algorithm: ClippingAlgorithm::default(),
            tolerance: constants::EPSILON * 10.0,
        }
    }
}

impl PolygonClipper {
    /// Creates a new `PolygonClipper` with the specified algorithm.
    ///
    /// # Arguments
    /// * `algorithm` - The `ClippingAlgorithm` to use.
    pub fn new(algorithm: ClippingAlgorithm) -> Self {
        Self {
            algorithm,
            ..Default::default()
        }
    }

    /// Sets the tolerance for floating-point comparisons used in the clipping process.
    ///
    /// # Arguments
    /// * `tolerance` - The desired tolerance value. Must be non-negative.
    pub fn with_tolerance(mut self, tolerance: f32) -> Self {
        self.tolerance = tolerance.max(0.0);
        self
    }

    /// Clips a subject polygon (represented by its vertices) against a clipper polygon.
    ///
    /// Both polygons are assumed to have their vertices ordered (e.g., counter-clockwise for
    /// Sutherland-Hodgman, where the clipper must be convex).
    ///
    /// # Arguments
    /// * `subject_polygon_points` - Vertices of the polygon to be clipped.
    /// * `clipper_polygon_points` - Vertices of the clipping polygon.
    ///
    /// # Returns
    /// A `MathResult` containing a `Vec<Vec<Vec2>>`. Each inner `Vec<Vec2>` is a
    /// resulting clipped polygon. Sutherland-Hodgman typically returns a single polygon or an empty list.
    /// Returns an error if input polygons have fewer than 3 points.
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
                    Ok(vec![result_points]) // Sutherland-Hodgman produces one polygon
                } else {
                    Ok(Vec::new()) // No valid polygon resulted from clipping
                }
            } // Other algorithms could be implemented and matched here.
        }
    }

    /// Clips a subject polygon (represented by its vertices) against an axis-aligned rectangle (`Bounds2D`).
    ///
    /// # Arguments
    /// * `subject_polygon_points` - Vertices of the polygon to be clipped.
    /// * `clip_bounds` - The axis-aligned bounding box to clip against.
    ///
    /// # Returns
    /// A `MathResult` similar to `clip_polygon_points`.
    /// Returns an error if the subject polygon has fewer than 3 points or if `clip_bounds` is invalid.
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

        // Convert Bounds2D to a polygon vertex list (CCW order for Sutherland-Hodgman)
        let clipper_rect_points = [
            clip_bounds.min,                                 // Bottom-left
            Vec2::new(clip_bounds.max.x, clip_bounds.min.y), // Bottom-right
            clip_bounds.max,                                 // Top-right
            Vec2::new(clip_bounds.min.x, clip_bounds.max.y), // Top-left
        ];

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

    // --- Sutherland-Hodgman Algorithm Implementation ---

    /// Performs polygon clipping using the Sutherland-Hodgman algorithm.
    /// Assumes `clipper_points` form a convex polygon with vertices in CCW order.
    /// The subject polygon is clipped against each edge of the clipper polygon.
    ///
    /// # Arguments
    /// * `subject_points` - Vertices of the polygon to be clipped.
    /// * `clipper_points` - Vertices of the convex clipper polygon (CCW order).
    ///
    /// # Returns
    /// A `MathResult<Vec<Vec2>>` representing the vertices of the clipped polygon.
    /// The list might be empty if the subject is entirely outside the clipper.
    fn sutherland_hodgman_clip(
        &self,
        subject_points: &[Vec2],
        clipper_points: &[Vec2], // Must be convex and CCW
    ) -> MathResult<Vec<Vec2>> {
        let mut current_clipped_points = subject_points.to_vec();

        for i in 0..clipper_points.len() {
            // If at any stage the polygon degenerates or vanishes, stop.
            if current_clipped_points.len() < 3 && !current_clipped_points.is_empty() {
                current_clipped_points.clear();
                break;
            }
            if current_clipped_points.is_empty() {
                break;
            }

            let clip_edge_p1 = clipper_points[i];
            let clip_edge_p2 = clipper_points[(i + 1) % clipper_points.len()]; // Next vertex, wraps around

            let input_for_this_edge = current_clipped_points.clone();
            current_clipped_points.clear();

            if input_for_this_edge.is_empty() {
                // Should be caught by above checks
                continue;
            }

            // `s` is the start point of the current subject edge being tested
            let mut s = *input_for_this_edge.last().unwrap();

            for &e in &input_for_this_edge {
                // `e` is the end point of the current subject edge
                let s_is_inside = self.is_point_inside_clip_edge(s, clip_edge_p1, clip_edge_p2);
                let e_is_inside = self.is_point_inside_clip_edge(e, clip_edge_p1, clip_edge_p2);

                if e_is_inside {
                    // Case: Current vertex `e` is inside the clip edge
                    if !s_is_inside {
                        // Case: Edge (s,e) enters the clipping region
                        if let Some(intersection) =
                            self.line_segment_intersection(s, e, clip_edge_p1, clip_edge_p2)
                        {
                            current_clipped_points.push(intersection);
                        }
                        // else: Error or edge case in intersection/inside check
                    }
                    current_clipped_points.push(e); // Add `e` as it's inside
                } else if s_is_inside {
                    // Case: Edge (s,e) exits the clipping region
                    if let Some(intersection) =
                        self.line_segment_intersection(s, e, clip_edge_p1, clip_edge_p2)
                    {
                        current_clipped_points.push(intersection);
                    }
                    // else: Error or edge case
                }
                // Case: Both s and e are outside - do nothing
                s = e; // Advance to the next subject edge
            }
        }
        Ok(current_clipped_points)
    }

    /// Checks if `p_test` is on the "inside" of the directed clipping edge (clip_p1 -> clip_p2).
    /// Assumes the clipper polygon is CCW, so "inside" is to the left of the edge.
    ///
    /// # Arguments
    /// * `p_test` - The point to test.
    /// * `clip_p1` - Start point of the clipping edge.
    /// * `clip_p2` - End point of the clipping edge.
    ///
    /// # Returns
    /// `true` if `p_test` is inside or on the edge, `false` otherwise.
    #[inline]
    fn is_point_inside_clip_edge(&self, p_test: Vec2, clip_p1: Vec2, clip_p2: Vec2) -> bool {
        // Using the 2D cross product (or perpendicular dot product) to determine orientation.
        // (clip_p2 - clip_p1).perp_dot(p_test - clip_p1)
        // = (clip_p2.x - clip_p1.x) * (p_test.y - clip_p1.y) - (clip_p2.y - clip_p1.y) * (p_test.x - clip_p1.x)
        // Result > 0: p_test is to the left (inside for CCW clipper).
        // Result = 0: p_test is collinear.
        // Result < 0: p_test is to the right (outside).
        let cross_product = (clip_p2.x - clip_p1.x) * (p_test.y - clip_p1.y)
            - (clip_p2.y - clip_p1.y) * (p_test.x - clip_p1.x);
        cross_product >= -self.tolerance // Points on the edge (or very close due to tolerance) are considered inside.
    }

    /// Computes the intersection point of two line segments (p1-p2) and (p3-p4).
    ///
    /// # Arguments
    /// * `p1`, `p2` - Endpoints of the first line segment.
    /// * `p3`, `p4` - Endpoints of the second line segment.
    ///
    /// # Returns
    /// `Some(intersection_point)` if the segments intersect, `None` otherwise.
    fn line_segment_intersection(&self, p1: Vec2, p2: Vec2, p3: Vec2, p4: Vec2) -> Option<Vec2> {
        let d1 = p2 - p1; // Vector for segment p1-p2
        let d2 = p4 - p3; // Vector for segment p3-p4

        // Denominator for parametric line intersection formula
        let denominator = d1.x * d2.y - d1.y * d2.x; // d1.perp_dot(d2)

        if denominator.abs() < self.tolerance * constants::EPSILON {
            // Lines are parallel or collinear.
            // For robust clipping, collinear overlapping segments might need special handling.
            // Here, we treat them as non-intersecting for simplicity unless they are identical.
            return None;
        }

        // Numerators for parameters t and u
        // Intersection P = p1 + t * d1 = p3 + u * d2
        let t_num = (p3.x - p1.x) * d2.y - (p3.y - p1.y) * d2.x; // (p3 - p1).perp_dot(d2)
        let u_num = (p3.x - p1.x) * d1.y - (p3.y - p1.y) * d1.x; // (p3 - p1).perp_dot(d1)
        // Note: u_num should be (p1 - p3).perp_dot(d1) for standard u parameter, or use -u_num.
        // Let's use the common form: u = ((p1-p3) x d1) / (d2 x d1)
        // which is u = ((p1.x-p3.x)*d1.y - (p1.y-p3.y)*d1.x) / (-denominator)
        // The current u_num is (p3-p1).perp_dot(d1). With denom = d1.perp_dot(d2),
        // u = u_num / denominator is parameter for line p3 + u*d2 if u_num was (p1-p3).perp_dot(d1).
        // It's simpler to use established formulas or vector algebra carefully.
        // Let's use the t parameter for p1 + t*d1.
        let t = t_num / denominator;
        // For u, using the formula u = ((p1-p3) x d1) / (d2 x d1)
        let u = u_num / (-denominator); // d2.perp_dot(d1) = -denominator

        // Check if the intersection point lies on both segments (0 <= t <= 1 and 0 <= u <= 1).
        // Include tolerance for endpoints.
        let t_on_segment = t >= -self.tolerance && t <= 1.0 + self.tolerance;
        let u_on_segment = u >= -self.tolerance && u <= 1.0 + self.tolerance;

        if t_on_segment && u_on_segment {
            // Clamp t to [0,1] to avoid slight over/undershoots due to tolerance when calculating intersection point.
            Some(p1 + d1 * t.clamp(0.0, 1.0))
        } else {
            None // Intersection point is outside one or both segments.
        }
    }
}

/// Provides utility functions for common clipping operations.
pub struct ClippingOperations;

impl ClippingOperations {
    /// Clips a list of points (representing a polygon) against an axis-aligned rectangle.
    ///
    /// # Arguments
    /// * `subject_points` - Vertices of the polygon to be clipped.
    /// * `clip_bounds` - The `Bounds2D` rectangle to clip against.
    /// * `algorithm` - The `ClippingAlgorithm` to use.
    /// * `tolerance` - Optional floating point tolerance for the operation. Uses default if `None`.
    ///
    /// # Returns
    /// A `MathResult<Vec<Vec<Vec2>>>` containing the clipped polygon(s).
    pub fn clip_points_to_rectangle(
        subject_points: &[Vec2],
        clip_bounds: &Bounds2D,
        algorithm: ClippingAlgorithm,
        tolerance: Option<f32>,
    ) -> MathResult<Vec<Vec<Vec2>>> {
        let mut clipper = PolygonClipper::new(algorithm);
        if let Some(tol) = tolerance {
            clipper = clipper.with_tolerance(tol);
        }
        clipper.clip_polygon_against_rectangle(subject_points, clip_bounds)
    }

    // Other specific clipping utility functions (e.g., clip_to_circle) could be added here.
    // These would require their own logic or approximation of the clipper shape by a polygon.
}
