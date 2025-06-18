// src/math/algorithms/convex_hull.rs

//! # Convex Hull Algorithms Module
//!
//! This module provides implementations of various algorithms to compute the
//! convex hull of a set of 2D points. The convex hull is the smallest convex polygon
//! that encloses all the given points.
//!
//! ## Supported Algorithms:
//! - Graham Scan
//! - Andrew's Monotone Chain
//! - Jarvis March (Gift Wrapping)
//! - QuickHull
//!
//! The resulting hull points are typically returned in a counter-clockwise (CCW) order.

use crate::math::{
    error::{MathError, MathResult},
    utils::constants,
};
use bevy::math::Vec2;

/// Enumerates the available algorithms for computing the convex hull.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvexHullAlgorithm {
    /// Graham Scan algorithm. Time complexity: O(n log n).
    GrahamScan,
    /// Andrew's Monotone Chain algorithm. Time complexity: O(n log n). Often performs well.
    AndrewMonotone,
    /// Jarvis March (Gift Wrapping) algorithm. Time complexity: O(nh), where h is the number of hull points.
    /// Efficient when h is small.
    JarvisMarch,
    /// QuickHull algorithm. Average time complexity: O(n log n), worst-case: O(n^2).
    /// Similar to QuickSort.
    QuickHull,
}

impl Default for ConvexHullAlgorithm {
    /// Returns the default convex hull algorithm.
    ///
    /// Default: `ConvexHullAlgorithm::AndrewMonotone` (a good general-purpose choice).
    fn default() -> Self {
        ConvexHullAlgorithm::AndrewMonotone
    }
}

/// Computes the convex hull of a set of 2D points using a specified algorithm.
pub struct ConvexHullComputer {
    algorithm: ConvexHullAlgorithm,
    /// Whether to include collinear points that lie on the boundary of the hull.
    include_collinear_points: bool,
    /// Tolerance for floating-point comparisons, e.g., for determining collinearity.
    tolerance: f32,
}

impl Default for ConvexHullComputer {
    /// Creates a default `ConvexHullComputer` instance.
    ///
    /// Default values:
    /// - `algorithm`: `ConvexHullAlgorithm::AndrewMonotone`
    /// - `include_collinear_points`: `false`
    /// - `tolerance`: `constants::EPSILON * 10.0`
    fn default() -> Self {
        Self {
            algorithm: ConvexHullAlgorithm::default(),
            include_collinear_points: false,
            tolerance: constants::EPSILON * 10.0,
        }
    }
}

impl ConvexHullComputer {
    /// Creates a new `ConvexHullComputer` with the specified algorithm.
    ///
    /// # Arguments
    /// * `algorithm` - The `ConvexHullAlgorithm` to use.
    pub fn new(algorithm: ConvexHullAlgorithm) -> Self {
        Self {
            algorithm,
            ..Default::default()
        }
    }

    /// Configures whether to include collinear points on the hull's boundary.
    ///
    /// # Arguments
    /// * `include` - If `true`, collinear points on the hull edge will be included.
    pub fn include_collinear(mut self, include: bool) -> Self {
        self.include_collinear_points = include;
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

    /// Computes the convex hull of the given set of input points.
    ///
    /// The method first removes duplicate points (within tolerance) and sorts them
    /// if required by the chosen algorithm.
    ///
    /// # Arguments
    /// * `input_points` - A slice of `Vec2` points for which to compute the convex hull.
    ///
    /// # Returns
    /// A `MathResult<Vec<Vec2>>` containing the vertices of the convex hull in
    /// counter-clockwise (CCW) order. The list does not repeat the first point at the end.
    /// Returns an error if fewer than 3 unique points are provided after deduplication,
    /// or if a geometric failure occurs.
    pub fn compute_hull_points(&self, input_points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        if input_points.is_empty() {
            return Ok(Vec::new()); // Empty input yields an empty hull.
        }

        // Remove duplicate points and sort for algorithms like Andrew's Monotone Chain.
        let mut points = input_points.to_vec();
        points.sort_by(|a, b| {
            a.x.partial_cmp(&b.x)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
        });
        points.dedup_by(|a, b| a.distance_squared(*b) < self.tolerance * self.tolerance);

        if points.len() < 3 {
            // A convex hull is typically defined for 3 or more points.
            // For 1 or 2 points, the "hull" is the points themselves.
            // However, to form a polygon, at least 3 points are needed.
            // We return an error consistent with polygon formation requirements.
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: points.len(),
            });
        }

        let hull_points = match self.algorithm {
            ConvexHullAlgorithm::GrahamScan => self.graham_scan(&points)?,
            ConvexHullAlgorithm::AndrewMonotone => self.andrew_monotone(&points)?,
            ConvexHullAlgorithm::JarvisMarch => self.jarvis_march(&points)?,
            ConvexHullAlgorithm::QuickHull => self.quick_hull(&points)?,
        };

        if hull_points.len() < 3 {
            // This case should ideally not be reached if the input `points` had >= 3 unique points
            // and the algorithm is correct. It's a safeguard.
            return Err(MathError::GeometricFailure {
                operation: "Convex hull computation resulted in fewer than 3 unique hull points."
                    .to_string(),
            });
        }
        Ok(hull_points)
    }

    /// Computes the 2D cross product of vectors (q-p) and (r-q).
    /// Also known as orientation test for triplet (p, q, r).
    /// - Result > 0: (p,q,r) makes a counter-clockwise (left) turn.
    /// - Result < 0: (p,q,r) makes a clockwise (right) turn.
    /// - Result ≈ 0: p, q, r are collinear.
    #[inline]
    fn cross_product_orientation(p: Vec2, q: Vec2, r: Vec2) -> f32 {
        (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
    }

    /// Implements the Graham Scan algorithm for convex hull.
    /// `points` must be pre-processed (e.g., at least 3 unique points).
    fn graham_scan(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        let mut local_points = points.to_vec();
        let n = local_points.len();

        // 1. Find p0: the point with the smallest y-coordinate. If ties, pick the one with smallest x.
        let mut min_idx = 0;
        for i in 1..n {
            if local_points[i].y < local_points[min_idx].y
                || (local_points[i].y == local_points[min_idx].y // Using direct float compare, consider tolerance if needed
                    && local_points[i].x < local_points[min_idx].x)
            {
                min_idx = i;
            }
        }
        local_points.swap(0, min_idx);
        let p0 = local_points[0];

        // 2. Sort remaining points by polar angle with p0. If angles are equal, use distance.
        local_points[1..].sort_unstable_by(|pa, pb| {
            let orientation = Self::cross_product_orientation(p0, *pa, *pb);
            if orientation.abs() < self.tolerance {
                // Collinear with p0
                let dist_sq_a = p0.distance_squared(*pa);
                let dist_sq_b = p0.distance_squared(*pb);
                // If including collinear, sort by distance (closer first).
                // If not, sort by distance (further first to discard closer ones later, or handle in loop).
                // Standard Graham often keeps the farthest if collinear and on hull.
                if self.include_collinear_points {
                    dist_sq_a
                        .partial_cmp(&dist_sq_b)
                        .unwrap_or(std::cmp::Ordering::Equal)
                } else {
                    // Keep farthest if not including (effectively, only one will make it to hull)
                    dist_sq_b
                        .partial_cmp(&dist_sq_a)
                        .unwrap_or(std::cmp::Ordering::Equal)
                }
            } else if orientation > 0.0 {
                // pa is CCW relative to pb (from p0)
                std::cmp::Ordering::Less
            } else {
                // pa is CW
                std::cmp::Ordering::Greater
            }
        });

        // 3. Build the hull using a stack.
        let mut hull: Vec<Vec2> = Vec::with_capacity(n);
        hull.push(p0);
        if n > 1 {
            hull.push(local_points[1]);
        } // Second point always on hull initially

        for i in 2..n {
            // Keep popping from hull while the turn (hull.second_last, hull.last, local_points[i]) is not CCW.
            while hull.len() >= 2 {
                let top = hull.pop().unwrap(); // Temporarily remove current top
                let orientation = Self::cross_product_orientation(
                    *hull.last().unwrap(), // New top (second to last)
                    top,                   // Old top
                    local_points[i],       // Candidate point
                );

                if orientation > self.tolerance {
                    // Strict left turn (CCW)
                    hull.push(top); // Put old top back
                    break; // This forms a valid part of the hull
                } else if self.include_collinear_points && orientation.abs() < self.tolerance {
                    // Collinear and we want to include them.
                    hull.push(top);
                    break;
                }
                // Otherwise (right turn, or collinear and not including), `top` is popped. Loop continues.
            }
            hull.push(local_points[i]);
        }
        Ok(hull)
    }

    /// Implements Andrew's Monotone Chain algorithm.
    /// `points` must be sorted lexicographically (by x, then y) and contain at least 3 unique points.
    fn andrew_monotone(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        let n = points.len();
        let mut hull: Vec<Vec2> = Vec::with_capacity(2 * n); // Max possible size

        // Build lower hull
        for &p_i in points {
            while hull.len() >= 2 {
                let orientation = Self::cross_product_orientation(
                    hull[hull.len() - 2],
                    hull[hull.len() - 1],
                    p_i,
                );
                // Pop if not a left turn (or collinear and not including)
                if orientation < -self.tolerance
                    || (!self.include_collinear_points && orientation.abs() < self.tolerance)
                {
                    hull.pop();
                } else {
                    break;
                }
            }
            hull.push(p_i);
        }

        // Build upper hull (iterate in reverse, excluding the last point added to lower hull)
        // `t` stores the size of the lower hull before adding upper hull points
        let t = hull.len();
        // Iterate from the second to last point in the original sorted list, down to the first.
        for i in (0..n - 1).rev() {
            // Exclude points[n-1] as it's already the end of lower hull
            let p_i = points[i];
            while hull.len() > t {
                // Ensure we don't pop beyond the lower hull's connection point
                let orientation = Self::cross_product_orientation(
                    hull[hull.len() - 2],
                    hull[hull.len() - 1],
                    p_i,
                );
                if orientation < -self.tolerance
                    || (!self.include_collinear_points && orientation.abs() < self.tolerance)
                {
                    hull.pop();
                } else {
                    break;
                }
            }
            hull.push(p_i);
        }

        // The last point pushed for the upper hull is points[0], which is the start of the lower hull.
        // Remove this duplicate if the hull has more than one point.
        if hull.len() > 1 {
            hull.pop();
        }

        Ok(hull)
    }

    /// Implements Jarvis March (Gift Wrapping) algorithm.
    /// `points` must contain at least 3 unique points.
    fn jarvis_march(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        let n = points.len();

        // 1. Find the starting point: the one with the smallest y-coordinate (leftmost if ties).
        let mut start_idx = 0;
        for i in 1..n {
            if points[i].y < points[start_idx].y
                || (points[i].y == points[start_idx].y && points[i].x < points[start_idx].x)
            {
                start_idx = i;
            }
        }

        let mut hull_indices = Vec::new();
        let mut current_hull_idx = start_idx;
        loop {
            hull_indices.push(current_hull_idx);
            // Find the next point that makes the most "rightward" turn (or is furthest if collinear).
            // Initialize next_candidate_idx to a point different from current_hull_idx.
            let mut next_candidate_idx = (current_hull_idx + 1) % n;

            for i in 0..n {
                // Iterate through all points to find the next hull point
                if i == current_hull_idx {
                    continue;
                }

                let orientation = Self::cross_product_orientation(
                    points[current_hull_idx],   // Current hull point
                    points[next_candidate_idx], // Current best candidate
                    points[i],                  // Point being tested
                );

                // If points[i] is to the left of the line (current_hull -> next_candidate),
                // then points[i] is a better candidate.
                if orientation > self.tolerance {
                    // points[i] is more CCW
                    next_candidate_idx = i;
                } else if orientation.abs() < self.tolerance {
                    // Collinear
                    // If collinear, pick the one further away from current_hull_idx.
                    let dist_sq_i = points[current_hull_idx].distance_squared(points[i]);
                    let dist_sq_cand =
                        points[current_hull_idx].distance_squared(points[next_candidate_idx]);
                    if dist_sq_i > dist_sq_cand {
                        next_candidate_idx = i;
                    }
                }
            }
            current_hull_idx = next_candidate_idx;
            if current_hull_idx == start_idx {
                break;
            } // Wrapped around to the start point

            if hull_indices.len() > n {
                // Safety break for non-termination
                return Err(MathError::GeometricFailure {
                    operation: "Jarvis March did not terminate correctly.".to_string(),
                });
            }
        }
        Ok(hull_indices.into_iter().map(|idx| points[idx]).collect())
    }

    /// Implements QuickHull algorithm.
    /// `points` must be sorted lexicographically and contain at least 3 unique points.
    fn quick_hull(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        let n = points.len();

        // Initial line: from the leftmost (points[0]) to the rightmost (points[n-1]) point.
        let p_min_x = points[0];
        let p_max_x = points[n - 1];

        let mut hull = Vec::new();

        // Points to the left of directed line p_min_x -> p_max_x (upper hull candidates)
        let mut s1 = Vec::new();
        // Points to the left of directed line p_max_x -> p_min_x (lower hull candidates)
        let mut s2 = Vec::new();

        for i in 1..(n - 1) {
            // Exclude the initial min/max x points
            let p = points[i];
            let orientation = Self::cross_product_orientation(p_min_x, p_max_x, p);
            if orientation > self.tolerance {
                // p is to the left of p_min_x -> p_max_x
                s1.push(p);
            } else if orientation < -self.tolerance {
                // p is to the right of p_min_x -> p_max_x
                s2.push(p);
            } else if self.include_collinear_points {
                // If collinear and including, add to both if it's between p_min_x and p_max_x
                // This simplifies by letting the recursion handle them, but ensure they are only added once.
                // A robust way is to handle them by checking if they are on the segment.
                // For now, QuickHull often discards collinear points implicitly unless handled carefully.
                if p.x > p_min_x.x && p.x < p_max_x.x { // Basic check if between
                    // Could add to s1 and s2 and let recursion discard if not furthest.
                }
            }
        }

        hull.push(p_min_x);
        self.quick_hull_recursive(&mut hull, &s1, p_min_x, p_max_x);
        hull.push(p_max_x);
        self.quick_hull_recursive(&mut hull, &s2, p_max_x, p_min_x);

        // The recursive calls might add p_min_x again at the end if hull was empty initially.
        // Remove last if it's same as first (and hull has >1 elements)
        if hull.len() > 1 && hull.first() == hull.last() {
            hull.pop();
        }
        Ok(hull)
    }

    /// Recursive helper for QuickHull. Finds points on one side of the line (p_a, p_b).
    /// Adds points to `hull_segment` in CCW order.
    fn quick_hull_recursive(
        &self,
        hull_segment: &mut Vec<Vec2>, // Accumulates hull points for this segment
        points_set: &[Vec2],          // Candidate points on one side of line (p_a, p_b)
        p_a: Vec2,                    // Start point of the dividing line segment
        p_b: Vec2,                    // End point of the dividing line segment
    ) {
        if points_set.is_empty() {
            return; // No points on this side
        }

        // Find point p_c in points_set furthest from line (p_a, p_b) on the "left" side.
        let mut max_dist_signed = 0.0; // Signed distance (from cross product)
        let mut p_c_idx = usize::MAX;

        for (idx, &p_candidate) in points_set.iter().enumerate() {
            let current_signed_dist = Self::cross_product_orientation(p_a, p_b, p_candidate);
            if current_signed_dist > max_dist_signed + self.tolerance {
                max_dist_signed = current_signed_dist;
                p_c_idx = idx;
            } else if self.include_collinear_points && current_signed_dist.abs() < self.tolerance {
                // If collinear and including, pick the one "projected" furthest along (p_a,p_b)
                // or simply the one with largest actual distance if multiple collinear points exist.
                // For simplicity, if current_signed_dist is near 0, and we allow collinear,
                // we might prefer the one that extends the hull further.
                // A common strategy is to pick the one with the largest perpendicular distance,
                // and if tied, the one with the largest projection onto the perpendicular bisector.
                // Here, we just take the first one encountered or one that increases max_dist_signed (even if slightly).
                if p_c_idx == usize::MAX || // if no point found yet
                    p_a.distance_squared(p_candidate) > p_a.distance_squared(points_set[p_c_idx])
                // simplistic tie-break
                {
                    max_dist_signed = current_signed_dist; // will be close to 0
                    p_c_idx = idx;
                }
            }
        }

        if p_c_idx == usize::MAX {
            // No point found to the left of the line
            return;
        }
        let p_c = points_set[p_c_idx];

        // Points to the left of directed line (p_a, p_c)
        let mut s_ac = Vec::new();
        for &p_s in points_set {
            if Self::cross_product_orientation(p_a, p_c, p_s) > self.tolerance {
                s_ac.push(p_s);
            } // Collinear points on (p_a, p_c) are handled if include_collinear_points (by not being > tolerance)
        }
        self.quick_hull_recursive(hull_segment, &s_ac, p_a, p_c);

        hull_segment.push(p_c); // p_c is part of the hull

        // Points to the left of directed line (p_c, p_b)
        let mut s_cb = Vec::new();
        for &p_s in points_set {
            if Self::cross_product_orientation(p_c, p_b, p_s) > self.tolerance {
                s_cb.push(p_s);
            }
        }
        self.quick_hull_recursive(hull_segment, &s_cb, p_c, p_b);
    }
}
