// src/math/geometry/polygon/operations/convex_hull.rs
use super::super::Polygon;
use crate::math::{error::*, utils::*};
use bevy::math::Vec2;

/// Verschiedene Algorithmen für die konvexe Hülle
#[derive(Debug, Clone, Copy)]
pub enum ConvexHullAlgorithm {
    /// Graham Scan (O(n log n))
    GrahamScan,
    /// Andrew's Monotone Chain (O(n log n))
    AndrewMonotone,
    /// Jarvis March / Gift Wrapping (O(nh), h = Anzahl Hull-Punkte)
    JarvisMarch,
    /// QuickHull (O(n log n) average, O(n²) worst case)
    QuickHull,
}

/// Konvexe Hülle Computer
pub struct ConvexHullComputer {
    algorithm: ConvexHullAlgorithm,
    include_collinear: bool,
    tolerance: f32,
}

impl ConvexHullComputer {
    /// Erstellt einen neuen ConvexHull-Computer
    pub fn new(algorithm: ConvexHullAlgorithm) -> Self {
        Self {
            algorithm,
            include_collinear: false,
            tolerance: constants::EPSILON * 1000.0,
        }
    }

    /// Setzt ob kollineare Punkte inkludiert werden sollen
    pub fn include_collinear(mut self, include: bool) -> Self {
        self.include_collinear = include;
        self
    }

    /// Setzt die Toleranz für numerische Vergleiche
    pub fn with_tolerance(mut self, tolerance: f32) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Berechnet die konvexe Hülle von Punkten
    pub fn compute_hull(&self, points: &[Vec2]) -> MathResult<Polygon> {
        if points.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: points.len(),
            });
        }

        let hull_points = match self.algorithm {
            ConvexHullAlgorithm::GrahamScan => self.graham_scan(points)?,
            ConvexHullAlgorithm::AndrewMonotone => self.andrew_monotone(points)?,
            ConvexHullAlgorithm::JarvisMarch => self.jarvis_march(points)?,
            ConvexHullAlgorithm::QuickHull => self.quickhull(points)?,
        };

        if hull_points.len() < 3 {
            return Err(MathError::GeometricFailure {
                operation: "Convex hull resulted in fewer than 3 points".to_string(),
            });
        }

        Polygon::closed(hull_points)
    }

    /// Berechnet die konvexe Hülle eines Polygons
    pub fn compute_polygon_hull(&self, polygon: &Polygon) -> MathResult<Polygon> {
        self.compute_hull(polygon.vertices())
    }

    // === Algorithmus-Implementierungen ===

    /// Graham Scan Algorithmus
    fn graham_scan(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        let mut points = points.to_vec();

        // 1. Finde den Punkt mit der niedrigsten Y-Koordinate (bei Gleichstand: niedrigste X)
        let start_idx = points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                a.y.partial_cmp(&b.y)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal))
            })
            .map(|(i, _)| i)
            .unwrap();

        points.swap(0, start_idx);
        let start_point = points[0];

        // 2. Sortiere Punkte nach Polarwinkel
        let mut remaining_points = points[1..].to_vec();
        remaining_points.sort_by(|a, b| {
            let cross = Self::cross_product(start_point, *a, *b);

            if cross.abs() < self.tolerance {
                // Kollineare Punkte: sortiere nach Abstand
                let dist_a = start_point.distance_squared(*a);
                let dist_b = start_point.distance_squared(*b);
                dist_a
                    .partial_cmp(&dist_b)
                    .unwrap_or(std::cmp::Ordering::Equal)
            } else if cross > 0.0 {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            }
        });

        // 3. Graham Scan
        let mut hull = vec![start_point];

        for &point in &remaining_points {
            // Entferne Punkte die eine Rechts-Kurve machen
            while hull.len() > 1 {
                let cross = Self::cross_product(hull[hull.len() - 2], hull[hull.len() - 1], point);

                if cross < -self.tolerance
                    || (!self.include_collinear && cross.abs() < self.tolerance)
                {
                    hull.pop();
                } else {
                    break;
                }
            }

            hull.push(point);
        }

        Ok(hull)
    }

    /// Andrew's Monotone Chain Algorithmus
    fn andrew_monotone(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        let mut points = points.to_vec();

        // Sortiere Punkte lexikographisch (erst X, dann Y)
        points.sort_by(|a, b| {
            a.x.partial_cmp(&b.x)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
        });

        // Entferne Duplikate
        points.dedup_by(|a, b| {
            (a.x - b.x).abs() < self.tolerance && (a.y - b.y).abs() < self.tolerance
        });

        if points.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: points.len(),
            });
        }

        // Baue unteren Hull
        let mut lower_hull = Vec::new();
        for &point in &points {
            while lower_hull.len() >= 2 {
                let cross = Self::cross_product(
                    lower_hull[lower_hull.len() - 2],
                    lower_hull[lower_hull.len() - 1],
                    point,
                );

                if cross <= 0.0 {
                    lower_hull.pop();
                } else {
                    break;
                }
            }
            lower_hull.push(point);
        }

        // Baue oberen Hull
        let mut upper_hull = Vec::new();
        for &point in points.iter().rev() {
            while upper_hull.len() >= 2 {
                let cross = Self::cross_product(
                    upper_hull[upper_hull.len() - 2],
                    upper_hull[upper_hull.len() - 1],
                    point,
                );

                if cross <= 0.0 {
                    upper_hull.pop();
                } else {
                    break;
                }
            }
            upper_hull.push(point);
        }

        // Entferne letzten Punkt von beiden Hulls (Duplikate)
        lower_hull.pop();
        upper_hull.pop();

        // Kombiniere beide Hulls
        lower_hull.extend(upper_hull);

        Ok(lower_hull)
    }

    /// Jarvis March (Gift Wrapping) Algorithmus
    fn jarvis_march(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        if points.is_empty() {
            return Ok(Vec::new());
        }

        // Finde den linkesten Punkt
        let start_idx = points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                a.x.partial_cmp(&b.x)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
            })
            .map(|(i, _)| i)
            .unwrap();

        let mut hull = Vec::new();
        let mut current_idx = start_idx;

        loop {
            hull.push(points[current_idx]);

            let mut next_idx = 0;
            for i in 1..points.len() {
                if i == current_idx {
                    continue;
                }

                if next_idx == current_idx {
                    next_idx = i;
                    continue;
                }

                let cross = Self::cross_product(points[current_idx], points[next_idx], points[i]);

                if cross > self.tolerance {
                    next_idx = i;
                } else if cross.abs() < self.tolerance {
                    // Kollinear: wähle den ferneren Punkt
                    let dist_next = points[current_idx].distance_squared(points[next_idx]);
                    let dist_i = points[current_idx].distance_squared(points[i]);

                    if dist_i > dist_next {
                        next_idx = i;
                    }
                }
            }

            current_idx = next_idx;

            // Wenn wir zum Startpunkt zurückkehren
            if current_idx == start_idx {
                break;
            }

            // Verhindere Endlosschleifen
            if hull.len() > points.len() {
                return Err(MathError::GeometricFailure {
                    operation: "Infinite loop in Jarvis march".to_string(),
                });
            }
        }

        Ok(hull)
    }

    /// QuickHull Algorithmus
    fn quickhull(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        if points.len() < 3 {
            return Ok(points.to_vec());
        }

        // Finde die extremen Punkte (linker und rechter)
        let (min_x_idx, max_x_idx) =
            points
                .iter()
                .enumerate()
                .fold((0, 0), |(min_idx, max_idx), (i, point)| {
                    let min_idx = if point.x < points[min_idx].x {
                        i
                    } else {
                        min_idx
                    };
                    let max_idx = if point.x > points[max_idx].x {
                        i
                    } else {
                        max_idx
                    };
                    (min_idx, max_idx)
                });

        let min_point = points[min_x_idx];
        let max_point = points[max_x_idx];

        // Teile Punkte in zwei Sets
        let mut upper_set = Vec::new();
        let mut lower_set = Vec::new();

        for (i, &point) in points.iter().enumerate() {
            if i == min_x_idx || i == max_x_idx {
                continue;
            }

            let cross = Self::cross_product(min_point, max_point, point);

            if cross > self.tolerance {
                upper_set.push(point);
            } else if cross < -self.tolerance {
                lower_set.push(point);
            }
        }

        let mut hull = vec![min_point];

        // Baue oberen Hull
        self.quickhull_recursive(&mut hull, &upper_set, min_point, max_point);

        hull.push(max_point);

        // Baue unteren Hull
        self.quickhull_recursive(&mut hull, &lower_set, max_point, min_point);

        Ok(hull)
    }

    fn quickhull_recursive(&self, hull: &mut Vec<Vec2>, points: &[Vec2], p1: Vec2, p2: Vec2) {
        if points.is_empty() {
            return;
        }

        // Finde den Punkt mit dem größten Abstand zur Linie p1-p2
        let farthest = points
            .iter()
            .max_by(|&&a, &&b| {
                let dist_a = Self::point_line_distance(a, p1, p2);
                let dist_b = Self::point_line_distance(b, p1, p2);
                dist_a
                    .partial_cmp(&dist_b)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
            .unwrap();

        // Teile die Punkte in zwei neue Sets
        let mut left_set = Vec::new();
        let mut right_set = Vec::new();

        for &point in points {
            if point == farthest {
                continue;
            }

            let cross1 = Self::cross_product(p1, farthest, point);
            let cross2 = Self::cross_product(farthest, p2, point);

            if cross1 > self.tolerance {
                left_set.push(point);
            } else if cross2 > self.tolerance {
                right_set.push(point);
            }
        }

        // Rekursive Aufrufe
        self.quickhull_recursive(hull, &left_set, p1, farthest);
        hull.push(farthest);
        self.quickhull_recursive(hull, &right_set, farthest, p2);
    }

    // === Hilfsfunktionen ===

    fn cross_product(a: Vec2, b: Vec2, c: Vec2) -> f32 {
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    }

    fn point_line_distance(point: Vec2, line_start: Vec2, line_end: Vec2) -> f32 {
        let line_vec = line_end - line_start;
        let point_vec = point - line_start;

        let cross = line_vec.x * point_vec.y - line_vec.y * point_vec.x;
        cross.abs() / line_vec.length().max(constants::EPSILON)
    }
}

/// Konvexitäts-Tests und Utilities
pub struct ConvexityUtils;

impl ConvexityUtils {
    /// Prüft ob ein Polygon konvex ist
    pub fn is_convex(polygon: &Polygon) -> bool {
        let vertices = polygon.vertices();
        let n = vertices.len();

        if n < 3 {
            return false;
        }

        let mut sign = None;

        for i in 0..n {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % n];
            let p3 = vertices[(i + 2) % n];

            let cross = ConvexHullComputer::cross_product(p1, p2, p3);

            if cross.abs() > constants::EPSILON {
                let current_sign = cross > 0.0;

                match sign {
                    None => sign = Some(current_sign),
                    Some(s) if s != current_sign => return false,
                    _ => {}
                }
            }
        }

        true
    }

    /// Berechnet die Konvexitäts-Defizite (wie weit ist das Polygon von konvex entfernt)
    pub fn convexity_defects(polygon: &Polygon) -> Vec<ConvexityDefect> {
        let vertices = polygon.vertices();
        let n = vertices.len();
        let mut defects = Vec::new();

        if n < 3 {
            return defects;
        }

        // Berechne die konvexe Hülle
        let hull_computer = ConvexHullComputer::new(ConvexHullAlgorithm::AndrewMonotone);
        let hull = match hull_computer.compute_hull(vertices) {
            Ok(h) => h,
            Err(_) => return defects,
        };

        let hull_vertices = hull.vertices();

        // Finde Punkte die nicht auf der konvexen Hülle liegen
        for (i, &vertex) in vertices.iter().enumerate() {
            let mut min_distance = f32::INFINITY;
            let mut closest_edge_start = 0;

            // Finde den nächsten Punkt zur konvexen Hülle
            for j in 0..hull_vertices.len() - 1 {
                let distance = ConvexHullComputer::point_line_distance(
                    vertex,
                    hull_vertices[j],
                    hull_vertices[j + 1],
                );

                if distance < min_distance {
                    min_distance = distance;
                    closest_edge_start = j;
                }
            }

            // Wenn der Punkt signifikant innerhalb der Hülle liegt
            if min_distance > constants::EPSILON * 100.0 {
                defects.push(ConvexityDefect {
                    vertex_index: i,
                    vertex: vertex,
                    distance_to_hull: min_distance,
                    closest_hull_edge: (
                        hull_vertices[closest_edge_start],
                        hull_vertices[closest_edge_start + 1],
                    ),
                });
            }
        }

        defects
    }

    /// Approximiert ein Polygon durch ein konvexes Polygon
    pub fn convex_approximation(polygon: &Polygon, tolerance: f32) -> MathResult<Polygon> {
        let hull_computer =
            ConvexHullComputer::new(ConvexHullAlgorithm::AndrewMonotone).with_tolerance(tolerance);

        hull_computer.compute_polygon_hull(polygon)
    }
}

/// Represents a convexity defect
#[derive(Debug, Clone)]
pub struct ConvexityDefect {
    pub vertex_index: usize,
    pub vertex: Vec2,
    pub distance_to_hull: f32,
    pub closest_hull_edge: (Vec2, Vec2),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::geometry::polygon::Polygon;

    #[test]
    fn test_convex_hull_square() {
        let points = vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(1.0, 0.0),
            Vec2::new(1.0, 1.0),
            Vec2::new(0.0, 1.0),
            Vec2::new(0.5, 0.5), // Interior point
        ];

        let computer = ConvexHullComputer::new(ConvexHullAlgorithm::GrahamScan);
        let hull = computer.compute_hull(&points).unwrap();

        // Hull should be the 4 corner points
        assert_eq!(hull.len(), 5); // 4 corners + closing point
        assert!(ConvexityUtils::is_convex(&hull));
    }

    #[test]
    fn test_all_algorithms_same_result() {
        let points = vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(3.0, 1.0),
            Vec2::new(2.0, 3.0),
            Vec2::new(-1.0, 2.0),
            Vec2::new(1.0, 1.0), // Interior
        ];

        let algorithms = [
            ConvexHullAlgorithm::GrahamScan,
            ConvexHullAlgorithm::AndrewMonotone,
            ConvexHullAlgorithm::JarvisMarch,
            ConvexHullAlgorithm::QuickHull,
        ];

        let mut hulls = Vec::new();
        for algorithm in &algorithms {
            let computer = ConvexHullComputer::new(*algorithm);
            let hull = computer.compute_hull(&points).unwrap();
            hulls.push(hull);
        }

        // All algorithms should produce hulls with the same number of vertices
        let first_len = hulls[0].len();
        for hull in &hulls {
            assert_eq!(hull.len(), first_len);
            assert!(ConvexityUtils::is_convex(hull));
        }
    }

    #[test]
    fn test_convexity_detection() {
        // Convex polygon
        let convex = Polygon::closed(vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(2.0, 0.0),
            Vec2::new(1.0, 2.0),
        ])
        .unwrap();

        assert!(ConvexityUtils::is_convex(&convex));

        // Non-convex polygon (arrow shape)
        let non_convex = Polygon::closed(vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(2.0, 1.0),
            Vec2::new(1.0, 1.0),
            Vec2::new(1.0, 2.0),
            Vec2::new(0.0, 1.0),
        ])
        .unwrap();

        assert!(!ConvexityUtils::is_convex(&non_convex));
    }

    #[test]
    fn test_convexity_defects() {
        let polygon = Polygon::closed(vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(2.0, 0.0),
            Vec2::new(2.0, 2.0),
            Vec2::new(1.0, 1.0), // Creates a defect
            Vec2::new(0.0, 2.0),
        ])
        .unwrap();

        let defects = ConvexityUtils::convexity_defects(&polygon);

        // Should find the interior point as a defect
        assert!(!defects.is_empty());
    }
}
