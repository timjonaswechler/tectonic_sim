// src/math/geometry/polygon/core/properties.rs

use crate::math::geometry::polygon::Polygon;
use crate::math::types::*;

/// Trait für Polygon-Eigenschaften
pub trait PolygonProperties {
    /// Berechnet die Fläche des Polygons (Shoelace-Formel)
    fn area(&self) -> f32;

    /// Berechnet den Umfang des Polygons
    fn perimeter(&self) -> f32;

    /// Prüft ob ein Punkt innerhalb des Polygons liegt (Ray-Casting)
    fn contains_point(&self, point: Point2D) -> bool;

    /// Prüft ob das Polygon konvex ist
    fn is_convex(&self) -> bool;

    /// Prüft die Orientierung (im Uhrzeigersinn oder gegen)
    fn orientation(&self) -> Orientation;

    /// Berechnet den Schwerpunkt (geometrisch korrekt für Polygone)
    fn geometric_centroid(&self) -> Option<Point2D>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Clockwise,
    CounterClockwise,
    Collinear,
}

impl PolygonProperties for Polygon {
    fn area(&self) -> f32 {
        if self.vertices.len() < 3 {
            return 0.0;
        }

        let mut area = 0.0;
        let n = if self.is_closed() {
            self.vertices.len() - 1 // Letzten Punkt ignorieren da er dem ersten entspricht
        } else {
            self.vertices.len()
        };

        for i in 0..n {
            let j = (i + 1) % n;
            area += self.vertices[i].x * self.vertices[j].y;
            area -= self.vertices[j].x * self.vertices[i].y;
        }

        (area * 0.5f32).abs()
    }

    fn perimeter(&self) -> f32 {
        if self.vertices.len() < 2 {
            return 0.0;
        }

        let mut perimeter = 0.0;
        let n = self.vertices.len();

        for i in 0..n {
            let j = if self.is_closed() {
                (i + 1) % (n - 1) // Bei geschlossenen Polygonen den duplizierten Punkt überspringen
            } else if i == n - 1 {
                break; // Bei offenen Polygonen beim letzten Punkt stoppen
            } else {
                i + 1
            };

            perimeter += self.vertices[i].distance(self.vertices[j]);
        }

        perimeter
    }

    fn contains_point(&self, point: Point2D) -> bool {
        if self.vertices.len() < 3 {
            return false;
        }

        let mut inside = false;
        let n = if self.is_closed() {
            self.vertices.len() - 1
        } else {
            self.vertices.len()
        };

        let mut j = n - 1;
        for i in 0..n {
            let vi = self.vertices[i];
            let vj = self.vertices[j];

            if ((vi.y > point.y) != (vj.y > point.y))
                && (point.x < (vj.x - vi.x) * (point.y - vi.y) / (vj.y - vi.y) + vi.x)
            {
                inside = !inside;
            }
            j = i;
        }

        inside
    }

    fn is_convex(&self) -> bool {
        if self.vertices.len() < 3 {
            return false;
        }

        let n = if self.is_closed() {
            self.vertices.len() - 1
        } else {
            self.vertices.len()
        };

        if n < 3 {
            return false;
        }

        let mut sign = None;

        for i in 0..n {
            let p1 = self.vertices[i];
            let p2 = self.vertices[(i + 1) % n];
            let p3 = self.vertices[(i + 2) % n];

            let cross_product = (p2.x - p1.x) * (p3.y - p2.y) - (p2.y - p1.y) * (p3.x - p2.x);

            if cross_product.abs() > 1e-10 {
                // Nicht kollinear
                let current_sign = cross_product > 0.0;

                match sign {
                    None => sign = Some(current_sign),
                    Some(s) if s != current_sign => return false,
                    _ => {}
                }
            }
        }

        true
    }

    fn orientation(&self) -> Orientation {
        if self.vertices.len() < 3 {
            return Orientation::Collinear;
        }

        let area_doubled = {
            let mut area = 0.0;
            let n = if self.is_closed() {
                self.vertices.len() - 1
            } else {
                self.vertices.len()
            };

            for i in 0..n {
                let j = (i + 1) % n;
                area += self.vertices[i].x * self.vertices[j].y;
                area -= self.vertices[j].x * self.vertices[i].y;
            }
            area
        };

        if (area_doubled as f32).abs() < 1e-10f32 {
            Orientation::Collinear
        } else if area_doubled > 0.0 {
            Orientation::CounterClockwise
        } else {
            Orientation::Clockwise
        }
    }

    fn geometric_centroid(&self) -> Option<Point2D> {
        if self.vertices.len() < 3 {
            return None;
        }

        let area = self.area();
        if area.abs() < 1e-10 {
            // Fallback auf arithmetisches Mittel
            return self.centroid();
        }

        let mut cx = 0.0;
        let mut cy = 0.0;

        let n = if self.is_closed() {
            self.vertices.len() - 1
        } else {
            self.vertices.len()
        };

        for i in 0..n {
            let j = (i + 1) % n;
            let factor =
                self.vertices[i].x * self.vertices[j].y - self.vertices[j].x * self.vertices[i].y;
            cx += (self.vertices[i].x + self.vertices[j].x) * factor;
            cy += (self.vertices[i].y + self.vertices[j].y) * factor;
        }

        let factor = 1.0 / (6.0 * area);
        Some(Point2D::new(cx * factor, cy * factor))
    }
}
