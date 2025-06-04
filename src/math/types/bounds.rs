// src/math/types/bounds.rs

use crate::math::{error::*, types::*};
use std::fmt;

/// 2D Bounding Box (Axis-Aligned Bounding Box)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds2D {
    pub min: Point2D,
    pub max: Point2D,
}

impl Bounds2D {
    /// Erstellt eine neue Bounding Box
    pub fn new(min: Point2D, max: Point2D) -> MathResult<Self> {
        if min.x > max.x || min.y > max.y {
            return Err(MathError::InvalidConfiguration {
                message: format!("Invalid bounds: min {:?} > max {:?}", min, max),
            });
        }

        Ok(Self { min, max })
    }

    /// Erstellt eine Bounding Box aus zwei beliebigen Punkten
    pub fn from_points(p1: Point2D, p2: Point2D) -> Self {
        Self {
            min: Point2D::new(p1.x.min(p2.x), p1.y.min(p2.y)),
            max: Point2D::new(p1.x.max(p2.x), p1.y.max(p2.y)),
        }
    }

    /// Erstellt eine Bounding Box aus Zentrum und Größe
    pub fn from_center_size(center: Point2D, size: Point2D) -> Self {
        let half_size = size * 0.5;
        Self {
            min: center - half_size,
            max: center + half_size,
        }
    }

    /// Erstellt eine Bounding Box die alle Punkte umschließt
    pub fn from_points_iter<I>(points: I) -> Option<Self>
    where
        I: IntoIterator<Item = Point2D>,
    {
        let mut points_iter = points.into_iter();
        let first_point = points_iter.next()?;

        let mut min = first_point;
        let mut max = first_point;

        for point in points_iter {
            min.x = min.x.min(point.x);
            min.y = min.y.min(point.y);
            max.x = max.x.max(point.x);
            max.y = max.y.max(point.y);
        }

        Some(Self { min, max })
    }

    /// Leere Bounding Box (ungültig)
    pub fn empty() -> Self {
        Self {
            min: Point2D::new(f32::INFINITY, f32::INFINITY),
            max: Point2D::new(f32::NEG_INFINITY, f32::NEG_INFINITY),
        }
    }

    /// Infinite Bounding Box (umschließt alles)
    pub fn infinite() -> Self {
        Self {
            min: Point2D::new(f32::NEG_INFINITY, f32::NEG_INFINITY),
            max: Point2D::new(f32::INFINITY, f32::INFINITY),
        }
    }

    /// Prüft ob die Bounding Box gültig ist
    pub fn is_valid(&self) -> bool {
        self.min.x <= self.max.x
            && self.min.y <= self.max.y
            && self.min.x.is_finite()
            && self.min.y.is_finite()
            && self.max.x.is_finite()
            && self.max.y.is_finite()
    }

    /// Prüft ob die Bounding Box leer ist
    pub fn is_empty(&self) -> bool {
        self.min.x > self.max.x || self.min.y > self.max.y
    }

    /// Breite der Bounding Box
    pub fn width(&self) -> f32 {
        (self.max.x - self.min.x).max(0.0)
    }

    /// Höhe der Bounding Box
    pub fn height(&self) -> f32 {
        (self.max.y - self.min.y).max(0.0)
    }

    /// Größe der Bounding Box
    pub fn size(&self) -> Point2D {
        Point2D::new(self.width(), self.height())
    }

    /// Zentrum der Bounding Box
    pub fn center(&self) -> Point2D {
        (self.min + self.max) * 0.5
    }

    /// Fläche der Bounding Box
    pub fn area(&self) -> f32 {
        if self.is_empty() {
            0.0
        } else {
            self.width() * self.height()
        }
    }

    /// Umfang der Bounding Box
    pub fn perimeter(&self) -> f32 {
        if self.is_empty() {
            0.0
        } else {
            2.0 * (self.width() + self.height())
        }
    }

    /// Prüft ob ein Punkt in der Bounding Box liegt
    pub fn contains_point(&self, point: Point2D) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
    }

    /// Prüft ob eine andere Bounding Box vollständig enthalten ist
    pub fn contains_bounds(&self, other: &Bounds2D) -> bool {
        if other.is_empty() {
            return true;
        }
        if self.is_empty() {
            return false;
        }

        self.min.x <= other.min.x
            && self.max.x >= other.max.x
            && self.min.y <= other.min.y
            && self.max.y >= other.max.y
    }

    /// Prüft ob sich zwei Bounding Boxes überschneiden
    pub fn intersects(&self, other: &Bounds2D) -> bool {
        if self.is_empty() || other.is_empty() {
            return false;
        }

        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
    }

    /// Berechnet die Überschneidung zweier Bounding Boxes
    pub fn intersection(&self, other: &Bounds2D) -> Self {
        if !self.intersects(other) {
            return Self::empty();
        }

        Self {
            min: Point2D::new(self.min.x.max(other.min.x), self.min.y.max(other.min.y)),
            max: Point2D::new(self.max.x.min(other.max.x), self.max.y.min(other.max.y)),
        }
    }

    /// Vereinigt zwei Bounding Boxes
    pub fn union(&self, other: &Bounds2D) -> Self {
        if self.is_empty() {
            return *other;
        }
        if other.is_empty() {
            return *self;
        }

        Self {
            min: Point2D::new(self.min.x.min(other.min.x), self.min.y.min(other.min.y)),
            max: Point2D::new(self.max.x.max(other.max.x), self.max.y.max(other.max.y)),
        }
    }

    /// Erweitert die Bounding Box um einen Punkt
    pub fn expand_to_include_point(&mut self, point: Point2D) {
        if self.is_empty() {
            self.min = point;
            self.max = point;
        } else {
            self.min.x = self.min.x.min(point.x);
            self.min.y = self.min.y.min(point.y);
            self.max.x = self.max.x.max(point.x);
            self.max.y = self.max.y.max(point.y);
        }
    }

    /// Erweitert die Bounding Box um eine andere Bounding Box
    pub fn expand_to_include_bounds(&mut self, other: &Bounds2D) {
        if !other.is_empty() {
            self.expand_to_include_point(other.min);
            self.expand_to_include_point(other.max);
        }
    }

    /// Erweitert die Bounding Box um einen Margin
    pub fn expand(&self, margin: f32) -> Self {
        if self.is_empty() {
            return *self;
        }

        Self {
            min: Point2D::new(self.min.x - margin, self.min.y - margin),
            max: Point2D::new(self.max.x + margin, self.max.y + margin),
        }
    }

    /// Erweitert die Bounding Box um verschiedene Margins
    pub fn expand_by(&self, left: f32, right: f32, top: f32, bottom: f32) -> Self {
        if self.is_empty() {
            return *self;
        }

        Self {
            min: Point2D::new(self.min.x - left, self.min.y - bottom),
            max: Point2D::new(self.max.x + right, self.max.y + top),
        }
    }

    /// Skaliert die Bounding Box um einen Faktor
    pub fn scale(&self, factor: f32) -> Self {
        if self.is_empty() {
            return *self;
        }

        let center = self.center();
        let half_size = self.size() * factor * 0.5;

        Self {
            min: center - half_size,
            max: center + half_size,
        }
    }

    /// Skaliert die Bounding Box um verschiedene Faktoren
    pub fn scale_xy(&self, factor_x: f32, factor_y: f32) -> Self {
        if self.is_empty() {
            return *self;
        }

        let center = self.center();
        let half_size = Point2D::new(
            self.width() * factor_x * 0.5,
            self.height() * factor_y * 0.5,
        );

        Self {
            min: center - half_size,
            max: center + half_size,
        }
    }

    /// Verschiebt die Bounding Box
    pub fn translate(&self, offset: Point2D) -> Self {
        Self {
            min: self.min + offset,
            max: self.max + offset,
        }
    }

    /// Berechnet den nächsten Punkt auf der Bounding Box zu einem gegebenen Punkt
    pub fn closest_point(&self, point: Point2D) -> Point2D {
        if self.is_empty() {
            return point;
        }

        Point2D::new(
            point.x.clamp(self.min.x, self.max.x),
            point.y.clamp(self.min.y, self.max.y),
        )
    }

    /// Berechnet den Abstand von einem Punkt zur Bounding Box
    pub fn distance_to_point(&self, point: Point2D) -> f32 {
        if self.is_empty() {
            return f32::INFINITY;
        }

        if self.contains_point(point) {
            return 0.0;
        }

        let closest = self.closest_point(point);
        point.distance(closest)
    }

    /// Erzeugt die vier Eckpunkte der Bounding Box
    pub fn corners(&self) -> [Point2D; 4] {
        [
            self.min,                             // unten links
            Point2D::new(self.max.x, self.min.y), // unten rechts
            self.max,                             // oben rechts
            Point2D::new(self.min.x, self.max.y), // oben links
        ]
    }
}

impl fmt::Display for Bounds2D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            write!(f, "Bounds2D(empty)")
        } else {
            write!(f, "Bounds2D({:?} to {:?})", self.min, self.max)
        }
    }
}

/// 3D Bounding Box
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds3D {
    pub min: Point3D,
    pub max: Point3D,
}

impl Bounds3D {
    pub fn new(min: Point3D, max: Point3D) -> MathResult<Self> {
        if min.x > max.x || min.y > max.y || min.z > max.z {
            return Err(MathError::InvalidConfiguration {
                message: format!("Invalid 3D bounds: min {:?} > max {:?}", min, max),
            });
        }

        Ok(Self { min, max })
    }

    pub fn from_points(p1: Point3D, p2: Point3D) -> Self {
        Self {
            min: Point3D::new(p1.x.min(p2.x), p1.y.min(p2.y), p1.z.min(p2.z)),
            max: Point3D::new(p1.x.max(p2.x), p1.y.max(p2.y), p1.z.max(p2.z)),
        }
    }

    pub fn center(&self) -> Point3D {
        (self.min + self.max) * 0.5
    }

    pub fn size(&self) -> Point3D {
        self.max - self.min
    }

    pub fn volume(&self) -> f32 {
        let size = self.size();
        size.x * size.y * size.z
    }

    pub fn contains_point(&self, point: Point3D) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
            && point.z >= self.min.z
            && point.z <= self.max.z
    }
}
