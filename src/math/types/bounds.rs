use crate::math::error::{MathError, MathResult};
use bevy::math::{Vec2, Vec3};
// --- Boolesche Hilfsstrukturen ---
// (Behalten, falls nützlich für komponentenweise Vergleiche)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Bool2 {
    pub x: bool,
    pub y: bool,
}
impl Bool2 {
    pub fn new(x: bool, y: bool) -> Self {
        Self { x, y }
    }
    pub fn all(self) -> bool {
        self.x && self.y
    }
    pub fn any(self) -> bool {
        self.x || self.y
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Bool3 {
    pub x: bool,
    pub y: bool,
    pub z: bool,
}
impl Bool3 {
    pub fn new(x: bool, y: bool, z: bool) -> Self {
        Self { x, y, z }
    }
    pub fn all(self) -> bool {
        self.x && self.y && self.z
    }
    pub fn any(self) -> bool {
        self.x || self.y || self.z
    }
}

// --- Bounding Box Strukturen ---

/// 2D Bounding Box (Axis-Aligned Bounding Box)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds2D {
    pub min: Vec2,
    pub max: Vec2,
}

impl Bounds2D {
    /// Erstellt eine neue Bounding Box.
    pub fn new(min: Vec2, max: Vec2) -> MathResult<Self> {
        if min.x > max.x || min.y > max.y {
            return Err(MathError::InvalidConfiguration {
                message: format!("Invalid bounds: min {:?} > max {:?}", min, max),
            });
        }
        Ok(Self { min, max })
    }

    /// Erstellt eine Bounding Box aus zwei beliebigen Punkten.
    pub fn from_points(p1: Vec2, p2: Vec2) -> Self {
        Self {
            min: Vec2::new(p1.x.min(p2.x), p1.y.min(p2.y)),
            max: Vec2::new(p1.x.max(p2.x), p1.y.max(p2.y)),
        }
    }

    /// Erstellt eine Bounding Box aus Zentrum und Größe.
    pub fn from_center_size(center: Vec2, size: Vec2) -> MathResult<Self> {
        if size.x < 0.0 || size.y < 0.0 {
            return Err(MathError::InvalidConfiguration {
                message: "Bounds size cannot be negative".to_string(),
            });
        }
        let half_size = size * 0.5;
        Self::new(center - half_size, center + half_size)
    }

    /// Erstellt eine Bounding Box die alle Punkte umschließt.
    pub fn from_points_iter<I>(points: I) -> Option<Self>
    where
        I: IntoIterator<Item = Vec2>,
    {
        let mut points_iter = points.into_iter();
        let first_point = points_iter.next()?;

        let mut min = first_point;
        let mut max = first_point;

        for point in points_iter {
            min = min.min(point); // Bevy's Vec2 hat min/max
            max = max.max(point);
        }
        Some(Self { min, max }) // new() wird nicht benötigt, da min/max korrekt sind
    }

    /// Leere Bounding Box (ungültig, aber nützlich als Startwert).
    pub fn empty() -> Self {
        Self {
            min: Vec2::new(f32::INFINITY, f32::INFINITY),
            max: Vec2::new(f32::NEG_INFINITY, f32::NEG_INFINITY),
        }
    }

    /// Infinite Bounding Box (umschließt alles).
    pub fn infinite() -> Self {
        Self {
            min: Vec2::new(f32::NEG_INFINITY, f32::NEG_INFINITY),
            max: Vec2::new(f32::INFINITY, f32::INFINITY),
        }
    }

    /// Prüft ob die Bounding Box gültig ist.
    pub fn is_valid(&self) -> bool {
        self.min.x <= self.max.x
            && self.min.y <= self.max.y
            && self.min.x.is_finite()
            && self.min.y.is_finite()
            && self.max.x.is_finite()
            && self.max.y.is_finite()
    }

    /// Prüft ob die Bounding Box leer ist (im Sinne von invertiert/ungültig).
    pub fn is_empty(&self) -> bool {
        self.min.x > self.max.x || self.min.y > self.max.y
    }

    pub fn width(&self) -> f32 {
        (self.max.x - self.min.x).max(0.0)
    }
    pub fn height(&self) -> f32 {
        (self.max.y - self.min.y).max(0.0)
    }
    pub fn size(&self) -> Vec2 {
        Vec2::new(self.width(), self.height())
    }
    pub fn center(&self) -> Vec2 {
        (self.min + self.max) * 0.5
    }

    pub fn area(&self) -> f32 {
        if self.is_empty() {
            0.0
        } else {
            self.width() * self.height()
        }
    }

    pub fn perimeter(&self) -> f32 {
        if self.is_empty() {
            0.0
        } else {
            2.0 * (self.width() + self.height())
        }
    }

    pub fn contains_point(&self, point: Vec2) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
    }

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

    pub fn intersects(&self, other: &Bounds2D) -> bool {
        if self.is_empty() || other.is_empty() {
            return false;
        }
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
    }

    pub fn intersection(&self, other: &Bounds2D) -> Self {
        if !self.intersects(other) {
            return Self::empty();
        }
        Self {
            min: self.min.max(other.min), // Vec2::max
            max: self.max.min(other.max), // Vec2::min
        }
    }

    pub fn union(&self, other: &Bounds2D) -> Self {
        if self.is_empty() {
            return *other;
        }
        if other.is_empty() {
            return *self;
        }
        Self {
            min: self.min.min(other.min), // Vec2::min
            max: self.max.max(other.max), // Vec2::max
        }
    }

    pub fn expand_to_include_point(&mut self, point: Vec2) {
        if self.is_empty() {
            self.min = point;
            self.max = point;
        } else {
            self.min = self.min.min(point);
            self.max = self.max.max(point);
        }
    }

    pub fn expand_to_include_bounds(&mut self, other: &Bounds2D) {
        if !other.is_empty() {
            self.expand_to_include_point(other.min);
            self.expand_to_include_point(other.max);
        }
    }

    pub fn expand(&self, margin: f32) -> Self {
        if self.is_empty() {
            return *self;
        }
        // new() kann fehlschlagen, wenn margin negativ und groß genug ist, um min > max zu machen
        // aber typischerweise ist margin positiv. Wir könnten hier ein unwrap_or_else verwenden
        // oder die Verantwortung dem Aufrufer überlassen. Fürs Erste:
        Self::new(
            self.min - Vec2::splat(margin),
            self.max + Vec2::splat(margin),
        )
        .unwrap_or_else(|_| *self) // Fallback auf self bei ungültigen Bounds
    }

    pub fn expand_by(&self, left: f32, right: f32, top: f32, bottom: f32) -> Self {
        if self.is_empty() {
            return *self;
        }
        Self::new(
            Vec2::new(self.min.x - left, self.min.y - bottom),
            Vec2::new(self.max.x + right, self.max.y + top),
        )
        .unwrap_or_else(|_| *self)
    }

    pub fn scale_around_center(&self, factor: f32) -> MathResult<Self> {
        if self.is_empty() {
            return Ok(*self);
        }
        let center = self.center();
        let half_size = self.size() * factor * 0.5;
        Self::from_center_size(center, half_size)
    }

    pub fn scale_xy_around_center(&self, factor_x: f32, factor_y: f32) -> MathResult<Self> {
        if self.is_empty() {
            return Ok(*self);
        }
        let center = self.center();
        let new_size = Vec2::new(self.width() * factor_x, self.height() * factor_y);
        Self::from_center_size(center, new_size)
    }

    pub fn translate(&self, offset: Vec2) -> Self {
        // new() kann hier nicht fehlschlagen, da es nur eine Verschiebung ist
        Self::new(self.min + offset, self.max + offset).expect("Translate should not fail")
    }

    pub fn closest_point(&self, point: Vec2) -> Vec2 {
        if self.is_empty() {
            return point;
        } // Oder ein anderer Fallback
        Vec2::new(
            point.x.clamp(self.min.x, self.max.x),
            point.y.clamp(self.min.y, self.max.y),
        )
    }

    pub fn distance_to_point_sq(&self, point: Vec2) -> f32 {
        if self.is_empty() {
            return f32::INFINITY;
        }
        if self.contains_point(point) {
            return 0.0;
        }
        let closest = self.closest_point(point);
        point.distance_squared(closest)
    }

    pub fn distance_to_point(&self, point: Vec2) -> f32 {
        self.distance_to_point_sq(point).sqrt()
    }

    pub fn corners(&self) -> [Vec2; 4] {
        [
            self.min,
            Vec2::new(self.max.x, self.min.y),
            self.max,
            Vec2::new(self.min.x, self.max.y),
        ]
    }
}

impl std::fmt::Display for Bounds2D {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
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
    pub min: Vec3,
    pub max: Vec3,
}

impl Bounds3D {
    pub fn new(min: Vec3, max: Vec3) -> MathResult<Self> {
        if min.x > max.x || min.y > max.y || min.z > max.z {
            return Err(MathError::InvalidConfiguration {
                message: format!("Invalid 3D bounds: min {:?} > max {:?}", min, max),
            });
        }
        Ok(Self { min, max })
    }

    pub fn from_points(p1: Vec3, p2: Vec3) -> Self {
        Self {
            min: p1.min(p2),
            max: p1.max(p2),
        }
    }

    pub fn center(&self) -> Vec3 {
        (self.min + self.max) * 0.5
    }
    pub fn size(&self) -> Vec3 {
        self.max - self.min
    }

    pub fn volume(&self) -> f32 {
        let s = self.size();
        if s.x < 0.0 || s.y < 0.0 || s.z < 0.0 {
            0.0
        } else {
            s.x * s.y * s.z
        }
    }

    pub fn contains_point(&self, point: Vec3) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
            && point.z >= self.min.z
            && point.z <= self.max.z
    }
}
