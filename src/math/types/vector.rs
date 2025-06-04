use crate::math::{error::*, types::*, utils::*};
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Erweiterte Vector-Operationen für Point2D
pub trait Vector2DExt {
    /// Länge des Vektors
    fn length(&self) -> f32;

    /// Quadrierte Länge (effizienter als length() wenn nur Vergleich nötig)
    fn length_squared(&self) -> f32;

    /// Normalisiert den Vektor (Länge = 1)
    fn normalize(&self) -> Self;

    /// Versucht den Vektor zu normalisieren, gibt None zurück wenn zu kurz
    fn try_normalize(&self, min_length: f32) -> Option<Self>
    where
        Self: Sized;

    /// Setzt die Länge des Vektors
    fn with_length(&self, length: f32) -> Self
    where
        Self: Sized;

    /// Begrenzt die Länge des Vektors
    fn clamp_length(&self, max_length: f32) -> Self
    where
        Self: Sized;

    /// Skalarprodukt
    fn dot(&self, other: Self) -> f32;

    /// Kreuzprodukt (2D gibt Skalar zurück)
    fn cross(&self, other: Self) -> f32;

    /// Winkel des Vektors
    fn angle(&self) -> f32;

    /// Winkel zwischen zwei Vektoren
    fn angle_to(&self, other: Self) -> f32;

    /// Rotiert den Vektor um einen Winkel
    fn rotate(&self, angle: f32) -> Self;

    /// Senkrechter Vektor (90° gedreht)
    fn perpendicular(&self) -> Self;

    /// Spiegelt den Vektor an einer Normalen
    fn reflect(&self, normal: Self) -> Self
    where
        Self: Sized;

    /// Projiziert den Vektor auf einen anderen
    fn project_onto(&self, other: Self) -> Self
    where
        Self: Sized;

    /// Rejection (Vektor minus Projektion)
    fn reject_from(&self, other: Self) -> Self
    where
        Self: Sized;

    /// Lineare Interpolation
    fn lerp(&self, other: Self, t: f32) -> Self
    where
        Self: Sized;

    /// Sphärische lineare Interpolation (für normalisierte Vektoren)
    fn slerp(&self, other: Self, t: f32) -> Self
    where
        Self: Sized;
}

impl Vector2DExt for Point2D {
    fn length(&self) -> f32 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    fn length_squared(&self) -> f32 {
        self.x * self.x + self.y * self.y
    }

    fn normalize(&self) -> Self {
        let len = self.length();
        if len > constants::EPSILON {
            *self / len
        } else {
            Point2D::ZERO
        }
    }

    fn try_normalize(&self, min_length: f32) -> Option<Self> {
        let len = self.length();
        if len >= min_length {
            Some(*self / len)
        } else {
            None
        }
    }

    fn with_length(&self, length: f32) -> Self {
        self.normalize() * length
    }

    fn clamp_length(&self, max_length: f32) -> Self {
        let len_sq = self.length_squared();
        if len_sq > max_length * max_length {
            self.with_length(max_length)
        } else {
            *self
        }
    }

    fn dot(&self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y
    }

    fn cross(&self, other: Self) -> f32 {
        self.x * other.y - self.y * other.x
    }

    fn angle(&self) -> f32 {
        self.y.atan2(self.x)
    }

    fn angle_to(&self, other: Self) -> f32 {
        let dot = self.dot(other);
        let cross = self.cross(other);
        cross.atan2(dot)
    }

    fn rotate(&self, angle: f32) -> Self {
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        Point2D::new(
            self.x * cos_a - self.y * sin_a,
            self.x * sin_a + self.y * cos_a,
        )
    }

    fn perpendicular(&self) -> Self {
        Point2D::new(-self.y, self.x)
    }

    fn reflect(&self, normal: Self) -> Self {
        *self - normal * (2.0 * self.dot(normal))
    }

    fn project_onto(&self, other: Self) -> Self {
        let other_len_sq = other.length_squared();
        if other_len_sq > constants::EPSILON {
            other * (self.dot(other) / other_len_sq)
        } else {
            Point2D::ZERO
        }
    }

    fn reject_from(&self, other: Self) -> Self {
        *self - self.project_onto(other)
    }

    fn lerp(&self, other: Self, t: f32) -> Self {
        *self + (other - *self) * t
    }

    fn slerp(&self, other: Self, t: f32) -> Self {
        let dot = self.dot(other).clamp(-1.0, 1.0);
        let angle = dot.acos() * t;

        if angle.abs() < constants::EPSILON {
            return self.lerp(other, t);
        }

        let other_normalized = (other - *self * dot).normalize();
        *self * angle.cos() + other_normalized * angle.sin()
    }
}

/// Vector Operationen für Point3D
pub trait Vector3DExt {
    fn length(&self) -> f32;
    fn length_squared(&self) -> f32;
    fn normalize(&self) -> Self;
    fn try_normalize(&self, min_length: f32) -> Option<Self>
    where
        Self: Sized;
    fn dot(&self, other: Self) -> f32;
    fn cross(&self, other: Self) -> Self;
    fn lerp(&self, other: Self, t: f32) -> Self
    where
        Self: Sized;
    fn project_onto(&self, other: Self) -> Self
    where
        Self: Sized;
    fn reflect(&self, normal: Self) -> Self
    where
        Self: Sized;
}

impl Vector3DExt for Point3D {
    fn length(&self) -> f32 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    fn length_squared(&self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    fn normalize(&self) -> Self {
        let len = self.length();
        if len > constants::EPSILON {
            *self / len
        } else {
            Point3D::ZERO
        }
    }

    fn try_normalize(&self, min_length: f32) -> Option<Self> {
        let len = self.length();
        if len >= min_length {
            Some(*self / len)
        } else {
            None
        }
    }

    fn dot(&self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn cross(&self, other: Self) -> Self {
        Point3D::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    fn lerp(&self, other: Self, t: f32) -> Self {
        *self + (other - *self) * t
    }

    fn project_onto(&self, other: Self) -> Self {
        let other_len_sq = other.length_squared();
        if other_len_sq > constants::EPSILON {
            other * (self.dot(other) / other_len_sq)
        } else {
            Point3D::ZERO
        }
    }

    fn reflect(&self, normal: Self) -> Self {
        *self - normal * (2.0 * self.dot(normal))
    }
}

/// Vector-Builder für komplexe Vector-Operationen
pub struct VectorBuilder2D {
    vector: Point2D,
}

impl VectorBuilder2D {
    pub fn new(x: f32, y: f32) -> Self {
        Self {
            vector: Point2D::new(x, y),
        }
    }

    pub fn from_angle(angle: f32) -> Self {
        Self {
            vector: Point2D::new(angle.cos(), angle.sin()),
        }
    }

    pub fn from_angle_length(angle: f32, length: f32) -> Self {
        Self {
            vector: Point2D::new(angle.cos() * length, angle.sin() * length),
        }
    }

    pub fn normalize(mut self) -> Self {
        self.vector = self.vector.normalize();
        self
    }

    pub fn with_length(mut self, length: f32) -> Self {
        self.vector = self.vector.with_length(length);
        self
    }

    pub fn rotate(mut self, angle: f32) -> Self {
        self.vector = Vector2DExt::rotate(&self.vector, angle);
        self
    }

    pub fn add(mut self, other: Point2D) -> Self {
        self.vector = self.vector + other;
        self
    }

    pub fn scale(mut self, factor: f32) -> Self {
        self.vector = self.vector * factor;
        self
    }

    pub fn build(self) -> Point2D {
        self.vector
    }
}
