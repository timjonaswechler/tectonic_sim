// src/math/types.rs

use crate::math::utils::constants; // utils::constants für EPSILON
use bevy::math::{Vec2, Vec3}; // Bevy's Vektor-Typen
use spade::Point2 as SpadeInternalPoint2; // Spade's Punkt-Typ, umbenannt zur Klarheit

// --- Kern-Typaliase ---
// Wir verwenden Vec2 und Vec3 direkt, keine zusätzlichen Aliase wie Point2D/Point3D hier,
// um Klarheit zu schaffen, dass wir Bevy-Typen verwenden.
// Anwender des math-Moduls können `use crate::math::types::{Vec2, Vec3};` verwenden,
// oder direkt `bevy::math::{Vec2, Vec3}` wenn sie Bevy eh schon importiert haben.
// Für die interne Konsistenz innerhalb des `math` Moduls werden wir Vec2/Vec3 nutzen.

/// Typ für 2D-Punkte, die mit Spade verwendet werden.
pub type SpadePoint = SpadeInternalPoint2<f32>;

// --- Konvertierungsfunktionen ---

/// Konvertiert einen Slice von Bevy Vec2 in einen Vec von Spade Points.
pub fn bevy_to_spade_points(points: &[Vec2]) -> Vec<SpadePoint> {
    points.iter().map(|p| SpadePoint::new(p.x, p.y)).collect()
}

/// Konvertiert einen Slice von Spade Points in einen Vec von Bevy Vec2.
pub fn spade_to_bevy_points(points: &[SpadePoint]) -> Vec<Vec2> {
    points.iter().map(|p| Vec2::new(p.x, p.y)).collect()
}

// --- Vektor Erweiterungen ---

/// Erweiterte Vektor-Operationen für Bevy's Vec2.
pub trait Vector2DExt {
    // fn length_squared(&self) -> f32; // length() ist schon in Vec2
    // fn try_normalize(&self) -> Option<Self>
    // where
    //     Self: Sized; // normalize() ist schon da, aber gibt ZERO bei Länge 0
    fn with_length(&self, length: f32) -> Self
    where
        Self: Sized;
    fn clamp_length(&self, max_length: f32) -> Self
    where
        Self: Sized;
    // dot() ist schon da
    fn cross_product(&self, other: Self) -> f32; // Bevy hat `perp_dot` was (a.perp()).dot(b) ist. cross ist a.x * b.y - a.y * b.x
    // fn angle_between(&self, other: Self) -> f32; // angle() ist schon da (Winkel zur X-Achse)
    // rotate() ist schon da (als Vec2::from_angle(self.angle() + angle_rad).normalize() * self.length())
    // oder man implementiert es direkt. Bevy's `Mat2::from_angle(rad).mul_vec2(*self)`
    fn perpendicular(&self) -> Self; // Vec2 hat `perp()`
    // fn reflect_around_normal(&self, normal: Self) -> Self
    // where
    //     Self: Sized; // Vec2 hat reflect
    // fn project_onto(&self, other: Self) -> Self
    // where
    //     Self: Sized; // Vec2 hat project_onto
    // fn reject_from(&self, other: Self) -> Self
    // where
    //     Self: Sized; // Vec2 hat reject_from
    // lerp() ist schon da
    fn slerp(&self, other: Self, t: f32) -> Self
    where
        Self: Sized;
}

impl Vector2DExt for Vec2 {
    // fn length_squared(&self) -> f32 {
    //     // Vec2::length_squared() existiert bereits
    //     self.length_squared()
    // }

    // fn try_normalize(&self) -> Option<Self> {
    //     // Vec2::try_normalize() existiert bereits
    //     self.try_normalize()
    // }

    fn with_length(&self, length: f32) -> Self {
        self.normalize_or_zero() * length // normalize_or_zero ist sicherer
    }

    fn clamp_length(&self, max_length: f32) -> Self {
        // Vec2::clamp_length_max() existiert bereits
        self.clamp_length_max(max_length)
    }

    fn cross_product(&self, other: Self) -> f32 {
        self.x * other.y - self.y * other.x
        // Alternative: self.perp_dot(other) -- aber Vorzeichen beachten!
        // perp_dot(v) = x * v.y - y * v.x (wenn self=(x,y)) => das ist es!
    }

    // fn angle_between(&self, other: Self) -> f32 {
    //     // Vec2::angle_between() existiert bereits
    //     self.angle_between(other)
    // }

    fn perpendicular(&self) -> Self {
        // Vec2::perp() existiert bereits
        self.perp()
    }

    // fn reflect_around_normal(&self, normal: Self) -> Self {
    //     // Vec2::reflect() existiert bereits
    //     self.reflect(normal)
    // }

    // fn project_onto(&self, other: Self) -> Self {
    //     // Vec2::project_onto() existiert bereits
    //     self.project_onto(other)
    // }

    // fn reject_from(&self, other: Self) -> Self {
    //     // Vec2::reject_from() existiert bereits
    //     self.reject_from(other)
    // }

    fn slerp(&self, other: Self, t: f32) -> Self {
        // Bevy's Vec2 hat keine direkte slerp-Methode.
        // Diese Implementierung ist eine gängige.
        let dot_product = self.dot(other).clamp(-1.0, 1.0);
        let angle = dot_product.acos();

        if angle.abs() < constants::EPSILON {
            // Nahezu kollinear
            return self.lerp(other, t);
        }

        let sin_angle = angle.sin();
        if sin_angle.abs() < constants::EPSILON {
            // Winkel ist 0 oder PI
            return self.lerp(other, t); // Fallback zu lerp
        }

        let a = ((1.0 - t) * angle).sin() / sin_angle;
        let b = (t * angle).sin() / sin_angle;
        *self * a + other * b
    }
}

/// Erweiterte Vektor-Operationen für Bevy's Vec3.
pub trait Vector3DExt {
    // length(), length_squared(), normalize(), try_normalize(), dot(), cross(), lerp(), project_onto(), reflect()
    // sind bereits in Bevy's Vec3 vorhanden.
    // Wir könnten hier spezifischere oder nicht vorhandene Methoden hinzufügen.
    // Für den Moment belassen wir es dabei, da Vec3 schon sehr mächtig ist.
    // Beispiel:
    // fn angle_between(&self, other: Self) -> f32;
}

impl Vector3DExt for Vec3 {
    // fn angle_between(&self, other: Self) -> f32 {
    //     self.normalize_or_zero().dot(other.normalize_or_zero()).acos()
    // }
}

/// Builder für komplexe Vec2-Operationen (Beispielhaft, da Vec2 schon viel kann)
pub struct Vec2Builder {
    vector: Vec2,
}

impl Vec2Builder {
    pub fn new(x: f32, y: f32) -> Self {
        Self {
            vector: Vec2::new(x, y),
        }
    }
    pub fn from_angle(angle_rad: f32) -> Self {
        Self {
            vector: Vec2::from_angle(angle_rad),
        }
    }

    pub fn from_angle_length(angle_rad: f32, length: f32) -> Self {
        Self {
            vector: Vec2::from_angle(angle_rad) * length,
        }
    }

    pub fn normalize(mut self) -> Self {
        self.vector = self.vector.normalize_or_zero();
        self
    }

    pub fn with_length(mut self, length: f32) -> Self {
        self.vector = self.vector.normalize_or_zero() * length;
        self
    }

    pub fn rotate(mut self, angle_rad: f32) -> Self {
        // Einfache Rotation, da Bevy's Vec2::rotate das nicht direkt so macht.
        let current_angle = self.vector.y.atan2(self.vector.x);
        let length = self.vector.length();
        self.vector = Vec2::from_angle(current_angle + angle_rad) * length;
        self
    }

    pub fn add(mut self, other: Vec2) -> Self {
        self.vector += other;
        self
    }

    pub fn scale(mut self, factor: f32) -> Self {
        self.vector *= factor;
        self
    }

    pub fn build(self) -> Vec2 {
        self.vector
    }
}
