// src/math/utils.rs

use crate::math::{error::*, types::*};
use std::f32::consts::{PI, TAU};

/// Mathematische Konstanten
pub mod constants {
    pub const EPSILON: f32 = 1e-6;
    pub const EPSILON_F64: f64 = 1e-10;
    pub const GOLDEN_RATIO: f32 = 1.618033988749895;
    pub const SQRT_2: f32 = 1.4142135623730951;
    pub const SQRT_3: f32 = 1.7320508075688772;
}

/// Vergleichsfunktionen mit Toleranz
pub mod comparison {
    use super::constants::EPSILON;

    /// Prüft ob zwei Floats (nahezu) gleich sind
    pub fn nearly_equal(a: f32, b: f32) -> bool {
        (a - b).abs() < EPSILON
    }

    /// Prüft ob zwei Floats mit custom Toleranz gleich sind
    pub fn nearly_equal_eps(a: f32, b: f32, epsilon: f32) -> bool {
        (a - b).abs() < epsilon
    }

    /// Prüft ob Float (nahezu) Null ist
    pub fn nearly_zero(a: f32) -> bool {
        a.abs() < EPSILON
    }

    /// Clamping für Floats
    pub fn clamp(value: f32, min: f32, max: f32) -> f32 {
        if value < min {
            min
        } else if value > max {
            max
        } else {
            value
        }
    }

    /// Lineare Interpolation
    pub fn lerp(a: f32, b: f32, t: f32) -> f32 {
        a + (b - a) * t
    }

    /// Inverse lineare Interpolation
    pub fn inverse_lerp(a: f32, b: f32, value: f32) -> f32 {
        if nearly_equal(a, b) {
            0.0
        } else {
            (value - a) / (b - a)
        }
    }

    /// Smoothstep interpolation (hermite interpolation)
    pub fn smoothstep(edge0: f32, edge1: f32, x: f32) -> f32 {
        let t = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
        t * t * (3.0 - 2.0 * t)
    }
}

/// Winkel-Hilfsfunktionen
pub mod angles {
    use super::*;

    /// Konvertiert Grad zu Radiant
    pub fn deg_to_rad(degrees: f32) -> f32 {
        degrees * PI / 180.0
    }

    /// Konvertiert Radiant zu Grad
    pub fn rad_to_deg(radians: f32) -> f32 {
        radians * 180.0 / PI
    }

    /// Normalisiert einen Winkel auf [0, 2π)
    pub fn normalize_angle(angle: f32) -> f32 {
        let mut result = angle % TAU;
        if result < 0.0 {
            result += TAU;
        }
        result
    }

    /// Normalisiert einen Winkel auf [-π, π)
    pub fn normalize_angle_signed(angle: f32) -> f32 {
        let mut result = angle % TAU;
        if result > PI {
            result -= TAU;
        } else if result < -PI {
            result += TAU;
        }
        result
    }

    /// Berechnet den kleinsten Winkelunterschied zwischen zwei Winkeln
    pub fn angle_difference(a: f32, b: f32) -> f32 {
        let diff = normalize_angle_signed(b - a);
        diff
    }

    /// Interpoliert zwischen zwei Winkeln (kürzester Weg)
    pub fn lerp_angle(a: f32, b: f32, t: f32) -> f32 {
        a + angle_difference(a, b) * t
    }
}

/// Geometrische Hilfsfunktionen
pub mod geometry {
    use super::*;

    /// Berechnet den Abstand zwischen zwei Punkten
    pub fn distance(p1: Point2D, p2: Point2D) -> f32 {
        ((p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)).sqrt()
    }

    /// Berechnet das Skalarprodukt zweier 2D-Vektoren
    pub fn dot_product(a: Point2D, b: Point2D) -> f32 {
        a.x * b.x + a.y * b.y
    }

    /// Berechnet das Kreuzprodukt zweier 2D-Vektoren (Skalar)
    pub fn cross_product_2d(a: Point2D, b: Point2D) -> f32 {
        a.x * b.y - a.y * b.x
    }

    /// Berechnet den Winkel eines 2D-Vektors
    pub fn vector_angle(v: Point2D) -> f32 {
        v.y.atan2(v.x)
    }

    /// Rotiert einen 2D-Vektor um einen Winkel
    pub fn rotate_vector_2d(v: Point2D, angle: f32) -> Point2D {
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        Point2D::new(v.x * cos_a - v.y * sin_a, v.x * sin_a + v.y * cos_a)
    }

    /// Projiziert einen Punkt auf eine Linie
    pub fn project_point_on_line(
        point: Point2D,
        line_start: Point2D,
        line_end: Point2D,
    ) -> Point2D {
        let line_vec = line_end - line_start;
        let point_vec = point - line_start;

        let line_length_sq = line_vec.length_squared();
        if line_length_sq < constants::EPSILON {
            return line_start; // Linie ist ein Punkt
        }

        let projection_length = dot_product(point_vec, line_vec) / line_length_sq;
        line_start + line_vec * projection_length
    }

    /// Berechnet den Abstand von einem Punkt zu einer Linie
    pub fn point_line_distance(point: Point2D, line_start: Point2D, line_end: Point2D) -> f32 {
        let projection = project_point_on_line(point, line_start, line_end);
        distance(point, projection)
    }

    /// Prüft ob ein Punkt auf einer Linie liegt
    pub fn point_on_line(
        point: Point2D,
        line_start: Point2D,
        line_end: Point2D,
        tolerance: f32,
    ) -> bool {
        point_line_distance(point, line_start, line_end) < tolerance
    }

    /// Berechnet den Flächeninhalt eines Dreiecks
    pub fn triangle_area(p1: Point2D, p2: Point2D, p3: Point2D) -> f32 {
        0.5 * (cross_product_2d(p2 - p1, p3 - p1)).abs()
    }

    /// Prüft ob ein Punkt in einem Dreieck liegt (Barycentric coordinates)
    pub fn point_in_triangle(point: Point2D, a: Point2D, b: Point2D, c: Point2D) -> bool {
        let v0 = c - a;
        let v1 = b - a;
        let v2 = point - a;

        let dot00 = dot_product(v0, v0);
        let dot01 = dot_product(v0, v1);
        let dot02 = dot_product(v0, v2);
        let dot11 = dot_product(v1, v1);
        let dot12 = dot_product(v1, v2);

        let inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        let u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
        let v = (dot00 * dot12 - dot01 * dot02) * inv_denom;

        u >= 0.0 && v >= 0.0 && u + v <= 1.0
    }
}

/// Numerische Algorithmen
pub mod numerical {
    use super::*;

    /// Newtons Methode für Wurzel-Approximation
    pub fn sqrt_newton(value: f32, initial_guess: f32, iterations: usize) -> f32 {
        let mut x = initial_guess;
        for _ in 0..iterations {
            x = 0.5 * (x + value / x);
        }
        x
    }

    /// Berechnet Fakultät (iterativ)
    pub fn factorial(n: u32) -> u64 {
        (1..=n as u64).product()
    }

    /// Berechnet Binomialkoeffizient "n choose k"
    pub fn binomial_coefficient(n: u32, k: u32) -> u64 {
        if k > n {
            return 0;
        }
        if k == 0 || k == n {
            return 1;
        }

        let k = k.min(n - k); // Nutze Symmetrie für Effizienz
        let mut result = 1u64;

        for i in 0..k {
            result = result * (n - i) as u64 / (i + 1) as u64;
        }

        result
    }

    /// Gaussian (Normal) Distribution
    pub fn gaussian(x: f32, mean: f32, std_dev: f32) -> f32 {
        let variance = std_dev * std_dev;
        let coefficient = 1.0 / (std_dev * (2.0 * PI).sqrt());
        let exponent = -0.5 * ((x - mean) / std_dev).powi(2);
        coefficient * exponent.exp()
    }

    /// Sigmoid function
    pub fn sigmoid(x: f32) -> f32 {
        1.0 / (1.0 + (-x).exp())
    }
}

/// Random utilities (erweitert vorhandene rand-Funktionalität)
pub mod random {
    use super::*;
    use rand::Rng;

    /// Generiert zufälligen Punkt in einem Rechteck
    pub fn random_point_in_rect(min: Point2D, max: Point2D) -> Point2D {
        let mut rng = rand::rng();
        Point2D::new(
            rng.random_range(min.x..=max.x),
            rng.random_range(min.y..=max.y),
        )
    }

    /// Generiert zufälligen Punkt in einem Kreis
    pub fn random_point_in_circle(center: Point2D, radius: f32) -> Point2D {
        let mut rng = rand::rng();
        let angle = rng.random_range(0.0..TAU);
        let r = radius * rng.random::<f32>().sqrt(); // Gleichmäßige Verteilung

        Point2D::new(center.x + r * angle.cos(), center.y + r * angle.sin())
    }

    /// Generiert zufälligen Punkt auf einem Kreis (nur Rand)
    pub fn random_point_on_circle(center: Point2D, radius: f32) -> Point2D {
        let mut rng = rand::rng();
        let angle = rng.random_range(0.0..TAU);

        Point2D::new(
            center.x + radius * angle.cos(),
            center.y + radius * angle.sin(),
        )
    }

    /// Gaussian distributed random number (Box-Muller transform)
    pub fn random_gaussian(mean: f32, std_dev: f32) -> f32 {
        use std::sync::Mutex;
        static CACHED: Mutex<Option<f32>> = Mutex::new(None);

        let mut cached = CACHED.lock().unwrap();
        if let Some(value) = cached.take() {
            return mean + std_dev * value;
        }

        let mut rng = rand::rng();
        let u1 = rng.random::<f32>();
        let u2 = rng.random::<f32>();

        let mag = std_dev * (-2.0 * u1.ln()).sqrt();
        let z0 = mag * (TAU * u2).cos();
        let z1 = mag * (TAU * u2).sin();

        *cached = Some(z1);
        mean + z0
    }
}

/// Easing-Funktionen für Animationen
pub mod easing {
    use super::*;

    /// Linear easing (keine Beschleunigung)
    pub fn linear(t: f32) -> f32 {
        t
    }

    /// Quadratic ease-in
    pub fn ease_in_quad(t: f32) -> f32 {
        t * t
    }

    /// Quadratic ease-out
    pub fn ease_out_quad(t: f32) -> f32 {
        1.0 - (1.0 - t) * (1.0 - t)
    }

    /// Quadratic ease-in-out
    pub fn ease_in_out_quad(t: f32) -> f32 {
        if t < 0.5 {
            2.0 * t * t
        } else {
            1.0 - 2.0 * (1.0 - t) * (1.0 - t)
        }
    }

    /// Cubic ease-in
    pub fn ease_in_cubic(t: f32) -> f32 {
        t * t * t
    }

    /// Cubic ease-out
    pub fn ease_out_cubic(t: f32) -> f32 {
        1.0 - (1.0 - t).powi(3)
    }

    /// Elastic ease-out
    pub fn ease_out_elastic(t: f32) -> f32 {
        if t == 0.0 || t == 1.0 {
            t
        } else {
            let c4 = TAU / 3.0;
            2.0_f32.powf(-10.0 * t) * ((t * 10.0 - 0.75) * c4).sin() + 1.0
        }
    }

    /// Bounce ease-out
    pub fn ease_out_bounce(t: f32) -> f32 {
        let n1 = 7.5625;
        let d1 = 2.75;

        if t < 1.0 / d1 {
            n1 * t * t
        } else if t < 2.0 / d1 {
            let t = t - 1.5 / d1;
            n1 * t * t + 0.75
        } else if t < 2.5 / d1 {
            let t = t - 2.25 / d1;
            n1 * t * t + 0.9375
        } else {
            let t = t - 2.625 / d1;
            n1 * t * t + 0.984375
        }
    }
}
