// src/math/utils.rs

/// Mathematische Konstanten
pub mod constants {
    pub const EPSILON: f32 = 1e-6;
    pub const EPSILON_F64: f64 = 1e-10; // Behalten für den Fall, dass wir es doch brauchen
    pub const EPSILON_SQUARED: f32 = EPSILON * EPSILON; // Für Vergleiche mit Längen
    pub const GOLDEN_RATIO: f32 = 1.618033988749895;
    pub const SQRT_2: f32 = 1.4142135623730951;
    pub const SQRT_3: f32 = 1.7320508075688772;
    pub const TAU: f32 = std::f32::consts::TAU; // Explizit TAU von std verwenden
    pub const PI: f32 = std::f32::consts::PI;
    pub const PI_OVER_2: f32 = std::f32::consts::PI / 2.0;
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
        value.clamp(min, max) // Standard-clamp verwenden
    }

    /// Lineare Interpolation
    pub fn lerp(a: f32, b: f32, t: f32) -> f32 {
        a + (b - a) * t
    }

    /// Inverse lineare Interpolation
    pub fn inverse_lerp(a: f32, b: f32, value: f32) -> f32 {
        if nearly_equal(a, b) {
            0.0 // Oder je nach Anwendungsfall NaN oder panic
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
    use super::constants::{PI, TAU}; // Importiere PI und TAU aus unserem constants Modul

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
        normalize_angle_signed(b - a)
    }

    /// Interpoliert zwischen zwei Winkeln (kürzester Weg)
    pub fn lerp_angle(a: f32, b: f32, t: f32) -> f32 {
        normalize_angle(a + angle_difference(a, b) * t) // normalize_angle, da das Ergebnis außerhalb [0, 2pi) sein kann
    }
}

/// Geometrische Hilfsfunktionen (einfach, ohne komplexe Strukturen)
pub mod simple_geometry {
    use crate::math::utils::constants;
    use bevy::math::Vec2;

    /// Berechnet den Abstand zwischen zwei Punkten
    pub fn distance_sq(p1: Vec2, p2: Vec2) -> f32 {
        (p2.x - p1.x).powi(2) + (p2.y - p1.y).powi(2)
    }

    pub fn distance(p1: Vec2, p2: Vec2) -> f32 {
        distance_sq(p1, p2).sqrt()
    }

    /// Berechnet das Skalarprodukt zweier 2D-Vektoren
    pub fn dot_product(a: Vec2, b: Vec2) -> f32 {
        a.x * b.x + a.y * b.y
    }

    /// Berechnet das Kreuzprodukt zweier 2D-Vektoren (Skalar)
    pub fn cross_product_2d(a: Vec2, b: Vec2) -> f32 {
        a.x * b.y - a.y * b.x
    }

    /// Berechnet den Winkel eines 2D-Vektors
    pub fn vector_angle(v: Vec2) -> f32 {
        v.y.atan2(v.x)
    }

    /// Rotiert einen 2D-Vektor um einen Winkel
    pub fn rotate_vector_2d(v: Vec2, angle_rad: f32) -> Vec2 {
        let cos_a = angle_rad.cos();
        let sin_a = angle_rad.sin();
        Vec2::new(v.x * cos_a - v.y * sin_a, v.x * sin_a + v.y * cos_a)
    }

    /// Projiziert einen Punkt auf eine Linie
    pub fn project_point_on_line(point: Vec2, line_start: Vec2, line_end: Vec2) -> Vec2 {
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
    pub fn point_line_distance(point: Vec2, line_start: Vec2, line_end: Vec2) -> f32 {
        let projection = project_point_on_line(point, line_start, line_end);
        distance(point, projection)
    }

    /// Prüft ob ein Punkt auf einer Linie liegt (innerhalb einer Toleranz)
    pub fn point_on_line_segment(
        point: Vec2,
        segment_start: Vec2,
        segment_end: Vec2,
        tolerance: f32,
    ) -> bool {
        let dist_to_line = point_line_distance(point, segment_start, segment_end);
        if dist_to_line >= tolerance {
            return false;
        }
        // Prüfe, ob der Punkt zwischen den Endpunkten des Segments liegt
        let dot = (point - segment_start).dot(segment_end - segment_start);
        if dot < 0.0 {
            return false;
        }
        let squared_length = (segment_end - segment_start).length_squared();
        if dot > squared_length {
            return false;
        }
        true
    }

    /// Berechnet den Flächeninhalt eines Dreiecks
    pub fn triangle_area(p1: Vec2, p2: Vec2, p3: Vec2) -> f32 {
        0.5 * (cross_product_2d(p2 - p1, p3 - p1)).abs()
    }

    /// Prüft ob ein Punkt in einem Dreieck liegt (Barycentric coordinates)
    pub fn point_in_triangle(point: Vec2, a: Vec2, b: Vec2, c: Vec2) -> bool {
        let v0 = c - a;
        let v1 = b - a;
        let v2 = point - a;

        let dot00 = dot_product(v0, v0);
        let dot01 = dot_product(v0, v1);
        let dot02 = dot_product(v0, v2);
        let dot11 = dot_product(v1, v1);
        let dot12 = dot_product(v1, v2);

        let inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01 + constants::EPSILON); // Add EPSILON to avoid div by zero
        let u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
        let v = (dot00 * dot12 - dot01 * dot02) * inv_denom;

        u >= -constants::EPSILON && v >= -constants::EPSILON && u + v <= 1.0 + constants::EPSILON
    }
}

/// Numerische Algorithmen
pub mod numerical {
    use super::constants::PI; // Stelle sicher, dass PI hier verfügbar ist

    /// Newtons Methode für Wurzel-Approximation
    pub fn sqrt_newton(value: f32, initial_guess: f32, iterations: usize) -> f32 {
        let mut x = initial_guess;
        if x <= 0.0 {
            x = 1.0;
        } // Sicherstellen, dass initial_guess positiv ist
        for _ in 0..iterations {
            if x.abs() < 1e-10 {
                return 0.0;
            } // Vermeide Division durch Null
            x = 0.5 * (x + value / x);
        }
        x
    }

    /// Berechnet Fakultät (iterativ)
    pub fn factorial(n: u32) -> u64 {
        if n > 20 {
            return u64::MAX;
        } // Vermeide Overflow für u64
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
        if k > n / 2 {
            // Nutze Symmetrie für Effizienz
            return binomial_coefficient(n, n - k);
        }

        let mut result = 1u64;
        for i in 0..k {
            match result.checked_mul((n - i) as u64) {
                Some(val) => result = val / (i + 1) as u64,
                None => return u64::MAX, // Overflow
            }
        }
        result
    }

    /// Gaussian (Normal) Distribution PDF
    pub fn gaussian_pdf(x: f32, mean: f32, std_dev: f32) -> f32 {
        if std_dev <= 0.0 {
            return 0.0;
        } // std_dev muss positiv sein
        let variance = std_dev * std_dev;
        let coefficient = 1.0 / (std_dev * (2.0 * PI).sqrt());
        let exponent = -0.5 * ((x - mean) * (x - mean)) / variance;
        coefficient * exponent.exp()
    }

    /// Sigmoid function
    pub fn sigmoid(x: f32) -> f32 {
        1.0 / (1.0 + (-x).exp())
    }
}

/// Random utilities (erweitert vorhandene rand-Funktionalität)
pub mod random {
    use crate::math::utils::constants::TAU;
    use bevy::math::Vec2;
    use rand::Rng; // Standard Rng Trait

    /// Generiert zufälligen Punkt in einem Rechteck
    pub fn random_point_in_rect(min: Vec2, max: Vec2, rng: &mut impl Rng) -> Vec2 {
        Vec2::new(
            rng.random_range(min.x..=max.x),
            rng.random_range(min.y..=max.y),
        )
    }

    /// Generiert zufälligen Punkt in einem Kreis
    pub fn random_point_in_circle(center: Vec2, radius: f32, rng: &mut impl Rng) -> Vec2 {
        let angle = rng.random_range(0.0..TAU);
        let r = radius * rng.random::<f32>().sqrt(); // Gleichmäßige Verteilung
        Vec2::new(center.x + r * angle.cos(), center.y + r * angle.sin())
    }

    /// Generiert zufälligen Punkt auf einem Kreis (nur Rand)
    pub fn random_point_on_circle(center: Vec2, radius: f32, rng: &mut impl Rng) -> Vec2 {
        let angle = rng.random_range(0.0..TAU);
        Vec2::new(
            center.x + radius * angle.cos(),
            center.y + radius * angle.sin(),
        )
    }

    /// Gaussian distributed random number (Box-Muller transform)
    pub fn random_gaussian(mean: f32, std_dev: f32, rng: &mut impl Rng) -> f32 {
        // Für Box-Muller braucht man zwei Zufallszahlen.
        // Die static CACHED Version ist problematisch, besonders mit externem RNG.
        // Einfachere Version ohne Caching:
        let u1: f32 = rng.random(); // Sollte (0, 1] sein, gen() gibt [0,1)
        let u2: f32 = rng.random();

        let mag = std_dev * (-2.0 * u1.ln()).sqrt();
        let z0 = mag * (TAU * u2).cos();
        // z1 = mag * (TAU * u2).sin(); // Zweiter Wert, wird hier nicht genutzt

        mean + z0
    }
}

/// Easing-Funktionen für Animationen
pub mod easing {
    use crate::math::utils::constants::TAU; // TAU für Elastic

    /// Linear easing (keine Beschleunigung)
    pub fn linear(t: f32) -> f32 {
        t.clamp(0.0, 1.0)
    }

    /// Quadratic ease-in
    pub fn ease_in_quad(t: f32) -> f32 {
        let t = t.clamp(0.0, 1.0);
        t * t
    }

    /// Quadratic ease-out
    pub fn ease_out_quad(t: f32) -> f32 {
        let t = t.clamp(0.0, 1.0);
        1.0 - (1.0 - t) * (1.0 - t)
    }

    /// Quadratic ease-in-out
    pub fn ease_in_out_quad(t: f32) -> f32 {
        let t = t.clamp(0.0, 1.0);
        if t < 0.5 {
            2.0 * t * t
        } else {
            1.0 - (-2.0 * t + 2.0).powi(2) / 2.0 // Korrigierte Formel
        }
    }

    /// Cubic ease-in
    pub fn ease_in_cubic(t: f32) -> f32 {
        let t = t.clamp(0.0, 1.0);
        t * t * t
    }

    /// Cubic ease-out
    pub fn ease_out_cubic(t: f32) -> f32 {
        let t = t.clamp(0.0, 1.0);
        1.0 - (1.0 - t).powi(3)
    }

    /// Elastic ease-out
    pub fn ease_out_elastic(t: f32) -> f32 {
        let t = t.clamp(0.0, 1.0);
        if t == 0.0 {
            return 0.0;
        }
        if t == 1.0 {
            return 1.0;
        }
        let c4 = TAU / 3.0;
        2.0_f32.powf(-10.0 * t) * ((t * 10.0 - 0.75) * c4).sin() + 1.0
    }

    /// Bounce ease-out
    pub fn ease_out_bounce(t: f32) -> f32 {
        let t = t.clamp(0.0, 1.0);
        let n1 = 7.5625;
        let d1 = 2.75;

        if t < 1.0 / d1 {
            n1 * t * t
        } else if t < 2.0 / d1 {
            let t_adj = t - 1.5 / d1;
            n1 * t_adj * t_adj + 0.75
        } else if t < 2.5 / d1 {
            let t_adj = t - 2.25 / d1;
            n1 * t_adj * t_adj + 0.9375
        } else {
            let t_adj = t - 2.625 / d1;
            n1 * t_adj * t_adj + 0.984375
        }
    }
}
