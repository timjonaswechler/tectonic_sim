/// Wie der Abstand angegeben wird
#[derive(Clone, Copy, Debug)]
pub enum Dist {
    /// δ bereits als Winkel in **Radiant**
    AngleRad(f64),
    /// δ als Winkel in **Grad**
    AngleDeg(f64),
    /// Bogenlänge s entlang der Kugeloberfläche
    Arc(f64),
    /// Gerade Schnurlänge L (Chord) durchs Kugelinnere
    Chord(f64),
}

/// Kreuzprodukt zweier Vektoren (3-D)
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Skalarprodukt
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Norm eines Vektors
fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}

/// Vektor auf Länge 1 normieren
fn normalize(mut v: [f64; 3]) -> [f64; 3] {
    let n = norm(v);
    v.iter_mut().for_each(|x| *x /= n);
    v
}

use rand::Rng;
use std::f64::consts::PI; // Added for random number generation

use bevy::math::Vec3;

/// Zufälliger Punkt **innerhalb** des Abstands `dist` zu `c`
pub fn random_point_within_distance(c: Vec3, dist: Dist, radius: f64) -> Vec3 {
    // ---------- 1) max. Winkel δmax bestimmen -----------------
    let delta_max = match dist {
        Dist::AngleRad(rad) => rad,
        Dist::AngleDeg(deg) => deg.to_radians(),
        Dist::Arc(s) => s / radius,
        Dist::Chord(L) => 2.0 * (L / (2.0 * radius)).asin(),
    };

    let mut rng = rand::rng();
    let c_arr = [c.x as f64, c.y as f64, c.z as f64];

    // ---------- 2) Orthonormale Basis im Tangentialraum -------
    // irgendein Vektor, der nicht (fast) parallel zu c_arr ist
    let mut v;
    loop {
        v = [
            rng.random_range(-1.0..=1.0),
            rng.random_range(-1.0..=1.0),
            rng.random_range(-1.0..=1.0),
        ];
        if dot(v, c_arr).abs() < 0.999 {
            break;
        }
    }
    let axis1 = normalize(cross(c_arr, v));
    let axis2 = cross(axis1, c_arr); // bereits normiert

    // ---------- 3) zufälligen Winkel (θ, φ) wählen ------------
    // gleichmäßige Fläche ⇒ cos θ gleichmäßig
    let cos_theta = rng.random_range(delta_max.cos()..=1.0);
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
    let phi = rng.random_range(0.0..(2.0 * PI));

    // ---------- 4) Punkt aus spherical-to-Cartesian -----------
    let p = [
        cos_theta * c_arr[0] + sin_theta * (axis1[0] * phi.cos() + axis2[0] * phi.sin()),
        cos_theta * c_arr[1] + sin_theta * (axis1[1] * phi.cos() + axis2[1] * phi.sin()),
        cos_theta * c_arr[2] + sin_theta * (axis1[2] * phi.cos() + axis2[2] * phi.sin()),
    ];

    Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32)
}
