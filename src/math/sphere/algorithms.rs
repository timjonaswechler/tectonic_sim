// src/math/util.rs
use bevy::prelude::Vec3;

/// Konvertiert geografische Koordinaten (Breitengrad, Längengrad) in kartesische 3D-Koordinaten.
/// Die Y-Achse ist die Rotationsachse der Kugel (zeigt zum Nordpol).
/// Die X-Achse geht durch den Punkt (lat=0, lon=0).
/// Die Z-Achse geht durch den Punkt (lat=0, lon=90 Grad Ost).
///
/// # Arguments
/// * `latitude_rad` - Breitengrad in Radiant (-PI/2 bis PI/2).
/// * `longitude_rad` - Längengrad in Radiant (-PI bis PI).
/// * `radius` - Der Radius der Kugel.
pub fn latlon_to_cartesian(latitude_rad: f32, longitude_rad: f32, radius: f32) -> Vec3 {
    let cos_lat = latitude_rad.cos();
    // Standardkonvention: x = r * cos(lat) * cos(lon), y = r * cos(lat) * sin(lon), z = r * sin(lat)
    // Hier verwenden wir y als "oben":
    let x = radius * cos_lat * longitude_rad.cos();
    let y = radius * latitude_rad.sin();
    let z = radius * cos_lat * longitude_rad.sin(); // Positive Z für östliche Longituden
    Vec3::new(x, y, z)
}

/// Konvertiert kartesische 3D-Koordinaten (auf einer Kugel) in geografische Koordinaten.
///
/// # Arguments
/// * `point` - Der `Vec3` Punkt auf der Kugeloberfläche.
///
/// # Returns
/// Ein Tupel `(latitude_rad: f32, longitude_rad: f32)`.
pub fn cartesian_to_latlon(point: Vec3) -> (f32, f32) {
    // Stelle sicher, dass der Punkt normalisiert ist, um den Radius-Faktor zu eliminieren,
    // da asin und atan2 auf einem Einheitsvektor operieren sollten.
    let p_norm = point.normalize();
    let latitude_rad = p_norm.y.asin(); // asin gibt Werte im Bereich [-PI/2, PI/2]
    let longitude_rad = p_norm.z.atan2(p_norm.x); // atan2(y, x) gibt den Winkel des Vektors (x,y)
    (latitude_rad, longitude_rad)
}

/// Erzeugt einen zufälligen Punkt auf der Oberfläche einer Einheitskugel.
/// Für eine Kugel mit spezifischem Radius: `random_point_on_unit_sphere() * radius`.
/// Nutzt die Methode von Marsaglia (1972) für gleichmäßige Verteilung, angepasst.
pub fn random_point_on_unit_sphere() -> Vec3 {
    use rand::Rng;
    let mut rng = rand::rng();
    // Alternativer, oft einfacher zu implementierender Algorithmus:
    // Generiere drei Normalverteilte Zufallszahlen (x,y,z) und normalisiere den Vektor (x,y,z).
    // Hier die Version mit gleichverteilten Zahlen und Zurückweisung:
    loop {
        let x1: f32 = rng.random_range(-1.0..=1.0);
        let x2: f32 = rng.random_range(-1.0..=1.0);
        let x3: f32 = rng.random_range(-1.0..=1.0); // Hinzugefügt für 3D Kugel
        let r_squared = x1 * x1 + x2 * x2 + x3 * x3;
        // Punkte innerhalb einer Kugel, die nicht der Nullvektor sind
        if r_squared < 1.0 && r_squared > 1e-6 {
            // Projiziere den Punkt auf die Oberfläche der Einheitskugel
            let inv_len = 1.0 / r_squared.sqrt();
            return Vec3::new(x1 * inv_len, x2 * inv_len, x3 * inv_len);
        }
    }
}

/// Erzeugt zwei orthogonale Einheitsvektoren, die in der Tangentialebene
/// an dem gegebenen `surface_normal` (ein normalisierter Vektor vom Kugelzentrum zum Punkt) liegen.
///
/// - `tangent1` ist grob "östlich" ausgerichtet (orthogonal zur globalen Y-Achse, wenn möglich).
/// - `tangent2` ist grob "nördlich" in der Tangentialebene ausgerichtet.
///
/// # Arguments
/// * `surface_normal` - Der normalisierte Vektor, der vom Kugelzentrum zum Punkt auf der Oberfläche zeigt.
///
/// # Returns
/// Ein Tupel `(tangent1: Vec3, tangent2: Vec3)`.
pub fn create_orthogonal_tangent_vectors(surface_normal: Vec3) -> (Vec3, Vec3) {
    // Wähle einen "up_vector", der nicht kollinear mit surface_normal ist.
    // Die globale Y-Achse ist meistens eine gute Wahl, außer wenn surface_normal selbst Y oder -Y ist.
    let global_up = if (surface_normal.y.abs() - 1.0).abs() < 1e-6 {
        // surface_normal ist (fast) (0,1,0) oder (0,-1,0), also Nord- oder Südpol.
        // Wähle X-Achse als Referenz, um eine stabile Tangente zu bekommen.
        Vec3::X
    } else {
        Vec3::Y
    };

    let tangent1 = global_up.cross(surface_normal).normalize();
    // Wenn global_up und surface_normal kollinear waren, ist tangent1 (0,0,0).
    // Dies sollte durch die obige if-Bedingung behandelt werden. Falls nicht, Fallback:
    let tangent1 = if tangent1.length_squared() < 1e-6 {
        // Fallback, falls die erste Wahl für global_up doch (fast) kollinear war.
        // Passiert, wenn surface_normal z.B. (0,1,0) und global_up (0,1,0)
        // In diesem Fall ist z.B. Vec3::X eine gute Wahl für die erste Tangente, wenn der Punkt der Nordpol ist.
        // Die Logik oben sollte das aber abfangen. Hier eine robustere Alternative für alle Fälle:
        let mut alt_axis = Vec3::X;
        if surface_normal.cross(alt_axis).length_squared() < 1e-6 {
            // Wenn surface_normal auch X ist
            alt_axis = Vec3::Y; // Dann nimm Y
        }
        surface_normal.cross(alt_axis).normalize() // Erzeugt eine Tangente, die senkrecht zu surface_normal ist
    } else {
        tangent1
    };

    let tangent2 = surface_normal.cross(tangent1).normalize(); // tangent1 ist bereits normalisiert

    (tangent1, tangent2)
}

/// Rotiert einen Vektor `point_to_rotate` um eine Achse `axis` (muss normalisiert sein)
/// um einen Winkel `angle_rad`.
pub fn rotate_vector_around_axis(point_to_rotate: Vec3, axis: Vec3, angle_rad: f32) -> Vec3 {
    use nalgebra::{Rotation3, Unit, Vector3 as NVec3}; // Alias für Klarheit

    // Stelle sicher, dass die Achse normalisiert ist, falls nicht schon geschehen.
    // Eine Achse mit Länge Null führt zu Problemen.
    if axis.length_squared() < 1e-9 {
        // Prüfe auf nahe Null
        bevy::log::warn!("Rotationsachse ist nahe Null, keine Rotation durchgeführt.");
        return point_to_rotate;
    }
    let unit_axis = Unit::new_normalize(NVec3::new(axis.x, axis.y, axis.z));
    let rotation = Rotation3::from_axis_angle(&unit_axis, angle_rad);

    let nalgebra_point = NVec3::new(point_to_rotate.x, point_to_rotate.y, point_to_rotate.z);
    let rotated_nalgebra_point = rotation * nalgebra_point;

    Vec3::new(
        rotated_nalgebra_point.x,
        rotated_nalgebra_point.y,
        rotated_nalgebra_point.z,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::FRAC_PI_2;
    use std::f32::consts::PI;
    #[test]
    fn test_latlon_cartesian_conversion() {
        let radius = 100.0;
        // Nordpol
        let (lat_np, lon_np) = (FRAC_PI_2, 0.0); // lat = 90, lon = 0
        let cart_np = latlon_to_cartesian(lat_np, lon_np, radius);
        assert!((cart_np.x - 0.0).abs() < 1e-5, "NP X: {}", cart_np.x);
        assert!((cart_np.y - radius).abs() < 1e-5, "NP Y: {}", cart_np.y);
        assert!((cart_np.z - 0.0).abs() < 1e-5, "NP Z: {}", cart_np.z);
        let (r_lat_np, r_lon_np) = cartesian_to_latlon(cart_np);
        assert!((r_lat_np - lat_np).abs() < 1e-5, "NP rLat: {}", r_lat_np);
        // Longitude am Pol ist instabil/mehrdeutig, kann 0 oder PI sein.
        // assert!((r_lon_np - lon_np).abs() < 1e-5 || (r_lon_np - lon_np + 2.0*PI).abs() < 1e-5 || (r_lon_np - lon_np - 2.0*PI).abs() < 1e-5 );

        // Punkt auf Äquator, an lon=0 (positiv X)
        let (lat_eq_x, lon_eq_x) = (0.0, 0.0);
        let cart_eq_x = latlon_to_cartesian(lat_eq_x, lon_eq_x, radius);
        assert!((cart_eq_x.x - radius).abs() < 1e-5);
        assert!((cart_eq_x.y - 0.0).abs() < 1e-5);
        assert!((cart_eq_x.z - 0.0).abs() < 1e-5);
        let (r_lat_eq_x, r_lon_eq_x) = cartesian_to_latlon(cart_eq_x);
        assert!((r_lat_eq_x - lat_eq_x).abs() < 1e-5);
        assert!((r_lon_eq_x - lon_eq_x).abs() < 1e-5);

        // Punkt auf Äquator, an lon=PI/2 (positiv Z)
        let (lat_eq_z, lon_eq_z) = (0.0, FRAC_PI_2);
        let cart_eq_z = latlon_to_cartesian(lat_eq_z, lon_eq_z, radius);
        assert!((cart_eq_z.x - 0.0).abs() < 1e-5, "EQ_Z X: {}", cart_eq_z.x);
        assert!((cart_eq_z.y - 0.0).abs() < 1e-5, "EQ_Z Y: {}", cart_eq_z.y);
        assert!(
            (cart_eq_z.z - radius).abs() < 1e-5,
            "EQ_Z Z: {}",
            cart_eq_z.z
        );
        let (r_lat_eq_z, r_lon_eq_z) = cartesian_to_latlon(cart_eq_z);
        assert!(
            (r_lat_eq_z - lat_eq_z).abs() < 1e-5,
            "EQ_Z rLat: {}",
            r_lat_eq_z
        );
        assert!(
            (r_lon_eq_z - lon_eq_z).abs() < 1e-5,
            "EQ_Z rLon: {}",
            r_lon_eq_z
        );
    }

    #[test]
    fn test_random_point_on_sphere() {
        for _ in 0..100 {
            let p = random_point_on_unit_sphere();
            assert!(
                (p.length() - 1.0).abs() < 1e-6,
                "Punkt nicht auf Einheitskugel: {:?}, Länge: {}",
                p,
                p.length()
            );
        }
    }

    #[test]
    fn test_rotate_vector() {
        let point = Vec3::X; // (1,0,0)
        let axis = Vec3::Y; // Rotiere um Y-Achse
        let angle = PI / 2.0; // 90 Grad

        let rotated = rotate_vector_around_axis(point, axis, angle);
        assert!((rotated.x - 0.0).abs() < 1e-6, "RX: {}", rotated.x);
        assert!((rotated.y - 0.0).abs() < 1e-6, "RY: {}", rotated.y);
        assert!((rotated.z - (-1.0)).abs() < 1e-6, "RZ: {}", rotated.z);

        let point2 = Vec3::new(1.0, 1.0, 0.0);
        let axis2 = Vec3::Z; // Rotiere um Z
        let angle2 = PI; // 180 Grad
        let rotated2 = rotate_vector_around_axis(point2, axis2, angle2);
        assert!((rotated2.x - (-1.0)).abs() < 1e-6, "R2X: {}", rotated2.x);
        assert!((rotated2.y - (-1.0)).abs() < 1e-6, "R2Y: {}", rotated2.y);
        assert!((rotated2.z - 0.0).abs() < 1e-6, "R2Z: {}", rotated2.z);
    }
}
