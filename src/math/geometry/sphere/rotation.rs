use crate::math::geometry::sphere::coordinates::{CoordinateConverter, GeographicCoordinates};
use crate::math::{
    error::*,
    utils::*, // constants, comparison
};

use bevy::math::{Quat as BevyQuat, Vec3}; // Bevy's Vec3 und Quat

/// Hilfsfunktionen für Rotationen auf einer Kugel, primär unter Verwendung von `bevy::math::Quat`.
pub struct SphereRotation;

impl SphereRotation {
    /// Rotiert einen Punkt auf der Oberfläche einer Kugel um eine gegebene Achse und einen Winkel.
    /// `point_on_sphere` sollte idealerweise bereits auf der Kugel liegen.
    /// `rotation_axis` muss nicht normalisiert sein.
    /// `angle_rad` ist der Rotationswinkel in Radiant.
    /// `sphere_radius` ist der Radius der Kugel.
    pub fn rotate_point_on_sphere(
        point_on_sphere: Vec3,
        rotation_axis: Vec3,
        angle_rad: f32,
        sphere_radius: f32,
    ) -> MathResult<Vec3> {
        let axis_norm = rotation_axis.normalize_or_zero();
        if axis_norm.length_squared() < constants::EPSILON * constants::EPSILON {
            return Err(MathError::InvalidConfiguration {
                message: "Rotation axis cannot be a zero vector.".to_string(),
            });
        }

        let rotation_quat = BevyQuat::from_axis_angle(axis_norm, angle_rad);
        let rotated_point = rotation_quat * point_on_sphere; // Bevy's Quat kann Vec3 direkt rotieren

        // Stelle sicher, dass der Punkt auf der Kugel bleibt (falls der ursprüngliche Punkt nicht exakt drauf war)
        Ok(CoordinateConverter::project_to_sphere(
            rotated_point,
            sphere_radius,
        ))
    }

    /// Führt eine sphärische lineare Interpolation (Slerp) eines Punktes auf einer Kugel
    /// von einem Startpunkt zu einem Endpunkt durch.
    /// `start_point_on_sphere` und `end_point_on_sphere` sollten auf der Kugeloberfläche liegen.
    /// `t` ist der Interpolationsfaktor [0, 1].
    pub fn slerp_point_on_sphere(
        start_point_on_sphere: Vec3,
        end_point_on_sphere: Vec3,
        t: f32,
        // sphere_radius: f32, // Radius ist implizit durch die Punkte gegeben
    ) -> MathResult<Vec3> {
        if t <= 0.0 {
            return Ok(start_point_on_sphere);
        }
        if t >= 1.0 {
            return Ok(end_point_on_sphere);
        }

        // Erstelle Quaternions, die vom Ursprung zu den Punkten zeigen (als reine Rotationen von einer Referenzachse).
        // Einfacher ist es, die direkte Rotation von start zu end zu finden.
        let rotation_from_start_to_end = BevyQuat::from_rotation_arc(
            start_point_on_sphere.normalize_or_zero(),
            end_point_on_sphere.normalize_or_zero(),
        );

        // Slerp zwischen Identitätsrotation und der vollen Rotation
        let interpolated_rotation = BevyQuat::IDENTITY.slerp(rotation_from_start_to_end, t);
        let interpolated_point = interpolated_rotation * start_point_on_sphere;

        // Radius beibehalten (falls die Punkte unterschiedliche Längen hatten, wird hier der von start_point genommen)
        Ok(CoordinateConverter::project_to_sphere(
            interpolated_point,
            start_point_on_sphere.length(),
        ))
    }

    /// Rotiert einen Punkt auf einer Kugel basierend auf Änderungen in geografischen Koordinaten (Delta Latitude/Longitude).
    /// `point_on_sphere` ist der zu rotierende Punkt.
    /// `delta_latitude_rad` und `delta_longitude_rad` sind die Änderungen in Radiant.
    /// `sphere_radius` ist der Radius der Kugel.
    pub fn rotate_by_geographic_deltas(
        point_on_sphere: Vec3,
        delta_latitude_rad: f32,
        delta_longitude_rad: f32,
        sphere_radius: f32,
    ) -> Vec3 {
        // Aktuelle geografische Koordinaten des Punktes (ohne Altitude)
        let _geo_current = GeographicCoordinates::from_cartesian(point_on_sphere, sphere_radius);

        // Neue geografische Koordinaten
        // Wichtig: Longitude-Rotation hängt von der aktuellen Latitude ab.
        // Einfacher Ansatz: Zwei separate Rotationen.
        // 1. Rotiere um die lokale "Ost"-Achse für die Latitude-Änderung.
        // Die Ost-Achse ist senkrecht zur Meridianebene des Punktes.
        // Meridianebene wird von Z-Achse und dem Punkt aufgespannt. Normale dazu ist Z.cross(P).
        let east_axis = Vec3::Z
            .cross(point_on_sphere.normalize_or_zero())
            .normalize_or_zero();
        let lat_rotation = BevyQuat::from_axis_angle(east_axis, delta_latitude_rad);
        let point_after_lat_rot = lat_rotation * point_on_sphere;

        // 2. Rotiere um die globale Z-Achse (Polachse) für die Longitude-Änderung.
        let lon_rotation = BevyQuat::from_axis_angle(Vec3::Z, delta_longitude_rad);
        let final_point = lon_rotation * point_after_lat_rot;

        CoordinateConverter::project_to_sphere(final_point, sphere_radius)
    }

    /// Generiert eine Anzahl von zufälligen, gleichverteilten Rotationsquaternions.
    /// Verwendet die Methode von Shoemake (basierend auf Arvo) für gleichverteilte Rotationen.
    pub fn generate_uniform_random_rotations(
        count: usize,
        rng: &mut impl rand::Rng, // Nimmt einen generischen RNG
    ) -> Vec<BevyQuat> {
        (0..count)
            .map(|_| {
                // Methode von Ken Shoemake, "Uniform Random Rotations", Graphics Gems III
                let u1: f32 = rng.random(); // [0,1)
                let u2: f32 = rng.random(); // [0,1)
                let u3: f32 = rng.random(); // [0,1)

                let q1_w = (1.0 - u1).sqrt() * (constants::TAU * u2).sin();
                let q1_x = (1.0 - u1).sqrt() * (constants::TAU * u2).cos();
                let q2_y = u1.sqrt() * (constants::TAU * u3).sin();
                let q2_z = u1.sqrt() * (constants::TAU * u3).cos();

                // BevyQuat erwartet (x, y, z, w)
                BevyQuat::from_xyzw(q1_x, q2_y, q2_z, q1_w).normalize() // Normalize für numerische Stabilität
            })
            .collect()
    }

    /// Berechnet die kürzeste Rotation, die `from_orientation` in `to_orientation` überführt.
    /// Das Ergebnis ist `to_orientation * from_orientation.inverse()`.
    pub fn minimal_rotation_between(
        from_orientation: BevyQuat,
        to_orientation: BevyQuat,
    ) -> BevyQuat {
        (to_orientation * from_orientation.inverse()).normalize()
    }

    /// Interpoliert sphärisch linear zwischen mehreren Quaternions.
    /// `t` ist der globale Interpolationsfaktor [0, 1] über die gesamte Sequenz.
    pub fn slerp_sequence(rotations: &[BevyQuat], t: f32) -> Option<BevyQuat> {
        if rotations.is_empty() {
            return None;
        }
        if rotations.len() == 1 {
            return Some(rotations[0]);
        }

        let t_clamped = t.clamp(0.0, 1.0);

        // Skaliere t auf den Bereich der Indizes [0, rotations.len() - 1)
        let max_idx = (rotations.len() - 1) as f32;
        let scaled_t = t_clamped * max_idx;

        let idx0 = scaled_t.floor() as usize;
        let local_t = scaled_t - idx0 as f32;

        // idx1 sollte nicht über das Array-Ende hinausgehen
        let idx1 = (idx0 + 1).min(rotations.len() - 1);

        if idx0 == idx1 {
            // Am Ende der Sequenz
            return Some(rotations[idx0]);
        }

        Some(rotations[idx0].slerp(rotations[idx1], local_t))
    }
}
