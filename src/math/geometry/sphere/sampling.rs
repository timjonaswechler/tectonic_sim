// src/math/geometry/sphere/sampling.rs

use crate::math::{
    MathError,
    MathResult,
    geometry::sphere::coordinates::{CoordinateConverter, SphericalCoordinates},

    prelude::SeedResource, // << NEU
    utils::*,
};
use crate::next_random_range;
use bevy::math::{Quat as BevyQuat, Vec2, Vec3}; // << NEU
// use rand::Rng; // Standard Rng Trait << ENTFERNEN
use std::f32::consts::{PI, TAU};

/// Verschiedene Sampling-Methoden für Punkte auf einer Kugeloberfläche.
#[derive(Debug, Clone, Copy)]
pub enum SphereSamplingMethod {
    // Umbenannt von SamplingMethod zur Klarheit
    /// Uniform verteilte Punkte (statistisch).
    Uniform,
    /// Punkte verteilt entlang einer Fibonacci-Spirale (sehr gleichmäßig).
    Fibonacci,
    /// Stratified Sampling (Grid-basiert mit Jitter).
    Stratified,
    /// Poisson Disk Sampling (garantiert minimalen Abstand zwischen Punkten).
    PoissonDisk { min_angular_distance_rad: f32 }, // Abstand als Winkel
}

/// Sampler zum Erzeugen von Punkten auf einer Kugeloberfläche.
pub struct SphereSampler {
    method: SphereSamplingMethod,
    radius: f32,
}

impl SphereSampler {
    /// Erstellt einen neuen Kugel-Sampler.
    pub fn new(method: SphereSamplingMethod, radius: f32) -> MathResult<Self> {
        if radius <= 0.0 {
            return Err(MathError::InvalidConfiguration {
                message: "Sphere radius must be positive for sampling.".to_string(),
            });
        }
        Ok(Self { method, radius })
    }

    /// Generiert eine bestimmte Anzahl von Punkten auf der Kugeloberfläche.
    /// `rng` ist die Quelle für Zufallszahlen.
    pub fn sample_points(
        &self,
        count: usize,
        seed_resource: &mut SeedResource, // << Geändert
    ) -> Vec<Vec3> {
        if count == 0 {
            return Vec::new();
        }

        match self.method {
            SphereSamplingMethod::Uniform => self.uniform_sampling(count, seed_resource), // << Geändert
            SphereSamplingMethod::Fibonacci => self.fibonacci_sampling(count),
            SphereSamplingMethod::Stratified => self.stratified_sampling(count, seed_resource), // << Geändert
            SphereSamplingMethod::PoissonDisk {
                min_angular_distance_rad,
            } => {
                let chord_length_unit_sphere = 2.0 * (min_angular_distance_rad / 2.0).sin();
                let min_cartesian_distance = chord_length_unit_sphere * self.radius;
                self.poisson_disk_sampling(count, min_cartesian_distance, seed_resource) // << Geändert
            }
        }
    }

    /// Generiert einen einzelnen zufälligen Punkt auf der Kugeloberfläche.
    pub fn sample_single_point(&self, seed_resource: &mut SeedResource) -> Vec3 {
        // << Geändert
        match self.method {
            SphereSamplingMethod::Uniform | SphereSamplingMethod::Stratified => {
                self.uniform_point_marsaglia(seed_resource) // << Geändert
            }
            SphereSamplingMethod::Fibonacci => self
                .fibonacci_sampling(1)
                .pop()
                .unwrap_or_else(|| self.uniform_point_marsaglia(seed_resource)), // << Geändert
            SphereSamplingMethod::PoissonDisk { .. } => {
                self.uniform_point_marsaglia(seed_resource) // << Geändert
            }
        }
    }

    /// Generiert uniform verteilte Punkte mittels Marsaglia-Methode.
    fn uniform_sampling(&self, count: usize, seed_resource: &mut SeedResource) -> Vec<Vec3> {
        // << Geändert
        (0..count)
            .map(|_| self.uniform_point_marsaglia(seed_resource)) // << Geändert
            .collect()
    }

    /// Generiert einen einzelnen uniform verteilten Punkt mittels Marsaglia-Methode.
    /// Gibt einen Punkt auf der Einheitskugel zurück, skaliert mit `self.radius`.
    fn uniform_point_marsaglia(&self, seed_resource: &mut SeedResource) -> Vec3 {
        // << Geändert
        loop {
            let x1: f32 = next_random_range!(seed_resource, -1.0, 1.0); // << Geändert
            let x2: f32 = next_random_range!(seed_resource, -1.0, 1.0); // << Geändert
            let sum_sq = x1 * x1 + x2 * x2;

            if sum_sq < 1.0 && sum_sq > constants::EPSILON {
                let factor_sqrt = (1.0 - sum_sq).sqrt();
                return Vec3::new(
                    2.0 * x1 * factor_sqrt,
                    2.0 * x2 * factor_sqrt,
                    1.0 - 2.0 * sum_sq,
                ) * self.radius;
            }
        }
    }

    /// Generiert Punkte mittels Fibonacci-Spirale (oder "Golden Spiral").
    /// Liefert eine sehr gleichmäßige Verteilung.
    fn fibonacci_sampling(&self, count: usize) -> Vec<Vec3> {
        if count == 0 {
            return Vec::new();
        }
        let mut points = Vec::with_capacity(count);
        let golden_angle = PI * (3.0 - 5.0_f32.sqrt()); // Näherung an Goldenen Winkel für Spirale

        for i in 0..count {
            let y = 1.0 - (i as f32 / (count - 1).max(1) as f32) * 2.0; // y geht von 1 bis -1
            let radius_at_y = (1.0 - y * y).sqrt(); // Radius der Scheibe bei Höhe y

            let theta = golden_angle * i as f32; // Golden Angle

            let x = radius_at_y * theta.cos();
            let z = radius_at_y * theta.sin(); // Beachte: übliche Spirale ist in XZ, Y ist "oben"

            points.push(Vec3::new(x, y, z) * self.radius);
        }
        points
    }

    /// Generiert Punkte mittels Stratified Sampling.
    /// Teilt die Kugeloberfläche in (annähernd) gleich große Flächen auf und platziert
    /// einen gejitterten Punkt in jeder Zelle.
    fn stratified_sampling(&self, count: usize, seed_resource: &mut SeedResource) -> Vec<Vec3> {
        // << Geändert
        let num_latitude_bands = ((count as f32 * PI).sqrt().round() as usize).max(1);
        let mut points = Vec::with_capacity(count);
        let mut points_generated = 0;

        for i in 0..num_latitude_bands {
            let u_min = i as f32 / num_latitude_bands as f32;
            let u_max = (i + 1) as f32 / num_latitude_bands as f32;
            let points_this_band = (count as f32 / num_latitude_bands as f32).ceil() as usize;

            for _ in 0..points_this_band {
                if points_generated >= count {
                    break;
                }

                let u_rand = next_random_range!(seed_resource, u_min, u_max); // << Geändert
                let v_rand: f32 = seed_resource.next_value(); // << Geändert (Annahme: rand_f32() existiert)
                let theta = TAU * v_rand;
                let phi = (2.0 * u_rand - 1.0).acos();
                let sin_phi = phi.sin();
                points.push(
                    Vec3::new(sin_phi * theta.cos(), sin_phi * theta.sin(), phi.cos())
                        * self.radius,
                );
                points_generated += 1;
            }
            if points_generated >= count {
                break;
            }
        }
        points.truncate(count);
        points
    }

    /// Generiert Punkte mittels Poisson Disk Sampling auf der Kugeloberfläche.
    /// `min_cartesian_distance` ist der minimale euklidische Abstand zwischen Punkten.
    fn poisson_disk_sampling(
        &self,
        target_count: usize,
        min_cartesian_distance: f32,
        seed_resource: &mut SeedResource, // << Geändert
    ) -> Vec<Vec3> {
        if min_cartesian_distance <= 0.0 {
            return self.uniform_sampling(target_count, seed_resource); // << Geändert
        }

        let mut points = Vec::new();
        let mut active_list = Vec::new();

        let first_point = self.uniform_point_marsaglia(seed_resource); // << Geändert
        points.push(first_point);
        active_list.push(0);

        let k_attempts = 30;

        while !active_list.is_empty() && points.len() < target_count {
            let active_idx_in_list =
                next_random_range!(seed_resource, 0, active_list.len() - 1) as usize;
            let current_point_idx = active_list[active_idx_in_list];
            let current_point = points[current_point_idx];

            let mut found_candidate_for_current = false;
            for _ in 0..k_attempts {
                let random_angle_dist_rad = next_random_range!(
                    // << Geändert
                    seed_resource,
                    (min_cartesian_distance / self.radius),
                    (2.0 * min_cartesian_distance / self.radius)
                );
                let random_azimuth_rad: f32 = next_random_range!(seed_resource, 0.0, TAU); // << Geändert

                let random_axis = self
                    .uniform_point_marsaglia(seed_resource)
                    .normalize_or_zero(); // << Geändert
                let _candidate_rotation =
                    BevyQuat::from_axis_angle(random_axis, random_angle_dist_rad);

                let (tangent_u, tangent_v) =
                    CoordinateConverter::sphere_tangents_at_point(current_point);
                let local_offset = Vec2::from_angle(random_azimuth_rad).normalize()
                    * min_cartesian_distance
                    * next_random_range!(seed_resource, 1.0, 2.0); // << Geändert
                let candidate_on_tangent_plane =
                    current_point + tangent_u * local_offset.x + tangent_v * local_offset.y;
                let candidate =
                    CoordinateConverter::project_to_sphere(candidate_on_tangent_plane, self.radius);

                if self.is_valid_poisson_candidate(&candidate, &points, min_cartesian_distance) {
                    points.push(candidate);
                    active_list.push(points.len() - 1);
                    found_candidate_for_current = true;
                    if points.len() >= target_count {
                        break;
                    }
                }
            }

            if !found_candidate_for_current {
                active_list.swap_remove(active_idx_in_list as usize);
            }
        }
        points
    }

    fn is_valid_poisson_candidate(
        &self,
        candidate: &Vec3,
        existing_points: &[Vec3],
        min_cartesian_distance: f32,
    ) -> bool {
        let min_dist_sq = min_cartesian_distance * min_cartesian_distance;
        for point in existing_points {
            if candidate.distance_squared(*point) < min_dist_sq {
                return false;
            }
        }
        true
    }
}

/// Erzeugt spezielle Punktmuster auf einer Kugeloberfläche.
pub struct SamplingPatterns;

impl SamplingPatterns {
    /// Erzeugt die 12 Eckpunkte eines Ikosaeders, normalisiert auf den gegebenen Radius.
    pub fn icosahedron_vertices(radius: f32) -> Vec<Vec3> {
        let phi = constants::GOLDEN_RATIO; // (1 + sqrt(5)) / 2
        let _norm_factor = (1.0 + phi * phi).sqrt(); // Für Normalisierung auf Einheitskugel
        // Punkte auf Einheitskugel
        let points_unit = vec![
            Vec3::new(-1.0, phi, 0.0),
            Vec3::new(1.0, phi, 0.0),
            Vec3::new(-1.0, -phi, 0.0),
            Vec3::new(1.0, -phi, 0.0),
            Vec3::new(0.0, -1.0, phi),
            Vec3::new(0.0, 1.0, phi),
            Vec3::new(0.0, -1.0, -phi),
            Vec3::new(0.0, 1.0, -phi),
            Vec3::new(phi, 0.0, -1.0),
            Vec3::new(phi, 0.0, 1.0),
            Vec3::new(-phi, 0.0, -1.0),
            Vec3::new(-phi, 0.0, 1.0),
        ];
        points_unit
            .into_iter()
            .map(|p| p.normalize_or_zero() * radius)
            .collect()
    }

    // Geodesic sphere ist komplexer, da es Triangulation und Mittelpunktberechnung erfordert.
    // Hier eine sehr vereinfachte Version, die nur existierende Punkte mittelt.
    // Eine korrekte Implementierung würde die Dreiecke des Ikosaeders unterteilen.
    // pub fn geodesic_sphere_simplified(radius: f32, subdivisions: usize) -> Vec<Vec3> { ... }

    /// Erzeugt Punkte entlang von Breitengradkreisen.
    pub fn latitude_circles(
        radius: f32,
        num_latitude_circles: usize, // Anzahl der Breitengradkreise (inkl. Pole, wenn count_per_circle > 0)
        points_per_circle: usize,    // Punkte pro Kreis (außer an Polen)
    ) -> Vec<Vec3> {
        if num_latitude_circles == 0 || points_per_circle == 0 {
            return Vec::new();
        }
        let mut points = Vec::new();

        // Inklusive Pole, falls num_latitude_circles >= 1
        if num_latitude_circles == 1 {
            // Nur Äquator (oder Pole, je nach Definition)
            if points_per_circle > 0 {
                points.push(Vec3::new(0.0, 0.0, radius)); // Nordpol
                if points_per_circle > 1 {
                    // Um Duplikat zu vermeiden, wenn nur 1 Punkt
                    points.push(Vec3::new(0.0, 0.0, -radius)); // Südpol
                }
            }
            return points;
        }

        // Schleife von Pol zu Pol
        for i in 0..=num_latitude_circles {
            // Inklusive Start- und End-Breitengrad
            let t_lat = i as f32 / num_latitude_circles as f32; // [0, 1]
            let latitude_rad = PI * t_lat - constants::PI_OVER_2; // Von -PI/2 (Süd) bis PI/2 (Nord)

            if latitude_rad.abs() >= constants::PI_OVER_2 - constants::EPSILON {
                // An den Polen
                points.push(Vec3::new(0.0, 0.0, latitude_rad.signum() * radius));
            } else {
                let cos_lat = latitude_rad.cos();
                for j in 0..points_per_circle {
                    let longitude_rad = TAU * (j as f32 / points_per_circle as f32);
                    points.push(Vec3::new(
                        radius * cos_lat * longitude_rad.cos(),
                        radius * cos_lat * longitude_rad.sin(),
                        radius * latitude_rad.sin(),
                    ));
                }
            }
        }
        // Dedupliziere Pole, falls sie mehrfach generiert wurden
        points.sort_by(|a: &Vec3, b: &Vec3| {
            a.x.partial_cmp(&b.x)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
                .then_with(|| a.z.partial_cmp(&b.z).unwrap_or(std::cmp::Ordering::Equal))
        });
        points.dedup_by(|a, b| a.distance_squared(*b) < constants::EPSILON * constants::EPSILON);
        points
    }

    /// Generiert zufällige Punkte in einem Kugel-Segment (sphärische Kappe).
    /// `cap_axis` ist die Achse der Kappe (zeigt vom Kugelzentrum zur Mitte der Kappe).
    /// `cap_angle_rad` ist der halbe Öffnungswinkel der Kappe in Radiant.
    pub fn spherical_cap_points(
        radius: f32,
        cap_axis: Vec3,
        cap_angle_rad: f32,
        count: usize,
        seed_resource: &mut SeedResource, // << Geändert
    ) -> Vec<Vec3> {
        if count == 0 || cap_angle_rad <= 0.0 {
            return Vec::new();
        }
        let cap_angle_clamped = cap_angle_rad.clamp(0.0, PI);
        let mut points = Vec::with_capacity(count);
        let cos_max_polar_offset = cap_angle_clamped.cos();
        let rotation_to_z = BevyQuat::from_rotation_arc(cap_axis.normalize_or_zero(), Vec3::Z);
        let rotation_from_z = rotation_to_z.inverse();

        for _ in 0..count {
            let cos_polar_offset = next_random_range!(seed_resource, cos_max_polar_offset, 1.0); // << Geändert
            let polar_offset = cos_polar_offset.acos();
            let azimuth = next_random_range!(seed_resource, 0.0, TAU); // << Geändert
            let point_on_z_cap = SphericalCoordinates::new(radius, polar_offset, azimuth)
                .expect("Spherical coord creation failed")
                .to_cartesian();
            points.push(rotation_from_z * point_on_z_cap);
        }
        points
    }
}
