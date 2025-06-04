// src/math/sphere/coordinates.rs

use crate::math::{error::*, types::*, utils::*};
use std::f32::consts::{PI, TAU};

/// Sphärische Koordinaten
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphericalCoordinates {
    /// Radius (Abstand vom Ursprung)
    pub radius: f32,
    /// Polar-Winkel (von der positiven Z-Achse, 0 bis π)
    pub polar: f32,
    /// Azimut-Winkel (um die Z-Achse, 0 bis 2π)
    pub azimuth: f32,
}

impl SphericalCoordinates {
    /// Erstellt neue sphärische Koordinaten
    pub fn new(radius: f32, polar: f32, azimuth: f32) -> MathResult<Self> {
        if radius < 0.0 {
            return Err(MathError::InvalidConfiguration {
                message: "Radius must be non-negative".to_string(),
            });
        }

        Ok(Self {
            radius,
            polar: polar.clamp(0.0, PI),
            azimuth: angles::normalize_angle(azimuth),
        })
    }

    /// Erstellt sphärische Koordinaten aus kartesischen
    pub fn from_cartesian(point: Point3D) -> Self {
        let radius = point.length();

        if radius < constants::EPSILON {
            return Self {
                radius: 0.0,
                polar: 0.0,
                azimuth: 0.0,
            };
        }

        let polar = (point.z / radius).acos();
        let azimuth = point.y.atan2(point.x);

        Self {
            radius,
            polar,
            azimuth: angles::normalize_angle(azimuth),
        }
    }

    /// Konvertiert zu kartesischen Koordinaten
    pub fn to_cartesian(&self) -> Point3D {
        let sin_polar = self.polar.sin();
        let cos_polar = self.polar.cos();
        let sin_azimuth = self.azimuth.sin();
        let cos_azimuth = self.azimuth.cos();

        Point3D::new(
            self.radius * sin_polar * cos_azimuth,
            self.radius * sin_polar * sin_azimuth,
            self.radius * cos_polar,
        )
    }

    /// Normalisiert auf Einheitskugel
    pub fn normalized(&self) -> Self {
        Self {
            radius: 1.0,
            polar: self.polar,
            azimuth: self.azimuth,
        }
    }
}

/// Geografische Koordinaten (Lat/Lon)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeographicCoordinates {
    /// Breitengrad in Radiant (-π/2 bis π/2)
    pub latitude: f32,
    /// Längengrad in Radiant (-π bis π)
    pub longitude: f32,
    /// Höhe über der Referenz-Sphäre
    pub altitude: f32,
}

impl GeographicCoordinates {
    /// Erstellt neue geografische Koordinaten
    pub fn new(latitude: f32, longitude: f32, altitude: f32) -> Self {
        Self {
            latitude: latitude.clamp(-PI * 0.5, PI * 0.5),
            longitude: angles::normalize_angle_signed(longitude),
            altitude,
        }
    }

    /// Erstellt aus Grad-Werten
    pub fn from_degrees(lat_deg: f32, lon_deg: f32, altitude: f32) -> Self {
        Self::new(
            angles::deg_to_rad(lat_deg),
            angles::deg_to_rad(lon_deg),
            altitude,
        )
    }

    /// Konvertiert zu kartesischen Koordinaten
    pub fn to_cartesian(&self, radius: f32) -> Point3D {
        let r = radius + self.altitude;
        let cos_lat = self.latitude.cos();

        Point3D::new(
            r * cos_lat * self.longitude.cos(),
            r * self.latitude.sin(),
            r * cos_lat * self.longitude.sin(),
        )
    }

    /// Erstellt aus kartesischen Koordinaten
    pub fn from_cartesian(point: Point3D, reference_radius: f32) -> Self {
        let radius = point.length();

        if radius < constants::EPSILON {
            return Self::new(0.0, 0.0, -reference_radius);
        }

        let normalized = point / radius;
        let latitude = normalized.y.asin();
        let longitude = normalized.z.atan2(normalized.x);
        let altitude = radius - reference_radius;

        Self::new(latitude, longitude, altitude)
    }

    /// Konvertiert zu Grad-Werten
    pub fn to_degrees(&self) -> (f32, f32, f32) {
        (
            angles::rad_to_deg(self.latitude),
            angles::rad_to_deg(self.longitude),
            self.altitude,
        )
    }

    /// Berechnet den Großkreis-Abstand zu einem anderen Punkt (Haversine-Formel)
    pub fn distance_to(&self, other: &GeographicCoordinates, radius: f32) -> f32 {
        let dlat = other.latitude - self.latitude;
        let dlon = other.longitude - self.longitude;

        let a = (dlat * 0.5).sin().powi(2)
            + self.latitude.cos() * other.latitude.cos() * (dlon * 0.5).sin().powi(2);

        let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());
        radius * c
    }

    /// Berechnet den Bärungs-Winkel zu einem anderen Punkt
    pub fn bearing_to(&self, other: &GeographicCoordinates) -> f32 {
        let dlon = other.longitude - self.longitude;
        let y = dlon.sin() * other.latitude.cos();
        let x = self.latitude.cos() * other.latitude.sin()
            - self.latitude.sin() * other.latitude.cos() * dlon.cos();

        angles::normalize_angle(y.atan2(x))
    }

    /// Interpoliert zwischen zwei geografischen Punkten
    pub fn interpolate_to(&self, other: &GeographicCoordinates, t: f32) -> Self {
        let d = self.distance_to(other, 1.0); // Einheitskugel

        if d < constants::EPSILON {
            return *self;
        }

        let a = (1.0 - t) * d / d.sin();
        let b = t * d / d.sin();

        let x = a * self.latitude.cos() * self.longitude.cos()
            + b * other.latitude.cos() * other.longitude.cos();
        let y = a * self.latitude.cos() * self.longitude.sin()
            + b * other.latitude.cos() * other.longitude.sin();
        let z = a * self.latitude.sin() + b * other.latitude.sin();

        let lat = z.atan2((x * x + y * y).sqrt());
        let lon = y.atan2(x);
        let alt = comparison::lerp(self.altitude, other.altitude, t);

        Self::new(lat, lon, alt)
    }
}

/// Koordinaten-Konvertierungs-Utilities
pub struct CoordinateConverter;

impl CoordinateConverter {
    /// Konvertiert kartesische zu sphärischen Koordinaten
    pub fn cartesian_to_spherical(point: Point3D) -> SphericalCoordinates {
        SphericalCoordinates::from_cartesian(point)
    }

    /// Konvertiert sphärische zu kartesischen Koordinaten
    pub fn spherical_to_cartesian(coords: SphericalCoordinates) -> Point3D {
        coords.to_cartesian()
    }

    /// Konvertiert kartesische zu geografischen Koordinaten
    pub fn cartesian_to_geographic(point: Point3D, radius: f32) -> GeographicCoordinates {
        GeographicCoordinates::from_cartesian(point, radius)
    }

    /// Konvertiert geografische zu kartesischen Koordinaten
    pub fn geographic_to_cartesian(coords: GeographicCoordinates, radius: f32) -> Point3D {
        coords.to_cartesian(radius)
    }

    /// Projiziert einen Punkt auf eine Kugel
    pub fn project_to_sphere(point: Point3D, radius: f32) -> Point3D {
        let distance = point.length();
        if distance < constants::EPSILON {
            return Point3D::new(radius, 0.0, 0.0);
        }
        point * (radius / distance)
    }

    /// Berechnet die Normale an einem Punkt auf der Kugel
    pub fn sphere_normal(point: Point3D) -> Point3D {
        point.normalize()
    }

    /// Berechnet Tangenten-Vektoren an einem Punkt auf der Kugel
    pub fn sphere_tangents(point: Point3D) -> (Point3D, Point3D) {
        let normal = Self::sphere_normal(point);

        // Wähle einen "up"-Vektor, der nicht parallel zur Normale ist
        let up = if normal.y.abs() < 0.9 {
            Point3D::new(0.0, 1.0, 0.0)
        } else {
            Point3D::new(1.0, 0.0, 0.0)
        };

        let tangent1 = normal.cross(up).normalize();
        let tangent2 = normal.cross(tangent1);

        (tangent1, tangent2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spherical_cartesian_conversion() {
        let point = Point3D::new(1.0, 0.0, 0.0);
        let spherical = SphericalCoordinates::from_cartesian(point);
        let back_to_cartesian = spherical.to_cartesian();

        assert!(comparison::nearly_equal(back_to_cartesian.x, 1.0));
        assert!(comparison::nearly_equal(back_to_cartesian.y, 0.0));
        assert!(comparison::nearly_equal(back_to_cartesian.z, 0.0));
    }

    #[test]
    fn test_geographic_coordinates() {
        let geographic = GeographicCoordinates::from_degrees(45.0, 90.0, 0.0);
        let cartesian = geographic.to_cartesian(1.0);
        let back_to_geographic = GeographicCoordinates::from_cartesian(cartesian, 1.0);

        let (lat_deg, lon_deg, _) = back_to_geographic.to_degrees();
        assert!(comparison::nearly_equal(lat_deg, 45.0));
        assert!(comparison::nearly_equal(lon_deg, 90.0));
    }

    #[test]
    fn test_great_circle_distance() {
        let coord1 = GeographicCoordinates::from_degrees(0.0, 0.0, 0.0);
        let coord2 = GeographicCoordinates::from_degrees(0.0, 90.0, 0.0);

        let distance = coord1.distance_to(&coord2, 1.0);
        assert!(comparison::nearly_equal(distance, PI * 0.5));
    }

    #[test]
    fn test_coordinate_converter() {
        let point = Point3D::new(2.0, 3.0, 4.0);
        let projected = CoordinateConverter::project_to_sphere(point, 5.0);

        assert!(comparison::nearly_equal(projected.length(), 5.0));

        let normal = CoordinateConverter::sphere_normal(projected);
        assert!(comparison::nearly_equal(normal.length(), 1.0));
    }
}
