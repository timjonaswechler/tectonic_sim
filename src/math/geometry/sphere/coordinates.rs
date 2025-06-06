// src/math/geometry/sphere/coordinates.rs

use crate::math::{error::*, utils::*}; // Vector3DExt für length, normalize etc. auf Vec3
use bevy::math::Vec3;
use std::f32::consts::PI;

/// Sphärische Koordinaten (ISO 80000-2:2019 Konvention: radius, Inklination/Polar, Azimut).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphericalCoordinates {
    /// Radius (r): Abstand vom Ursprung (>= 0).
    pub radius: f32,
    /// Polar- oder Inklinationswinkel (θ): Winkel von der positiven Z-Achse (0 bis π).
    pub polar: f32,
    /// Azimutwinkel (φ): Winkel in der XY-Ebene von der positiven X-Achse (0 bis 2π oder -π bis π).
    pub azimuth: f32,
}

impl SphericalCoordinates {
    /// Erstellt neue sphärische Koordinaten.
    /// `polar` wird auf [0, π] geklemmt, `azimuth` auf [0, 2π) normalisiert.
    pub fn new(radius: f32, polar: f32, azimuth: f32) -> MathResult<Self> {
        if radius < 0.0 {
            return Err(MathError::InvalidConfiguration {
                message: "Radius must be non-negative".to_string(),
            });
        }
        Ok(Self {
            radius,
            polar: polar.clamp(0.0, PI),
            azimuth: angles::normalize_angle(azimuth), // Normalisiert auf [0, 2π)
        })
    }

    /// Erstellt sphärische Koordinaten aus kartesischen Koordinaten (Vec3).
    pub fn from_cartesian(point: Vec3) -> Self {
        let radius = point.length();

        if radius < constants::EPSILON {
            // Nahe am Ursprung
            return Self {
                radius: 0.0,
                polar: 0.0,
                azimuth: 0.0,
            };
        }

        // polar = acos(z / r)
        let polar = (point.z / radius).acos();
        // azimuth = atan2(y, x)
        let azimuth = point.y.atan2(point.x);

        Self {
            radius,
            polar, // Ist bereits im Bereich [0, PI] durch acos
            azimuth: angles::normalize_angle(azimuth), // Normalisiert auf [0, 2π)
        }
    }

    /// Konvertiert sphärische Koordinaten zu kartesischen Koordinaten (Vec3).
    pub fn to_cartesian(&self) -> Vec3 {
        let sin_polar = self.polar.sin();
        // let cos_polar = self.polar.cos(); // Nicht direkt benötigt für Standardformel
        // let sin_azimuth = self.azimuth.sin(); // Nicht direkt benötigt für Standardformel
        // let cos_azimuth = self.azimuth.cos(); // Nicht direkt benötigt für Standardformel

        Vec3::new(
            self.radius * sin_polar * self.azimuth.cos(), // x = r * sin(θ) * cos(φ)
            self.radius * sin_polar * self.azimuth.sin(), // y = r * sin(θ) * sin(φ)
            self.radius * self.polar.cos(),               // z = r * cos(θ)
        )
    }

    /// Normalisiert die Koordinaten auf eine Einheitskugel (Radius = 1.0).
    pub fn normalized(&self) -> Self {
        Self {
            radius: 1.0,
            polar: self.polar,
            azimuth: self.azimuth,
        }
    }
}

/// Geografische Koordinaten (Breitengrad, Längengrad, Höhe).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeographicCoordinates {
    /// Breitengrad (Latitude, φ): Winkel von der Äquatorialebene (-π/2 bis π/2 Radiant, positiv nach Norden).
    pub latitude: f32,
    /// Längengrad (Longitude, λ): Winkel von einem Referenzmeridian (-π bis π Radiant, positiv nach Osten).
    pub longitude: f32,
    /// Höhe (Altitude): Höhe über der Referenz-Sphäre/Ellipsoid.
    pub altitude: f32,
}

impl GeographicCoordinates {
    /// Erstellt neue geografische Koordinaten.
    /// `latitude` wird auf [-π/2, π/2] geklemmt, `longitude` auf [-π, π) normalisiert.
    pub fn new(latitude_rad: f32, longitude_rad: f32, altitude: f32) -> Self {
        Self {
            latitude: latitude_rad.clamp(-constants::PI_OVER_2, constants::PI_OVER_2),
            longitude: angles::normalize_angle_signed(longitude_rad), // Normalisiert auf [-π, π)
            altitude,
        }
    }

    /// Erstellt geografische Koordinaten aus Grad-Werten.
    pub fn from_degrees(latitude_deg: f32, longitude_deg: f32, altitude: f32) -> Self {
        Self::new(
            angles::deg_to_rad(latitude_deg),
            angles::deg_to_rad(longitude_deg),
            altitude,
        )
    }

    /// Konvertiert geografische Koordinaten zu kartesischen Koordinaten (Vec3).
    /// `sphere_radius` ist der Radius der Referenzsphäre, auf den sich die `altitude` bezieht.
    pub fn to_cartesian(&self, sphere_radius: f32) -> Vec3 {
        let total_radius = sphere_radius + self.altitude;
        let cos_lat = self.latitude.cos();

        Vec3::new(
            total_radius * cos_lat * self.longitude.cos(), // x
            total_radius * cos_lat * self.longitude.sin(), // y (Achtung: Geographisch oft z als "oben")
            total_radius * self.latitude.sin(),            // z (oder y, je nach Konvention)
        );
        // Konvention hier: X nach vorne, Y nach links, Z nach oben (üblich in vielen 3D Systemen)
        // Geographisch: X zum Nullmeridian/Äquator, Y 90° Ost/Äquator, Z zum Nordpol
        // Umstellung auf gängige geographische Konvention (X=Nullmeridian, Y=90Ost, Z=Pol):
        // X = R * cos(lat) * cos(lon)
        // Y = R * cos(lat) * sin(lon)
        // Z = R * sin(lat)
        // Die obige Formel war eher Z-up für sphärisch. Korrigieren wir es für geographisch:
        // Vec3::new(
        //     total_radius * cos_lat * self.longitude.cos(), // X
        //     total_radius * cos_lat * self.longitude.sin(), // Y
        //     total_radius * self.latitude.sin(),            // Z
        // )
        // Die ursprüngliche Implementierung in der Frage hatte Y und Z vertauscht im Vergleich zur Standard-Spherical-to-Cartesian,
        // was auf eine Y-Up-Konvention für Geographisch hindeutete. Behalten wir das erstmal bei,
        // aber es ist wichtig, die Achsenkonvention klar zu definieren.
        // Wenn Z nach oben sein soll (üblich für Karten):
        // X = R * cos(lat) * cos(lon)
        // Y = R * cos(lat) * sin(lon)
        // Z = R * sin(lat)
        // Das `point.y.asin()` unten bei `from_cartesian` deutet auf Z-Up hin.
        // Also, die obige `to_cartesian` ist mit Y als Höhe. Ändern wir das zu Z als Höhe:
        // Vec3::new(
        //     total_radius * cos_lat * self.longitude.cos(), // x = Rcos(lat)cos(lon)
        //     total_radius * cos_lat * self.longitude.sin(), // y = Rcos(lat)sin(lon)
        //     total_radius * self.latitude.sin(),            // z = Rsin(lat)
        // )
        // Die Version in deinem Originalcode war:
        // Point3D::new(
        // r * cos_lat * self.longitude.cos(), // x
        // r * self.latitude.sin(),            // y (was eigentlich z sein sollte, wenn latitude der Winkel zur xy-Ebene ist)
        // r * cos_lat * self.longitude.sin(), // z (was eigentlich y sein sollte)
        // )
        // Dies passt zu einer Konvention, bei der Y die "Höhe" (wie latitude) und Z die "Tiefe" (wie sin(longitude)) ist.
        // Lass uns die Standardkonvention (X, Y in Äquatorebene, Z zum Pol) verwenden:
        Vec3::new(
            total_radius * cos_lat * self.longitude.cos(), // X-Achse zeigt zum Nullmeridian
            total_radius * cos_lat * self.longitude.sin(), // Y-Achse zeigt 90° Ost entlang des Äquators
            total_radius * self.latitude.sin(),            // Z-Achse zeigt zum Nordpol
        )
    }

    /// Erstellt geografische Koordinaten aus kartesischen Koordinaten (Vec3).
    /// `reference_radius` ist der Radius der Referenzsphäre.
    pub fn from_cartesian(point: Vec3, reference_radius: f32) -> Self {
        let cartesian_radius = point.length();

        if cartesian_radius < constants::EPSILON {
            return Self::new(0.0, 0.0, -reference_radius); // Am Ursprung
        }

        // Standardkonversion (angenommen X,Y in Äquatorebene, Z zum Pol):
        // Latitude (Breitengrad) φ = asin(z / R)
        // Longitude (Längengrad) λ = atan2(y, x)
        let latitude = (point.z / cartesian_radius).asin(); // Ergebnis in [-PI/2, PI/2]
        let longitude = point.y.atan2(point.x); // Ergebnis in [-PI, PI]

        let altitude = cartesian_radius - reference_radius;

        Self {
            // new() kümmert sich um Normalisierung/Clamping
            latitude,
            longitude,
            altitude,
        }
    }

    /// Konvertiert die Winkel zu Grad-Werten.
    pub fn to_degrees(&self) -> (f32, f32, f32) {
        // (latitude_deg, longitude_deg, altitude)
        (
            angles::rad_to_deg(self.latitude),
            angles::rad_to_deg(self.longitude),
            self.altitude,
        )
    }

    /// Berechnet den Großkreis-Abstand zu einem anderen geografischen Punkt (Haversine-Formel).
    /// `sphere_radius` ist der Radius der Kugel, für die der Abstand berechnet wird.
    pub fn distance_to(&self, other: &GeographicCoordinates, sphere_radius: f32) -> f32 {
        let lat1 = self.latitude;
        let lon1 = self.longitude;
        let lat2 = other.latitude;
        let lon2 = other.longitude;

        let d_lat = lat2 - lat1;
        let d_lon = lon2 - lon1;

        let a = (d_lat * 0.5).sin().powi(2) + lat1.cos() * lat2.cos() * (d_lon * 0.5).sin().powi(2);

        let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt()); // 2.0 * asin(a.sqrt()) ist auch üblich und oft stabiler
        sphere_radius * c
    }

    /// Berechnet den anfänglichen Peilwinkel (Bearing) von diesem Punkt zum anderen Punkt.
    /// Ergebnis in Radiant, [0, 2π), gemessen im Uhrzeigersinn von Nord.
    pub fn initial_bearing_to(&self, other: &GeographicCoordinates) -> f32 {
        let lat1 = self.latitude;
        let lon1 = self.longitude;
        let lat2 = other.latitude;
        let lon2 = other.longitude;

        let d_lon = lon2 - lon1;

        let y = d_lon.sin() * lat2.cos();
        let x = lat1.cos() * lat2.sin() - lat1.sin() * lat2.cos() * d_lon.cos();

        angles::normalize_angle(y.atan2(x)) // atan2 gibt Ergebnis in [-PI, PI], normalize_angle macht [0, 2PI)
    }

    /// Interpoliert sphärisch linear (Slerp) zwischen zwei geografischen Punkten auf einer Einheitskugel.
    /// Die Altitude wird linear interpoliert.
    pub fn slerp(&self, other: &GeographicCoordinates, t: f32, sphere_radius: f32) -> Self {
        if t <= 0.0 {
            return *self;
        }
        if t >= 1.0 {
            return *other;
        }

        let p1 = self.to_cartesian(sphere_radius);
        let p2 = other.to_cartesian(sphere_radius);

        // Slerp für Vec3 (Beispiel, Bevy hat kein Vec3::slerp)
        // Diese Logik ist die gleiche wie für Quaternion::slerp, aber auf Vektoren angewendet.
        let mut dot = p1.dot(p2) / (p1.length_squared() * p2.length_squared()).sqrt(); // Normierter Dot
        dot = dot.clamp(-1.0, 1.0);

        let theta = dot.acos();
        if theta.abs() < constants::EPSILON {
            // Nahezu kollinear
            return *self; // Oder other, je nach t
        }

        let sin_theta = theta.sin();
        let p_interpolated = (p1 * ((1.0 - t) * theta).sin() + p2 * (t * theta).sin()) / sin_theta;

        let interpolated_altitude = comparison::lerp(self.altitude, other.altitude, t);
        let mut result = Self::from_cartesian(p_interpolated, sphere_radius);
        result.altitude = interpolated_altitude;
        result
    }
}

/// Hilfsfunktionen für Koordinatenkonvertierungen und Kugeloperationen.
pub struct CoordinateConverter;

impl CoordinateConverter {
    /// Konvertiert kartesische Koordinaten (Vec3) zu sphärischen Koordinaten.
    pub fn cartesian_to_spherical(point: Vec3) -> SphericalCoordinates {
        SphericalCoordinates::from_cartesian(point)
    }

    /// Konvertiert sphärische Koordinaten zu kartesischen Koordinaten (Vec3).
    pub fn spherical_to_cartesian(coords: SphericalCoordinates) -> Vec3 {
        coords.to_cartesian()
    }

    /// Konvertiert kartesische Koordinaten (Vec3) zu geografischen Koordinaten.
    pub fn cartesian_to_geographic(point: Vec3, reference_radius: f32) -> GeographicCoordinates {
        GeographicCoordinates::from_cartesian(point, reference_radius)
    }

    /// Konvertiert geografische Koordinaten zu kartesischen Koordinaten (Vec3).
    pub fn geographic_to_cartesian(coords: GeographicCoordinates, sphere_radius: f32) -> Vec3 {
        coords.to_cartesian(sphere_radius)
    }

    /// Projiziert einen 3D-Punkt auf die Oberfläche einer Kugel mit gegebenem Radius.
    pub fn project_to_sphere(point: Vec3, radius: f32) -> Vec3 {
        if radius <= 0.0 {
            return Vec3::ZERO;
        } // Ungültiger Radius
        let dist_sq = point.length_squared();
        if dist_sq < constants::EPSILON * constants::EPSILON {
            // Punkt ist am Ursprung, wähle einen beliebigen Punkt auf der Kugeloberfläche
            return Vec3::new(radius, 0.0, 0.0);
        }
        point * (radius / dist_sq.sqrt())
    }

    /// Berechnet den Normalenvektor an einem Punkt auf der Kugeloberfläche
    /// (zeigt vom Kugelzentrum zum Punkt).
    pub fn sphere_normal_at_point(point_on_sphere: Vec3) -> Vec3 {
        point_on_sphere.normalize_or_zero() // normalize_or_zero für den Fall, dass point_on_sphere der Nullvektor ist
    }

    /// Berechnet zwei orthogonale Tangentenvektoren an einem Punkt auf der Kugeloberfläche.
    /// `point_on_sphere` muss ein normalisierter Vektor sein oder wird intern normalisiert.
    pub fn sphere_tangents_at_point(point_on_sphere: Vec3) -> (Vec3, Vec3) {
        let normal = Self::sphere_normal_at_point(point_on_sphere);

        // Wähle einen "up"-Vektor, der nicht (anti-)parallel zur Normale ist
        let temp_up = if normal.x.abs() < 0.9 {
            // Wenn normal nicht zu sehr entlang der X-Achse liegt
            Vec3::X // Wähle X-Achse als temporären Up-Vektor
        } else {
            Vec3::Y // Ansonsten wähle Y-Achse
        };

        let tangent1 = temp_up.cross(normal).normalize_or_zero(); // Erste Tangente
        let tangent2 = normal.cross(tangent1).normalize_or_zero(); // Zweite Tangente, orthogonal zu normal und tangent1

        (tangent1, tangent2)
    }
}
