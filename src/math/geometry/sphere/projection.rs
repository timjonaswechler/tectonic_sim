// src/math/geometry/sphere/projection.rs

use crate::math::{
    error::*,
    geometry::sphere::coordinates::{CoordinateConverter, GeographicCoordinates},
    utils::*, // constants, angles, comparison
};
use bevy::math::{Quat as BevyQuat, Vec2, Vec3}; // Bevy's Vec3 und Quat
use std::f32::consts::{PI, TAU};

/// Verschiedene Projektionstypen, um Punkte von einer Kugeloberfläche auf eine 2D-Ebene abzubilden.
#[derive(Debug, Clone, Copy)]
pub enum SphereProjectionType {
    /// Stereographische Projektion von einem Pol aus.
    Stereographic { pole_axis: Vec3 }, // Achse, die zum Projektionspol zeigt (z.B. Vec3::Z für Nordpol)
    /// Orthographische Projektion (Parallelprojektion).
    Orthographic { view_direction: Vec3 }, // Richtung, aus der die Kugel betrachtet wird
    /// Mercator-Projektion.
    Mercator,
    /// Lambert Azimuthal Equal-Area Projektion.
    LambertAzimuthal { center_on_sphere: Vec3 }, // Zentrum der Projektion auf der Kugeloberfläche
    /// Gnomonische Projektion (Zentralprojektion vom Kugelmittelpunkt).
    Gnomonic { center_on_sphere: Vec3 }, // Zentrum der Projektion auf der Kugeloberfläche
    /// Mollweide-Projektion (flächengetreu).
    Mollweide,
    /// Equirectangular-Projektion (auch Plate Carrée genannt).
    Equirectangular,
}

/// Projektions-Engine für Transformationen von Kugeloberflächenpunkten (Vec3) zu 2D-Punkten (Vec2).
pub struct SphereProjector {
    projection_type: SphereProjectionType,
    sphere_radius: f32,
    projection_scale: f32, // Skalierungsfaktor für die Ausgabekoordinaten der Projektion
}

impl SphereProjector {
    /// Erstellt einen neuen Kugel-Projektor.
    /// `sphere_radius` ist der Radius der Kugel, von der projiziert wird.
    pub fn new(projection_type: SphereProjectionType, sphere_radius: f32) -> MathResult<Self> {
        if sphere_radius <= 0.0 {
            return Err(MathError::InvalidConfiguration {
                message: "Sphere radius must be positive for projection.".to_string(),
            });
        }
        Ok(Self {
            projection_type,
            sphere_radius,
            projection_scale: 1.0, // Standardmäßig keine zusätzliche Skalierung
        })
    }

    /// Setzt einen zusätzlichen Skalierungsfaktor für die projizierten 2D-Koordinaten.
    pub fn with_projection_scale(mut self, scale: f32) -> Self {
        self.projection_scale = scale.max(constants::EPSILON); // Skalierung sollte positiv sein
        self
    }

    /// Projiziert einen einzelnen 3D-Punkt von der Kugeloberfläche auf die 2D-Ebene.
    /// Der Eingabepunkt wird intern auf den `sphere_radius` normalisiert.
    pub fn project_point(&self, point_on_sphere_surface: Vec3) -> MathResult<Vec2> {
        // Stelle sicher, dass der Punkt auf der Kugeloberfläche liegt (oder normalisiere ihn dorthin)
        let point_normalized_to_radius = if point_on_sphere_surface.length_squared()
            < constants::EPSILON * constants::EPSILON
        {
            // Nahe am Ursprung, kann nicht sinnvoll projiziert werden
            return Err(MathError::InvalidConfiguration {
                message: "Cannot project zero vector or point very close to origin.".to_string(),
            });
        } else {
            CoordinateConverter::project_to_sphere(point_on_sphere_surface, self.sphere_radius)
        };

        let projected_2d = match self.projection_type {
            SphereProjectionType::Stereographic { pole_axis } => self
                .project_stereographic(point_normalized_to_radius, pole_axis.normalize_or_zero())?,
            SphereProjectionType::Orthographic { view_direction } => self.project_orthographic(
                point_normalized_to_radius,
                view_direction.normalize_or_zero(),
            )?,
            SphereProjectionType::Mercator => self.project_mercator(point_normalized_to_radius)?,
            SphereProjectionType::LambertAzimuthal { center_on_sphere } => self
                .project_lambert_azimuthal(
                    point_normalized_to_radius,
                    center_on_sphere.normalize_or_zero(),
                )?,
            SphereProjectionType::Gnomonic { center_on_sphere } => self.project_gnomonic(
                point_normalized_to_radius,
                center_on_sphere.normalize_or_zero(),
            )?,
            SphereProjectionType::Mollweide => {
                self.project_mollweide(point_normalized_to_radius)?
            }
            SphereProjectionType::Equirectangular => {
                self.project_equirectangular(point_normalized_to_radius)?
            }
        };

        Ok(projected_2d * self.projection_scale)
    }

    /// Projiziert eine Liste von 3D-Punkten auf die 2D-Ebene.
    pub fn project_points(&self, points_on_sphere_surface: &[Vec3]) -> MathResult<Vec<Vec2>> {
        points_on_sphere_surface
            .iter()
            .map(|&point| self.project_point(point))
            .collect()
    }

    /// Inverse Projektion: Konvertiert einen 2D-Punkt zurück zu einem 3D-Punkt auf der Kugeloberfläche.
    /// Nicht für alle Projektionstypen trivial oder eindeutig implementierbar.
    pub fn unproject_point(&self, point_2d: Vec2) -> MathResult<Vec3> {
        let scaled_point_2d = point_2d / self.projection_scale;

        match self.projection_type {
            SphereProjectionType::Stereographic { pole_axis } => {
                self.unproject_stereographic(scaled_point_2d, pole_axis.normalize_or_zero())
            }
            SphereProjectionType::Orthographic { view_direction } => {
                self.unproject_orthographic(scaled_point_2d, view_direction.normalize_or_zero())
            }
            SphereProjectionType::Equirectangular => {
                self.unproject_equirectangular(scaled_point_2d)
            }
            SphereProjectionType::Mercator => self.unproject_mercator(scaled_point_2d),
            SphereProjectionType::LambertAzimuthal { center_on_sphere } => self
                .unproject_lambert_azimuthal(scaled_point_2d, center_on_sphere.normalize_or_zero()),
            // TODO: Implement inverse for other projections if needed
            _ => Err(MathError::InvalidConfiguration {
                message: "Inverse projection not implemented for this projection type".to_string(),
            }),
        }
    }

    // --- Projektions-Implementierungen ---

    fn project_stereographic(&self, point: Vec3, projection_pole: Vec3) -> MathResult<Vec2> {
        // Der Punkt wird vom Gegenpol des `projection_pole` projiziert.
        // Wir rotieren das System so, dass der `projection_pole` zum Nordpol (0,0,R) wird.
        // Dann projizieren wir vom Südpol (0,0,-R).
        let rotation = BevyQuat::from_rotation_arc(
            projection_pole, // projection_pole sollte bereits normalisiert sein oder hier normalisieren
            Vec3::Z,
        );
        let point_rotated = rotation * point;

        // Standard Stereographische Projektion vom Südpol (0,0,-R) auf die Ebene z=0
        // x' = R * x / (R + z)
        // y' = R * y / (R + z)
        // Wenn point_rotated.z sehr nahe an -self.sphere_radius ist (also der Projektionspol selbst), dann Singularität.
        if (point_rotated.z + self.sphere_radius).abs() < constants::EPSILON {
            return Err(MathError::GeometricFailure {
                operation: "Stereographic projection of the projection pole itself".to_string(),
            });
        }
        let k = self.sphere_radius / (self.sphere_radius + point_rotated.z); // Beachte: Vorzeichen von z, wenn von Südpol projiziert
        Ok(Vec2::new(k * point_rotated.x, k * point_rotated.y))
    }

    fn unproject_stereographic(&self, point_2d: Vec2, projection_pole: Vec3) -> MathResult<Vec3> {
        // Inverse der Standard Stereographischen Projektion (vom Südpol auf z=0)
        // x = 2R * x' / (R^2 + x'^2 + y'^2)  <- Fehler in Formel, R ist nicht quadriert im Nenner für Standard
        // x = 2R * x' / (1 + x'^2/R^2 + y'^2/R^2) wenn x',y' die unskalierten sind.
        // Oder:
        // x = 2R * x' / (k^2 + x'^2 + y'^2) wenn k der Skalierungsfaktor war.
        // Hier: x', y' sind bereits die projizierten Koordinaten.

        // Formel für Projektion vom Südpol (0,0,-R) auf Ebene z=0:
        // x' = R * x / (R + z)
        // y' = R * y / (R + z)
        // Inverse:
        // Sei d_sq = x'^2 + y'^2
        // z = R * (R^2 - d_sq) / (R^2 + d_sq)
        // x = x' * (R + z) / R
        // y = y' * (R + z) / R
        // Wichtig: Dies gilt für eine Projektion, bei der der Projektionspol (von dem aus geschaut wird) der Südpol ist
        // und der Punkt, von dem projiziert wird, der Nordpol ist.
        // In unserem `project_stereographic` rotieren wir so, dass `projection_pole` zu `Vec3::Z` (Nordpol) wird
        // und projizieren vom Südpol. Also ist die obige Inverse korrekt für `point_rotated`.

        let r = self.sphere_radius;
        let d_sq = point_2d.length_squared();

        // Vermeide Division durch Null, wenn d_sq sehr groß ist (was nicht passieren sollte, wenn point_2d aus einer gültigen Projektion stammt)
        if (r * r + d_sq).abs() < constants::EPSILON {
            return Err(MathError::GeometricFailure {
                operation: "Stereographic unprojection, denominator near zero".to_string(),
            });
        }

        let z_rotated = r * (r * r - d_sq) / (r * r + d_sq);
        // Wenn z_rotated sehr nahe an -r ist, kann (r + z_rotated) / r problematisch werden.
        // Dies entspricht dem Fall, dass der unprojizierte Punkt der Projektionspol selbst ist.
        if (r + z_rotated).abs() < constants::EPSILON && r.abs() > constants::EPSILON {
            // point_2d ist unendlich, der Punkt ist der Projektionspol (Südpol im rotierten System)
            let point_on_sphere_at_projection_pole = Vec3::new(0.0, 0.0, -r);
            let rotation_to_original_pole = BevyQuat::from_rotation_arc(
                Vec3::Z,         // Im rotierten System war der Nordpol Vec3::Z
                projection_pole, // Ziel ist der ursprüngliche Projektionspol
            );
            return Ok(rotation_to_original_pole * point_on_sphere_at_projection_pole);
        }

        let factor = (r + z_rotated) / r; // Kann problematisch sein, wenn r = 0, aber r > 0 wird im Konstruktor geprüft.
        // Wenn r + z_rotated = 0, dann ist z_rotated = -r.
        // Das bedeutet r * (r^2 - d_sq) / (r^2 + d_sq) = -r
        // r^2 - d_sq = -(r^2 + d_sq)
        // r^2 - d_sq = -r^2 - d_sq
        // 2 * r^2 = 0 => r = 0, was ausgeschlossen ist.
        // Also ist factor nur dann unendlich, wenn r = 0.

        let x_rotated = point_2d.x * factor;
        let y_rotated = point_2d.y * factor;

        let point_on_sphere_rotated_system =
            Vec3::new(x_rotated, y_rotated, z_rotated).normalize() * r;

        // Rotiere zurück in das ursprüngliche Koordinatensystem.
        // Die ursprüngliche Rotation war von `projection_pole` zu `Vec3::Z`.
        // Die inverse Rotation ist von `Vec3::Z` zu `projection_pole`.
        let rotation_to_original_pole = BevyQuat::from_rotation_arc(
            Vec3::Z,
            projection_pole, // projection_pole sollte bereits normalisiert sein
        );
        let point_on_sphere_original_system =
            rotation_to_original_pole * point_on_sphere_rotated_system;

        Ok(point_on_sphere_original_system)
    }

    fn project_orthographic(&self, point: Vec3, view_direction: Vec3) -> MathResult<Vec2> {
        // Rotiere das System so, dass die `view_direction` entlang der Z-Achse liegt.
        // Die Projektion ist dann einfach die X- und Y-Koordinate.
        let rotation = BevyQuat::from_rotation_arc(
            view_direction, // view_direction sollte bereits normalisiert sein
            Vec3::Z,
        );
        let point_rotated = rotation * point;

        // Nur Punkte auf der "Vorderseite" der Kugel sind sichtbar.
        if point_rotated.z < 0.0 {
            // TODO: Entscheiden, wie mit unsichtbaren Punkten umgegangen wird.
            // Option 1: Fehler zurückgeben
            // Option 2: Speziellen Wert zurückgeben (z.B. NaN, oder außerhalb eines definierten Bereichs)
            // Option 3: Den Punkt am Rand projizieren (kann zu Artefakten führen)
            // Fürs Erste geben wir einen Fehler zurück, da dies oft unerwartet ist.
            return Err(MathError::GeometricFailure {
                operation: "Orthographic projection of a point on the back hemisphere".to_string(),
            });
        }

        Ok(Vec2::new(point_rotated.x, point_rotated.y))
    }

    fn unproject_orthographic(&self, point_2d: Vec2, view_direction: Vec3) -> MathResult<Vec3> {
        // Gegeben x', y' in der Projektionsebene.
        // z' muss berechnet werden, da der Punkt auf der Kugeloberfläche liegen soll.
        // x'^2 + y'^2 + z'^2 = R^2  => z' = sqrt(R^2 - (x'^2 + y'^2))
        // Wir nehmen die positive Wurzel, da wir auf die Vorderseite projizieren.
        let r_sq = self.sphere_radius * self.sphere_radius;
        let xy_sq_sum = point_2d.length_squared();

        if xy_sq_sum > r_sq + constants::EPSILON {
            // Der 2D-Punkt liegt außerhalb des Kugelradius in der Projektionsebene.
            return Err(MathError::GeometricFailure {
                operation: "Orthographic unprojection: 2D point is outside the sphere's disk"
                    .to_string(),
            });
        }

        // um numerische Probleme bei xy_sq_sum sehr nahe an r_sq zu vermeiden:
        let z_rotated_sq = r_sq - xy_sq_sum;
        let z_rotated = if z_rotated_sq < 0.0 {
            0.0
        } else {
            z_rotated_sq.sqrt()
        };

        let point_on_sphere_view_aligned = Vec3::new(point_2d.x, point_2d.y, z_rotated);

        // Rotiere zurück in das ursprüngliche Koordinatensystem.
        // Die ursprüngliche Rotation war von `view_direction` zu `Vec3::Z`.
        // Die inverse Rotation ist von `Vec3::Z` zu `view_direction`.
        let rotation_to_original_view = BevyQuat::from_rotation_arc(
            Vec3::Z,
            view_direction, // view_direction sollte bereits normalisiert sein
        );
        let point_on_sphere_rotated_back = rotation_to_original_view * point_on_sphere_view_aligned;

        Ok(point_on_sphere_rotated_back.normalize_or_zero() * self.sphere_radius) // Sicherstellen, dass es auf der Kugel liegt
    }

    fn project_mercator(&self, point: Vec3) -> MathResult<Vec2> {
        let geo = GeographicCoordinates::from_cartesian(point, self.sphere_radius);
        // Mercator kann die Pole nicht darstellen (log(tan(PI/2)))
        if geo.latitude.abs() >= constants::PI_OVER_2 - constants::EPSILON * 10.0 {
            return Err(MathError::GeometricFailure {
                operation: "Mercator projection cannot represent poles.".to_string(),
            });
        }
        let x = self.sphere_radius * geo.longitude;
        let y = self.sphere_radius * ((constants::PI_OVER_2 + geo.latitude) * 0.5).tan().ln();
        Ok(Vec2::new(x, y))
    }

    fn unproject_mercator(&self, point_2d: Vec2) -> MathResult<Vec3> {
        let longitude = point_2d.x / self.sphere_radius;
        let latitude = 2.0 * (point_2d.y / self.sphere_radius).exp().atan() - constants::PI_OVER_2;

        let geo = GeographicCoordinates::new(latitude, longitude, 0.0); // Altitude 0
        Ok(geo.to_cartesian(self.sphere_radius))
    }

    fn project_lambert_azimuthal(&self, point: Vec3, center_on_sphere: Vec3) -> MathResult<Vec2> {
        let cos_c = point.normalize().dot(center_on_sphere); // Winkel zwischen Punkt und Zentrum
        if cos_c < -1.0 + constants::EPSILON {
            // Punkt ist antipodal zum Zentrum
            return Err(MathError::GeometricFailure {
                operation: "Lambert Azimuthal: point is antipodal to projection center."
                    .to_string(),
            });
        }

        let k_prime = (2.0 / (1.0 + cos_c)).sqrt();

        // Rotiere System, sodass center_on_sphere = (0,0,R) ist, dann x'=k*x, y'=k*y
        // Einfacher: Projektion auf Tangentialebene am center_on_sphere
        let (u_axis, v_axis) = CoordinateConverter::sphere_tangents_at_point(center_on_sphere);

        // Vektor vom Kugelzentrum zum Punkt
        let _r_vec = point;

        // Die Formel ist x = R * k' * cos(lat_p) * sin(lon_p - lon_0)
        // y = R * k' * (cos(lat_0)*sin(lat_p) - sin(lat_0)*cos(lat_p)*cos(lon_p-lon_0))
        // Wobei lat_0, lon_0 das Zentrum sind.
        // Alternative:
        Ok(Vec2::new(
            self.sphere_radius * k_prime * point.dot(u_axis.normalize_or_zero()),
            self.sphere_radius * k_prime * point.dot(v_axis.normalize_or_zero()),
        ))
    }
    fn unproject_lambert_azimuthal(
        &self,
        point_2d: Vec2,
        center_on_sphere: Vec3,
    ) -> MathResult<Vec3> {
        let r_sphere = self.sphere_radius;
        let x_proj = point_2d.x;
        let y_proj = point_2d.y;

        let rho_sq = x_proj * x_proj + y_proj * y_proj;
        if rho_sq > (2.0 * r_sphere).powi(2) + constants::EPSILON {
            // Punkte außerhalb des max. proj. Radius (2R)
            return Err(MathError::InvalidConfiguration {
                message: "Point is outside the valid domain for Lambert Azimuthal unprojection"
                    .to_string(),
            });
        }

        let rho = rho_sq.sqrt();
        if rho < constants::EPSILON {
            // Punkt ist im Projektionszentrum
            return Ok(center_on_sphere.normalize_or_zero() * r_sphere);
        }

        // c ist der Winkelabstand vom Projektionszentrum zum Punkt auf der Kugel
        // Für Lambert Azimuthal Equal Area: rho = 2 * R * sin(c/2)
        // Also sin(c/2) = rho / (2 * R)
        // c/2 = asin(rho / (2 * R))
        // c = 2 * asin(rho / (2 * R))
        let c_half = (rho / (2.0 * r_sphere)).clamp(-1.0, 1.0); // Clamp für asin
        let c = 2.0 * c_half.asin();

        // Sinus und Kosinus von c
        let sin_c = c.sin();
        let cos_c = c.cos();

        // Geographische Koordinaten des Projektionszentrums (lat0, lon0)
        // center_on_sphere muss auf geografische Koordinaten umgerechnet werden.
        // Annahme: SphereProjector speichert center_on_sphere als Vec3 auf der Kugeloberfläche.
        // Konvertiere center_on_sphere zu sphärischen Koordinaten, dann zu geographischen.
        let geo_center = GeographicCoordinates::from_cartesian(center_on_sphere, r_sphere);
        let lat0 = geo_center.latitude;
        let lon0 = geo_center.longitude;
        let cos_lat0 = lat0.cos();
        let sin_lat0 = lat0.sin();

        // Inverse Formeln (aus https://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html oder ähnlichen Quellen)
        // lat = asin(cos(c) * sin(lat0) + (y_proj * sin(c) * cos(lat0)) / rho)
        // lon = lon0 + atan2(x_proj * sin(c), rho * cos(lat0) * cos(c) - y_proj * sin(lat0) * sin(c))

        let latitude = (cos_c * sin_lat0 + (y_proj * sin_c * cos_lat0) / rho)
            .clamp(-1.0, 1.0)
            .asin();

        let lon_denom = rho * cos_lat0 * cos_c - y_proj * sin_lat0 * sin_c;
        let longitude = if lon_denom.abs() < constants::EPSILON
            && (x_proj * sin_c).abs() < constants::EPSILON
        {
            lon0 // Wenn sowohl Zähler als auch Nenner für atan2 klein sind, ist der Winkel unbestimmt -> nimm lon0
        } else {
            lon0 + (x_proj * sin_c).atan2(lon_denom)
        };

        let result_geo = GeographicCoordinates::new(latitude, longitude, 0.0); // Altitude 0
        Ok(result_geo.to_cartesian(r_sphere))
    }
    fn project_gnomonic(&self, point: Vec3, center_on_sphere: Vec3) -> MathResult<Vec2> {
        let cos_c = point.normalize().dot(center_on_sphere);
        if cos_c <= constants::EPSILON {
            // Punkt auf oder hinter dem Horizont des Projektionszentrums
            return Err(MathError::GeometricFailure {
                operation: "Gnomonic projection: point is on or behind the horizon.".to_string(),
            });
        }
        let (u_axis, v_axis) = CoordinateConverter::sphere_tangents_at_point(center_on_sphere);
        // x = R * tan(theta_x), y = R * tan(theta_y) in der Tangentialebene
        // oder x = R * x_rotated / z_rotated
        Ok(Vec2::new(
            self.sphere_radius * point.dot(u_axis.normalize_or_zero()) / cos_c,
            self.sphere_radius * point.dot(v_axis.normalize_or_zero()) / cos_c,
        ))
    }

    fn project_mollweide(&self, point: Vec3) -> MathResult<Vec2> {
        let geo = GeographicCoordinates::from_cartesian(point, self.sphere_radius);
        let mut theta_aux = geo.latitude; // Hilfswinkel

        // Iterative Lösung für den Hilfswinkel theta_aux
        for _ in 0..10 {
            // Max 10 Iterationen (typischerweise schnell konvergent)
            let delta_theta = -(2.0 * theta_aux + (2.0 * theta_aux).sin()
                - PI * geo.latitude.sin())
                / (2.0 + 2.0 * (2.0 * theta_aux).cos());
            theta_aux += delta_theta;
            if delta_theta.abs() < constants::EPSILON * 10.0 {
                break;
            }
        }
        let x =
            self.sphere_radius * (2.0 * 2.0_f32.sqrt() / PI) * geo.longitude * (theta_aux).cos();
        let y = self.sphere_radius * 2.0_f32.sqrt() * (theta_aux).sin();
        Ok(Vec2::new(x, y))
    }

    fn project_equirectangular(&self, point: Vec3) -> MathResult<Vec2> {
        let geo = GeographicCoordinates::from_cartesian(point, self.sphere_radius);
        // x = R * lambda
        // y = R * phi
        Ok(Vec2::new(
            self.sphere_radius * geo.longitude,
            self.sphere_radius * geo.latitude,
        ))
    }

    fn unproject_equirectangular(&self, point_2d: Vec2) -> MathResult<Vec3> {
        let longitude = point_2d.x / self.sphere_radius;
        let latitude = point_2d.y / self.sphere_radius;

        // Überprüfe Gültigkeit der Winkel
        if latitude.abs() > constants::PI_OVER_2 + constants::EPSILON {
            return Err(MathError::InvalidConfiguration {
                message: "Latitude out of bounds for inverse equirectangular projection."
                    .to_string(),
            });
        }
        if longitude.abs() > PI + constants::EPSILON {
            return Err(MathError::InvalidConfiguration {
                message: "Longitude out of bounds for inverse equirectangular projection."
                    .to_string(),
            });
        }

        let geo = GeographicCoordinates::new(latitude, longitude, 0.0); // Altitude 0
        Ok(geo.to_cartesian(self.sphere_radius))
    }
}

/// UV-Mapping für Texturen auf Kugeln.
pub struct SphereUVMapper;

impl SphereUVMapper {
    /// Standard sphärische UV-Koordinaten (basierend auf Längen- und Breitengrad).
    /// `point_on_sphere` sollte ein Punkt auf der Kugeloberfläche sein.
    /// Gibt (u, v) im Bereich [0, 1] zurück.
    pub fn spherical_uv(point_on_sphere: Vec3) -> Vec2 {
        // Nimmt an, dass point_on_sphere bereits auf der Kugel ist, Radius spielt keine Rolle für Winkel.
        let geo = GeographicCoordinates::from_cartesian(point_on_sphere, 1.0); // Referenzradius 1.0 für Winkel
        Vec2::new(
            (geo.longitude + PI) / TAU, // Longitude von [-PI, PI] nach [0, 1]
            (geo.latitude + constants::PI_OVER_2) / PI, // Latitude von [-PI/2, PI/2] nach [0, 1]
        )
    }

    /// Kubische UV-Koordinaten (für Skybox-ähnliche Projektionen).
    /// `point_from_center` ist ein Vektor vom Kugelzentrum zum Punkt.
    /// Gibt (face_index, uv_on_face) zurück. face_index: 0:+X, 1:-X, 2:+Y, 3:-Y, 4:+Z, 5:-Z.
    pub fn cubic_uv(point_from_center: Vec3) -> (usize, Vec2) {
        let abs_point = point_from_center.abs();
        let point = point_from_center; // Kürzerer Alias

        let (face_index, u_raw, v_raw) = if abs_point.x >= abs_point.y && abs_point.x >= abs_point.z
        {
            // X-dominante Achse
            if point.x > 0.0 {
                (0, -point.z / point.x, -point.y / point.x)
            }
            // +X
            else {
                (1, point.z / point.x, -point.y / point.x)
            } // -X
        } else if abs_point.y >= abs_point.z {
            // Y-dominante Achse
            if point.y > 0.0 {
                (2, point.x / point.y, point.z / point.y)
            }
            // +Y
            else {
                (3, point.x / point.y, -point.z / point.y)
            } // -Y
        } else {
            // Z-dominante Achse
            if point.z > 0.0 {
                (4, point.x / point.z, -point.y / point.z)
            }
            // +Z
            else {
                (5, -point.x / point.z, -point.y / point.z)
            } // -Z
        };

        let _ = Vec2::new((u_raw + 1.0) * 0.5, (v_raw + 1.0) * 0.5); // Normalisiere auf [0,1]
        (
            face_index,
            Vec2::new((u_raw + 1.0) * 0.5, (v_raw + 1.0) * 0.5),
        )
    }

    /// Zylindrische UV-Koordinaten.
    /// `point_on_cylinder` ist ein Punkt, dessen UVs berechnet werden sollen.
    /// `height_range` ist (min_y, max_y) des Zylinders.
    pub fn cylindrical_uv(point_on_cylinder: Vec3, height_range: (f32, f32)) -> Vec2 {
        let angle_rad = point_on_cylinder.z.atan2(point_on_cylinder.x); // Winkel in der XZ-Ebene
        let u = angles::normalize_angle(angle_rad) / TAU; // u von [0, 1]

        let v = if (height_range.1 - height_range.0).abs() < constants::EPSILON {
            0.5 // Fallback, wenn Höhe null ist
        } else {
            ((point_on_cylinder.y - height_range.0) / (height_range.1 - height_range.0))
                .clamp(0.0, 1.0)
        };
        Vec2::new(u, v)
    }
}

/// Hilfsfunktionen für Projektionen.
pub struct SphereProjectionUtils;

/// Hilfsfunktionen für Projektionen.
impl SphereProjectionUtils {
    /// Berechnet eine Näherung der Flächenverzerrung einer Projektion an einem gegebenen Punkt auf der Kugel.
    /// `projector` ist der verwendete Kugel-Projektor.
    /// `point_on_sphere` ist der Punkt, an dem die Verzerrung berechnet wird.
    /// `angular_delta` ist ein kleiner Winkel-Offset für die finite Differenzenberechnung (in Radiant).
    pub fn calculate_area_distortion(
        projector: &SphereProjector,
        point_on_sphere: Vec3,
        angular_delta: f32,
    ) -> MathResult<f32> {
        let p0_2d = projector.project_point(point_on_sphere)?;

        // Erzeuge zwei orthogonale Tangentenvektoren am Punkt auf der Kugel
        // point_on_sphere sollte hier idealerweise schon normalisiert sein oder zumindest nicht der Nullvektor.
        let (tangent_u, tangent_v) = CoordinateConverter::sphere_tangents_at_point(point_on_sphere);

        // Kleine Rotationen, um benachbarte Punkte auf der Kugeloberfläche zu erhalten
        let rotation_around_v_axis = BevyQuat::from_axis_angle(tangent_v, angular_delta);
        let point_offset_u_3d = rotation_around_v_axis * point_on_sphere;

        let rotation_around_u_axis = BevyQuat::from_axis_angle(tangent_u, -angular_delta); // Negative Drehung für die andere Richtung
        let point_offset_v_3d = rotation_around_u_axis * point_on_sphere;

        // Projiziere die verschobenen Punkte
        // project_point normalisiert intern auf den sphere_radius des Projektors.
        let p_u_2d = projector.project_point(point_offset_u_3d)?;
        let p_v_2d = projector.project_point(point_offset_v_3d)?;

        // Vektoren in der Projektionsebene
        let du_vec = p_u_2d - p0_2d;
        let dv_vec = p_v_2d - p0_2d;

        // Fläche des Parallelogramms, aufgespannt durch du_vec, dv_vec (Betrag des 2D-Kreuzprodukts)
        let projected_area_element_approx = (du_vec.x * dv_vec.y - du_vec.y * dv_vec.x).abs();

        // Fläche des ursprünglichen "quadratischen" Elements auf der Kugel (angenähert)
        // Bogenlänge = Winkel * Radius
        let side_length_on_sphere_approx = angular_delta * projector.sphere_radius;
        let original_area_element_approx =
            side_length_on_sphere_approx * side_length_on_sphere_approx;

        if original_area_element_approx < constants::EPSILON * constants::EPSILON {
            // Vermeide Division durch Null, wenn Originalfläche winzig ist.
            // In diesem Fall ist die Verzerrung schwer zu definieren oder sollte als 1 (keine Verzerrung) betrachtet werden.
            return Ok(1.0);
        }
        Ok(projected_area_element_approx / original_area_element_approx)
    }

    /// Findet den optimalen Projektionszentrum-Vektor für eine gegebene Menge von Punkten.
    /// Gibt den normalisierten Schwerpunkt der Punkte zurück.
    pub fn optimize_projection_center(points: &[Vec3]) -> Option<Vec3> {
        if points.is_empty() {
            return None;
        }
        let sum = points.iter().fold(Vec3::ZERO, |acc, &p| acc + p);
        let centroid = sum / points.len() as f32;
        Some(centroid.normalize_or_zero()) // normalize_or_zero für den Fall, dass der Centroid (0,0,0) ist
    }

    /// Testet verschiedene Projektionstypen und gibt diejenige mit der im Durchschnitt geringsten Flächenverzerrung zurück.
    pub fn find_best_projection_type(
        points_on_sphere: &[Vec3],
        sphere_radius: f32,
        angular_delta_for_distortion: f32, // z.B. 0.01 Radiant
    ) -> MathResult<(SphereProjectionType, f32)> {
        // (Beste Projektion, Durchschnittliche Verzerrung)
        if points_on_sphere.is_empty() {
            return Err(MathError::InsufficientPoints {
                expected: 1,
                actual: 0,
            });
        }

        let optimized_center =
            Self::optimize_projection_center(points_on_sphere).unwrap_or(Vec3::Z); // Fallback auf Z-Achse (Nordpol)

        let projection_candidates = [
            // Array statt Vec für feste Größe
            SphereProjectionType::Stereographic {
                pole_axis: optimized_center,
            },
            SphereProjectionType::Orthographic {
                view_direction: optimized_center,
            },
            SphereProjectionType::LambertAzimuthal {
                center_on_sphere: optimized_center,
            },
            SphereProjectionType::Gnomonic {
                center_on_sphere: optimized_center,
            },
            SphereProjectionType::Equirectangular, // Hat kein Zentrum als Parameter
            SphereProjectionType::Mollweide,       // Hat kein Zentrum als Parameter
                                                   // Mercator ist oft problematisch wegen Polen, kann aber hinzugefügt werden
                                                   // SphereProjectionType::Mercator,
        ];

        let mut best_projection_type = projection_candidates[0];
        let mut min_avg_distortion = f32::INFINITY;

        for &current_projection_type in &projection_candidates {
            let projector = SphereProjector::new(current_projection_type, sphere_radius)?; // Kann fehlschlagen, wenn Radius ungültig

            let mut total_distortion_metric = 0.0;
            let mut successful_calculations = 0;

            for &point in points_on_sphere {
                // Stelle sicher, dass der Punkt für die Verzerrungsberechnung nicht der Nullvektor ist
                if point.length_squared() < constants::EPSILON * constants::EPSILON {
                    continue;
                }
                match Self::calculate_area_distortion(
                    &projector,
                    point,
                    angular_delta_for_distortion,
                ) {
                    Ok(distortion) => {
                        // Wir wollen Verzerrung minimieren. Ein Wert von 1 ist ideal (keine Verzerrung).
                        // Daher ist (distortion - 1.0)^2 ein guter Metrikwert, den wir minimieren wollen.
                        total_distortion_metric += (distortion - 1.0).powi(2);
                        successful_calculations += 1;
                    }
                    Err(_) => {
                        // Punkt konnte nicht projiziert oder Verzerrung nicht berechnet werden (z.B. an Singularität)
                        // Dies könnte die Projektion "bestrafen", indem wir einen hohen Wert addieren.
                        total_distortion_metric += 100.0; // Hoher Strafwert
                        successful_calculations += 1; // Zähle es trotzdem, um nicht durch 0 zu teilen
                    }
                }
            }

            if successful_calculations > 0 {
                let avg_distortion_metric =
                    total_distortion_metric / successful_calculations as f32;
                if avg_distortion_metric < min_avg_distortion {
                    min_avg_distortion = avg_distortion_metric;
                    best_projection_type = current_projection_type;
                }
            }
        }

        if min_avg_distortion == f32::INFINITY {
            return Err(MathError::GeometricFailure {
                operation:
                    "Could not find a suitable projection or calculate distortion for any point."
                        .to_string(),
            });
        }

        Ok((best_projection_type, min_avg_distortion.sqrt())) // sqrt, um zur ursprünglichen Skala der Verzerrung zurückzukehren
    }
}
