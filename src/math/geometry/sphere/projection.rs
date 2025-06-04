// src/math/sphere/projection.rs

use crate::math::{error::*, geometry::sphere::coordinates::*, types::*, utils::*};
use std::f32::consts::{PI, TAU};

/// Verschiedene Projektionstypen von Kugel zu Ebene
#[derive(Debug, Clone, Copy)]
pub enum ProjectionType {
    /// Stereographische Projektion
    Stereographic { pole: Point3D },
    /// Orthographische Projektion (parallel)
    Orthographic { view_direction: Point3D },
    /// Merkator-Projektion
    Mercator,
    /// Lambert Azimuthal Equal-Area
    LambertAzimuthal { center: Point3D },
    /// Gnomonic Projektion (Zentral-Projektion)
    Gnomonic { center: Point3D },
    /// Mollweide-Projektion (Equal-Area)
    Mollweide,
    /// Equirectangular (Plate Carrée)
    Equirectangular,
}

/// Projektions-Engine für Sphere-zu-Plane Transformationen
pub struct SphereProjector {
    projection_type: ProjectionType,
    sphere_radius: f32,
    output_scale: f32,
}

impl SphereProjector {
    /// Erstellt einen neuen Projektor
    pub fn new(projection_type: ProjectionType, sphere_radius: f32) -> Self {
        Self {
            projection_type,
            sphere_radius,
            output_scale: 1.0,
        }
    }

    /// Setzt den Ausgabe-Skalierungsfaktor
    pub fn with_scale(mut self, scale: f32) -> Self {
        self.output_scale = scale;
        self
    }

    /// Projiziert einen 3D-Punkt auf die 2D-Ebene
    pub fn project_point(&self, point: Point3D) -> MathResult<Point2D> {
        // Normiere auf Kugel-Radius
        let normalized_point = if point.length() > constants::EPSILON {
            point.normalize() * self.sphere_radius
        } else {
            return Err(MathError::InvalidConfiguration {
                message: "Cannot project zero vector".to_string(),
            });
        };

        let projected = match self.projection_type {
            ProjectionType::Stereographic { pole } => {
                self.stereographic_projection(normalized_point, pole)?
            }
            ProjectionType::Orthographic { view_direction } => {
                self.orthographic_projection(normalized_point, view_direction)?
            }
            ProjectionType::Mercator => self.mercator_projection(normalized_point)?,
            ProjectionType::LambertAzimuthal { center } => {
                self.lambert_azimuthal_projection(normalized_point, center)?
            }
            ProjectionType::Gnomonic { center } => {
                self.gnomonic_projection(normalized_point, center)?
            }
            ProjectionType::Mollweide => self.mollweide_projection(normalized_point)?,
            ProjectionType::Equirectangular => self.equirectangular_projection(normalized_point)?,
        };

        Ok(projected * self.output_scale)
    }

    /// Projiziert eine Liste von Punkten
    pub fn project_points(&self, points: &[Point3D]) -> MathResult<Vec<Point2D>> {
        points
            .iter()
            .map(|&point| self.project_point(point))
            .collect()
    }

    /// Inverse Projektion (2D zu 3D) - nicht für alle Projektionen verfügbar
    pub fn unproject_point(&self, point: Point2D) -> MathResult<Point3D> {
        let scaled_point = point / self.output_scale;

        match self.projection_type {
            ProjectionType::Stereographic { pole } => {
                self.inverse_stereographic_projection(scaled_point, pole)
            }
            ProjectionType::Orthographic { view_direction } => {
                self.inverse_orthographic_projection(scaled_point, view_direction)
            }
            ProjectionType::Equirectangular => {
                self.inverse_equirectangular_projection(scaled_point)
            }
            _ => Err(MathError::InvalidConfiguration {
                message: "Inverse projection not implemented for this projection type".to_string(),
            }),
        }
    }

    // === Projektions-Implementierungen ===

    fn stereographic_projection(&self, point: Point3D, pole: Point3D) -> MathResult<Point2D> {
        let pole_normalized = pole.normalize();

        // Prüfe ob Punkt am Pol liegt
        if (point.normalize().dot(pole_normalized) - 1.0).abs() < constants::EPSILON {
            return Err(MathError::GeometricFailure {
                operation: "Stereographic projection at pole".to_string(),
            });
        }

        // Stereographische Projektion von Nordpol aus
        let z_offset = if pole_normalized.z > 0.0 { 1.0 } else { -1.0 };
        let denominator = self.sphere_radius * (z_offset + point.z / self.sphere_radius);

        if denominator.abs() < constants::EPSILON {
            return Err(MathError::GeometricFailure {
                operation: "Division by zero in stereographic projection".to_string(),
            });
        }

        Ok(Point2D::new(point.x / denominator, point.y / denominator))
    }

    fn inverse_stereographic_projection(
        &self,
        point: Point2D,
        pole: Point3D,
    ) -> MathResult<Point3D> {
        let r_sq = point.x * point.x + point.y * point.y;
        let factor = 2.0 * self.sphere_radius / (1.0 + r_sq);

        Ok(Point3D::new(
            factor * point.x,
            factor * point.y,
            self.sphere_radius * (r_sq - 1.0) / (r_sq + 1.0),
        ))
    }

    fn orthographic_projection(
        &self,
        point: Point3D,
        view_direction: Point3D,
    ) -> MathResult<Point2D> {
        let view_dir = view_direction.normalize();

        // Prüfe ob Punkt auf der Rückseite liegt
        if point.dot(view_dir) < 0.0 {
            return Err(MathError::GeometricFailure {
                operation: "Point on back side of sphere".to_string(),
            });
        }

        // Orthogonale Basis im Tangentialraum
        let (tangent1, tangent2) = CoordinateConverter::sphere_tangents(view_dir);

        Ok(Point2D::new(point.dot(tangent1), point.dot(tangent2)))
    }

    fn inverse_orthographic_projection(
        &self,
        point: Point2D,
        view_direction: Point3D,
    ) -> MathResult<Point3D> {
        let view_dir = view_direction.normalize();
        let (tangent1, tangent2) = CoordinateConverter::sphere_tangents(view_dir);

        let r_sq = point.x * point.x + point.y * point.y;
        if r_sq > self.sphere_radius * self.sphere_radius {
            return Err(MathError::GeometricFailure {
                operation: "Point outside projection circle".to_string(),
            });
        }

        let z = (self.sphere_radius * self.sphere_radius - r_sq).sqrt();
        let surface_point = tangent1 * point.x + tangent2 * point.y + view_dir * z;

        Ok(surface_point)
    }

    fn mercator_projection(&self, point: Point3D) -> MathResult<Point2D> {
        let geographic = GeographicCoordinates::from_cartesian(point, self.sphere_radius);

        // Mercator kann nicht die Pole projizieren
        if geographic.latitude.abs() > PI * 0.4999 {
            return Err(MathError::GeometricFailure {
                operation: "Mercator projection near poles".to_string(),
            });
        }

        let x = geographic.longitude;
        let y = (PI * 0.25 + geographic.latitude * 0.5).tan().ln();

        Ok(Point2D::new(x, y))
    }

    fn lambert_azimuthal_projection(&self, point: Point3D, center: Point3D) -> MathResult<Point2D> {
        let center_normalized = center.normalize();
        let point_normalized = point.normalize();

        let cos_c = center_normalized.dot(point_normalized);

        // Punkt liegt auf gegenüberliegender Seite
        if cos_c < -0.999 {
            return Err(MathError::GeometricFailure {
                operation: "Lambert projection at antipodal point".to_string(),
            });
        }

        let k = (2.0 / (1.0 + cos_c)).sqrt();

        // Berechne lokale Koordinaten
        let (tangent1, tangent2) = CoordinateConverter::sphere_tangents(center_normalized);

        let local_x = point_normalized.dot(tangent1);
        let local_y = point_normalized.dot(tangent2);

        Ok(Point2D::new(k * local_x, k * local_y))
    }

    fn gnomonic_projection(&self, point: Point3D, center: Point3D) -> MathResult<Point2D> {
        let center_normalized = center.normalize();
        let point_normalized = point.normalize();

        let cos_c = center_normalized.dot(point_normalized);

        // Punkt liegt auf der Rückseite oder am Horizont
        if cos_c <= constants::EPSILON {
            return Err(MathError::GeometricFailure {
                operation: "Gnomonic projection at horizon or behind".to_string(),
            });
        }

        let (tangent1, tangent2) = CoordinateConverter::sphere_tangents(center_normalized);

        let local_x = point_normalized.dot(tangent1) / cos_c;
        let local_y = point_normalized.dot(tangent2) / cos_c;

        Ok(Point2D::new(local_x, local_y))
    }

    fn mollweide_projection(&self, point: Point3D) -> MathResult<Point2D> {
        let geographic = GeographicCoordinates::from_cartesian(point, self.sphere_radius);

        // Newton-Iteration für Mollweide
        let mut theta = geographic.latitude;
        for _ in 0..10 {
            let sin_theta = theta.sin();
            let cos_theta = theta.cos();
            let delta = -(theta + sin_theta - PI * geographic.latitude.sin()) / (1.0 + cos_theta);
            theta += delta;
            if delta.abs() < constants::EPSILON {
                break;
            }
        }

        let sqrt_2 = 2.0_f32.sqrt();
        let x = 2.0 * sqrt_2 * geographic.longitude * theta.cos() / PI;
        let y = sqrt_2 * theta.sin();

        Ok(Point2D::new(x, y))
    }

    fn equirectangular_projection(&self, point: Point3D) -> MathResult<Point2D> {
        let geographic = GeographicCoordinates::from_cartesian(point, self.sphere_radius);

        Ok(Point2D::new(geographic.longitude, geographic.latitude))
    }

    fn inverse_equirectangular_projection(&self, point: Point2D) -> MathResult<Point3D> {
        let geographic = GeographicCoordinates::new(point.y, point.x, 0.0);
        Ok(geographic.to_cartesian(self.sphere_radius))
    }
}

/// UV-Mapping für Texturen auf Kugeln
pub struct SphereUVMapper;

impl SphereUVMapper {
    /// Standard sphärische UV-Koordinaten
    pub fn spherical_uv(point: Point3D, radius: f32) -> Point2D {
        let geographic = GeographicCoordinates::from_cartesian(point, radius);

        Point2D::new(
            (geographic.longitude / TAU + 0.5) % 1.0,
            geographic.latitude / PI + 0.5,
        )
    }

    /// Kubische UV-Koordinaten (Skybox-Style)
    pub fn cubic_uv(point: Point3D) -> (usize, Point2D) {
        let abs_point = Point3D::new(point.x.abs(), point.y.abs(), point.z.abs());

        let (face, u, v) = if abs_point.x >= abs_point.y && abs_point.x >= abs_point.z {
            // X-Faces
            if point.x > 0.0 {
                (0, -point.z / point.x, -point.y / point.x) // +X
            } else {
                (1, point.z / point.x, -point.y / point.x) // -X
            }
        } else if abs_point.y >= abs_point.z {
            // Y-Faces
            if point.y > 0.0 {
                (2, point.x / point.y, point.z / point.y) // +Y
            } else {
                (3, point.x / point.y, -point.z / point.y) // -Y
            }
        } else {
            // Z-Faces
            if point.z > 0.0 {
                (4, point.x / point.z, -point.y / point.z) // +Z
            } else {
                (5, -point.x / point.z, -point.y / point.z) // -Z
            }
        };

        // Konvertiere zu [0,1] Bereich
        let u_normalized = (u + 1.0) * 0.5;
        let v_normalized = (v + 1.0) * 0.5;

        (face, Point2D::new(u_normalized, v_normalized))
    }

    /// Zylindrische UV-Koordinaten
    pub fn cylindrical_uv(point: Point3D, height_range: (f32, f32)) -> Point2D {
        let angle = point.z.atan2(point.x);
        let u = (angle / TAU + 0.5) % 1.0;

        let v = (point.y - height_range.0) / (height_range.1 - height_range.0);
        let v_clamped = v.clamp(0.0, 1.0);

        Point2D::new(u, v_clamped)
    }
}

/// Spezielle Projektions-Utilities
pub struct ProjectionUtils;

impl ProjectionUtils {
    /// Berechnet die Verzerrung einer Projektion an einem Punkt
    pub fn calculate_distortion(
        projector: &SphereProjector,
        point: Point3D,
        delta: f32,
    ) -> MathResult<f32> {
        let projected_center = projector.project_point(point)?;

        // Berechne Jacobian durch finite Differenzen
        let (tangent1, tangent2) = CoordinateConverter::sphere_tangents(point);

        let point_u = point + tangent1 * delta;
        let point_v = point + tangent2 * delta;

        let projected_u = projector.project_point(point_u.normalize() * projector.sphere_radius)?;
        let projected_v = projector.project_point(point_v.normalize() * projector.sphere_radius)?;

        let du_dx = (projected_u - projected_center) / delta;
        let dv_dy = (projected_v - projected_center) / delta;

        // Determinante der Jacobian-Matrix
        let determinant = du_dx.x * dv_dy.y - du_dx.y * dv_dy.x;

        Ok(determinant.abs())
    }

    /// Findet optimale Projektions-Parameter für gegebene Punkte
    pub fn optimize_projection_center(points: &[Point3D]) -> Option<Point3D> {
        if points.is_empty() {
            return None;
        }

        // Berechne den Schwerpunkt und projiziere auf Kugel
        let sum = points.iter().fold(Point3D::ZERO, |acc, &p| acc + p);
        let centroid = sum / points.len() as f32;

        Some(centroid.normalize())
    }

    /// Testet verschiedene Projektionen und gibt die beste zurück
    pub fn find_best_projection(
        points: &[Point3D],
        radius: f32,
    ) -> MathResult<(ProjectionType, f32)> {
        let center =
            Self::optimize_projection_center(points).unwrap_or(Point3D::new(0.0, 0.0, 1.0));

        let projections = vec![
            ProjectionType::Stereographic { pole: center },
            ProjectionType::Orthographic {
                view_direction: center,
            },
            ProjectionType::LambertAzimuthal { center },
            ProjectionType::Gnomonic { center },
            ProjectionType::Equirectangular,
        ];

        let mut best_projection = projections[0];
        let mut best_score = f32::INFINITY;

        for projection_type in projections {
            let projector = SphereProjector::new(projection_type, radius);

            let mut total_distortion = 0.0;
            let mut success_count = 0;

            for &point in points {
                if let Ok(distortion) = Self::calculate_distortion(&projector, point, 0.01) {
                    total_distortion += distortion;
                    success_count += 1;
                }
            }

            if success_count > 0 {
                let avg_distortion = total_distortion / success_count as f32;
                if avg_distortion < best_score {
                    best_score = avg_distortion;
                    best_projection = projection_type;
                }
            }
        }

        Ok((best_projection, best_score))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stereographic_projection() {
        let projector = SphereProjector::new(
            ProjectionType::Stereographic {
                pole: Point3D::new(0.0, 0.0, 1.0),
            },
            1.0,
        );

        let point = Point3D::new(1.0, 0.0, 0.0);
        let projected = projector.project_point(point).unwrap();

        // Stereographic sollte Punkte am Äquator korrekt projizieren
        assert!(projected.x.abs() > 0.0);
        assert!(comparison::nearly_equal(projected.y, 0.0));
    }

    #[test]
    fn test_orthographic_projection() {
        let projector = SphereProjector::new(
            ProjectionType::Orthographic {
                view_direction: Point3D::new(0.0, 0.0, 1.0),
            },
            1.0,
        );

        let point = Point3D::new(0.5, 0.0, (0.75_f32).sqrt());
        let projected = projector.project_point(point).unwrap();

        assert!(comparison::nearly_equal(projected.x, 0.5));
        assert!(comparison::nearly_equal(projected.y, 0.0));
    }

    #[test]
    fn test_equirectangular_projection() {
        let projector = SphereProjector::new(ProjectionType::Equirectangular, 1.0);

        let point = Point3D::new(1.0, 0.0, 0.0);
        let projected = projector.project_point(point).unwrap();
        let unprojected = projector.unproject_point(projected).unwrap();

        assert!(comparison::nearly_equal(unprojected.distance(point), 0.0));
    }

    #[test]
    fn test_sphere_uv_mapping() {
        let point = Point3D::new(1.0, 0.0, 0.0);
        let uv = SphereUVMapper::spherical_uv(point, 1.0);

        // Punkt bei (1,0,0) sollte UV (0.5, 0.5) haben
        assert!(comparison::nearly_equal(uv.x, 0.5));
        assert!(comparison::nearly_equal(uv.y, 0.5));
    }

    #[test]
    fn test_cubic_uv_mapping() {
        let point = Point3D::new(1.0, 0.5, 0.0);
        let (face, uv) = SphereUVMapper::cubic_uv(point);

        // Sollte auf +X Face landen
        assert_eq!(face, 0);
        assert!(uv.x >= 0.0 && uv.x <= 1.0);
        assert!(uv.y >= 0.0 && uv.y <= 1.0);
    }

    #[test]
    fn test_projection_optimization() {
        let points = vec![
            Point3D::new(1.0, 0.0, 0.0),
            Point3D::new(0.0, 1.0, 0.0),
            Point3D::new(0.0, 0.0, 1.0),
        ];

        let center = ProjectionUtils::optimize_projection_center(&points).unwrap();
        assert!(center.length() > 0.5); // Sollte normalisiert und sinnvoll sein

        let (best_projection, score) = ProjectionUtils::find_best_projection(&points, 1.0).unwrap();
        assert!(score < f32::INFINITY);
    }
}
