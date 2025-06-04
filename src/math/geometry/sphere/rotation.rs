
use crate::math::{error::*, types::*, utils::*};
use std::f32::consts::{PI, TAU};

/// Quaternion für 3D-Rotationen
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quaternion {
    pub w: f32,
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Quaternion {
    /// Erstellt ein neues Quaternion
    pub fn new(w: f32, x: f32, y: f32, z: f32) -> Self {
        Self { w, x, y, z }
    }
    
    /// Identitäts-Quaternion (keine Rotation)
    pub fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 0.0)
    }
    
    /// Erstellt Quaternion aus Achse und Winkel
    pub fn from_axis_angle(axis: Point3D, angle: f32) -> MathResult<Self> {
        let axis_length = axis.length();
        if axis_length < constants::EPSILON {
            return Err(MathError::InvalidConfiguration {
                message: "Rotation axis cannot be zero vector".to_string(),
            });
        }
        
        let normalized_axis = axis / axis_length;
        let half_angle = angle * 0.5;
        let sin_half = half_angle.sin();
        
        Ok(Self::new(
            half_angle.cos(),
            normalized_axis.x * sin_half,
            normalized_axis.y * sin_half,
            normalized_axis.z * sin_half,
        ))
    }
    
    /// Erstellt Quaternion aus Euler-Winkeln (Reihenfolge: Z-Y-X)
    pub fn from_euler(yaw: f32, pitch: f32, roll: f32) -> Self {
        let cy = (yaw * 0.5).cos();
        let sy = (yaw * 0.5).sin();
        let cp = (pitch * 0.5).cos();
        let sp = (pitch * 0.5).sin();
        let cr = (roll * 0.5).cos();
        let sr = (roll * 0.5).sin();
        
        Self::new(
            cr * cp * cy + sr * sp * sy,
            sr * cp * cy - cr * sp * sy,
            cr * sp * cy + sr * cp * sy,
            cr * cp * sy - sr * sp * cy,
        )
    }
    
    /// Erstellt Quaternion für Rotation von einem Vektor zu einem anderen
    pub fn from_vectors(from: Point3D, to: Point3D) -> MathResult<Self> {
        let from_normalized = from.normalize();
        let to_normalized = to.normalize();
        
        let dot = from_normalized.dot(to_normalized);
        
        // Vektoren sind bereits gleich
        if dot > 0.999999 {
            return Ok(Self::identity());
        }
        
        // Vektoren sind entgegengesetzt
        if dot < -0.999999 {
            // Finde orthogonalen Vektor
            let orthogonal = if from_normalized.x.abs() < 0.9 {
                Point3D::new(1.0, 0.0, 0.0)
            } else {
                Point3D::new(0.0, 1.0, 0.0)
            };
            
            let axis = from_normalized.cross(orthogonal).normalize();
            return Self::from_axis_angle(axis, PI);
        }
        
        let cross = from_normalized.cross(to_normalized);
        let w = 1.0 + dot;
        
        Ok(Self::new(w, cross.x, cross.y, cross.z).normalize())
    }
    
    /// Länge des Quaternions
    pub fn length(&self) -> f32 {
        (self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }
    
    /// Normalisiert das Quaternion
    pub fn normalize(&self) -> Self {
        let length = self.length();
        if length < constants::EPSILON {
            return Self::identity();
        }
        Self::new(
            self.w / length,
            self.x / length,
            self.y / length,
            self.z / length,
        )
    }
    
    /// Konjugat des Quaternions
    pub fn conjugate(&self) -> Self {
        Self::new(self.w, -self.x, -self.y, -self.z)
    }
    
    /// Inverses Quaternion
    pub fn inverse(&self) -> Self {
        let length_sq = self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z;
        if length_sq < constants::EPSILON {
            return Self::identity();
        }
        
        let conjugate = self.conjugate();
        Self::new(
            conjugate.w / length_sq,
            conjugate.x / length_sq,
            conjugate.y / length_sq,
            conjugate.z / length_sq,
        )
    }
    
    /// Multipliziert zwei Quaternions
    pub fn multiply(&self, other: &Quaternion) -> Self {
        Self::new(
            self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
            self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
            self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x,
            self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w,
        )
    }
    
    /// Rotiert einen Vektor
    pub fn rotate_vector(&self, vector: Point3D) -> Point3D {
        let quat_vector = Quaternion::new(0.0, vector.x, vector.y, vector.z);
        let result = self.multiply(&quat_vector).multiply(&self.conjugate());
        Point3D::new(result.x, result.y, result.z)
    }
    
    /// Konvertiert zu Achse und Winkel
    pub fn to_axis_angle(&self) -> (Point3D, f32) {
        let normalized = self.normalize();
        let angle = 2.0 * normalized.w.acos();
        
        let sin_half_angle = (1.0 - normalized.w * normalized.w).sqrt();
        if sin_half_angle < constants::EPSILON {
            return (Point3D::new(1.0, 0.0, 0.0), 0.0);
        }
        
        let axis = Point3D::new(
            normalized.x / sin_half_angle,
            normalized.y / sin_half_angle,
            normalized.z / sin_half_angle,
        );
        
        (axis, angle)
    }
    
    /// Konvertiert zu Euler-Winkeln (Z-Y-X Reihenfolge)
    pub fn to_euler(&self) -> (f32, f32, f32) {
        let normalized = self.normalize();
        
        // Roll (x-axis rotation)
        let sin_r_cp = 2.0 * (normalized.w * normalized.x + normalized.y * normalized.z);
        let cos_r_cp = 1.0 - 2.0 * (normalized.x * normalized.x + normalized.y * normalized.y);
        let roll = sin_r_cp.atan2(cos_r_cp);
        
        // Pitch (y-axis rotation)
        let sin_p = 2.0 * (normalized.w * normalized.y - normalized.z * normalized.x);
        let pitch = if sin_p.abs() >= 1.0 {
            sin_p.signum() * PI * 0.5 // Gimbal lock
        } else {
            sin_p.asin()
        };
        
        // Yaw (z-axis rotation)
        let sin_y_cp = 2.0 * (normalized.w * normalized.z + normalized.x * normalized.y);
        let cos_y_cp = 1.0 - 2.0 * (normalized.y * normalized.y + normalized.z * normalized.z);
        let yaw = sin_y_cp.atan2(cos_y_cp);
        
        (yaw, pitch, roll)
    }
    
    /// Spherical Linear Interpolation (SLERP)
    pub fn slerp(&self, other: &Quaternion, t: f32) -> Self {
        let mut dot = self.w * other.w + self.x * other.x + self.y * other.y + self.z * other.z;
        
        // Kürzester Pfad
        let other_corrected = if dot < 0.0 {
            dot = -dot;
            Quaternion::new(-other.w, -other.x, -other.y, -other.z)
        } else {
            *other
        };
        
        // Lineare Interpolation für nahe Quaternions
        if dot > 0.9995 {
            let result = Quaternion::new(
                self.w + t * (other_corrected.w - self.w),
                self.x + t * (other_corrected.x - self.x),
                self.y + t * (other_corrected.y - self.y),
                self.z + t * (other_corrected.z - self.z),
            );
            return result.normalize();
        }
        
        let angle = dot.acos();
        let sin_angle = angle.sin();
        
        let scale_a = ((1.0 - t) * angle).sin() / sin_angle;
        let scale_b = (t * angle).sin() / sin_angle;
        
        Quaternion::new(
            scale_a * self.w + scale_b * other_corrected.w,
            scale_a * self.x + scale_b * other_corrected.x,
            scale_a * self.y + scale_b * other_corrected.y,
            scale_a * self.z + scale_b * other_corrected.z,
        )
    }
}

/// Sphere-Rotation Utilities
pub struct SphereRotation;

impl SphereRotation {
    /// Rotiert einen Punkt auf der Kugel um eine Achse
    pub fn rotate_point_on_sphere(
        point: Point3D,
        axis: Point3D,
        angle: f32,
        radius: f32,
    ) -> MathResult<Point3D> {
        let rotation = Quaternion::from_axis_angle(axis, angle)?;
        let rotated = rotation.rotate_vector(point);
        
        // Stelle sicher, dass der Punkt auf der Kugel bleibt
        Ok(rotated.normalize() * radius)
    }
    
    /// Berechnet Rotation um einen großen Kreis
    pub fn great_circle_rotation(
        start: Point3D,
        end: Point3D,
        t: f32,
    ) -> MathResult<Point3D> {
        let rotation = Quaternion::from_vectors(start, end)?;
        Ok(rotation.slerp(&Quaternion::identity(), 1.0 - t).rotate_vector(start))
    }
    
    /// Rotiert um geografische Achsen
    pub fn rotate_by_geographic(
        point: Point3D,
        delta_latitude: f32,
        delta_longitude: f32,
        radius: f32,
    ) -> Point3D {
        // Latitude rotation (um X-Achse)
        let lat_rotation = Quaternion::from_axis_angle(Point3D::new(1.0, 0.0, 0.0), delta_latitude)
            .unwrap_or(Quaternion::identity());
        
        // Longitude rotation (um Y-Achse)
        let lon_rotation = Quaternion::from_axis_angle(Point3D::new(0.0, 1.0, 0.0), delta_longitude)
            .unwrap_or(Quaternion::identity());
        
        let total_rotation = lon_rotation.multiply(&lat_rotation);
        let rotated = total_rotation.rotate_vector(point);
        
        rotated.normalize() * radius
    }
    
    /// Generiert gleichmäßige Rotationen um alle Achsen
    pub fn uniform_rotations(count: usize, seed: Option<u64>) -> Vec<Quaternion> {
        use rand::{Rng, SeedableRng};
        
        let mut rng = match seed {
            Some(s) => rand::rngs::StdRng::seed_from_u64(s),
            None => rand::rngs::StdRng::from_entropy(),
        };
        
        (0..count)
            .map(|_| {
                // Marsaglia-Methode für uniform verteilte Quaternions
                let u1 = rng.gen::<f32>();
                let u2 = rng.gen_range(0.0..TAU);
                let u3 = rng.gen_range(0.0..TAU);
                
                let sqrt1_u1 = (1.0 - u1).sqrt();
                let sqrt_u1 = u1.sqrt();
                
                Quaternion::new(
                    sqrt1_u1 * u2.sin(),
                    sqrt1_u1 * u2.cos(),
                    sqrt_u1 * u3.sin(),
                    sqrt_u1 * u3.cos(),
                )
            })
            .collect()
    }
    
    /// Berechnet minimale Rotation zwischen zwei Orientierungen
    pub fn minimal_rotation(from: Quaternion, to: Quaternion) -> Quaternion {
        let relative = to.multiply(&from.inverse());
        relative.normalize()
    }
    
    /// Interpoliert zwischen mehreren Rotationen
    pub fn interpolate_rotations(rotations: &[Quaternion], t: f32) -> Option<Quaternion> {
        if rotations.is_empty() {
            return None;
        }
        
        if rotations.len() == 1 {
            return Some(rotations[0]);
        }
        
        let scaled_t = t * (rotations.len() - 1) as f32;
        let index = scaled_t.floor() as usize;
        let local_t = scaled_t - index as f32;
        
        if index >= rotations.len() - 1 {
            return Some(rotations[rotations.len() - 1]);
        }
        
        Some(rotations[index].slerp(&rotations[index + 1], local_t))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_quaternion_creation() {
        let q = Quaternion::from_axis_angle(Point3D::new(0.0, 0.0, 1.0), PI * 0.5).unwrap();
        let (axis, angle) = q.to_axis_angle();
        
        assert!(comparison::nearly_equal(angle, PI * 0.5));
        assert!(comparison::nearly_equal(axis.z, 1.0));
    }
    
    #[test]
    fn test_quaternion_rotation() {
        let q = Quaternion::from_axis_angle(Point3D::new(0.0, 0.0, 1.0), PI * 0.5).unwrap();
        let point = Point3D::new(1.0, 0.0, 0.0);
        let rotated = q.rotate_vector(point);
        
        assert!(comparison::nearly_equal(rotated.x, 0.0));
        assert!(comparison::nearly_equal(rotated.y, 1.0));
        assert!(comparison::nearly_equal(rotated.z, 0.0));
    }
    
    #[test]
    fn test_quaternion_slerp() {
        let q1 = Quaternion::identity();
        let q2 = Quaternion::from_axis_angle(Point3D::new(0.0, 0.0, 1.0), PI).unwrap();
        
        let interpolated = q1.slerp(&q2, 0.5);
        let (_, angle) = interpolated.to_axis_angle();
        
        assert!(comparison::nearly_equal(angle, PI * 0.5));
    }
    
    #[test]
    fn test_sphere_rotation() {
        let point = Point3D::new(1.0, 0.0, 0.0);
        let axis = Point3D::new(0.0, 0.0, 1.0);
        let rotated = SphereRotation::rotate_point_on_sphere(point, axis, PI * 0.5, 1.0).unwrap();
        
        assert!(comparison::nearly_equal(rotated.length(), 1.0));
        assert!(comparison::nearly_equal(rotated.x, 0.0));
        assert!(comparison::nearly_equal(rotated.y, 1.0));
    }
}