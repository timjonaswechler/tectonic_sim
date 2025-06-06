// src/math/geometry/polygon/transformations/projection.rs

use super::super::Polygon;
use crate::math::{error::*, geometry::sphere::projection::*, types::*, utils::*};
use bevy::math::{Vec2, Vec3};
use std::f32::consts::PI;
/// 2D zu 2D Projektionstypen
#[derive(Debug, Clone, Copy)]
pub enum PolygonProjection2DType {
    /// Orthogonale Projektion auf eine Linie
    OrthogonalLine { line_start: Vec2, line_end: Vec2 },
    /// Perspektivische Projektion von einem Punkt
    Perspective { center: Vec2, distance: f32 },
    /// Projektion auf einen Kreis
    CircularProjection { center: Vec2, radius: f32 },
    /// Logarithmische Spirale Projektion
    LogSpiral { center: Vec2, growth_factor: f32 },
    /// Konforme Mapping (komplexe Zahlen)
    Conformal { transform_type: ConformalType },
}

/// Konforme Transformationstypen
#[derive(Debug, Clone, Copy)]
pub enum ConformalType {
    /// Möbius-Transformation
    Mobius {
        a: Complex,
        b: Complex,
        c: Complex,
        d: Complex,
    },
    /// Joukowsky-Transformation
    Joukowsky,
    /// Logarithmische Transformation
    Logarithmic,
    /// Exponential-Transformation
    Exponential,
}

/// Komplexe Zahlen für konforme Mappings
#[derive(Debug, Clone, Copy)]
pub struct Complex {
    pub real: f32,
    pub imag: f32,
}

impl Complex {
    pub fn new(real: f32, imag: f32) -> Self {
        Self { real, imag }
    }

    pub fn from_point(point: Vec2) -> Self {
        Self {
            real: point.x,
            imag: point.y,
        }
    }

    pub fn to_point(&self) -> Vec2 {
        Vec2::new(self.real, self.imag)
    }

    pub fn multiply(&self, other: &Complex) -> Complex {
        Complex {
            real: self.real * other.real - self.imag * other.imag,
            imag: self.real * other.imag + self.imag * other.real,
        }
    }

    pub fn divide(&self, other: &Complex) -> Option<Complex> {
        let denom = other.real * other.real + other.imag * other.imag;
        if denom < constants::EPSILON {
            return None;
        }

        Some(Complex {
            real: (self.real * other.real + self.imag * other.imag) / denom,
            imag: (self.imag * other.real - self.real * other.imag) / denom,
        })
    }

    pub fn magnitude(&self) -> f32 {
        (self.real * self.real + self.imag * self.imag).sqrt()
    }

    pub fn argument(&self) -> f32 {
        self.imag.atan2(self.real)
    }

    pub fn exp(&self) -> Complex {
        let mag = self.real.exp();
        Complex {
            real: mag * self.imag.cos(),
            imag: mag * self.imag.sin(),
        }
    }

    pub fn log(&self) -> Option<Complex> {
        let mag = self.magnitude();
        if mag < constants::EPSILON {
            return None;
        }

        Some(Complex {
            real: mag.ln(),
            imag: self.argument(),
        })
    }
}

/// 2D Projektions-Engine
pub struct Polygon2DProjector {
    projection_type: PolygonProjection2DType,
    output_bounds: Option<Bounds2D>,
}

impl Polygon2DProjector {
    /// Erstellt einen neuen 2D Projektor
    pub fn new(projection_type: PolygonProjection2DType) -> Self {
        Self {
            projection_type,
            output_bounds: None,
        }
    }

    /// Setzt die Ausgabe-Bounds (für Clipping)
    pub fn with_bounds(mut self, bounds: Bounds2D) -> Self {
        self.output_bounds = Some(bounds);
        self
    }

    /// Projiziert ein Polygon
    pub fn project_polygon(&self, polygon: &Polygon) -> MathResult<Polygon> {
        let projected_vertices = self.project_points(polygon.vertices())?;

        if polygon.is_closed() {
            Polygon::closed(projected_vertices)
        } else {
            Polygon::new(projected_vertices)
        }
    }

    /// Projiziert eine Liste von Punkten
    pub fn project_points(&self, points: &[Vec2]) -> MathResult<Vec<Vec2>> {
        let mut projected_points = Vec::with_capacity(points.len());

        for &point in points {
            let projected = self.project_point(point)?;

            // Optional: Clipping an Ausgabe-Bounds
            if let Some(bounds) = &self.output_bounds {
                if bounds.contains_point(projected) {
                    projected_points.push(projected);
                }
            } else {
                projected_points.push(projected);
            }
        }

        Ok(projected_points)
    }

    /// Projiziert einen einzelnen Punkt
    pub fn project_point(&self, point: Vec2) -> MathResult<Vec2> {
        match self.projection_type {
            PolygonProjection2DType::OrthogonalLine {
                line_start,
                line_end,
            } => Ok(self.orthogonal_line_projection(point, line_start, line_end)),
            PolygonProjection2DType::Perspective { center, distance } => {
                self.perspective_projection(point, center, distance)
            }
            PolygonProjection2DType::CircularProjection { center, radius } => {
                Ok(self.circular_projection(point, center, radius))
            }
            PolygonProjection2DType::LogSpiral {
                center,
                growth_factor,
            } => Ok(self.log_spiral_projection(point, center, growth_factor)),
            PolygonProjection2DType::Conformal { transform_type } => {
                self.conformal_projection(point, transform_type)
            }
        }
    }

    // === Projektions-Implementierungen ===

    fn orthogonal_line_projection(&self, point: Vec2, line_start: Vec2, line_end: Vec2) -> Vec2 {
        let line_vec = line_end - line_start;
        let point_vec = point - line_start;

        let line_length_sq = line_vec.length_squared();
        if line_length_sq < constants::EPSILON {
            return line_start; // Degenerate line
        }

        let projection_factor = point_vec.dot(line_vec) / line_length_sq;
        line_start + line_vec * projection_factor
    }

    fn perspective_projection(&self, point: Vec2, center: Vec2, distance: f32) -> MathResult<Vec2> {
        if distance < constants::EPSILON {
            return Err(MathError::InvalidConfiguration {
                message: "Perspective distance must be positive".to_string(),
            });
        }

        let relative_point = point - center;

        // Perspektivische Projektion: projiziere auf Ebene bei z = distance
        // Annahme: Eingabe-Punkte sind bei z = 0, Projektor bei (center.x, center.y, -distance)
        let perspective_factor = distance / (distance + 0.0); // Vereinfacht für 2D

        Ok(center + relative_point * perspective_factor)
    }

    fn circular_projection(&self, point: Vec2, center: Vec2, radius: f32) -> Vec2 {
        let relative_point = point - center;
        let distance = relative_point.length();

        if distance < constants::EPSILON {
            return center + Vec2::new(radius, 0.0); // Fallback
        }

        // Projiziere auf Kreis
        center + relative_point * (radius / distance)
    }

    fn log_spiral_projection(&self, point: Vec2, center: Vec2, growth_factor: f32) -> Vec2 {
        let relative_point = point - center;
        let distance = relative_point.length();
        let angle = relative_point.y.atan2(relative_point.x);

        if distance < constants::EPSILON {
            return center;
        }

        // Logarithmische Spirale: r = a * e^(b * θ)
        let new_radius = distance * (growth_factor * angle).exp();

        center + Vec2::new(new_radius * angle.cos(), new_radius * angle.sin())
    }

    fn conformal_projection(&self, point: Vec2, transform_type: ConformalType) -> MathResult<Vec2> {
        let z = Complex::from_point(point);

        let transformed = match transform_type {
            ConformalType::Mobius { a, b, c, d } => self.mobius_transform(z, a, b, c, d)?,
            ConformalType::Joukowsky => self.joukowsky_transform(z),
            ConformalType::Logarithmic => z.log().ok_or_else(|| MathError::GeometricFailure {
                operation: "Logarithmic transformation at origin".to_string(),
            })?,
            ConformalType::Exponential => z.exp(),
        };

        Ok(transformed.to_point())
    }

    fn mobius_transform(
        &self,
        z: Complex,
        a: Complex,
        b: Complex,
        c: Complex,
        d: Complex,
    ) -> MathResult<Complex> {
        let numerator = a.multiply(&z).real + b.real;
        let denominator_z = c.multiply(&z);
        let denominator = Complex::new(denominator_z.real + d.real, denominator_z.imag + d.imag);

        let numerator_complex = Complex::new(numerator, a.multiply(&z).imag + b.imag);

        numerator_complex
            .divide(&denominator)
            .ok_or_else(|| MathError::GeometricFailure {
                operation: "Division by zero in Möbius transformation".to_string(),
            })
    }

    fn joukowsky_transform(&self, z: Complex) -> Complex {
        // Joukowsky: f(z) = (z + 1/z) / 2
        let one = Complex::new(1.0, 0.0);
        let z_inv = one.divide(&z).unwrap_or(Complex::new(0.0, 0.0));

        Complex::new((z.real + z_inv.real) * 0.5, (z.imag + z_inv.imag) * 0.5)
    }
}

/// 3D zu 2D Projektions-Wrapper
pub struct Polygon3DProjector {
    sphere_projector: SphereProjector,
}

impl Polygon3DProjector {
    /// Erstellt einen 3D zu 2D Projektor
    pub fn new(projection_type: SphereProjectionType, sphere_radius: f32) -> MathResult<Self> {
        Ok(Self {
            sphere_projector: SphereProjector::new(projection_type, sphere_radius)?,
        })
    }

    /// Projiziert 3D-Punkte auf 2D-Ebene
    pub fn project_3d_points(&self, points: &[Vec3]) -> MathResult<Vec<Vec2>> {
        self.sphere_projector.project_points(points)
    }

    /// Projiziert 3D-Polygon auf 2D-Ebene  
    pub fn project_3d_polygon_to_2d(&self, points_3d: &[Vec3]) -> MathResult<Polygon> {
        let points_2d = self.project_3d_points(points_3d)?;

        // Bestimme ob das Polygon geschlossen werden sollte
        let should_close = points_3d.len() > 2
            && points_3d
                .first()
                .map(|p| p.distance(*points_3d.last().unwrap()))
                < Some(0.01);

        if should_close {
            Polygon::closed(points_2d)
        } else {
            Polygon::new(points_2d)
        }
    }
}

/// Batch-Projektions-Utilities
pub struct ProjectionBatch;

impl ProjectionBatch {
    /// Projiziert mehrere Polygone mit derselben Projektion
    pub fn project_polygons_2d(
        polygons: &[Polygon],
        projector: &Polygon2DProjector,
    ) -> MathResult<Vec<Polygon>> {
        polygons
            .iter()
            .map(|polygon| projector.project_polygon(polygon))
            .collect()
    }

    /// Adaptive Projektion basierend auf Polygon-Eigenschaften
    pub fn adaptive_projection(polygon: &Polygon) -> MathResult<Polygon> {
        use crate::math::geometry::polygon::PolygonProperties;

        let area = polygon.area();
        let perimeter = polygon.perimeter();

        // Wähle Projektion basierend auf Eigenschaften
        let projection_type = if area < 100.0 {
            // Kleine Polygone: Perspektivische Projektion
            PolygonProjection2DType::Perspective {
                center: polygon.centroid().unwrap_or(Vec2::ZERO),
                distance: 50.0,
            }
        } else if perimeter / area > 10.0 {
            // Langgestreckte Polygone: Lineare Projektion
            if let Some(bounds) = polygon.bounds() {
                let center = (bounds.min + bounds.max) * 0.5;
                PolygonProjection2DType::OrthogonalLine {
                    line_start: Vec2::new(bounds.min.x, center.y),
                    line_end: Vec2::new(bounds.max.x, center.y),
                }
            } else {
                PolygonProjection2DType::Perspective {
                    center: Vec2::ZERO,
                    distance: 100.0,
                }
            }
        } else {
            // Standard-Polygone: Kreisprojektion
            PolygonProjection2DType::CircularProjection {
                center: polygon.centroid().unwrap_or(Vec2::ZERO),
                radius: (area / PI).sqrt(),
            }
        };

        let projector = Polygon2DProjector::new(projection_type);
        projector.project_polygon(polygon)
    }

    /// Multi-Scale Projektion für verschiedene Detail-Level
    pub fn multi_scale_projection(polygon: &Polygon, scales: &[f32]) -> MathResult<Vec<Polygon>> {
        let center = polygon.centroid().unwrap_or(Vec2::ZERO);
        let mut results = Vec::with_capacity(scales.len());

        for &scale in scales {
            let projection_type = PolygonProjection2DType::Perspective {
                center,
                distance: 100.0 * scale,
            };

            let projector = Polygon2DProjector::new(projection_type);
            results.push(projector.project_polygon(polygon)?);
        }

        Ok(results)
    }
}

/// Spezielle Projektions-Effects
pub struct ProjectionEffects;

impl ProjectionEffects {
    /// Fisheye-Effekt
    pub fn fisheye(polygon: &Polygon, center: Vec2, strength: f32) -> MathResult<Polygon> {
        let vertices = polygon.vertices();
        let mut projected_vertices = Vec::with_capacity(vertices.len());

        for &vertex in vertices {
            let relative = vertex - center;
            let distance = relative.length();

            if distance < constants::EPSILON {
                projected_vertices.push(vertex);
                continue;
            }

            // Fisheye-Transformation
            let angle = distance * strength;
            let fisheye_distance = angle.tan() / strength;
            let factor = fisheye_distance / distance;

            projected_vertices.push(center + relative * factor);
        }

        if polygon.is_closed() {
            Polygon::closed(projected_vertices)
        } else {
            Polygon::new(projected_vertices)
        }
    }

    /// Tunnel-Effekt (inverse radiale Projektion)
    pub fn tunnel_effect(
        polygon: &Polygon,
        center: Vec2,
        tunnel_radius: f32,
    ) -> MathResult<Polygon> {
        let vertices = polygon.vertices();
        let mut projected_vertices = Vec::with_capacity(vertices.len());

        for &vertex in vertices {
            let relative = vertex - center;
            let distance = relative.length();

            if distance < constants::EPSILON {
                projected_vertices.push(center + Vec2::new(tunnel_radius, 0.0));
                continue;
            }

            // Tunnel-Transformation: r' = tunnel_radius / r
            let new_distance = tunnel_radius * tunnel_radius / distance;
            let factor = new_distance / distance;

            projected_vertices.push(center + relative * factor);
        }

        if polygon.is_closed() {
            Polygon::closed(projected_vertices)
        } else {
            Polygon::new(projected_vertices)
        }
    }

    /// Spiral-Verzerrung
    pub fn spiral_distortion(
        polygon: &Polygon,
        center: Vec2,
        twist_factor: f32,
    ) -> MathResult<Polygon> {
        let vertices = polygon.vertices();
        let mut projected_vertices = Vec::with_capacity(vertices.len());

        for &vertex in vertices {
            let relative = vertex - center;
            let distance = relative.length();
            let angle = relative.y.atan2(relative.x);

            // Spiral-Transformation: θ' = θ + twist_factor * r
            let new_angle = angle + twist_factor * distance;

            projected_vertices
                .push(center + Vec2::new(distance * new_angle.cos(), distance * new_angle.sin()));
        }

        if polygon.is_closed() {
            Polygon::closed(projected_vertices)
        } else {
            Polygon::new(projected_vertices)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::geometry::polygon::Polygon;

    #[test]
    fn test_orthogonal_projection() {
        let square = Polygon::closed(vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(1.0, 0.0),
            Vec2::new(1.0, 1.0),
            Vec2::new(0.0, 1.0),
        ])
        .unwrap();

        let projector = Polygon2DProjector::new(PolygonProjection2DType::OrthogonalLine {
            line_start: Vec2::new(0.0, 0.0),
            line_end: Vec2::new(1.0, 0.0),
        });

        let projected = projector.project_polygon(&square).unwrap();

        // All points should be on the x-axis
        for vertex in projected.vertices() {
            assert!(vertex.y.abs() < constants::EPSILON);
        }
    }

    #[test]
    fn test_circular_projection() {
        let triangle = Polygon::closed(vec![
            Vec2::new(1.0, 0.0),
            Vec2::new(0.0, 1.0),
            Vec2::new(-1.0, 0.0),
        ])
        .unwrap();

        let projector = Polygon2DProjector::new(PolygonProjection2DType::CircularProjection {
            center: Vec2::ZERO,
            radius: 2.0,
        });

        let projected = projector.project_polygon(&triangle).unwrap();

        // All points should be at distance 2.0 from origin
        for vertex in projected.vertices() {
            let distance = vertex.length();
            assert!(comparison::nearly_equal(distance, 2.0));
        }
    }

    #[test]
    fn test_conformal_joukowsky() {
        let circle = Polygon::closed(vec![
            Vec2::new(1.0, 0.0),
            Vec2::new(0.0, 1.0),
            Vec2::new(-1.0, 0.0),
            Vec2::new(0.0, -1.0),
        ])
        .unwrap();

        let projector = Polygon2DProjector::new(PolygonProjection2DType::Conformal {
            transform_type: ConformalType::Joukowsky,
        });

        let projected = projector.project_polygon(&circle).unwrap();

        // Joukowsky transform should produce valid polygon
        assert_eq!(projected.len(), circle.len());
        assert!(projected.is_closed());
    }

    #[test]
    fn test_projection_effects() {
        let square = Polygon::closed(vec![
            Vec2::new(-1.0, -1.0),
            Vec2::new(1.0, -1.0),
            Vec2::new(1.0, 1.0),
            Vec2::new(-1.0, 1.0),
        ])
        .unwrap();

        let fisheye = ProjectionEffects::fisheye(&square, Vec2::ZERO, 0.5).unwrap();
        let tunnel = ProjectionEffects::tunnel_effect(&square, Vec2::ZERO, 2.0).unwrap();
        let spiral = ProjectionEffects::spiral_distortion(&square, Vec2::ZERO, 0.1).unwrap();

        // All effects should produce valid polygons
        assert!(fisheye.is_closed());
        assert!(tunnel.is_closed());
        assert!(spiral.is_closed());

        // Each effect should change the original shape
        assert_ne!(fisheye.vertices(), square.vertices());
        assert_ne!(tunnel.vertices(), square.vertices());
        assert_ne!(spiral.vertices(), square.vertices());
    }

    #[test]
    fn test_complex_number_operations() {
        let z1 = Complex::new(1.0, 2.0);
        let z2 = Complex::new(3.0, 4.0);

        let product = z1.multiply(&z2);
        assert!(comparison::nearly_equal(product.real, -5.0));
        assert!(comparison::nearly_equal(product.imag, 10.0));

        let quotient = z1.divide(&z2).unwrap();
        assert!(comparison::nearly_equal(quotient.real, 0.44));
        assert!(comparison::nearly_equal(quotient.imag, 0.08));

        let magnitude = z1.magnitude();
        assert!(comparison::nearly_equal(magnitude, 5.0_f32.sqrt()));
    }
}
