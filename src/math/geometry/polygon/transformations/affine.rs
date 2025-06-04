use super::super::Polygon;
use crate::math::{error::MathResult, types::*};

/// Affine Transformations-Matrix (3x3 für 2D)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AffineTransform {
    // Matrix in der Form: [a c tx]
    //                    [b d ty]
    //                    [0 0  1]
    pub a: f32,  // x-Skalierung
    pub b: f32,  // y-Scherung in x
    pub c: f32,  // x-Scherung in y
    pub d: f32,  // y-Skalierung
    pub tx: f32, // x-Translation
    pub ty: f32, // y-Translation
}

impl AffineTransform {
    /// Identitäts-Transformation
    pub fn identity() -> Self {
        Self {
            a: 1.0,
            b: 0.0,
            c: 0.0,
            d: 1.0,
            tx: 0.0,
            ty: 0.0,
        }
    }

    /// Translation
    pub fn translation(tx: f32, ty: f32) -> Self {
        Self {
            a: 1.0,
            b: 0.0,
            c: 0.0,
            d: 1.0,
            tx,
            ty,
        }
    }

    /// Skalierung
    pub fn scale(sx: f32, sy: f32) -> Self {
        Self {
            a: sx,
            b: 0.0,
            c: 0.0,
            d: sy,
            tx: 0.0,
            ty: 0.0,
        }
    }

    /// Uniform Skalierung
    pub fn uniform_scale(scale: f32) -> Self {
        Self::scale(scale, scale)
    }

    /// Rotation um den Ursprung
    pub fn rotation(angle_rad: f32) -> Self {
        let cos_a = angle_rad.cos();
        let sin_a = angle_rad.sin();

        Self {
            a: cos_a,
            b: sin_a,
            c: -sin_a,
            d: cos_a,
            tx: 0.0,
            ty: 0.0,
        }
    }

    /// Rotation um einen Punkt
    pub fn rotation_around(angle_rad: f32, center: Point2D) -> Self {
        let rotation = Self::rotation(angle_rad);
        let translate_back = Self::translation(center.x, center.y);
        let translate_to_origin = Self::translation(-center.x, -center.y);

        translate_back
            .compose(&rotation)
            .compose(&translate_to_origin)
    }

    /// Scherung
    pub fn shear(shx: f32, shy: f32) -> Self {
        Self {
            a: 1.0,
            b: shy,
            c: shx,
            d: 1.0,
            tx: 0.0,
            ty: 0.0,
        }
    }

    /// Spiegelung an der X-Achse
    pub fn reflect_x() -> Self {
        Self::scale(1.0, -1.0)
    }

    /// Spiegelung an der Y-Achse
    pub fn reflect_y() -> Self {
        Self::scale(-1.0, 1.0)
    }

    /// Spiegelung an beliebiger Linie durch den Ursprung
    pub fn reflect_line(angle_rad: f32) -> Self {
        let cos_2a = (2.0 * angle_rad).cos();
        let sin_2a = (2.0 * angle_rad).sin();

        Self {
            a: cos_2a,
            b: sin_2a,
            c: sin_2a,
            d: -cos_2a,
            tx: 0.0,
            ty: 0.0,
        }
    }

    /// Transformations-Komposition
    pub fn compose(&self, other: &AffineTransform) -> Self {
        Self {
            a: self.a * other.a + self.c * other.b,
            b: self.b * other.a + self.d * other.b,
            c: self.a * other.c + self.c * other.d,
            d: self.b * other.c + self.d * other.d,
            tx: self.a * other.tx + self.c * other.ty + self.tx,
            ty: self.b * other.tx + self.d * other.ty + self.ty,
        }
    }

    /// Inverse Transformation
    pub fn inverse(&self) -> Option<Self> {
        let det = self.a * self.d - self.b * self.c;

        if det.abs() < 1e-10 {
            return None; // Nicht invertierbar
        }

        let inv_det = 1.0 / det;

        Some(Self {
            a: self.d * inv_det,
            b: -self.b * inv_det,
            c: -self.c * inv_det,
            d: self.a * inv_det,
            tx: (self.c * self.ty - self.d * self.tx) * inv_det,
            ty: (self.b * self.tx - self.a * self.ty) * inv_det,
        })
    }

    /// Transformiert einen Punkt
    pub fn transform_point(&self, point: Point2D) -> Point2D {
        Point2D::new(
            self.a * point.x + self.c * point.y + self.tx,
            self.b * point.x + self.d * point.y + self.ty,
        )
    }

    /// Transformiert nur die Richtung (ohne Translation)
    pub fn transform_direction(&self, direction: Point2D) -> Point2D {
        Point2D::new(
            self.a * direction.x + self.c * direction.y,
            self.b * direction.x + self.d * direction.y,
        )
    }
}

/// Trait für affine Transformationen
pub trait AffineTransformable {
    fn transform(&self, transform: &AffineTransform) -> MathResult<Self>
    where
        Self: Sized;

    fn transform_mut(&mut self, transform: &AffineTransform) -> MathResult<()>;

    // Convenience-Methoden
    fn translate(&self, tx: f32, ty: f32) -> MathResult<Self>
    where
        Self: Sized,
    {
        self.transform(&AffineTransform::translation(tx, ty))
    }

    fn scale(&self, sx: f32, sy: f32) -> MathResult<Self>
    where
        Self: Sized,
    {
        self.transform(&AffineTransform::scale(sx, sy))
    }

    fn rotate(&self, angle_rad: f32) -> MathResult<Self>
    where
        Self: Sized,
    {
        self.transform(&AffineTransform::rotation(angle_rad))
    }

    fn rotate_around(&self, angle_rad: f32, center: Point2D) -> MathResult<Self>
    where
        Self: Sized,
    {
        self.transform(&AffineTransform::rotation_around(angle_rad, center))
    }
}

impl AffineTransformable for Polygon {
    fn transform(&self, transform: &AffineTransform) -> MathResult<Self> {
        let transformed_vertices: Vec<Point2D> = self
            .vertices()
            .iter()
            .map(|&vertex| transform.transform_point(vertex))
            .collect();

        if self.is_closed() {
            Polygon::closed(transformed_vertices)
        } else {
            Polygon::new(transformed_vertices)
        }
    }

    fn transform_mut(&mut self, transform: &AffineTransform) -> MathResult<()> {
        for vertex in self.vertices_mut() {
            *vertex = transform.transform_point(*vertex);
        }
        Ok(())
    }
}

/// Builder für komplexe Transformationen
pub struct TransformBuilder {
    transform: AffineTransform,
}

impl TransformBuilder {
    pub fn new() -> Self {
        Self {
            transform: AffineTransform::identity(),
        }
    }

    pub fn translate(mut self, tx: f32, ty: f32) -> Self {
        self.transform = self
            .transform
            .compose(&AffineTransform::translation(tx, ty));
        self
    }

    pub fn scale(mut self, sx: f32, sy: f32) -> Self {
        self.transform = self.transform.compose(&AffineTransform::scale(sx, sy));
        self
    }

    pub fn rotate(mut self, angle_rad: f32) -> Self {
        self.transform = self
            .transform
            .compose(&AffineTransform::rotation(angle_rad));
        self
    }

    pub fn rotate_around(mut self, angle_rad: f32, center: Point2D) -> Self {
        self.transform = self
            .transform
            .compose(&AffineTransform::rotation_around(angle_rad, center));
        self
    }

    pub fn shear(mut self, shx: f32, shy: f32) -> Self {
        self.transform = self.transform.compose(&AffineTransform::shear(shx, shy));
        self
    }

    pub fn reflect_x(mut self) -> Self {
        self.transform = self.transform.compose(&AffineTransform::reflect_x());
        self
    }

    pub fn reflect_y(mut self) -> Self {
        self.transform = self.transform.compose(&AffineTransform::reflect_y());
        self
    }

    pub fn build(self) -> AffineTransform {
        self.transform
    }

    pub fn apply_to<T: AffineTransformable>(&self, target: &T) -> MathResult<T> {
        target.transform(&self.transform)
    }
}

impl Default for TransformBuilder {
    fn default() -> Self {
        Self::new()
    }
}
