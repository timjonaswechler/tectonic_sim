// src/math/geometry/polygon/transformations/distortion.rs

use super::super::Polygon;
use crate::math::error::MathResult;
use bevy::math::Vec2;

/// Trait für benutzerdefinierte Verzerrungsfunktionen
pub trait DistortionFn: Send + Sync + std::fmt::Debug {
    fn transform(&self, point: Vec2, center: Vec2) -> Vec2;
    fn clone_box(&self) -> Box<dyn DistortionFn + Send + Sync>;
}

impl Clone for Box<dyn DistortionFn + Send + Sync> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

/// Verschiedene Verzerrungstypen
#[derive(Debug, Clone)]
pub enum DistortionType {
    /// Linsenverzerrung: factor * r²
    Barrel(f32),
    /// Kissenverzerrung: -factor * r²
    Pincushion(f32),
    /// Wellenverzerrung
    Wave {
        amplitude: f32,
        frequency: f32,
        phase: f32,
    },
    /// Spiral-Verzerrung
    Spiral { strength: f32, center: Vec2 },
    /// Benutzerdefinierte Funktion
    Custom(Box<dyn DistortionFn + Send + Sync + 'static>),
}

/// Trait für Verzerrungsoperationen
pub trait Distortion {
    /// Verzerrt einen einzelnen Punkt
    fn distort_point(&self, point: Vec2, center: Vec2, params: DistortionType) -> Vec2;

    /// Verzerrt das gesamte Polygon
    fn distort(&self, center: Vec2, params: DistortionType) -> MathResult<Polygon>;

    /// Verzerrt das Polygon in-place
    fn distort_mut(&mut self, center: Vec2, params: DistortionType) -> MathResult<()>;
}

impl Distortion for Polygon {
    fn distort_point(&self, point: Vec2, center: Vec2, params: DistortionType) -> Vec2 {
        match params {
            DistortionType::Barrel(strength) | DistortionType::Pincushion(strength) => {
                let dx = point.x - center.x;
                let dy = point.y - center.y;
                let r2 = dx * dx + dy * dy;

                let factor = match params {
                    DistortionType::Barrel(_) => 1.0 + strength * r2,
                    DistortionType::Pincushion(_) => 1.0 - strength * r2,
                    _ => unreachable!(),
                };

                Vec2::new(center.x + dx * factor, center.y + dy * factor)
            }

            DistortionType::Wave {
                amplitude,
                frequency,
                phase,
            } => {
                let dx = point.x - center.x;
                let dy = point.y - center.y;
                let distance = (dx * dx + dy * dy).sqrt();

                let wave_offset = amplitude * (frequency * distance + phase).sin();
                let angle = dy.atan2(dx);

                Vec2::new(
                    point.x + wave_offset * angle.cos(),
                    point.y + wave_offset * angle.sin(),
                )
            }

            DistortionType::Spiral {
                strength,
                center: spiral_center,
            } => {
                let dx = point.x - spiral_center.x;
                let dy = point.y - spiral_center.y;
                let r = (dx * dx + dy * dy).sqrt();
                let angle = dy.atan2(dx) + strength * r;

                Vec2::new(
                    spiral_center.x + r * angle.cos(),
                    spiral_center.y + r * angle.sin(),
                )
            }
            DistortionType::Custom(func) => func.transform(point, center),
        }
    }

    fn distort(&self, center: Vec2, params: DistortionType) -> MathResult<Polygon> {
        let distorted_vertices: Vec<Vec2> = self
            .vertices()
            .iter()
            .map(|&point| self.distort_point(point, center, params.clone()))
            .collect();

        if self.is_closed() {
            Polygon::closed(distorted_vertices)
        } else {
            Polygon::new(distorted_vertices)
        }
    }

    fn distort_mut(&mut self, center: Vec2, params: DistortionType) -> MathResult<()> {
        let new_points: Vec<Vec2> = self
            .vertices()
            .iter()
            .map(|&point| self.distort_point(point, center, params.clone()))
            .collect();

        for (vertex, new_point) in self.vertices_mut().iter_mut().zip(new_points.into_iter()) {
            *vertex = new_point;
        }
        Ok(())
    }
}

/// Erweiterte Verzerrungsfunktionen
pub struct DistortionEffects;

impl DistortionEffects {
    /// Sinusförmige Verzerrung entlang einer Achse
    pub fn sine_wave(amplitude: f32, frequency: f32, axis: SineAxis) -> DistortionType {
        #[derive(Debug)]
        struct SineWaveDistortion {
            amplitude: f32,
            frequency: f32,
            axis: SineAxis,
        }

        impl DistortionFn for SineWaveDistortion {
            fn transform(&self, point: Vec2, _center: Vec2) -> Vec2 {
                match self.axis {
                    SineAxis::X => Vec2::new(
                        point.x,
                        point.y + self.amplitude * (self.frequency * point.x).sin(),
                    ),
                    SineAxis::Y => Vec2::new(
                        point.x + self.amplitude * (self.frequency * point.y).sin(),
                        point.y,
                    ),
                }
            }

            fn clone_box(&self) -> Box<dyn DistortionFn + Send + Sync> {
                Box::new(Self { ..*self })
            }
        }

        DistortionType::Custom(Box::new(SineWaveDistortion {
            amplitude,
            frequency,
            axis,
        }))
    }

    /// Turbulenz-ähnliche Verzerrung
    pub fn turbulence(strength: f32, scale: f32) -> DistortionType {
        #[derive(Debug)]
        struct TurbulenceDistortion {
            strength: f32,
            scale: f32,
        }

        impl DistortionFn for TurbulenceDistortion {
            fn transform(&self, point: Vec2, _center: Vec2) -> Vec2 {
                let noise_x = ((point.x * self.scale).sin() * (point.y * self.scale * 1.3).cos())
                    * self.strength;
                let noise_y = ((point.y * self.scale).sin() * (point.x * self.scale * 0.7).cos())
                    * self.strength;
                Vec2::new(point.x + noise_x, point.y + noise_y)
            }
            fn clone_box(&self) -> Box<dyn DistortionFn + Send + Sync> {
                Box::new(Self { ..*self })
            }
        }

        DistortionType::Custom(Box::new(TurbulenceDistortion { strength, scale }))
    }

    /// Magnetfeld-ähnliche Verzerrung
    pub fn magnetic_field(poles: &[(Vec2, f32)]) -> DistortionType {
        #[derive(Debug, Clone)]
        struct MagneticFieldDistortion {
            poles: Vec<(Vec2, f32)>,
        }

        impl DistortionFn for MagneticFieldDistortion {
            fn transform(&self, point: Vec2, _center: Vec2) -> Vec2 {
                let mut total_force = Vec2::ZERO;

                for (pole_pos, strength) in &self.poles {
                    let dx = point.x - pole_pos.x;
                    let dy = point.y - pole_pos.y;
                    let dist_sq = dx * dx + dy * dy;

                    if dist_sq > 1e-6 {
                        let force_magnitude = *strength / dist_sq;
                        total_force.x += dx * force_magnitude;
                        total_force.y += dy * force_magnitude;
                    }
                }

                point + total_force
            }
            fn clone_box(&self) -> Box<dyn DistortionFn + Send + Sync> {
                Box::new(self.clone())
            }
        }

        DistortionType::Custom(Box::new(MagneticFieldDistortion {
            poles: poles.to_vec(),
        }))
    }
}

#[derive(Debug, Clone, Copy)]
pub enum SineAxis {
    X,
    Y,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_barrel_distortion() {
        let square = Polygon::closed(vec![
            Vec2::new(-1.0, -1.0),
            Vec2::new(1.0, -1.0),
            Vec2::new(1.0, 1.0),
            Vec2::new(-1.0, 1.0),
        ])
        .unwrap();

        let distorted = square
            .distort(Vec2::ZERO, DistortionType::Barrel(0.1))
            .unwrap();

        // Barrel-Verzerrung sollte Punkte nach außen verschieben
        assert!(
            distorted.vertices()[0].distance(Vec2::ZERO)
                > square.vertices()[0].distance(Vec2::ZERO)
        );
    }

    #[test]
    fn test_wave_distortion() {
        let line = Polygon::new(vec![
            Vec2::new(0.0, 0.0),
            Vec2::new(1.0, 0.0),
            Vec2::new(2.0, 0.0),
        ])
        .unwrap();

        let distorted = line
            .distort(
                Vec2::ZERO,
                DistortionType::Wave {
                    amplitude: 0.5,
                    frequency: 1.0,
                    phase: 0.0,
                },
            )
            .unwrap();

        // Wellenverzerrung sollte Y-Koordinaten ändern
        assert_ne!(distorted.vertices()[1].y, 0.0);
    }
}
