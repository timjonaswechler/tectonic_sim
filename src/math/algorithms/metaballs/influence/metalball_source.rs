// src/math/algorithms/metaballs/influence/metalball_source.rs

use crate::math::algorithms::metaballs::influence::FieldInfluence;
use crate::math::types::Bounds2D;
use crate::math::utils::constants;
use bevy::math::Vec2;

/// Repräsentiert einen einzelnen klassischen Metaball (oder Meta-Kreis in 2D) als Einflussquelle für das Feld.
/// Dieser Typ von Quelle zeichnet sich durch eine Position, einen Radius und eine Stärke aus,
/// deren Einfluss mit der Distanz gemäß einem wählbaren Falloff-Typ abnimmt.
#[derive(Debug, Clone)]
pub struct MetaballSource {
    pub position: Vec2,
    /// Effektiver Radius, bis zu dem der Einfluss primär wirkt.
    /// Je nach Falloff-Typ kann der Einfluss auch darüber hinausgehen oder davor Null werden.
    pub radius: f32,
    /// Maximale Stärke des Einflusses (typischerweise im Zentrum oder nahe dem Zentrum).
    pub strength: f32,
    pub falloff_type: MetaballFalloff,
}

/// Definiert, wie der Einfluss eines `MetaballSource` mit der Distanz abnimmt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MetaballFalloff {
    /// Algebraischer Falloff, oft für den klassischen Metaball-Look verwendet.
    /// Formel: `strength / ( (distance / (radius * scale_factor))^2 + 1.0 )`
    /// `scale_factor` (z.B. 0.5) bestimmt, bei welcher relativen Distanz zum Radius die Stärke halbiert wird.
    Algebraic {
        /// Skalierungsfaktor für den Radius in der Formel.
        /// Ein Wert von 0.5 bedeutet, dass bei `distance = 0.5 * radius` die Stärke ca. `strength / 2` ist.
        radius_scale_factor: f32,
    },
    /// Exponentieller Falloff: `strength * exp(-steepness * normalized_distance)`.
    Exponential {
        /// Steuert, wie schnell der Einfluss abfällt. Größere Werte = schnellerer Abfall.
        steepness: f32,
    },
    /// Gausscher Falloff: `strength * exp(-(normalized_distance / (2.0 * variance_factor))^2)`.
    /// Alternativ: `strength * exp(-distance_sq / (2.0 * sigma_sq))`.
    Gaussian {
        /// Beeinflusst die "Breite" der Gausskurve. Ein kleinerer Wert führt zu einer schmaleren Kurve.
        /// `variance_factor` wird verwendet, um Sigma relativ zum Radius zu bestimmen: `sigma = radius * variance_factor`.
        variance_factor: f32,
    },
    /// Linearer Falloff von `strength` (bei Distanz 0) zu `0` (bei `distance = radius`).
    Linear,
    /// Weicher Übergang (kubische Hermite-Interpolation) von `strength` zu `0` über den `radius`.
    SmoothStep,
    /// Wyvill-Falloff-Familie, bekannt für "blobby" Look.
    /// Formel: `strength * (1 - (distance/radius)^p1)^p2`. Hier mit p1=2, p2=3.
    WyvillCubic, // (1 - (d/R)^2)^3
}

impl MetaballSource {
    /// Erstellt einen neuen `MetaballSource`.
    /// Beachte die Parameterreihenfolge: `position`, `radius`, `strength`.
    pub fn new(position: Vec2, radius: f32, strength: f32) -> Self {
        Self {
            position,
            radius: radius.max(constants::EPSILON), // Radius muss positiv sein
            strength,
            // Standard-Falloff: Ein algebraischer Falloff, der bei halbem Radius etwa halbe Stärke hat.
            falloff_type: MetaballFalloff::Algebraic {
                radius_scale_factor: 0.5,
            },
        }
    }

    /// Setzt den Falloff-Typ für diesen `MetaballSource`.
    pub fn with_falloff(mut self, falloff_type: MetaballFalloff) -> Self {
        self.falloff_type = falloff_type;
        self
    }

    /// Berechnet den Einflusswert basierend auf dem Abstand und dem Falloff-Typ.
    fn calculate_influence_value(&self, point_distance_sq: f32) -> f32 {
        if point_distance_sq > self.radius * self.radius + constants::EPSILON
            && matches!(
                self.falloff_type,
                MetaballFalloff::Linear
                    | MetaballFalloff::SmoothStep
                    | MetaballFalloff::WyvillCubic
            )
        {
            return 0.0; // Außerhalb des definierten Einflussradius für diese Typen
        }

        let distance = point_distance_sq.sqrt();

        let normalized_distance = if self.radius > constants::EPSILON {
            (distance / self.radius).clamp(0.0, 1.0)
        } else {
            if distance < constants::EPSILON {
                0.0
            } else {
                return 0.0;
            }
        };

        match self.falloff_type {
            MetaballFalloff::Algebraic {
                radius_scale_factor,
            } => {
                let scaled_denominator_dist =
                    self.radius * radius_scale_factor.max(constants::EPSILON);
                if scaled_denominator_dist < constants::EPSILON {
                    return if distance < constants::EPSILON {
                        self.strength
                    } else {
                        0.0
                    };
                }
                let term_sq = (distance / scaled_denominator_dist).powi(2);
                self.strength / (term_sq + 1.0)
            }
            MetaballFalloff::Exponential { steepness } => {
                self.strength * (-normalized_distance * steepness.max(0.01)).exp()
            }
            MetaballFalloff::Gaussian { variance_factor } => {
                if variance_factor < constants::EPSILON {
                    return if normalized_distance < constants::EPSILON {
                        self.strength
                    } else {
                        0.0
                    };
                }
                let exponent = -normalized_distance.powi(2) / (2.0 * variance_factor.powi(2));
                self.strength * exponent.exp()
            }
            MetaballFalloff::Linear => self.strength * (1.0 - normalized_distance),
            MetaballFalloff::SmoothStep => {
                let t = 1.0 - normalized_distance;
                self.strength * t * t * (3.0 - 2.0 * t)
            }
            MetaballFalloff::WyvillCubic => {
                if normalized_distance >= 1.0 {
                    return 0.0;
                }
                let term = 1.0 - normalized_distance.powi(2);
                self.strength * term.powi(3)
            }
        }
    }
}

impl FieldInfluence for MetaballSource {
    fn influence_at(&self, point: Vec2) -> f32 {
        let distance_sq = self.position.distance_squared(point);
        self.calculate_influence_value(distance_sq)
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        let min_coord = self.position - Vec2::splat(self.radius);
        let max_coord = self.position + Vec2::splat(self.radius);
        Some(
            Bounds2D::new(min_coord, max_coord)
                .expect("MetaballSource bounds should be valid due to positive radius"),
        )
    }

    fn max_influence_distance(&self) -> Option<f32> {
        Some(self.radius)
    }

    fn influence_type(&self) -> &'static str {
        "MetaballSource"
    }
}
