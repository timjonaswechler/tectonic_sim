// src/math/geometry/metalballs/influence/metalball_source.rs
// (oder wie auch immer du die Datei nach der Refaktorierung nennst)

use crate::math::geometry::metalballs::influence::FieldInfluence; // Der Trait
use crate::math::types::Bounds2D;
use crate::math::utils::constants;
use bevy::math::Vec2;

/// Repräsentiert einen einzelnen Metaball (oder Meta-Kreis in 2D) als Einflussquelle für das Metaballs-System.
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
        // Die meisten Falloffs sind so definiert, dass sie am oder vor dem Radius Null werden.
        // Diese Prüfung stellt sicher, dass Punkte weit außerhalb keinen Einfluss haben.
        // Für Falloffs, die theoretisch unendlich reichen (wie manche algebraische oder Gaußsche),
        // dient der `radius` als primärer Bereich von Interesse.
        if point_distance_sq >= self.radius * self.radius &&
           // Ausnahmen für Falloffs, die exakt bei `radius` Null werden (oder davor)
           !matches!(self.falloff_type, MetaballFalloff::Linear | MetaballFalloff::SmoothStep | MetaballFalloff::WyvillCubic)
        {
            // Für Falloffs, die über den Radius hinausgehen könnten, aber hier abgeschnitten werden sollen.
            // Wenn der Falloff selbst bei self.radius null wird, ist dieser Check redundant, aber nicht schädlich.
            // Bei Algebraic oder Gaussian könnte dies zu einem harten Cutoff führen.
            // Alternativ: Den Check entfernen und sich auf die Falloff-Formel verlassen, dass sie klein genug wird.
            // Für jetzt: Expliziter Cutoff bei `radius` für die meisten Typen.
        }

        // Strikter Check: Wenn Distanz >= Radius, ist der Einfluss 0 für die meisten gängigen Definitionen,
        // die innerhalb des Radius wirken.
        // Die Falloff-Funktionen selbst sollten dies auch berücksichtigen.
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

        // normalized_distance ist nur sinnvoll, wenn radius > EPSILON
        let normalized_distance = if self.radius > constants::EPSILON {
            (distance / self.radius).clamp(0.0, 1.0) // Clamp auf 1.0, da Falloffs meist bis Radius definiert sind
        } else {
            // Wenn Radius 0 ist, ist der Punkt entweder genau auf der Position (Distanz 0) oder nicht.
            if distance < constants::EPSILON {
                0.0
            } else {
                return 0.0;
            } // Einfluss nur am exakten Punkt
        };

        match self.falloff_type {
            MetaballFalloff::Algebraic {
                radius_scale_factor,
            } => {
                // Nenner wird nie 0, da distance_sq >= 0 und 1.0 addiert wird.
                // Skalierung der Distanz: (distance / (self.radius * radius_scale_factor))
                // Wenn radius_scale_factor = 0.5, und distance = 0.5 * radius, dann ist (dist_scaled)^2 = 1.
                // strength / (1+1) = strength / 2.
                let scaled_denominator_dist =
                    self.radius * radius_scale_factor.max(constants::EPSILON);
                if scaled_denominator_dist < constants::EPSILON {
                    // Vermeide Division durch Null wenn radius oder factor zu klein
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
                // normalized_distance ist hier [0,1]
                self.strength * (-normalized_distance * steepness.max(0.01)).exp()
            }
            MetaballFalloff::Gaussian { variance_factor } => {
                // sigma = radius * variance_factor.
                // normalized_distance = distance / radius.
                // -(distance^2) / (2 * sigma^2) = -( (normalized_distance * radius)^2 ) / (2 * (radius * variance_factor)^2 )
                // = -(normalized_distance^2) / (2 * variance_factor^2)
                if variance_factor < constants::EPSILON {
                    // Vermeide Division durch Null
                    return if normalized_distance < constants::EPSILON {
                        self.strength
                    } else {
                        0.0
                    };
                }
                let exponent = -normalized_distance.powi(2) / (2.0 * variance_factor.powi(2));
                self.strength * exponent.exp()
            }
            MetaballFalloff::Linear => {
                // normalized_distance ist [0,1]
                self.strength * (1.0 - normalized_distance) // .max(0.0) ist durch clamp oben schon abgedeckt
            }
            MetaballFalloff::SmoothStep => {
                // t geht von 1 (Mitte, normalized_distance=0) zu 0 (Rand, normalized_distance=1)
                let t = 1.0 - normalized_distance;
                self.strength * t * t * (3.0 - 2.0 * t)
            }
            MetaballFalloff::WyvillCubic => {
                // t_sq = (distance/radius)^2 = normalized_distance^2
                // Formel: strength * (1 - t_sq)^3
                // normalized_distance ist bereits auf [0,1] geklemmt.
                if normalized_distance >= 1.0 {
                    return 0.0;
                } // Exakt am Rand oder darüber ist 0
                let term = 1.0 - normalized_distance.powi(2);
                self.strength * term.powi(3)
            }
        }
    }
}

impl FieldInfluence for MetaballSource {
    fn influence_at(&self, point: Vec2) -> f32 {
        let distance_sq = self.position.distance_squared(point);
        // Die Logik, ob der Punkt außerhalb des Radius liegt und 0 zurückgeben soll,
        // ist nun primär in calculate_influence_value und den Falloff-Definitionen.
        self.calculate_influence_value(distance_sq)
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        let min_coord = self.position - Vec2::splat(self.radius);
        let max_coord = self.position + Vec2::splat(self.radius);
        // .expect ist hier okay, da der Konstruktor sicherstellt, dass radius >= EPSILON.
        Some(
            Bounds2D::new(min_coord, max_coord)
                .expect("MetaballSource bounds should be valid due to positive radius"),
        )
    }

    fn max_influence_distance(&self) -> Option<f32> {
        // Der nominelle Radius. Einige Falloffs (Gaussian, Algebraic) könnten theoretisch weiter reichen,
        // aber `radius` ist der primäre Parameter für den Einflussbereich.
        Some(self.radius)
    }

    fn influence_type(&self) -> &'static str {
        "MetaballSource"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::utils::constants::EPSILON; // Für Vergleiche

    const TEST_STRENGTH: f32 = 10.0;
    const TEST_RADIUS: f32 = 5.0;

    fn create_test_ball(falloff_type: MetaballFalloff) -> MetaballSource {
        MetaballSource {
            position: Vec2::ZERO,
            radius: TEST_RADIUS,
            strength: TEST_STRENGTH,
            falloff_type,
        }
    }

    #[test]
    fn test_metaball_at_center() {
        // Alle Falloffs sollten im Zentrum (Distanz 0) die volle Stärke haben.
        let falloffs = [
            MetaballFalloff::Algebraic {
                radius_scale_factor: 0.5,
            },
            MetaballFalloff::Exponential { steepness: 2.0 },
            MetaballFalloff::Gaussian {
                variance_factor: 0.5,
            },
            MetaballFalloff::Linear,
            MetaballFalloff::SmoothStep,
            MetaballFalloff::WyvillCubic,
        ];

        for falloff in falloffs.iter() {
            let ball = create_test_ball(*falloff);
            let influence = ball.influence_at(Vec2::ZERO);
            assert!(
                (influence - TEST_STRENGTH).abs() < EPSILON,
                "Center influence mismatch for {:?}: expected {}, got {}",
                falloff,
                TEST_STRENGTH,
                influence
            );
        }
    }

    #[test]
    fn test_metaball_at_radius() {
        // Linear, SmoothStep, Wyvill sollten am Radius (Distanz = TEST_RADIUS) genau 0 sein.
        let exact_zero_falloffs = [
            MetaballFalloff::Linear,
            MetaballFalloff::SmoothStep,
            MetaballFalloff::WyvillCubic,
        ];
        for falloff in exact_zero_falloffs.iter() {
            let ball = create_test_ball(*falloff);
            let influence = ball.influence_at(Vec2::new(TEST_RADIUS, 0.0));
            assert!(
                influence.abs() < EPSILON,
                "Radius influence mismatch for {:?}: expected 0, got {}",
                falloff,
                influence
            );
        }

        // Andere Falloffs (Algebraic, Exponential, Gaussian) sind am Radius nicht zwingend 0,
        // aber die `calculate_influence_value` schneidet sie evtl. ab, oder sie sind sehr klein.
        // Die `normalized_distance` wird auf 1.0 geklemmt.

        // Algebraic an R: strength / ( (R / (R*factor))^2 + 1 ) = strength / ( (1/factor)^2 + 1 )
        let ball_alg = create_test_ball(MetaballFalloff::Algebraic {
            radius_scale_factor: 0.5,
        });
        let expected_alg_at_r = TEST_STRENGTH / ((1.0 / 0.5_f32).powi(2) + 1.0); // 10 / (2^2+1) = 10 / 5 = 2.0
        let influence_alg_at_r = ball_alg.influence_at(Vec2::new(TEST_RADIUS, 0.0));
        assert!(
            (influence_alg_at_r - expected_alg_at_r).abs() < EPSILON,
            "Algebraic at Radius: expected {}, got {}",
            expected_alg_at_r,
            influence_alg_at_r
        );

        // Exponential an R: strength * exp(-steepness * 1.0)
        let steepness = 2.0;
        let ball_exp = create_test_ball(MetaballFalloff::Exponential { steepness });
        let expected_exp_at_r = TEST_STRENGTH * (-steepness).exp();
        let influence_exp_at_r = ball_exp.influence_at(Vec2::new(TEST_RADIUS, 0.0));
        assert!(
            (influence_exp_at_r - expected_exp_at_r).abs() < EPSILON,
            "Exponential at Radius: expected {}, got {}",
            expected_exp_at_r,
            influence_exp_at_r
        );

        // Gaussian an R: strength * exp( -(1.0)^2 / (2*var_factor^2) )
        let variance_factor = 0.5;
        let ball_gauss = create_test_ball(MetaballFalloff::Gaussian { variance_factor });
        let expected_gauss_at_r =
            TEST_STRENGTH * (-1.0_f32.powi(2) / (2.0 * variance_factor.powi(2))).exp();
        let influence_gauss_at_r = ball_gauss.influence_at(Vec2::new(TEST_RADIUS, 0.0));
        assert!(
            (influence_gauss_at_r - expected_gauss_at_r).abs() < EPSILON,
            "Gaussian at Radius: expected {}, got {}",
            expected_gauss_at_r,
            influence_gauss_at_r
        );
    }

    #[test]
    fn test_metaball_beyond_radius() {
        // Alle Falloffs sollten jenseits des Radius 0 oder sehr kleinen Einfluss haben,
        // je nach Cutoff-Logik. Die `calculate_influence_value` hat einen Cutoff.
        // Für Linear, SmoothStep, Wyvill ist es strikt 0.
        let falloffs = [
            MetaballFalloff::Algebraic {
                radius_scale_factor: 0.5,
            },
            MetaballFalloff::Exponential { steepness: 2.0 },
            MetaballFalloff::Gaussian {
                variance_factor: 0.5,
            },
            MetaballFalloff::Linear,
            MetaballFalloff::SmoothStep,
            MetaballFalloff::WyvillCubic,
        ];
        for falloff in falloffs.iter() {
            let ball = create_test_ball(*falloff);
            // Die `calculate_influence_value` schneidet bei `distance_sq > radius*radius + EPSILON` für die meisten Typen ab.
            // Die genaue Erwartung hängt davon ab, ob der Falloff selbst Null wird oder der Cutoff greift.
            // Der strikte Cutoff `distance_sq > self.radius * self.radius + constants::EPSILON` in `calculate_influence_value`
            // ist für Linear, SmoothStep, WyvillCubic relevant.
            // Für Algebraic, Exponential, Gaussian, die theoretisch weitergehen, wird ihr Wert durch die Formel bestimmt,
            // aber `normalized_distance` ist auf 1.0 geklemmt.
            // Die Prüfung `if point_distance_sq >= self.radius * self.radius && ...` ist etwas komplex.
            // Einfacher wäre, wenn die Falloffs bei norm_dist=1 alle sehr klein sind, und dann ein globaler Cutoff.

            // Testpunkt leicht außerhalb des Radius
            let point_beyond = Vec2::new(TEST_RADIUS + EPSILON * 10.0, 0.0);
            let influence = ball.influence_at(point_beyond);

            if matches!(
                *falloff,
                MetaballFalloff::Linear
                    | MetaballFalloff::SmoothStep
                    | MetaballFalloff::WyvillCubic
            ) {
                assert!(
                    influence.abs() < EPSILON,
                    "Beyond radius influence strictly zero for {:?}: expected 0, got {}",
                    falloff,
                    influence
                );
            } else {
                // Für die anderen wird der Wert bei norm_dist=1.0 berechnet.
                let influence_at_radius_edge = ball.influence_at(Vec2::new(TEST_RADIUS, 0.0));
                assert!(
                    (influence - influence_at_radius_edge).abs() < EPSILON,
                    "Beyond radius influence for {:?} should be same as at radius due to clamp: expected ~{}, got {}",
                    falloff,
                    influence_at_radius_edge,
                    influence
                );
            }

            // Testpunkt weit außerhalb des Radius
            let point_far_beyond = Vec2::new(TEST_RADIUS * 2.0, 0.0);
            let influence_far = ball.influence_at(point_far_beyond);
            if matches!(
                *falloff,
                MetaballFalloff::Linear
                    | MetaballFalloff::SmoothStep
                    | MetaballFalloff::WyvillCubic
            ) {
                assert!(
                    influence_far.abs() < EPSILON,
                    "Far beyond radius influence strictly zero for {:?}: expected 0, got {}",
                    falloff,
                    influence_far
                );
            } else {
                let influence_at_radius_edge = ball.influence_at(Vec2::new(TEST_RADIUS, 0.0));
                assert!(
                    (influence_far - influence_at_radius_edge).abs() < EPSILON,
                    "Far beyond radius influence for {:?} should be same as at radius due to clamp: expected ~{}, got {}",
                    falloff,
                    influence_at_radius_edge,
                    influence_far
                );
            }
        }
    }

    #[test]
    fn test_linear_falloff_halfway() {
        let ball = create_test_ball(MetaballFalloff::Linear);
        let halfway_point = Vec2::new(TEST_RADIUS * 0.5, 0.0); // Distanz = TEST_RADIUS * 0.5
        // normalized_distance = 0.5
        // Expected = strength * (1.0 - 0.5) = TEST_STRENGTH * 0.5
        let expected_influence = TEST_STRENGTH * 0.5;
        let influence = ball.influence_at(halfway_point);
        assert!(
            (influence - expected_influence).abs() < EPSILON,
            "Linear halfway: expected {}, got {}",
            expected_influence,
            influence
        );
    }

    #[test]
    fn test_wyvill_cubic_properties() {
        let ball = create_test_ball(MetaballFalloff::WyvillCubic);
        // Test an einem Punkt, wo 1 - (d/R)^2 = 0.5 ist.
        // (d/R)^2 = 0.5 => d/R = sqrt(0.5) => d = R * sqrt(0.5)
        let dist_for_half_term = TEST_RADIUS * (0.5_f32).sqrt();
        let point = Vec2::new(dist_for_half_term, 0.0);
        // term = 0.5. influence = strength * (0.5)^3 = strength * 0.125
        let expected_influence = TEST_STRENGTH * 0.125;
        let influence = ball.influence_at(point);
        assert!(
            (influence - expected_influence).abs() < EPSILON,
            "Wyvill specific point: expected {}, got {}",
            expected_influence,
            influence
        );
    }

    #[test]
    fn test_small_radius() {
        let ball = MetaballSource {
            position: Vec2::ZERO,
            radius: EPSILON * 0.5, // Sehr kleiner, aber positiver Radius
            strength: TEST_STRENGTH,
            falloff_type: MetaballFalloff::Linear,
        };

        // Am Zentrum: Volle Stärke
        let influence_center = ball.influence_at(Vec2::ZERO);
        assert!(
            (influence_center - TEST_STRENGTH).abs() < EPSILON,
            "Small radius, center"
        );

        // Knapp außerhalb des Mikroradius: Sollte 0 sein
        let influence_outside = ball.influence_at(Vec2::new(EPSILON * 2.0, 0.0));
        assert!(influence_outside.abs() < EPSILON, "Small radius, outside");
    }

    #[test]
    fn test_zero_radius() {
        // Radius wird im Konstruktor auf mind. EPSILON gesetzt.
        // Testen wir den Fall, wenn er manuell auf 0 gesetzt würde (was der Konstruktor verhindert).
        // Hier simulieren wir es, indem wir es fast 0 machen.
        let ball = MetaballSource {
            position: Vec2::ZERO,
            radius: constants::EPSILON, // Durch Konstruktor auf EPSILON gesetzt
            strength: TEST_STRENGTH,
            falloff_type: MetaballFalloff::Linear,
        };

        // Am Zentrum: Volle Stärke
        let influence_center = ball.influence_at(Vec2::ZERO);
        assert!(
            (influence_center - TEST_STRENGTH).abs() < EPSILON,
            "Zero-like radius, center: expected {}, got {}",
            TEST_STRENGTH,
            influence_center
        );

        // Knapp außerhalb des EPSILON-Radius: Sollte 0 sein, da norm_dist >= 1.0 wird
        let influence_outside = ball.influence_at(Vec2::new(constants::EPSILON * 2.0, 0.0));
        assert!(
            influence_outside.abs() < EPSILON,
            "Zero-like radius, outside: expected 0, got {}",
            influence_outside
        );
    }
}
