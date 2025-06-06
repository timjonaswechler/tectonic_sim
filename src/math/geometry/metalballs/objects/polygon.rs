// src/math/geometry/metalballs/objects/polygon.rs

use crate::math::geometry::metalballs::traits::FieldInfluence; // Der Trait
use crate::math::{
    geometry::{
        polygon::{core::Polygon, properties::PolygonProperties}, // Unser Polygon und seine Eigenschaften
    },
    types::Bounds2D,  // Unser Bounds2D
    utils::constants, // Für EPSILON
};
use bevy::math::Vec2;

/// Definiert, wie der Einfluss eines Polygons mit dem Abstand abnimmt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PolygonFalloff {
    Linear,
    /// Quadratischer Abfall von der Stärke zum Rand (strength * (1-norm_dist)^2)
    QuadraticSmooth,
    /// Algebraischer Abfall, ähnlich wie Metaballs (strength / (dist_sq_scaled + 1))
    AlgebraicInverseSq,
    Exponential {
        steepness: f32,
    },
    Gaussian {
        variance_factor: f32,
    }, // variance_factor < 1 -> schmaler, > 1 -> breiter als falloff_distance
    Step, // Harter Übergang
    Sigmoid {
        steepness: f32,
        midpoint_factor: f32,
    }, // midpoint_factor von 0 bis 1
}

impl Default for PolygonFalloff {
    fn default() -> Self {
        PolygonFalloff::Linear // Ein einfacher, berechenbarer Standard
    }
}

/// Definiert einen Einflussbereich, der von der Form eines Polygons ausgeht.
#[derive(Clone)] // Debug wird manuell implementiert
pub struct PolygonInfluence {
    pub polygon: Polygon, // Unser Polygon, das Vec<Vec2> verwendet
    pub strength: f32,
    pub falloff_distance: f32, // Abstand vom Polygonrand, über den der Einfluss abfällt
    pub falloff_type: PolygonFalloff,
    /// Wenn true, wirkt der Einfluss nur innerhalb des Polygons (ggf. mit Falloff zum Rand).
    /// Wenn false, wirkt der Einfluss auch außerhalb, abfallend vom Rand.
    pub interior_only: bool,
    /// Ein konstanter Offset, der zum berechneten Einfluss addiert wird (nach Falloff).
    pub offset: f32,
    // scale_factor wurde entfernt, da die Falloff-Typen jetzt eigene Parameter haben können.
    cached_bounds_expanded: Option<Bounds2D>, // Bounding Box des Polygons + falloff_distance
}

impl PolygonInfluence {
    /// Erstellt eine neue PolygonInfluence-Instanz.
    /// Der Builder wird für eine komfortablere Erstellung empfohlen.
    #[allow(clippy::too_many_arguments)] // Akzeptabel für einen internen Konstruktor
    fn new(
        polygon: Polygon,
        strength: f32,
        falloff_distance: f32,
        falloff_type: PolygonFalloff,
        interior_only: bool,
        offset: f32,
    ) -> Self {
        let mut new_self = Self {
            polygon,
            strength,
            falloff_distance: falloff_distance.max(constants::EPSILON), // Muss positiv sein
            falloff_type,
            interior_only,
            offset,
            cached_bounds_expanded: None,
        };
        new_self.recalculate_cached_bounds(); // Cache initial füllen
        new_self
    }

    /// Berechnet die erweiterte Bounding Box neu und speichert sie im Cache.
    fn recalculate_cached_bounds(&mut self) {
        if let Some(poly_bounds) = self.polygon.bounds() {
            // Polygon::bounds() gibt Option<Bounds2D>
            // Erweitere um Falloff-Distanz
            self.cached_bounds_expanded = Some(poly_bounds.expand(self.falloff_distance));
        } else {
            // Fallback für leere Polygone oder Polygone ohne gültige Bounds
            self.cached_bounds_expanded = None;
        }
    }

    /// Berechnet den signierten Abstand von einem Punkt zum Polygon.
    /// - Negativ: Punkt ist innerhalb des Polygons (Abstand zur nächsten Kante).
    /// - Positiv: Punkt ist außerhalb des Polygons (Abstand zur nächsten Kante).
    /// - Null: Punkt liegt auf der Kante.
    /// Gibt f32::INFINITY zurück, wenn das Polygon keine Kanten hat.
    fn signed_distance_to_polygon(&self, point: Vec2) -> f32 {
        let distance_to_boundary = self.distance_to_boundary(point);

        if self.polygon.contains_point(point) {
            // PolygonProperties Trait wird benötigt
            -distance_to_boundary
        } else {
            distance_to_boundary
        }
    }

    /// Berechnet den kürzesten (positiven) Abstand von einem Punkt zur Grenze des Polygons.
    /// Gibt f32::INFINITY zurück, wenn das Polygon keine Kanten hat.
    fn distance_to_boundary(&self, point: Vec2) -> f32 {
        let vertices = self.polygon.vertices();
        let n = vertices.len();

        if n == 0 {
            return f32::INFINITY;
        }
        if n == 1 {
            return point.distance(vertices[0]);
        }

        let mut min_dist_sq = f32::INFINITY;

        let num_segments = if self.polygon.is_closed() {
            if n > 0 && vertices.first() == vertices.last() {
                n - 1
            } else {
                n
            }
        } else {
            n - 1
        };
        if num_segments == 0 {
            return if n == 1 {
                point.distance(vertices[0])
            } else {
                f32::INFINITY
            };
        }

        for i in 0..num_segments {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % n]; // %n für korrekten Umlauf, auch wenn num_segments < n
            min_dist_sq = min_dist_sq.min(point_to_line_segment_distance_sq(point, p1, p2));
        }
        min_dist_sq.sqrt()
    }

    /// Wendet die konfigurierte Falloff-Funktion auf den Abstand an.
    /// `distance_abs` ist der absolute Abstand zum relevanten Teil des Polygons (Rand oder Inneres).
    fn apply_falloff(&self, distance_abs: f32) -> f32 {
        if distance_abs >= self.falloff_distance - constants::EPSILON {
            return 0.0;
        }
        if self.falloff_distance < constants::EPSILON {
            return if distance_abs < constants::EPSILON {
                self.strength
            } else {
                0.0
            };
        }

        let normalized_distance = (distance_abs / self.falloff_distance).clamp(0.0, 1.0);

        match self.falloff_type {
            PolygonFalloff::Linear => self.strength * (1.0 - normalized_distance),
            PolygonFalloff::QuadraticSmooth => self.strength * (1.0 - normalized_distance).powi(2),
            PolygonFalloff::AlgebraicInverseSq => {
                // Skaliert, sodass bei dist=0 -> strength; bei dist=falloff_distance -> strength / (k_scaled^2+1)
                // Wir wollen, dass es bei falloff_distance nahe 0 ist.
                // Eine übliche Form, die bei 0 max ist und abfällt: S / ( (d/D_half)^2 + 1 )
                // Wobei D_half z.B. falloff_distance / 2 ist.
                let effective_dist =
                    distance_abs / (self.falloff_distance * 0.3).max(constants::EPSILON); // Skalierungsfaktor
                self.strength / (effective_dist * effective_dist + 1.0)
            }
            PolygonFalloff::Exponential { steepness } => {
                self.strength * (-normalized_distance * steepness.max(0.1)).exp()
            }
            PolygonFalloff::Gaussian { variance_factor } => {
                // variance_factor: <1 schmal, >1 breit. Bei 1.0 erreicht es ca. 0.6*strength am Rand.
                let variance = (self.falloff_distance / variance_factor.max(0.1) / 2.0).powi(2); // sigma^2
                self.strength * (-distance_abs.powi(2) / (2.0 * variance)).exp()
            }
            PolygonFalloff::Step => {
                // Harter Übergang in der Mitte des Falloff-Bereichs
                if normalized_distance <= 0.5 {
                    self.strength
                } else {
                    0.0
                }
            }
            PolygonFalloff::Sigmoid {
                steepness,
                midpoint_factor,
            } => {
                let midpoint = self.falloff_distance * midpoint_factor.clamp(0.0, 1.0);
                // Die Normalisierung von steepness ist hier wichtig.
                // x für Sigmoid: (distance - midpoint) * (steepness / relevant_width)
                // relevant_width könnte falloff_distance sein.
                let x_val =
                    (distance_abs - midpoint) * (steepness.max(0.1) / self.falloff_distance);
                self.strength * (1.0 - crate::math::utils::numerical::sigmoid(x_val)) // 1-sigmoid für Abfall
            }
        }
    }
}

/// Hilfsfunktion: Quadratischer Abstand von Punkt zu Liniensegment.
fn point_to_line_segment_distance_sq(point: Vec2, seg_a: Vec2, seg_b: Vec2) -> f32 {
    let l2 = seg_a.distance_squared(seg_b);
    if l2 < constants::EPSILON * constants::EPSILON {
        // Segment ist ein Punkt
        return point.distance_squared(seg_a);
    }
    // Betrachte die Linie als p(t) = seg_a + t * (seg_b - seg_a)
    // Projektion des Punktes auf die Linie
    let t = ((point - seg_a).dot(seg_b - seg_a) / l2).clamp(0.0, 1.0);
    let projection = seg_a + t * (seg_b - seg_a);
    point.distance_squared(projection)
}

impl std::fmt::Debug for PolygonInfluence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PolygonInfluence")
            .field("polygon_vertices", &self.polygon.vertices().len())
            .field("strength", &self.strength)
            .field("falloff_distance", &self.falloff_distance)
            .field("falloff_type", &self.falloff_type)
            .field("interior_only", &self.interior_only)
            .field("offset", &self.offset)
            .field(
                "cached_bounds_expanded",
                &self.cached_bounds_expanded.is_some(),
            )
            .finish()
    }
}

impl FieldInfluence for PolygonInfluence {
    fn influence_at(&self, point: Vec2) -> f32 {
        // Früher Ausstieg, wenn Punkt außerhalb der erweiterten Bounding Box liegt.
        if let Some(bounds) = self.cached_bounds_expanded {
            if !bounds.contains_point(point) {
                return 0.0;
            }
        } else {
            // Sollte nicht passieren, wenn Konstruktor/Setter den Cache füllen,
            // aber als Fallback: keine Optimierung.
        }

        let signed_dist = self.signed_distance_to_polygon(point);

        if self.interior_only && signed_dist > constants::EPSILON {
            // Punkt ist außerhalb
            return 0.0;
        }

        let mut influence_val = if self.interior_only {
            // Im interior_only Modus ist der "Abstand" für den Falloff
            // der Abstand zur *nächsten* Grenze, auch wenn der Punkt tief im Inneren ist.
            // Oder man definiert, dass im Inneren immer volle Stärke herrscht bis zum Falloff-Bereich am Rand.
            // Hier nehmen wir an: Volle Stärke im Inneren, Falloff beginnt `falloff_distance` vom Rand entfernt nach innen.
            // `signed_dist` ist negativ im Inneren. `abs(signed_dist)` ist der Abstand zum Rand.
            let dist_from_boundary_abs = signed_dist.abs();
            if dist_from_boundary_abs <= self.falloff_distance {
                self.apply_falloff(dist_from_boundary_abs) // Falloff vom Rand nach innen
            } else {
                self.strength // Tief im Inneren, außerhalb des Falloff-Bereichs vom Rand
            }
        } else {
            // Nicht interior_only: Einfluss basiert auf Abstand zum Rand (positiv oder negativ)
            self.apply_falloff(signed_dist.abs())
        };

        // Wende Offset an, wenn ein relevanter Einfluss besteht
        if influence_val > constants::EPSILON
            || (self.interior_only && signed_dist <= constants::EPSILON)
        {
            influence_val += self.offset;
        }

        influence_val.max(0.0) // Endgültiger Einfluss darf nicht negativ sein
    }

    fn bounding_box(&self) -> Option<Bounds2D> {
        self.cached_bounds_expanded // Gibt die *erweiterte* Bounding Box zurück
    }

    fn max_influence_distance(&self) -> Option<f32> {
        // Die tatsächliche Reichweite hängt vom Polygon ab, aber falloff_distance
        // ist die Distanz vom Rand des Polygons aus.
        // Wenn das Polygon selbst eine Ausdehnung hat, ist die max. Distanz größer.
        // Für eine konservative Schätzung nehmen wir die Diagonale der BBox + falloff_distance.
        if let Some(bounds) = self.cached_bounds_expanded {
            Some(bounds.size().length() / 2.0 + self.falloff_distance) // Grobe Schätzung
        } else {
            Some(self.falloff_distance) // Fallback
        }
    }

    fn influence_type(&self) -> &'static str {
        "PolygonInfluence"
    }
}

/// Builder für `PolygonInfluence`-Instanzen.
pub struct PolygonInfluenceBuilder {
    polygon: Option<Polygon>,
    strength: f32,
    falloff_distance: f32,
    falloff_type: PolygonFalloff,
    interior_only: bool,
    offset: f32,
}

impl Default for PolygonInfluenceBuilder {
    fn default() -> Self {
        Self {
            polygon: None,
            strength: 1.0,
            falloff_distance: 10.0,
            falloff_type: PolygonFalloff::default(),
            interior_only: false,
            offset: 0.0,
        }
    }
}

impl PolygonInfluenceBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn polygon(mut self, polygon: Polygon) -> Self {
        self.polygon = Some(polygon);
        self
    }

    /// Erstellt das Polygon aus Vertices und schließt es.
    pub fn vertices(mut self, vertices: Vec<Vec2>) -> Self {
        match Polygon::closed(vertices) {
            // Verwendet Polygon::closed für eine garantierte geschlossene Form
            Ok(polygon) => self.polygon = Some(polygon),
            Err(e) => {
                // Fehler loggen oder behandeln. Hier setzen wir Polygon auf None.
                eprintln!(
                    "Failed to create polygon for PolygonInfluenceBuilder: {:?}",
                    e
                );
                self.polygon = None;
            }
        }
        self
    }

    pub fn strength(mut self, strength: f32) -> Self {
        self.strength = strength;
        self
    }

    pub fn falloff_distance(mut self, distance: f32) -> Self {
        self.falloff_distance = distance.max(constants::EPSILON);
        self
    }

    pub fn falloff_type(mut self, falloff_type: PolygonFalloff) -> Self {
        self.falloff_type = falloff_type;
        self
    }

    pub fn interior_only(mut self, interior_only: bool) -> Self {
        self.interior_only = interior_only;
        self
    }

    pub fn offset(mut self, offset: f32) -> Self {
        self.offset = offset;
        self
    }

    pub fn build(self) -> Option<PolygonInfluence> {
        self.polygon.map(|p| {
            PolygonInfluence::new(
                p,
                self.strength,
                self.falloff_distance,
                self.falloff_type,
                self.interior_only,
                self.offset,
            )
        })
    }
}
