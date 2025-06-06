// src/math/geometry/metalballs/influence/polygon_source.rs
// (oder wie auch immer du die Datei nach der Refaktorierung nennst)

use crate::math::geometry::metalballs::influence::FieldInfluence; // Der Trait
use crate::math::{
    geometry::polygon::{core::Polygon, properties::PolygonProperties}, // Polygon und seine Eigenschaften
    types::Bounds2D,
    utils::{constants, numerical}, // Für EPSILON und ggf. sigmoid
};
use bevy::math::Vec2;

/// Definiert, wie der Einfluss eines Polygons mit dem Abstand abnimmt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PolygonFalloff {
    Linear,
    /// Normalisierter quadratischer Abfall: `strength * (1 - normalized_distance)^2`. Erreicht 0 am `falloff_distance`.
    QuadraticSmooth,
    /// Algebraischer Abfall, ähnlich Metaballs: `strength / ( (scaled_distance)^2 + 1)`.
    /// Muss sorgfältig skaliert werden, um bei `falloff_distance` nahe Null zu sein.
    AlgebraicInverseSq {
        /// Faktor, durch den `distance_abs` geteilt wird, bevor es quadriert wird.
        /// Ein kleinerer `distance_scale` führt zu einem schnelleren Abfall.
        /// z.B. `falloff_distance * 0.3`
        distance_scale: f32,
    },
    Exponential {
        /// Steuert, wie schnell der Einfluss abfällt. Größere Werte = schnellerer Abfall.
        steepness: f32,
    },
    Gaussian {
        /// Beeinflusst die "Breite" der Gausskurve relativ zur `falloff_distance`.
        /// `variance_factor` < 1.0 führt zu einer schmaleren Kurve (schnellerer Abfall).
        /// `variance_factor` = 1.0: Bei `distance_abs = falloff_distance` ist der Wert ca. `strength * exp(-2.0)`.
        /// Um bei `falloff_distance` nahezu 0 zu sein, sollte `variance_factor` eher klein sein (z.B. 0.25-0.5).
        variance_factor: f32,
    },
    /// Harter Übergang. Wert ist `strength` bis `midpoint_factor * falloff_distance`, danach 0.
    Step {
        /// Relativer Punkt (0.0 bis 1.0) innerhalb der `falloff_distance`, an dem der Sprung stattfindet.
        #[deprecated(
            note = "Step Falloff hat einen festen Midpoint bei 0.5. Für variable Midpoints Sigmoid verwenden oder Custom."
        )]
        midpoint_factor: f32, // Typischerweise 0.5 für einen Sprung in der Mitte
    },
    Sigmoid {
        /// Steuert die Steilheit des Übergangs. Größere Werte = steilerer Übergang.
        steepness: f32,
        /// Relativer Punkt (0.0 bis 1.0) innerhalb der `falloff_distance`, um den der Übergang zentriert ist.
        midpoint_factor: f32,
    },
}

impl Default for PolygonFalloff {
    fn default() -> Self {
        PolygonFalloff::Linear
    }
}

/// Definiert einen Einflussbereich, der von der Form eines Polygons ausgeht.
#[derive(Clone)] // Debug wird manuell implementiert
pub struct PolygonInfluence {
    pub polygon: Polygon,
    pub strength: f32,
    pub falloff_distance: f32,
    pub falloff_type: PolygonFalloff,
    pub interior_only: bool,
    pub offset: f32,
    cached_bounds_expanded: Option<Bounds2D>,
}

impl PolygonInfluence {
    /// Erstellt eine neue PolygonInfluence-Instanz.
    /// Der Builder wird für eine komfortablere Erstellung empfohlen.
    #[allow(clippy::too_many_arguments)]
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
            falloff_distance: falloff_distance.max(constants::EPSILON),
            falloff_type,
            interior_only,
            offset,
            cached_bounds_expanded: None,
        };
        new_self.recalculate_cached_bounds(); // Cache initial füllen
        new_self
    }

    /// Berechnet die erweiterte Bounding Box neu und speichert sie im Cache.
    /// Wird intern vom Konstruktor aufgerufen. Kann auch extern verwendet werden,
    /// falls sich Polygon-Eigenschaften nach der Erstellung ändern könnten (was hier nicht der Fall ist).
    pub fn recalculate_cached_bounds(&mut self) {
        if let Some(poly_bounds) = self.polygon.bounds() {
            // Annahme: Polygon::bounds() -> Option<Bounds2D>
            self.cached_bounds_expanded = Some(poly_bounds.expand(self.falloff_distance));
        } else {
            self.cached_bounds_expanded = None;
        }
    }

    /// Berechnet den signierten Abstand von einem Punkt zum Polygon.
    /// - Negativ: Punkt ist innerhalb des Polygons (Abstand zur nächsten Kante).
    /// - Positiv: Punkt ist außerhalb des Polygons (Abstand zur nächsten Kante).
    /// - Null: Punkt liegt auf der Kante.
    /// Gibt `f32::INFINITY` zurück, wenn das Polygon keine Kanten hat oder leer ist.
    fn signed_distance_to_polygon(&self, point: Vec2) -> f32 {
        if self.polygon.vertices().is_empty() {
            return f32::INFINITY;
        }
        let distance_to_boundary = self.distance_to_boundary(point);

        if self.polygon.contains_point(point) {
            -distance_to_boundary
        } else {
            distance_to_boundary
        }
    }

    /// Berechnet den kürzesten (positiven) Abstand von einem Punkt zur Grenze des Polygons.
    /// Gibt `f32::INFINITY` zurück, wenn das Polygon keine Kanten hat.
    fn distance_to_boundary(&self, point: Vec2) -> f32 {
        let vertices = self.polygon.vertices();
        let n = vertices.len();

        if n == 0 {
            return f32::INFINITY;
        }
        if n == 1 {
            // Ein einzelner Punkt als Polygon
            return point.distance(vertices[0]);
        }

        let mut min_dist_sq = f32::INFINITY;

        // Die Schleife muss korrekt über die Segmente des Polygons iterieren.
        // Polygon::closed stellt sicher, dass der erste und letzte Vertex bei Bedarf identisch sind
        // oder die Segmentlogik das Schließen behandelt.
        // Wir gehen hier davon aus, dass `vertices()` die Sequenz für ein geschlossenes Polygon liefert.
        let num_segments = if n > 1 && vertices.first() == vertices.last() {
            n - 1 // Explizit geschlossenes Polygon (Duplicate End Vertex)
        } else {
            n // Vertices bilden implizit ein geschlossenes Polygon (letzter zu erstem)
        };

        if num_segments == 0 {
            // Sollte durch n==1 oben abgedeckt sein, aber als Sicherheit
            return if n == 1 {
                point.distance(vertices[0])
            } else {
                f32::INFINITY
            };
        }

        for i in 0..num_segments {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % n]; // %n sorgt für den Umlauf zum ersten Vertex für das letzte Segment
            min_dist_sq = min_dist_sq.min(point_to_line_segment_distance_sq(point, p1, p2));
        }
        min_dist_sq.sqrt()
    }

    /// Wendet die konfigurierte Falloff-Funktion auf den Abstand an.
    /// `distance_abs` ist der absolute Abstand zum relevanten Teil des Polygons.
    fn apply_falloff(&self, distance_abs: f32) -> f32 {
        if distance_abs >= self.falloff_distance - constants::EPSILON {
            return 0.0;
        }
        // Vermeide Division durch Null, wenn falloff_distance extrem klein ist
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
            PolygonFalloff::AlgebraicInverseSq { distance_scale } => {
                let effective_dist = distance_abs / distance_scale.max(constants::EPSILON);
                self.strength / (effective_dist * effective_dist + 1.0)
            }
            PolygonFalloff::Exponential { steepness } => {
                self.strength * (-normalized_distance * steepness.max(0.01)).exp()
            }
            PolygonFalloff::Gaussian { variance_factor } => {
                // sigma^2 = (falloff_distance / (2 * variance_factor_sqrt))^2 (damit bei D der Wert klein ist)
                // Einfacher: variance = (falloff_distance * variance_factor_coeff)^2
                // variance_factor von 0.01 bis ~1.0. Kleiner = schmaler.
                // (variance_factor * falloff_distance) kann als eine Art Standardabweichung (sigma) gesehen werden.
                let sigma = self.falloff_distance * variance_factor.max(0.01);
                self.strength * (-distance_abs.powi(2) / (2.0 * sigma.powi(2))).exp()
            }
            #[allow(deprecated)] // Temporär, bis die API final ist
            PolygonFalloff::Step { midpoint_factor: _ } => {
                // midpoint_factor wird ignoriert, fester Midpoint
                if normalized_distance <= 0.5 {
                    // Harter Übergang in der Mitte
                    self.strength
                } else {
                    0.0
                }
            }
            PolygonFalloff::Sigmoid {
                steepness,
                midpoint_factor,
            } => {
                let midpoint_abs = self.falloff_distance * midpoint_factor.clamp(0.0, 1.0);
                // x_val normalisiert die Distanz relativ zum Midpoint und skaliert mit Steepness
                // Ein positiver x_val bedeutet, dass distance_abs > midpoint_abs
                let x_val = (distance_abs - midpoint_abs)
                    * (steepness.max(0.01) / self.falloff_distance.max(constants::EPSILON));
                self.strength * (1.0 - numerical::sigmoid(x_val)) // 1.0 - sigmoid für abfallende Kurve
            }
        }
    }
}

/// Hilfsfunktion: Quadratischer Abstand von Punkt zu Liniensegment.
/// Außerhalb von PolygonInfluence platziert, da es eine generische geometrische Funktion ist.
fn point_to_line_segment_distance_sq(point: Vec2, seg_a: Vec2, seg_b: Vec2) -> f32 {
    let l2 = seg_a.distance_squared(seg_b);
    if l2 < constants::EPSILON_SQUARED {
        // EPSILON_SQUARED für Performance und Konsistenz
        return point.distance_squared(seg_a); // Segment ist quasi ein Punkt
    }
    // Projektion des Punktes `point` auf die Linie, die durch `seg_a` und `seg_b` definiert ist.
    // t = dot(point - seg_a, seg_b - seg_a) / |seg_b - seg_a|^2
    let t = ((point - seg_a).dot(seg_b - seg_a) / l2).clamp(0.0, 1.0);
    let projection = seg_a + t * (seg_b - seg_a);
    point.distance_squared(projection)
}

// constants::EPSILON_SQUARED müsste in deinem constants Modul definiert sein:
// pub const EPSILON_SQUARED: f32 = EPSILON * EPSILON;

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
        if let Some(bounds) = self.cached_bounds_expanded {
            if !bounds.contains_point(point) {
                return 0.0;
            }
        }
        // Kein `else` Block nötig, da `cached_bounds_expanded` im Konstruktor gefüllt wird.
        // Wenn es `None` ist (z.B. leeres Polygon), ist das Verhalten ohnehin problematisch.

        let signed_dist = self.signed_distance_to_polygon(point);

        if signed_dist == f32::INFINITY {
            // Z.B. bei leerem Polygon
            return 0.0;
        }

        if self.interior_only && signed_dist > constants::EPSILON {
            return 0.0;
        }

        let mut influence_val = if self.interior_only {
            // Im `interior_only` Modus:
            // `signed_dist` ist negativ im Inneren. `abs(signed_dist)` ist der Abstand zum Rand.
            // Falloff wirkt vom Rand nach innen über `falloff_distance`.
            let dist_from_boundary_abs = signed_dist.abs();
            if dist_from_boundary_abs <= self.falloff_distance {
                // Innerhalb des Falloff-Bereichs (vom Rand nach innen)
                self.apply_falloff(dist_from_boundary_abs)
            } else {
                // Tief im Inneren, außerhalb des Falloff-Bereichs vom Rand -> volle Stärke
                self.strength
            }
        } else {
            // Nicht `interior_only`: Einfluss basiert auf dem absoluten Abstand zum Polygonrand.
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
        self.cached_bounds_expanded
    }

    fn max_influence_distance(&self) -> Option<f32> {
        // Dies ist eine grobe Schätzung der maximalen Ausdehnung.
        // Es ist die Distanz vom Ursprung (oder einem Punkt) zum am weitesten entfernten Punkt der
        // *erweiterten* Bounding Box.
        if let Some(bounds) = self.cached_bounds_expanded {
            // Max Distanz vom Zentrum der Bounds zu einer Ecke der Bounds
            let half_size = bounds.size() * 0.5;
            Some(half_size.length() + bounds.center().length()) // Konservative obere Schranke
        } else {
            Some(self.falloff_distance) // Fallback, wenn keine Bounds (z.B. leeres Polygon)
        }
    }

    fn influence_type(&self) -> &'static str {
        "PolygonInfluence" // Oder "PolygonSourceInfluence", je nach Umbenennung
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

    pub fn vertices(mut self, vertices: Vec<Vec2>) -> Self {
        match Polygon::closed(vertices) {
            Ok(polygon) => self.polygon = Some(polygon),
            Err(e) => {
                eprintln!(
                    "Failed to create polygon from vertices for PolygonInfluenceBuilder: {:?}. Polygon will be None.",
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

// Die Tests müssten hier angepasst werden, um die neuen Falloff-Typen
// und ihre Parameter korrekt zu verwenden.
// z.B.
// .falloff_type(PolygonFalloff::Exponential { steepness: 3.0 })
// .falloff_type(PolygonFalloff::Gaussian { variance_factor: 0.5 })
// etc.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::geometry::metalballs::influence::FieldInfluence; // Korrekter Pfad nach Refactoring
    use crate::math::utils::constants; // Für EPSILON in Tests

    // Hilfsfunktionen für Test-Polygone (TestShapes) bleiben hier nützlich,
    // müssen aber ggf. angepasst werden, um die neuen Falloff-Typen zu nutzen.

    struct TestShapes;

    impl TestShapes {
        // ... (rectangle_influence_builder, circle_influence_builder etc. hier einfügen) ...
        // Beispiel für Anpassung:
        // pub fn rectangle_builder(center: Vec2, width: f32, height: f32) -> PolygonInfluenceBuilder {
        //     let half_w = width * 0.5;
        //     let half_h = height * 0.5;
        //     let vertices = vec![
        //         Vec2::new(center.x - half_w, center.y - half_h),
        //         Vec2::new(center.x + half_w, center.y - half_h),
        //         Vec2::new(center.x + half_w, center.y + half_h),
        //         Vec2::new(center.x - half_w, center.y + half_h),
        //     ];
        //     PolygonInfluenceBuilder::new().vertices(vertices)
        // }
        fn circle_vertices(center: Vec2, radius: f32, segments: usize) -> Vec<Vec2> {
            let mut vertices = Vec::with_capacity(segments);
            for i in 0..segments {
                let angle = (i as f32 / segments as f32) * std::f32::consts::TAU;
                let x = center.x + radius * angle.cos();
                let y = center.y + radius * angle.sin();
                vertices.push(Vec2::new(x, y));
            }
            vertices
        }

        pub fn circle_influence_builder(
            center: Vec2,
            radius: f32,
            strength: f32,
            segments: usize,
        ) -> PolygonInfluenceBuilder {
            let vertices = Self::circle_vertices(center, radius, segments);
            PolygonInfluenceBuilder::new()
                .vertices(vertices)
                .strength(strength)
        }

        pub fn rectangle_influence_builder(
            center: Vec2,
            width: f32,
            height: f32,
            strength: f32,
        ) -> PolygonInfluenceBuilder {
            let half_w = width * 0.5;
            let half_h = height * 0.5;
            let vertices = vec![
                Vec2::new(center.x - half_w, center.y - half_h),
                Vec2::new(center.x + half_w, center.y - half_h),
                Vec2::new(center.x + half_w, center.y + half_h),
                Vec2::new(center.x - half_w, center.y + half_h),
            ];
            PolygonInfluenceBuilder::new()
                .vertices(vertices)
                .strength(strength)
        }
    }

    #[test]
    fn test_polygon_influence_linear_falloff() {
        let square = TestShapes::rectangle_influence_builder(Vec2::ZERO, 2.0, 2.0, 1.0)
            .falloff_distance(1.0) // Falloff über 1 Einheit Distanz
            .falloff_type(PolygonFalloff::Linear)
            .build()
            .unwrap();
        // square.recalculate_cached_bounds(); // Nicht mehr nötig, da new() es macht

        // Punkt im Zentrum des Quadrats (-1,-1) bis (1,1).
        // Abstand zur Grenze ist 1.0. `signed_distance_to_polygon` ist -1.0.
        // `abs(signed_dist)` ist 1.0.
        // `normalized_distance` = 1.0 / 1.0 = 1.0.
        // Linear Falloff: strength * (1.0 - 1.0) = 1.0 * 0.0 = 0.0.
        let center_influence = square.influence_at(Vec2::ZERO);
        assert!(
            (center_influence - 0.0).abs() < constants::EPSILON,
            "Center influence (Linear, d=boundary): expected 0.0, got {}",
            center_influence
        );

        // Punkt auf halbem Weg zur Grenze vom Zentrum aus (z.B. (0.5, 0.0))
        // Abstand zur Grenze ist 0.5. `abs(signed_dist)` = 0.5.
        // `normalized_distance` = 0.5 / 1.0 = 0.5.
        // Linear Falloff: strength * (1.0 - 0.5) = 1.0 * 0.5 = 0.5.
        let halfway_influence = square.influence_at(Vec2::new(0.5, 0.0));
        assert!(
            (halfway_influence - 0.5).abs() < constants::EPSILON,
            "Halfway influence (Linear): expected 0.5, got {}",
            halfway_influence
        );

        // Punkt direkt auf der Grenze (z.B. (1.0, 0.0))
        // Abstand zur Grenze ist 0.0. `abs(signed_dist)` = 0.0.
        // `normalized_distance` = 0.0 / 1.0 = 0.0.
        // Linear Falloff: strength * (1.0 - 0.0) = 1.0 * 1.0 = 1.0.
        let boundary_influence = square.influence_at(Vec2::new(1.0, 0.0));
        assert!(
            (boundary_influence - 1.0).abs() < constants::EPSILON,
            "Boundary influence (Linear): expected 1.0, got {}",
            boundary_influence
        );

        // Punkt außerhalb, aber innerhalb des Falloffs (z.B. (1.5, 0.0))
        // Polygon ist von x=-1 bis x=1. Nächste Grenze bei x=1. Abstand ist 0.5.
        // `abs(signed_dist)` = 0.5.
        // `normalized_distance` = 0.5 / 1.0 = 0.5.
        // Linear Falloff: strength * (1.0 - 0.5) = 1.0 * 0.5 = 0.5.
        let outside_falloff_influence = square.influence_at(Vec2::new(1.5, 0.0));
        assert!(
            (outside_falloff_influence - 0.5).abs() < constants::EPSILON,
            "Outside falloff influence (Linear): expected 0.5, got {}",
            outside_falloff_influence
        );

        // Punkt weit außerhalb (z.B. (3.0, 0.0))
        // Abstand zur Grenze ist 2.0. Das ist >= falloff_distance (1.0).
        // Sollte 0.0 sein.
        let far_influence = square.influence_at(Vec2::new(3.0, 0.0));
        assert!(
            (far_influence - 0.0).abs() < constants::EPSILON,
            "Far influence (Linear): expected 0.0, got {}",
            far_influence
        );
    }

    #[test]
    fn test_polygon_influence_exponential_falloff() {
        let steepness = 4.0; // Relativ starker Abfall
        let circle = TestShapes::circle_influence_builder(Vec2::ZERO, 1.0, 1.0, 32)
            .falloff_distance(1.0) // Falloff über 1 Einheit Distanz vom Kreisrand
            .falloff_type(PolygonFalloff::Exponential { steepness })
            .build()
            .unwrap();

        // Punkt auf dem Rand des Kreises (Radius 1.0), z.B. (1.0, 0.0)
        // `distance_abs` = 0.0. `normalized_distance` = 0.0.
        // Exp Falloff: strength * exp(-0.0 * steepness) = 1.0 * exp(0) = 1.0 * 1.0 = 1.0.
        let boundary_influence = circle.influence_at(Vec2::new(1.0, 0.0));
        assert!(
            (boundary_influence - 1.0).abs() < constants::EPSILON,
            "Boundary influence (Exp): expected 1.0, got {}",
            boundary_influence
        );

        // Punkt außerhalb des Kreises, auf halbem Weg durch den Falloff-Bereich.
        // Kreisrand bei r=1.0. Falloff bis r=2.0. Punkt bei r=1.5 (z.B. (1.5, 0.0)).
        // `distance_abs` = 0.5 (Distanz vom Kreisrand bei r=1.0).
        // `normalized_distance` = 0.5 / 1.0 = 0.5.
        // Exp Falloff: strength * exp(-0.5 * steepness) = 1.0 * exp(-0.5 * 4.0) = 1.0 * exp(-2.0).
        let expected_halfway = 1.0 * (-2.0_f32).exp(); // ca. 0.1353
        let halfway_influence = circle.influence_at(Vec2::new(1.5, 0.0));
        assert!(
            (halfway_influence - expected_halfway).abs() < constants::EPSILON,
            "Halfway falloff influence (Exp): expected {}, got {}",
            expected_halfway,
            halfway_influence
        );

        // Punkt am Ende des Falloff-Bereichs (z.B. (2.0, 0.0))
        // `distance_abs` = 1.0. `normalized_distance` = 1.0.
        // Exp Falloff: strength * exp(-1.0 * steepness) = 1.0 * exp(-4.0).
        // Dieser Wert ist nicht exakt 0, aber apply_falloff prüft `distance_abs >= self.falloff_distance`
        // und sollte 0 zurückgeben.
        let end_of_falloff_influence = circle.influence_at(Vec2::new(2.0, 0.0));
        assert!(
            (end_of_falloff_influence - 0.0).abs() < constants::EPSILON,
            "End of falloff influence (Exp): expected 0.0, got {}",
            end_of_falloff_influence
        );
    }

    #[test]
    fn test_interior_only_mode() {
        let square = TestShapes::rectangle_influence_builder(Vec2::ZERO, 2.0, 2.0, 1.0)
            .interior_only(true)
            .falloff_distance(0.5) // Falloff wirkt vom Rand 0.5 Einheiten nach innen
            .falloff_type(PolygonFalloff::Linear)
            .build()
            .unwrap();

        // Punkt tief im Inneren, außerhalb des inneren Falloff-Bereichs (z.B. (0.0, 0.0))
        // Quadratgrenzen bei +/-1. `signed_dist` = -1.0. `dist_from_boundary_abs` = 1.0.
        // Da `dist_from_boundary_abs` (1.0) > `falloff_distance` (0.5), sollte es `strength` (1.0) sein.
        let deep_inside_influence = square.influence_at(Vec2::ZERO);
        assert!(
            (deep_inside_influence - 1.0).abs() < constants::EPSILON,
            "Deep inside influence (InteriorOnly): expected 1.0, got {}",
            deep_inside_influence
        );

        // Punkt im inneren Falloff-Bereich (z.B. (0.75, 0.0))
        // `signed_dist` = -0.25 (Distanz zur Grenze x=1.0 ist 0.25). `dist_from_boundary_abs` = 0.25.
        // `dist_from_boundary_abs` (0.25) <= `falloff_distance` (0.5).
        // `normalized_distance` = 0.25 / 0.5 = 0.5.
        // Linear Falloff: strength * (1.0 - 0.5) = 1.0 * 0.5 = 0.5.
        let inner_falloff_influence = square.influence_at(Vec2::new(0.75, 0.0));
        assert!(
            (inner_falloff_influence - 0.5).abs() < constants::EPSILON,
            "Inner falloff influence (InteriorOnly): expected 0.5, got {}",
            inner_falloff_influence
        );

        // Punkt genau auf der Grenze (z.B. (1.0, 0.0))
        // `signed_dist` = 0.0. `dist_from_boundary_abs` = 0.0.
        // `normalized_distance` = 0.0 / 0.5 = 0.0.
        // Linear Falloff: strength * (1.0 - 0.0) = 1.0 * 1.0 = 1.0.
        let boundary_influence = square.influence_at(Vec2::new(1.0, 0.0));
        assert!(
            (boundary_influence - 1.0).abs() < constants::EPSILON,
            "Boundary influence (InteriorOnly): expected 1.0, got {}",
            boundary_influence
        );

        // Punkt außerhalb sollte keinen Einfluss haben
        let outside_influence = square.influence_at(Vec2::new(2.0, 2.0));
        assert!(
            (outside_influence - 0.0).abs() < constants::EPSILON,
            "Outside influence (InteriorOnly): expected 0.0, got {}",
            outside_influence
        );
    }

    #[test]
    fn test_empty_polygon_influence_build_none() {
        let influence = PolygonInfluenceBuilder::new()
            .vertices(vec![]) // Leeres Polygon
            .build();
        assert!(
            influence.is_none(),
            "Building with empty vertices should result in None"
        );
    }

    #[test]
    fn test_empty_polygon_influence_at_zero() {
        // Erstelle ein Polygon-Objekt direkt mit leeren Vertices, um `build()` zu umgehen
        let empty_poly = Polygon::new(vec![]).unwrap(); // Safe unwrap as empty polygon is valid
        let influence_with_empty_poly =
            PolygonInfluence::new(empty_poly, 1.0, 1.0, PolygonFalloff::Linear, false, 0.0);
        // recalculate_cached_bounds sollte None setzen
        assert!(influence_with_empty_poly.cached_bounds_expanded.is_none());

        // influence_at sollte 0.0 zurückgeben, wenn signed_distance_to_polygon f32::INFINITY ist
        let influence_val = influence_with_empty_poly.influence_at(Vec2::ZERO);
        assert!(
            (influence_val - 0.0).abs() < constants::EPSILON,
            "influence_at for empty polygon should be 0, got {}",
            influence_val
        );
    }

    #[test]
    fn test_single_point_polygon_influence() {
        let single_point_poly = Polygon::closed(vec![Vec2::new(1.0, 1.0)]).unwrap();
        let influence = PolygonInfluenceBuilder::new()
            .polygon(single_point_poly)
            .strength(1.0)
            .falloff_distance(2.0)
            .falloff_type(PolygonFalloff::Linear)
            .build()
            .unwrap();

        // Punkt am Polygon-Punkt
        let val_at_point = influence.influence_at(Vec2::new(1.0, 1.0));
        // dist=0, norm_dist=0, influence = 1.0
        assert!(
            (val_at_point - 1.0).abs() < constants::EPSILON,
            "Val at single point: {}",
            val_at_point
        );

        // Punkt 1 Einheit entfernt
        let val_at_dist_1 = influence.influence_at(Vec2::new(2.0, 1.0));
        // dist=1, falloff_dist=2, norm_dist=0.5, influence = 1.0 * (1-0.5) = 0.5
        assert!(
            (val_at_dist_1 - 0.5).abs() < constants::EPSILON,
            "Val at dist 1 from single point: {}",
            val_at_dist_1
        );

        // Punkt am Ende des Falloffs
        let val_at_falloff_end = influence.influence_at(Vec2::new(3.0, 1.0));
        // dist=2, falloff_dist=2, norm_dist=1.0, influence = 1.0 * (1-1) = 0.0
        // apply_falloff sollte 0 geben, da dist >= falloff_distance
        assert!(
            (val_at_falloff_end - 0.0).abs() < constants::EPSILON,
            "Val at falloff end from single point: {}",
            val_at_falloff_end
        );
    }
}
