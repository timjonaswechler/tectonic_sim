// src/math/algorithms/boolean.rs

use crate::math::{
    algorithms::clipping::{ClippingAlgorithm, PolygonClipper},
    error::{MathError, MathResult},
    types::Bounds2D, // Unser Bounds2D mit Vec2
    utils::constants,
};
use bevy::math::Vec2;
// Wir benötigen eine Polygon-Struktur, um mit Polygonen zu arbeiten und deren Eigenschaften zu nutzen.
// Da dies ein Algorithmen-Modul ist, sollten wir idealerweise nicht direkt von `crate::math::geometry::polygon::Polygon`
// abhängig sein, sondern mit Punktlisten arbeiten. Die Konvertierung zu/von Polygon erfolgt dann außerhalb.
// Für Methoden wie `polygon_contains_polygon_points` benötigen wir aber eine Art "Polygon-Checker".
// Wir könnten einen Trait `PolygonLike` definieren oder Hilfsfunktionen, die mit `&[Vec2]` arbeiten.

// Fürs Erste verwenden wir Hilfsfunktionen, die `&[Vec2]` als Polygone interpretieren.

/// Definiert die Art der durchzuführenden Boolean-Operation auf Polygonen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOpType {
    Union,        // Vereinigung (A ∪ B)
    Intersection, // Schnittmenge (A ∩ B)
    Difference,   // Differenz (A - B)
    Xor,          // Symmetrische Differenz (A ⊕ B)
}

/// Führt Boolean-Operationen (Vereinigung, Schnitt, Differenz) auf Polygonen durch,
/// die als Listen von 2D-Punkten repräsentiert werden.
pub struct PolygonBooleanOps {
    clipping_algorithm: ClippingAlgorithm, // Algorithmus für Schnittmengen
    tolerance: f32,
    // simplify_result: bool, // Vereinfachung sollte außerhalb erfolgen, z.B. durch einen Post-Processing-Schritt
}

impl Default for PolygonBooleanOps {
    fn default() -> Self {
        Self {
            clipping_algorithm: ClippingAlgorithm::SutherlandHodgman,
            tolerance: constants::EPSILON * 10.0,
            // simplify_result: true,
        }
    }
}

impl PolygonBooleanOps {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_clipping_algorithm(mut self, algo: ClippingAlgorithm) -> Self {
        self.clipping_algorithm = algo;
        self
    }

    pub fn with_tolerance(mut self, tolerance: f32) -> Self {
        self.tolerance = tolerance.max(0.0);
        self
    }

    /// Führt die spezifizierte Boolean-Operation zwischen zwei Polygonen durch.
    /// `poly_a_points` und `poly_b_points` sind die Eckpunkte der Polygone.
    /// Gibt eine Liste von Punktlisten zurück, die die resultierenden Polygone darstellen.
    pub fn execute(
        &self,
        poly_a_points: &[Vec2],
        poly_b_points: &[Vec2],
        operation: BooleanOpType,
    ) -> MathResult<Vec<Vec<Vec2>>> {
        if poly_a_points.len() < 3 || poly_b_points.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: poly_a_points.len().min(poly_b_points.len()),
            });
        }
        // TODO: Validierung auf einfache Polygone (keine Selbstüberschneidungen) wäre hier gut.

        match operation {
            BooleanOpType::Intersection => self.intersection(poly_a_points, poly_b_points),
            BooleanOpType::Union => self.union(poly_a_points, poly_b_points),
            BooleanOpType::Difference => self.difference(poly_a_points, poly_b_points),
            BooleanOpType::Xor => self.xor(poly_a_points, poly_b_points),
        }
    }

    /// Berechnet die Schnittmenge (Intersection) von Polygon A und Polygon B.
    /// A ∩ B
    fn intersection(
        &self,
        poly_a_points: &[Vec2],
        poly_b_points: &[Vec2],
    ) -> MathResult<Vec<Vec<Vec2>>> {
        // Intersection kann oft als Clipping von A gegen B (oder B gegen A) realisiert werden.
        // Sutherland-Hodgman benötigt einen konvexen Clipper. Wenn poly_b nicht konvex ist,
        // müsste man poly_b in konvexe Teile zerlegen oder einen komplexeren Algorithmus verwenden.
        // Für diese Implementierung nehmen wir an, der Clipper (poly_b) ist konvex,
        // oder wir akzeptieren, dass das Ergebnis für nicht-konvexe Clipper möglicherweise nicht korrekt ist.
        let clipper = PolygonClipper::new(self.clipping_algorithm).with_tolerance(self.tolerance);

        // Clip A gegen B
        let result_a_clipped_by_b = clipper.clip_polygon_points(poly_a_points, poly_b_points)?;

        // Für eine symmetrischere Intersection könnte man auch B gegen A clippen und die Ergebnisse vereinen/mitteln,
        // aber das ist komplexer. Sutherland-Hodgman A gegen B ist der Standardansatz.
        Ok(result_a_clipped_by_b)
    }

    /// Berechnet die Vereinigung (Union) von Polygon A und Polygon B.
    /// A ∪ B
    /// Dies ist eine sehr vereinfachte Implementierung. Robuste Union ist komplex.
    fn union(&self, poly_a_points: &[Vec2], poly_b_points: &[Vec2]) -> MathResult<Vec<Vec<Vec2>>> {
        // Eine wirklich korrekte Union ist sehr komplex (z.B. Greiner-Hormann oder Vatti).
        // Eine vereinfachte Näherung für überlappende Polygone, die nicht ineinander liegen:
        // 1. Prüfe auf triviale Fälle: Kein Überlapp, A in B, B in A.
        // 2. Wenn sie überlappen:
        //    - Berechne A - B (Punkte von A, die nicht in B sind)
        //    - Berechne B - A (Punkte von B, die nicht in A sind)
        //    - Berechne A ∩ B
        //    - Die Union ist dann (A-B) ∪ (B-A) ∪ (A∩B), aber das ist nicht einfach zusammenzusetzen.
        //
        // Alternative vereinfachte Logik (oft ungenau für konkave Formen):
        // Wenn A B enthält, ist Union A. Wenn B A enthält, ist Union B.
        // Wenn sie sich schneiden, ist eine robuste Lösung nötig.
        // Wenn sie disjunkt sind, ist Union {A, B}.

        let bounds_a = Bounds2D::from_points_iter(poly_a_points.iter().copied());
        let bounds_b = Bounds2D::from_points_iter(poly_b_points.iter().copied());

        if let (Some(ba), Some(bb)) = (bounds_a, bounds_b) {
            if !ba.intersects(&bb) {
                // Schneller Test: Bounding Boxes überlappen nicht
                // Grobe Prüfung, ob Schnittpunkte existieren, um sicherzustellen, dass sie nicht nur an den Bounds scheitern
                if find_line_segment_intersections(poly_a_points, poly_b_points, self.tolerance)
                    .is_empty()
                {
                    return Ok(vec![poly_a_points.to_vec(), poly_b_points.to_vec()]); // Disjunkt
                }
            }
        }

        // Enthält A B vollständig?
        if polygon_contains_all_points(poly_a_points, poly_b_points, self.tolerance) {
            return Ok(vec![poly_a_points.to_vec()]);
        }
        // Enthält B A vollständig?
        if polygon_contains_all_points(poly_b_points, poly_a_points, self.tolerance) {
            return Ok(vec![poly_b_points.to_vec()]);
        }

        // Hier wäre ein robuster Algorithmus (z.B. Greiner-Hormann) nötig.
        // Als sehr grobe Näherung für diesen Platzhalter:
        // Berechne die konvexe Hülle aller Punkte. Das ist nur für konvexe Eingaben halbwegs sinnvoll.
        let mut all_points = poly_a_points.to_vec();
        all_points.extend_from_slice(poly_b_points);

        let hull_computer = super::convex_hull::ConvexHullComputer::new(
            super::convex_hull::ConvexHullAlgorithm::AndrewMonotone,
        )
        .with_tolerance(self.tolerance);
        let hull = hull_computer.compute_hull_points(&all_points)?;
        Ok(vec![hull])

        // TODO: Ersetze dies durch einen korrekten Union-Algorithmus.
        // Ein Ansatz mit Clipping: Union(A,B) = A + B - Intersection(A,B) (symbolisch)
        // Oder: A + (B - A)
        // let diff_b_minus_a = self.difference(poly_b_points, poly_a_points)?;
        // let mut result_polygons = vec![poly_a_points.to_vec()];
        // result_polygons.extend(diff_b_minus_a);
        // Diese müssten dann ggf. noch zu einem einzelnen Polygon verschmolzen werden, falls sie sich berühren.
        // Das ist das Kernproblem der Union/Differenz.
    }

    /// Berechnet die Differenz von Polygon A und Polygon B.
    /// A - B (Teile von A, die nicht in B liegen).
    fn difference(
        &self,
        poly_a_points: &[Vec2], // Das Polygon, von dem abgezogen wird
        poly_b_points: &[Vec2], // Das Polygon, das abgezogen wird
    ) -> MathResult<Vec<Vec<Vec2>>> {
        // A - B kann als Clipping von A gegen das Komplement von B betrachtet werden.
        // Einfacher mit Sutherland-Hodgman: Clippe A gegen jede Kante von B, aber invertiere die "Innen"-Bedingung.
        // S-H clippt A *innerhalb* des Clippers. Für A - B wollen wir A *außerhalb* von B.
        // Dies erfordert, dass B konvex ist.
        // Das Ergebnis kann mehrere disjunkte Polygone sein.

        // Invertiere die Orientierung von B (wenn CCW -> CW, wenn CW -> CCW)
        // und clippe A gegen dieses "invertierte" B. Das ist aber nicht ganz korrekt für S-H.
        // S-H mit invertierter Innen/Außen-Prüfung:
        let clipper = PolygonClipper::new(self.clipping_algorithm).with_tolerance(self.tolerance);

        // Die Standard-Clipping-Implementierung für S-H gibt den Teil von A *innerhalb* B zurück.
        // Für A - B:
        // 1. Finde Schnittpunkte zwischen A und B.
        // 2. Konstruiere die resultierenden Polygone. Dies ist der komplexe Teil.

        // Als vereinfachte Näherung, die nur funktioniert, wenn B konvex ist und A von B "durchstochen" wird:
        // Clippe A gegen die "Außenseite" jeder Kante von B.
        // Für diese Demo verwenden wir einen Trick, wenn `poly_b` als konvex angenommen werden kann:
        // Man kann A gegen ein sehr großes Polygon clippen, aus dem B "ausgeschnitten" wurde. Das ist aber aufwändig.

        // Fallback zu einer sehr einfachen Logik:
        // Wenn A B enthält, ist das Ergebnis leer (oder A mit einem Loch, was wir hier nicht unterstützen).
        if polygon_contains_all_points(poly_b_points, poly_a_points, self.tolerance) {
            // B enthält A: A - B = {}
            return Ok(Vec::new());
        }

        // Hier wäre ein robuster Algorithmus (z.B. Weiler-Atherton oder Greiner-Hormann) nötig.
        // Placeholder: Gibt A zurück, wenn keine offensichtliche Subtraktion möglich ist
        // oder wenn B disjunkt ist.
        // Dies ist FALSCH für die meisten Fälle von A-B!

        // Ein besserer, aber immer noch nicht perfekter Ansatz mit S-H:
        // Teile B in Kanten. Für jede Kante von B, betrachte die Halbebene *außerhalb* von B.
        // Clippe A gegen diese Halbebenen.
        // Dies ist äquivalent zum Clippen von A gegen ein Polygon, das die "Außenseite" von B darstellt.
        // Wenn B CCW ist, zeigt "links" nach innen. Wir wollen also rechts von jeder Kante clippen.
        // Dies bedeutet, wir müssen die `is_point_inside_clip_edge` Bedingung umkehren
        // oder die Kanten von B in umgekehrter Reihenfolge (CW) verwenden.

        let mut reversed_b_points = poly_b_points.to_vec();
        reversed_b_points.reverse(); // Macht aus CCW -> CW (oder umgekehrt)

        // Jetzt clippe A gegen dieses "orientierungsinvertierte" B.
        // Das Ergebnis ist der Teil von A, der im "Inneren" des orientierungsinvertierten B liegt.
        // Dies ist konzeptionell A ∩ (Komplement von B) wenn B konvex.
        let result = clipper.clip_polygon_points(poly_a_points, &reversed_b_points)?;
        Ok(result)

        // WICHTIG: Diese Implementierung der Differenz mit S-H und Umkehrung des Clippers
        // ist nur für konvexe Clipper `poly_b` annähernd korrekt und kann zu unerwarteten Ergebnissen führen.
        // Eine robuste Differenz ist komplex.
    }

    /// Berechnet die symmetrische Differenz (XOR) von Polygon A und Polygon B.
    /// A ⊕ B = (A - B) ∪ (B - A)
    fn xor(&self, poly_a_points: &[Vec2], poly_b_points: &[Vec2]) -> MathResult<Vec<Vec<Vec2>>> {
        let diff_a_minus_b = self.difference(poly_a_points, poly_b_points)?;
        let diff_b_minus_a = self.difference(poly_b_points, poly_a_points)?;

        let mut result_polygons = diff_a_minus_b;
        result_polygons.extend(diff_b_minus_a);

        // Hier könnte eine zusätzliche Vereinigung/Bereinigung der result_polygons nötig sein,
        // wenn sie sich überlappen oder berühren.
        Ok(result_polygons)
    }
}

// --- Hilfsfunktionen für Polygon-Punktlisten ---

/// Prüft, ob ein Polygon (definiert durch `poly_points`) einen Punkt enthält.
/// Verwendet den Ray-Casting-Algorithmus. `poly_points` sollte ein einfaches Polygon sein.
fn point_in_polygon(point: Vec2, poly_points: &[Vec2], tolerance: f32) -> bool {
    let n = poly_points.len();
    if n < 3 {
        return false;
    }

    let mut inside = false;
    let mut p1 = poly_points[n - 1]; // Letzter Punkt als Start für erste Kante

    for i in 0..n {
        let p2 = poly_points[i];
        // Ray (von point nach +X) schneidet Kante (p1, p2)?
        if (p1.y > point.y) != (p2.y > point.y) {
            // Kante kreuzt die horizontale Linie des Punktes
            // Berechne x-Koordinate des Schnittpunkts der Kante mit der horizontalen Linie
            // Handle horizontale Kanten und Kanten, die genau auf der Höhe des Punktes enden, sorgfältig.
            if (p2.y - p1.y).abs() > tolerance {
                // Nicht-horizontale Kante
                let x_intersection = (point.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;
                if x_intersection > point.x {
                    // Schnittpunkt liegt rechts vom Punkt
                    inside = !inside;
                }
            } else if p1.y == point.y
                && ((p1.x >= point.x && p2.x < point.x - tolerance)
                    || (p2.x >= point.x && p1.x < point.x - tolerance))
            {
                // Horizontale Kante, die rechts vom Punkt beginnt/endet und den Punkt überquert
                // Dieser Fall wird oft speziell behandelt oder ignoriert, um Robustheit zu gewährleisten.
                // Hier für Einfachheit: Zähle nur klare Kreuzungen.
            }
        }
        p1 = p2;
    }
    inside
}

/// Prüft, ob alle Punkte von `inner_poly_points` innerhalb von `outer_poly_points` liegen.
fn polygon_contains_all_points(
    outer_poly_points: &[Vec2],
    inner_poly_points: &[Vec2],
    tolerance: f32,
) -> bool {
    if outer_poly_points.len() < 3 {
        return false;
    } // Äußeres Polygon ist degeneriert
    if inner_poly_points.is_empty() {
        return true;
    } // Leeres inneres Polygon ist immer enthalten

    for &inner_point in inner_poly_points {
        if !point_in_polygon(inner_point, outer_poly_points, tolerance) {
            return false;
        }
    }
    true
}

/// Findet alle Schnittpunkte zwischen den Kanten zweier Polygone.
fn find_line_segment_intersections(poly_a: &[Vec2], poly_b: &[Vec2], tolerance: f32) -> Vec<Vec2> {
    let mut intersections = Vec::new();
    let n_a = poly_a.len();
    let n_b = poly_b.len();

    if n_a < 2 || n_b < 2 {
        return intersections;
    }

    for i in 0..n_a {
        let p1a = poly_a[i];
        let p2a = poly_a[(i + 1) % n_a];
        for j in 0..n_b {
            let p1b = poly_b[j];
            let p2b = poly_b[(j + 1) % n_b];

            // Verwende die Schnittpunktfunktion aus PolygonClipper oder eine äquivalente lokale.
            // Hier nehmen wir an, PolygonClipper::line_segment_intersection ist zugänglich
            // oder wir implementieren eine lokale Version.
            // Für dieses Beispiel eine lokale Kopie der Logik:
            let d1 = p2a - p1a;
            let d2 = p2b - p1b;
            let denominator = d1.x * d2.y - d1.y * d2.x;
            if denominator.abs() > tolerance * constants::EPSILON {
                let t_num = (p1b.x - p1a.x) * d2.y - (p1b.y - p1a.y) * d2.x;
                let u_num = (p1b.x - p1a.x) * d1.y - (p1b.y - p1a.y) * d1.x;
                let t = t_num / denominator;
                let u = u_num / denominator;
                if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
                    intersections.push(p1a + d1 * t);
                }
            }
        }
    }
    // Duplikate entfernen
    intersections.sort_by(|a, b| {
        a.x.partial_cmp(&b.x)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
    });
    intersections.dedup_by(|a, b| a.distance_squared(*b) < tolerance * tolerance);
    intersections
}
