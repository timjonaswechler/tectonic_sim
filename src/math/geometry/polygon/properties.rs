// src/math/geometry/polygon/properties.rs

use crate::math::geometry::polygon::core::Polygon; // Zugriff auf die Polygon-Struktur
use crate::math::utils::constants; // Für EPSILON
use bevy::math::Vec2;

/// Trait für fortgeschrittene geometrische Eigenschaften von Polygonen.
pub trait PolygonProperties {
    /// Berechnet die Fläche des Polygons unter Verwendung der Shoelace-Formel.
    /// Berücksichtigt die `is_closed()`-Eigenschaft des Polygons.
    fn area(&self) -> f32;

    /// Berechnet den Umfang des Polygons.
    /// Berücksichtigt die `is_closed()`-Eigenschaft.
    fn perimeter(&self) -> f32;

    /// Prüft, ob ein Punkt innerhalb des Polygons liegt (Ray-Casting-Algorithmus).
    /// Funktioniert für einfache Polygone (nicht selbst-überschneidend).
    fn contains_point(&self, point: Vec2) -> bool;

    /// Prüft, ob das Polygon konvex ist.
    /// Ein Polygon ist konvex, wenn alle inneren Winkel kleiner als 180 Grad sind.
    fn is_convex(&self) -> bool;

    /// Bestimmt die Orientierung des Polygons (im Uhrzeigersinn oder gegen den Uhrzeigersinn).
    /// Basiert auf dem Vorzeichen der (doppelten) Fläche.
    fn orientation(&self) -> Orientation;

    /// Berechnet den geometrischen Schwerpunkt (Zentroid) des Polygons.
    /// Dies ist der "Massenmittelpunkt", wenn das Polygon eine homogene Dichte hätte.
    /// Gibt `None` zurück, wenn das Polygon weniger als 3 Vertices hat oder die Fläche null ist.
    fn geometric_centroid(&self) -> Option<Vec2>;
}

/// Gibt die Orientierung eines Polygons an.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Clockwise,
    CounterClockwise,
    Collinear, // Alle Punkte liegen auf einer Linie
}

impl PolygonProperties for Polygon {
    fn area(&self) -> f32 {
        let vertices = self.vertices();
        let n = vertices.len();

        if n < 3 {
            return 0.0;
        }

        let mut area_sum = 0.0;

        // Wenn das Polygon als geschlossen markiert ist und der erste und letzte Punkt
        // identisch sind, müssen wir den letzten Punkt für die Shoelace-Formel nicht
        // explizit doppelt nehmen, da die Schleife dies korrekt handhabt.
        // Wenn es offen ist, oder geschlossen aber mit explizit dupliziertem Endpunkt,
        // ist die Schleife bis n-1 korrekt, wenn der letzte Punkt der Duplikat ist.
        // Die sauberste Methode ist, immer über die "effektiven" Kanten zu iterieren.

        let num_edges = if self.is_closed() {
            // Wenn explizit geschlossen und erster==letzter, dann n-1 Kanten.
            // Wenn geschlossen, aber erster!=letzter (was `close()` verhindern sollte, aber zur Sicherheit), dann n Kanten.
            if n > 0 && vertices.first() == vertices.last() {
                n - 1
            } else {
                n
            }
        } else {
            n - 1 // Offenes Polygon hat n-1 Kanten
        };
        if num_edges < 2 && n >= 3 {
            // Sonderfall: Polygon hat 3 Punkte, aber ist offen, nur 2 Kanten, Shoelace geht nicht direkt
            if n == 3 && !self.is_closed() {
                // Dreieck als offene Linienkette
                return 0.0; // Fläche ist null für eine offene Linienkette
            }
        }
        if num_edges < 1 {
            return 0.0;
        }

        for i in 0..num_edges {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % n]; // %n stellt sicher, dass wir für die letzte Kante zum ersten Punkt zurückkehren, falls nötig
            area_sum += (p1.x * p2.y) - (p2.x * p1.y);
        }

        // Für ein offenes Polygon wollen wir typischerweise keine Fläche, es sei denn,
        // es wird implizit als geschlossen für die Flächenberechnung angenommen.
        // Die Standard-Shoelace-Formel gilt für geschlossene Polygone.
        // Wenn das Polygon offen ist und trotzdem eine Fläche berechnet werden soll,
        // müsste es als temporär geschlossen behandelt werden.
        // Hier nehmen wir an, area() für ein offenes Polygon ist 0, es sei denn, es sind
        // genug Punkte da und es WURDE als geschlossen markiert.
        if !self.is_closed() && n < 3 {
            // Ein offenes Polygon mit <3 Punkten hat keine Fläche
            return 0.0;
        }
        // Wenn es offen ist, aber 3+ Punkte hat, könnten wir die Fläche des umschlossenen Bereichs berechnen,
        // indem wir es temporär schließen. Das aktuelle `is_closed()` Flag ist hier entscheidend.
        // Die obige Logik für `num_edges` und die Schleife `0..num_edges` mit `(i+1)%n`
        // sollte die Shoelace-Formel für geschlossene Polygone (auch mit dupliziertem Endpunkt) korrekt handhaben.
        // Für offene Polygone mit n Kanten (also n+1 Punkten) ist das Ergebnis der Shoelace nicht die "Fläche".

        (area_sum * 0.5).abs()
    }

    fn perimeter(&self) -> f32 {
        let vertices = self.vertices();
        let n = vertices.len();

        if n < 2 {
            return 0.0;
        }

        let mut perimeter_sum = 0.0;
        let num_segments = if self.is_closed() {
            // Wenn explizit geschlossen und erster==letzter, dann n-1 Segmente.
            // Ansonsten n Segmente (impliziert, dass das Polygon sich selbst schließt).
            if n > 0 && vertices.first() == vertices.last() {
                n - 1
            } else {
                n
            }
        } else {
            n - 1 // Offenes Polygon hat n-1 Segmente
        };

        if num_segments == 0 && n == 1 {
            return 0.0;
        } // Einzelner Punkt

        for i in 0..num_segments {
            perimeter_sum += vertices[i].distance(vertices[(i + 1) % n]);
        }
        perimeter_sum
    }

    fn contains_point(&self, point: Vec2) -> bool {
        let vertices = self.vertices();
        let n = vertices.len();

        if n < 3 {
            return false; // Ein Polygon mit weniger als 3 Punkten kann keinen Punkt enthalten (im Sinne einer Fläche)
        }

        let mut inside = false;

        // Nutze die effektive Anzahl von Kanten/Punkten für die Schleife
        let mut j = if self.is_closed() && n > 0 && vertices.first() == vertices.last() {
            n - 2 // vorletzter Index, da der letzte ein Duplikat ist
        } else {
            n - 1 // letzter Index
        };

        let num_unique_points = if self.is_closed() && n > 0 && vertices.first() == vertices.last()
        {
            n - 1
        } else {
            n
        };
        if num_unique_points < 3 {
            return false;
        }

        for i in 0..num_unique_points {
            let vi = vertices[i];
            let vj = vertices[j];

            // Ray-Casting: Prüfe, ob der Strahl von `point` nach rechts die Kante (vi, vj) schneidet
            let intersect = ((vi.y > point.y) != (vj.y > point.y))
                && (point.x
                    < (vj.x - vi.x) * (point.y - vi.y) / (vj.y - vi.y + constants::EPSILON) + vi.x);
            // EPSILON im Nenner zur Vermeidung von DivByZero bei horizontalen Kanten

            if intersect {
                inside = !inside;
            }
            j = i; // j ist der vorherige Vertex in der nächsten Iteration
        }
        inside
    }

    fn is_convex(&self) -> bool {
        let vertices = self.vertices();
        let n = vertices.len();

        if n < 3 {
            return false; // Nicht genug Punkte für Konvexitätsprüfung (oder sogar für ein Polygon)
        }

        // Anzahl der "einzigartigen" Punkte für die Prüfung
        let num_unique_points = if self.is_closed() && n > 0 && vertices.first() == vertices.last()
        {
            n - 1
        } else {
            n
        };

        if num_unique_points < 3 {
            return false;
        }

        let mut sign_of_cross_product: Option<bool> = None;

        for i in 0..num_unique_points {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % num_unique_points]; // % num_unique_points für korrekten Umlauf
            let p3 = vertices[(i + 2) % num_unique_points];

            // Kreuzprodukt (p2-p1) x (p3-p2)
            // (x2-x1)*(y3-y2) - (y2-y1)*(x3-x2)
            let cross_product_z = (p2.x - p1.x) * (p3.y - p2.y) - (p2.y - p1.y) * (p3.x - p2.x);

            if cross_product_z.abs() > constants::EPSILON {
                // Nur nicht-kollineare Punkte betrachten
                let current_sign = cross_product_z > 0.0;
                if let Some(expected_sign) = sign_of_cross_product {
                    if current_sign != expected_sign {
                        return false; // Zeichenwechsel deutet auf Nicht-Konvexität hin
                    }
                } else {
                    sign_of_cross_product = Some(current_sign);
                }
            }
        }
        // Wenn alle Kreuzprodukte dasselbe Vorzeichen haben (oder null sind), ist das Polygon konvex.
        // Wenn keine nicht-kollinearen Tripel gefunden wurden (z.B. Linie), ist es auch nicht konvex im Sinne einer Fläche.
        sign_of_cross_product.is_some()
    }

    fn orientation(&self) -> Orientation {
        let vertices = self.vertices();
        let n = vertices.len();

        if n < 3 {
            return Orientation::Collinear;
        }

        // Die `area()` Methode gibt den Absolutwert. Wir brauchen hier das vorzeichenbehaftete Ergebnis der Shoelace-Formel.
        let mut signed_area_doubled = 0.0;
        let num_edges = if self.is_closed() && n > 0 && vertices.first() == vertices.last() {
            n - 1
        } else {
            n
        };
        if num_edges < 2 && n >= 3 {
            if n == 3 && !self.is_closed() {
                return Orientation::Collinear;
            }
        }
        if num_edges < 1 {
            return Orientation::Collinear;
        }

        for i in 0..num_edges {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % n];
            signed_area_doubled += (p1.x * p2.y) - (p2.x * p1.y);
        }

        if !self.is_closed() && n < 3 {
            return Orientation::Collinear;
        }

        if signed_area_doubled.abs() < constants::EPSILON * n as f32 {
            // Toleranz basierend auf Anzahl der Punkte
            Orientation::Collinear
        } else if signed_area_doubled > 0.0 {
            Orientation::CounterClockwise // Standardorientierung für positive Fläche
        } else {
            Orientation::Clockwise
        }
    }

    fn geometric_centroid(&self) -> Option<Vec2> {
        let vertices = self.vertices();
        let n = vertices.len();

        if n < 3 {
            return None;
        }

        // Die Berechnung des geometrischen Schwerpunkts benötigt die vorzeichenbehaftete Fläche.
        let mut signed_area_doubled = 0.0;
        let num_edges = if self.is_closed() && n > 0 && vertices.first() == vertices.last() {
            n - 1
        } else {
            n
        };
        if num_edges < 2 && n >= 3 {
            if n == 3 && !self.is_closed() {
                return self.centroid(); // Fallback für offene Linienkette
            }
        }
        if num_edges < 1 {
            return self.centroid();
        }

        for i in 0..num_edges {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % n];
            signed_area_doubled += (p1.x * p2.y) - (p2.x * p1.y);
        }

        if signed_area_doubled.abs() < constants::EPSILON * n as f32 {
            return self.centroid(); // Fallback bei keiner Fläche (Linie)
        }

        if !self.is_closed() && n < 3 {
            return self.centroid();
        }

        let mut centroid_x = 0.0;
        let mut centroid_y = 0.0;

        for i in 0..num_edges {
            let p1 = vertices[i];
            let p2 = vertices[(i + 1) % n];
            let factor = (p1.x * p2.y) - (p2.x * p1.y);
            centroid_x += (p1.x + p2.x) * factor;
            centroid_y += (p1.y + p2.y) * factor;
        }

        let inv_area_factor = 1.0 / (3.0 * signed_area_doubled); // Eigentlich 1.0 / (6.0 * Area)
        Some(Vec2::new(
            centroid_x * inv_area_factor,
            centroid_y * inv_area_factor,
        ))
    }
}
