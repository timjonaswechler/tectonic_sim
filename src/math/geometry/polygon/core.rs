// src/math/geometry/polygon/core.rs

use crate::math::{error::*, types::Bounds2D}; // Bounds2D wird für die bounds() Methode benötigt
use bevy::math::Vec2; // Wir verwenden Vec2 direkt
use std::fmt;

/// Polygon-Struktur, die eine Sequenz von 2D-Punkten (Vertices) darstellt.
#[derive(Debug, Clone, PartialEq)]
pub struct Polygon {
    /// Die Vertices, die das Polygon definieren.
    /// Für geschlossene Polygone: Der erste und letzte Vertex können identisch sein,
    /// aber die `is_closed`-Logik behandelt dies.
    pub vertices: Vec<Vec2>,
    /// Gibt an, ob das Polygon als geschlossen betrachtet wird (letzter Vertex mit erstem verbunden).
    is_closed_flag: bool, // Umbenannt, um Verwirrung mit der Methode zu vermeiden
}

impl Polygon {
    /// Erstellt ein neues Polygon aus einer Liste von Vertices.
    /// Das Polygon wird standardmäßig als offen betrachtet, es sei denn, der erste und letzte
    /// Vertex sind identisch (innerhalb einer kleinen Toleranz).
    pub fn new(vertices: Vec<Vec2>) -> MathResult<Self> {
        if vertices.len() < 2 && !vertices.is_empty() {
            // Ein einzelner Punkt ist kein Polygon, Linie braucht 2
            return Err(MathError::InsufficientPoints {
                expected: 2, // Mindestens 2 für eine offene Linie, 3 für eine Fläche
                actual: vertices.len(),
            });
        }
        // Ein Polygon, das eine Fläche umschließt, braucht mind. 3 unterschiedliche Punkte
        // Diese Prüfung erfolgt später bei Operationen wie area() oder wenn explizit `closed` aufgerufen wird.

        let mut is_closed_flag = false;
        if vertices.len() >= 3 {
            if let (Some(first), Some(last)) = (vertices.first(), vertices.last()) {
                if first.distance_squared(*last)
                    < crate::math::utils::constants::EPSILON
                        * crate::math::utils::constants::EPSILON
                {
                    is_closed_flag = true;
                }
            }
        }

        Ok(Self {
            vertices,
            is_closed_flag,
        })
    }

    /// Erstellt ein Polygon und markiert es explizit als geschlossen.
    /// Wenn der erste und letzte Vertex nicht identisch sind, wird der erste Vertex am Ende angefügt.
    /// Benötigt mindestens 3 eindeutige Vertices, um eine geschlossene Fläche zu bilden.
    pub fn closed(mut vertices: Vec<Vec2>) -> MathResult<Self> {
        if vertices.len() < 3 {
            // Überprüfe auf eindeutige Punkte, falls nur Duplikate vorhanden sind
            let mut unique_points = vertices.clone();
            unique_points.sort_by(|a, b| {
                a.x.partial_cmp(&b.x)
                    .unwrap()
                    .then(a.y.partial_cmp(&b.y).unwrap())
            });
            unique_points.dedup_by(|a, b| a.distance_squared(*b) < 1e-9);

            if unique_points.len() < 3 {
                return Err(MathError::InsufficientPoints {
                    expected: 3, // Mindestens 3 für eine geschlossene Fläche
                    actual: unique_points.len(),
                });
            }
        }

        if let (Some(&first), Some(&last)) = (vertices.first(), vertices.last()) {
            if first.distance_squared(last)
                > crate::math::utils::constants::EPSILON * crate::math::utils::constants::EPSILON
            {
                vertices.push(first);
            }
        } else if vertices.is_empty() {
            // Sollte durch die obige Längenprüfung nicht passieren, aber sicher ist sicher
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: 0,
            });
        }

        Ok(Self {
            vertices,
            is_closed_flag: true,
        })
    }

    /// Gibt einen Slice der Vertices zurück.
    pub fn vertices(&self) -> &[Vec2] {
        &self.vertices
    }

    /// Gibt einen mutable Slice der Vertices zurück.
    /// Vorsicht: Direkte Modifikation kann die `is_closed_flag`-Konsistenz stören.
    /// Es wird empfohlen, `add_vertex`, `remove_vertex`, `close`, `open` zu verwenden.
    pub fn vertices_mut(&mut self) -> &mut Vec<Vec2> {
        // Wenn Vertices geändert werden, ist der is_closed_flag möglicherweise nicht mehr korrekt.
        // Man könnte hier den Flag auf false setzen oder eine Neuberechnung erzwingen.
        // Fürs Erste belassen wir es dabei, dass der Nutzer vorsichtig sein muss.
        &mut self.vertices
    }

    /// Anzahl der Vertices (inklusive des duplizierten Endpunktes bei geschlossenen Polygonen,
    /// wenn er explizit gespeichert ist).
    pub fn len(&self) -> usize {
        self.vertices.len()
    }

    /// Prüft, ob das Polygon keine Vertices hat.
    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }

    /// Prüft, ob das Polygon als geschlossen markiert ist.
    pub fn is_closed(&self) -> bool {
        self.is_closed_flag
    }

    /// Schließt das Polygon.
    /// Wenn es bereits als geschlossen markiert ist oder weniger als 3 Vertices hat, passiert nichts.
    /// Fügt den ersten Vertex am Ende an, falls noch nicht geschehen.
    pub fn close(&mut self) -> MathResult<()> {
        if self.is_closed_flag || self.vertices.len() < 3 {
            if self.vertices.len() < 3 && !self.vertices.is_empty() {
                return Err(MathError::InsufficientPoints {
                    expected: 3,
                    actual: self.vertices.len(),
                });
            }
            return Ok(()); // Bereits geschlossen oder nicht genug Punkte für eine Fläche
        }

        if let (Some(&first), Some(&last)) = (self.vertices.first(), self.vertices.last()) {
            if first.distance_squared(last)
                > crate::math::utils::constants::EPSILON * crate::math::utils::constants::EPSILON
            {
                self.vertices.push(first);
            }
        } else if self.vertices.is_empty() {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: 0,
            });
        }
        self.is_closed_flag = true;
        Ok(())
    }

    /// Öffnet das Polygon.
    /// Wenn es als geschlossen markiert war und der letzte Vertex dem ersten gleicht, wird der letzte entfernt.
    pub fn open(&mut self) {
        if self.is_closed_flag && !self.vertices.is_empty() {
            if let (Some(&first), Some(&last)) = (self.vertices.first(), self.vertices.last()) {
                if first.distance_squared(last)
                    < crate::math::utils::constants::EPSILON
                        * crate::math::utils::constants::EPSILON
                    && self.vertices.len() > 1
                {
                    // Nur entfernen, wenn es tatsächlich ein Duplikat ist und mehr als ein Punkt existiert
                    // (um zu vermeiden, dass ein Polygon aus einem einzelnen Punkt leer wird)
                    if self.vertices.len() > 1 {
                        // Sicherstellen, dass wir nicht den einzigen Punkt entfernen
                        self.vertices.pop();
                    }
                }
            }
        }
        self.is_closed_flag = false;
    }

    /// Berechnet die Bounding Box des Polygons.
    /// Gibt `None` zurück, wenn das Polygon leer ist.
    pub fn bounds(&self) -> Option<Bounds2D> {
        Bounds2D::from_points_iter(self.vertices.iter().copied())
    }

    /// Berechnet den arithmetischen Mittelpunkt der Vertices (nicht den geometrischen Schwerpunkt).
    pub fn centroid(&self) -> Option<Vec2> {
        if self.vertices.is_empty() {
            return None;
        }
        let sum = self.vertices.iter().fold(Vec2::ZERO, |acc, &v| acc + v);
        Some(sum / self.vertices.len() as f32)
    }

    /// Kehrt die Reihenfolge der Vertices um.
    /// Bei geschlossenen Polygonen wird der (ggf. duplizierte) Endpunkt korrekt behandelt.
    pub fn reverse(&mut self) {
        if self.vertices.is_empty() {
            return;
        }

        if self.is_closed_flag && self.vertices.len() > 1 {
            let first_is_last = self
                .vertices
                .first()
                .unwrap()
                .distance_squared(*self.vertices.last().unwrap())
                < 1e-9;
            if first_is_last {
                let first = self.vertices.remove(0); // Ersten entfernen
                self.vertices.reverse(); // Rest umkehren
                self.vertices.insert(0, first); // Ersten wieder an den Anfang
                // Der neue letzte Punkt ist jetzt der alte vorletzte.
                // Wir müssen sicherstellen, dass der neue erste Punkt auch der neue letzte ist.
                let new_first = self.vertices[0];
                *self.vertices.last_mut().unwrap() = new_first;
            } else {
                // Geschlossen, aber erster und letzter Punkt nicht identisch gespeichert
                self.vertices.reverse();
            }
        } else {
            self.vertices.reverse();
        }
    }

    /// Erstellt eine Kopie des Polygons mit umgekehrter Vertex-Reihenfolge.
    pub fn reversed(&self) -> Self {
        let mut copy = self.clone();
        copy.reverse();
        copy
    }

    /// Fügt einen Vertex hinzu.
    /// Bei geschlossenen Polygonen wird der Vertex vor dem (ggf. duplizierten) Endpunkt eingefügt.
    pub fn add_vertex(&mut self, vertex: Vec2) {
        if self.is_closed_flag && !self.vertices.is_empty() {
            // Wenn der erste und letzte Punkt identisch sind, füge vor dem letzten ein.
            let insert_pos =
                if self.vertices.first() == self.vertices.last() && self.vertices.len() > 1 {
                    self.vertices.len() - 1
                } else {
                    self.vertices.len()
                };
            self.vertices.insert(insert_pos, vertex);
        } else {
            self.vertices.push(vertex);
        }
    }

    /// Entfernt einen Vertex am gegebenen Index.
    /// Gibt den entfernten Vertex zurück oder einen Fehler, wenn der Index ungültig ist
    /// oder das Entfernen zu weniger als 3 Vertices für ein geschlossenes Polygon führen würde.
    pub fn remove_vertex(&mut self, index: usize) -> MathResult<Vec2> {
        if index >= self.vertices.len() {
            return Err(MathError::InvalidConfiguration {
                message: format!(
                    "Index {} out of bounds for polygon with {} vertices",
                    index,
                    self.vertices.len()
                ),
            });
        }

        // Mindestanzahl an Punkten für ein Polygon (offen: 2, geschlossen: 3 einzigartige)
        let min_points_for_shape = if self.is_closed_flag { 3 } else { 2 };
        let mut _unique_count_after_remove = self.vertices.len() - 1;
        if self.is_closed_flag
            && self.vertices.first() == self.vertices.last()
            && self.vertices.len() > 1
        {
            _unique_count_after_remove -= 1; // Wenn erster und letzter identisch sind, zählt das Entfernen eines "echten" Punkts doppelt
        }

        if self.vertices.len() - 1 < min_points_for_shape {
            // Einfache Längenprüfung
            return Err(MathError::InsufficientPoints {
                expected: min_points_for_shape,
                actual: self.vertices.len() - 1,
            });
        }

        let removed = self.vertices.remove(index);

        // Konsistenz für geschlossene Polygone wiederherstellen
        if self.is_closed_flag && !self.vertices.is_empty() {
            // Wenn der entfernte Punkt der erste oder letzte war (und sie identisch waren)
            // oder wenn der erste Punkt entfernt wurde, muss der letzte ggf. aktualisiert werden.
            if self.vertices.len() >= 1 && self.vertices.first() != self.vertices.last() {
                if self.vertices.len() >= min_points_for_shape {
                    // Nur schließen, wenn genug Punkte da sind
                    let first = self.vertices[0];
                    *self.vertices.last_mut().unwrap() = first;
                } else {
                    self.is_closed_flag = false; // Kann nicht mehr geschlossen sein
                }
            } else if self.vertices.is_empty() || self.vertices.len() < min_points_for_shape {
                self.is_closed_flag = false; // Nicht mehr genug Punkte für ein geschlossenes Polygon
            }
        }
        Ok(removed)
    }
}

/// Display-Implementierung für Debugging.
impl fmt::Display for Polygon {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Polygon({} vertices", self.vertices.len())?;
        if self.is_closed_flag {
            write!(f, ", closed")?;
        } else {
            write!(f, ", open")?;
        }
        write!(f, ")")?;
        Ok(())
    }
}

/// Konvertierung von `Vec<Vec2>` zu `Polygon`.
impl TryFrom<Vec<Vec2>> for Polygon {
    type Error = MathError;

    fn try_from(vertices: Vec<Vec2>) -> Result<Self, Self::Error> {
        Self::new(vertices)
    }
}

/// Konvertierung von `Polygon` zu `Vec<Vec2>`.
impl From<Polygon> for Vec<Vec2> {
    fn from(polygon: Polygon) -> Self {
        polygon.vertices
    }
}

/// Iterator-Implementierung für das Polygon (über seine Vertices).
impl IntoIterator for Polygon {
    type Item = Vec2;
    type IntoIter = std::vec::IntoIter<Vec2>;

    fn into_iter(self) -> Self::IntoIter {
        self.vertices.into_iter()
    }
}

impl<'a> IntoIterator for &'a Polygon {
    type Item = &'a Vec2;
    type IntoIter = std::slice::Iter<'a, Vec2>;

    fn into_iter(self) -> Self::IntoIter {
        self.vertices.iter()
    }
}
