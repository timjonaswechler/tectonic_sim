// src/math/geometry/polygon/core/polygon.rs

use crate::math::{error::*, types::*};
use std::fmt;

/// Hauptpolygon-Typ mit umfassender Funktionalität
#[derive(Debug, Clone, PartialEq)]
pub struct Polygon {
    pub vertices: Vec<Point2D>,
    pub is_closed: bool,
}

impl Polygon {
    /// Erstellt ein neues Polygon aus Vertices
    pub fn new(vertices: Vec<Point2D>) -> MathResult<Self> {
        Self::from_vertices(vertices, false)
    }

    /// Erstellt ein geschlossenes Polygon
    pub fn closed(vertices: Vec<Point2D>) -> MathResult<Self> {
        Self::from_vertices(vertices, true)
    }

    /// Erstellt Polygon mit Validierung
    fn from_vertices(mut vertices: Vec<Point2D>, force_closed: bool) -> MathResult<Self> {
        if vertices.len() < 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: vertices.len(),
            });
        }

        // Automatisch schließen wenn erwünscht und nicht bereits geschlossen
        let is_closed = if force_closed {
            if vertices.first() != vertices.last() {
                vertices.push(vertices[0]);
            }
            true
        } else {
            vertices.first() == vertices.last()
        };

        Ok(Self {
            vertices,
            is_closed,
        })
    }

    /// Zugriff auf Vertices
    pub fn vertices(&self) -> &[Point2D] {
        &self.vertices
    }

    /// Mutable Zugriff auf Vertices
    pub fn vertices_mut(&mut self) -> &mut Vec<Point2D> {
        &mut self.vertices
    }

    /// Anzahl der Vertices
    pub fn len(&self) -> usize {
        self.vertices.len()
    }

    /// Ist das Polygon leer?
    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }

    /// Ist das Polygon geschlossen?
    pub fn is_closed(&self) -> bool {
        self.is_closed
    }

    /// Polygon schließen
    pub fn close(&mut self) {
        if !self.is_closed && !self.vertices.is_empty() {
            if self.vertices.first() != self.vertices.last() {
                self.vertices.push(self.vertices[0]);
            }
            self.is_closed = true;
        }
    }

    /// Polygon öffnen (letzten Punkt entfernen wenn er dem ersten entspricht)
    pub fn open(&mut self) {
        if self.is_closed && self.vertices.len() > 1 {
            if self.vertices.first() == self.vertices.last() {
                self.vertices.pop();
            }
            self.is_closed = false;
        }
    }

    /// Bounding Box berechnen
    pub fn bounds(&self) -> Option<(Point2D, Point2D)> {
        if self.vertices.is_empty() {
            return None;
        }

        let mut min = self.vertices[0];
        let mut max = self.vertices[0];

        for vertex in &self.vertices[1..] {
            min.x = min.x.min(vertex.x);
            min.y = min.y.min(vertex.y);
            max.x = max.x.max(vertex.x);
            max.y = max.y.max(vertex.y);
        }

        Some((min, max))
    }

    /// Zentroid berechnen
    pub fn centroid(&self) -> Option<Point2D> {
        if self.vertices.is_empty() {
            return None;
        }

        let sum = self.vertices.iter().fold(Point2D::ZERO, |acc, v| acc + *v);
        Some(sum / self.vertices.len() as f32)
    }

    /// Polygon umkehren (Vertices in umgekehrter Reihenfolge)
    pub fn reverse(&mut self) {
        if self.is_closed && self.vertices.len() > 1 {
            // Bei geschlossenen Polygonen: ersten Punkt beibehalten, Rest umkehren
            let first = self.vertices[0];
            self.vertices[1..].reverse();
            if self.vertices.last() == Some(&first) {
                self.vertices.pop();
                self.vertices.push(first);
            }
        } else {
            self.vertices.reverse();
        }
    }

    /// Erstellt eine Kopie mit umgekehrten Vertices
    pub fn reversed(&self) -> Self {
        let mut copy = self.clone();
        copy.reverse();
        copy
    }

    /// Fügt einen Vertex hinzu
    pub fn add_vertex(&mut self, vertex: Point2D) -> MathResult<()> {
        if self.is_closed {
            // Bei geschlossenen Polygonen: vor dem letzten Punkt einfügen
            let last_idx = self.vertices.len() - 1;
            self.vertices.insert(last_idx, vertex);
        } else {
            self.vertices.push(vertex);
        }
        Ok(())
    }

    /// Entfernt einen Vertex an Index
    pub fn remove_vertex(&mut self, index: usize) -> MathResult<Point2D> {
        if index >= self.vertices.len() {
            return Err(MathError::InvalidConfiguration {
                message: format!("Index {} out of bounds", index),
            });
        }

        if self.vertices.len() <= 3 {
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: self.vertices.len() - 1,
            });
        }

        let removed = self.vertices.remove(index);

        // Bei geschlossenen Polygonen: Schließung beibehalten
        if self.is_closed && index == 0 && !self.vertices.is_empty() {
            *self.vertices.last_mut().unwrap() = self.vertices[0];
        }

        Ok(removed)
    }
}

/// Display-Implementierung für Debugging
impl fmt::Display for Polygon {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Polygon({} vertices", self.vertices.len())?;
        if self.is_closed {
            write!(f, ", closed")?;
        }
        write!(f, ")")?;
        Ok(())
    }
}

/// Konvertierung von Vec<Point2D>
impl TryFrom<Vec<Point2D>> for Polygon {
    type Error = MathError;

    fn try_from(vertices: Vec<Point2D>) -> Result<Self, Self::Error> {
        Self::new(vertices)
    }
}

/// Konvertierung zu Vec<Point2D>
impl From<Polygon> for Vec<Point2D> {
    fn from(polygon: Polygon) -> Self {
        polygon.vertices
    }
}

/// Iterator-Implementierung
impl IntoIterator for Polygon {
    type Item = Point2D;
    type IntoIter = std::vec::IntoIter<Point2D>;

    fn into_iter(self) -> Self::IntoIter {
        self.vertices.into_iter()
    }
}

impl<'a> IntoIterator for &'a Polygon {
    type Item = &'a Point2D;
    type IntoIter = std::slice::Iter<'a, Point2D>;

    fn into_iter(self) -> Self::IntoIter {
        self.vertices.iter()
    }
}
