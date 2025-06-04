// src/math/geometry/metalballs/marching_squares.rs

use super::MetaballField;
use bevy::prelude::*;

/// Eine Kontur-Linie extrahiert aus dem Skalarfeld
#[derive(Debug, Clone)]
pub struct Contour {
    pub vertices: Vec<Vec2>,
    pub is_closed: bool,
}

impl Contour {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            is_closed: false,
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            vertices: Vec::with_capacity(capacity),
            is_closed: false,
        }
    }

    pub fn add_vertex(&mut self, vertex: Vec2) {
        self.vertices.push(vertex);
    }

    pub fn close(&mut self) {
        self.is_closed = true;
    }

    pub fn len(&self) -> usize {
        self.vertices.len()
    }

    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }
}

/// Marching Squares Algorithmus Iterator
pub struct MarchingSquaresIterator<'a> {
    field: &'a MetaballField,
    threshold: f32,
    current_x: usize,
    current_y: usize,
    visited: Vec<Vec<bool>>,
    edge_table: [Vec<(f32, f32)>; 16],
}

impl<'a> MarchingSquaresIterator<'a> {
    pub fn new(field: &'a MetaballField, threshold: f32) -> Self {
        let visited =
            vec![vec![false; field.height.saturating_sub(1)]; field.width.saturating_sub(1)];

        let mut iterator = Self {
            field,
            threshold,
            current_x: 0,
            current_y: 0,
            visited,
            edge_table: Default::default(),
        };

        iterator.initialize_edge_table();
        iterator
    }

    /// Initialisiert die Lookup-Tabelle für Marching Squares
    fn initialize_edge_table(&mut self) {
        // Die 16 möglichen Konfigurationen für Marching Squares
        // Jede Konfiguration beschreibt, welche Kanten des Quadrats durchschnitten werden

        self.edge_table[0] = vec![]; // 0000 - keine Kanten
        self.edge_table[1] = vec![(0.0, 0.5), (0.5, 1.0)]; // 0001 - unten links
        self.edge_table[2] = vec![(0.5, 1.0), (1.0, 0.5)]; // 0010 - unten rechts
        self.edge_table[3] = vec![(0.0, 0.5), (1.0, 0.5)]; // 0011 - unten
        self.edge_table[4] = vec![(1.0, 0.5), (0.5, 0.0)]; // 0100 - oben rechts
        self.edge_table[5] = vec![(0.0, 0.5), (0.5, 0.0), (1.0, 0.5), (0.5, 1.0)]; // 0101 - diagonal (ambiguous)
        self.edge_table[6] = vec![(0.5, 1.0), (0.5, 0.0)]; // 0110 - rechts
        self.edge_table[7] = vec![(0.0, 0.5), (0.5, 0.0)]; // 0111 - nicht unten links
        self.edge_table[8] = vec![(0.5, 0.0), (0.0, 0.5)]; // 1000 - oben links
        self.edge_table[9] = vec![(0.5, 0.0), (0.5, 1.0)]; // 1001 - links
        self.edge_table[10] = vec![(0.5, 0.0), (0.0, 0.5), (0.5, 1.0), (1.0, 0.5)]; // 1010 - diagonal (ambiguous)
        self.edge_table[11] = vec![(0.5, 0.0), (1.0, 0.5)]; // 1011 - nicht unten rechts
        self.edge_table[12] = vec![(1.0, 0.5), (0.0, 0.5)]; // 1100 - oben
        self.edge_table[13] = vec![(0.5, 1.0), (1.0, 0.5)]; // 1101 - nicht oben links
        self.edge_table[14] = vec![(0.5, 1.0), (0.0, 0.5)]; // 1110 - nicht oben rechts
        self.edge_table[15] = vec![]; // 1111 - keine Kanten
    }

    /// Berechnet den Konfigurations-Index für eine Zelle
    fn get_configuration(&self, x: usize, y: usize) -> usize {
        let mut config = 0;

        // Prüfe die vier Eckpunkte der Zelle
        if x < self.field.width && y < self.field.height && self.field.get(x, y) >= self.threshold {
            config |= 8; // oben links
        }
        if x + 1 < self.field.width
            && y < self.field.height
            && self.field.get(x + 1, y) >= self.threshold
        {
            config |= 4; // oben rechts
        }
        if x + 1 < self.field.width
            && y + 1 < self.field.height
            && self.field.get(x + 1, y + 1) >= self.threshold
        {
            config |= 2; // unten rechts
        }
        if x < self.field.width
            && y + 1 < self.field.height
            && self.field.get(x, y + 1) >= self.threshold
        {
            config |= 1; // unten links
        }

        config
    }

    /// Konvertiert relative Koordinaten zu absoluten Weltkoordinaten
    fn to_world_coords(&self, x: usize, y: usize, rel_x: f32, rel_y: f32) -> Vec2 {
        Vec2::new(
            (x as f32 + rel_x) * self.field.cell_size,
            (y as f32 + rel_y) * self.field.cell_size,
        )
    }

    /// Interpoliert zwischen zwei Werten basierend auf dem Threshold
    fn interpolate_edge(&self, val1: f32, val2: f32) -> f32 {
        if (val2 - val1).abs() < 1e-6 {
            return 0.5;
        }

        (self.threshold - val1) / (val2 - val1)
    }

    /// Berechnet die tatsächlichen Schnittpunkte für eine Konfiguration
    fn get_edge_intersections(&self, x: usize, y: usize, config: usize) -> Vec<Vec2> {
        let mut intersections = Vec::new();

        if config == 0 || config == 15 {
            return intersections;
        }

        // Hole die Werte an den Eckpunkten
        let val_tl = if x < self.field.width && y < self.field.height {
            self.field.get(x, y)
        } else {
            0.0
        };
        let val_tr = if x + 1 < self.field.width && y < self.field.height {
            self.field.get(x + 1, y)
        } else {
            0.0
        };
        let val_br = if x + 1 < self.field.width && y + 1 < self.field.height {
            self.field.get(x + 1, y + 1)
        } else {
            0.0
        };
        let val_bl = if x < self.field.width && y + 1 < self.field.height {
            self.field.get(x, y + 1)
        } else {
            0.0
        };

        // Berechne Schnittpunkte für jede Kante
        let edges = &self.edge_table[config];

        for i in (0..edges.len()).step_by(2) {
            if i + 1 >= edges.len() {
                break;
            }

            let (start_x, start_y) = edges[i];
            let (end_x, end_y) = edges[i + 1];

            // Bestimme welche Kante geschnitten wird und interpoliere
            let start_pos = self
                .interpolate_edge_position(x, y, start_x, start_y, val_tl, val_tr, val_br, val_bl);
            let end_pos =
                self.interpolate_edge_position(x, y, end_x, end_y, val_tl, val_tr, val_br, val_bl);

            intersections.push(start_pos);
            intersections.push(end_pos);
        }

        intersections
    }

    /// Interpoliert die Position auf einer Kante
    fn interpolate_edge_position(
        &self,
        x: usize,
        y: usize,
        rel_x: f32,
        rel_y: f32,
        val_tl: f32,
        val_tr: f32,
        val_br: f32,
        val_bl: f32,
    ) -> Vec2 {
        let mut world_x = (x as f32 + rel_x) * self.field.cell_size;
        let mut world_y = (y as f32 + rel_y) * self.field.cell_size;

        // Lineare Interpolation basierend auf der Position
        if rel_y == 0.0 {
            // Obere Kante
            let t = self.interpolate_edge(val_tl, val_tr);
            world_x = (x as f32 + t) * self.field.cell_size;
        } else if rel_x == 1.0 {
            // Rechte Kante
            let t = self.interpolate_edge(val_tr, val_br);
            world_y = (y as f32 + t) * self.field.cell_size;
        } else if rel_y == 1.0 {
            // Untere Kante
            let t = self.interpolate_edge(val_bl, val_br);
            world_x = (x as f32 + t) * self.field.cell_size;
        } else if rel_x == 0.0 {
            // Linke Kante
            let t = self.interpolate_edge(val_tl, val_bl);
            world_y = (y as f32 + t) * self.field.cell_size;
        }

        Vec2::new(world_x, world_y)
    }

    /// Findet die nächste unbesuchte Zelle mit Kontur
    fn find_next_contour_start(&mut self) -> Option<(usize, usize)> {
        for y in self.current_y..self.field.height.saturating_sub(1) {
            let start_x = if y == self.current_y {
                self.current_x
            } else {
                0
            };

            for x in start_x..self.field.width.saturating_sub(1) {
                if !self.visited[x][y] {
                    let config = self.get_configuration(x, y);
                    if config != 0 && config != 15 {
                        self.current_x = x;
                        self.current_y = y;
                        return Some((x, y));
                    }
                    self.visited[x][y] = true;
                }
            }
            self.current_x = 0;
        }
        None
    }

    /// Verfolgt eine Kontur von einem Startpunkt aus
    fn trace_contour(&mut self, start_x: usize, start_y: usize) -> Contour {
        let mut contour = Contour::new();
        let mut current_x = start_x;
        let mut current_y = start_y;

        loop {
            if current_x >= self.field.width.saturating_sub(1)
                || current_y >= self.field.height.saturating_sub(1)
            {
                break;
            }

            if self.visited[current_x][current_y] {
                break;
            }

            self.visited[current_x][current_y] = true;

            let config = self.get_configuration(current_x, current_y);
            if config == 0 || config == 15 {
                break;
            }

            let intersections = self.get_edge_intersections(current_x, current_y, config);

            for vertex in intersections {
                contour.add_vertex(vertex);
            }

            // Finde nächste Zelle (vereinfachte Implementierung)
            let mut found_next = false;
            for dy in -1i32..=1i32 {
                for dx in -1i32..=1i32 {
                    if dx == 0 && dy == 0 {
                        continue;
                    }

                    let next_x = current_x as i32 + dx;
                    let next_y = current_y as i32 + dy;

                    if next_x >= 0 && next_y >= 0 {
                        let next_x = next_x as usize;
                        let next_y = next_y as usize;

                        if next_x < self.field.width.saturating_sub(1)
                            && next_y < self.field.height.saturating_sub(1)
                            && !self.visited[next_x][next_y]
                        {
                            let next_config = self.get_configuration(next_x, next_y);
                            if next_config != 0 && next_config != 15 {
                                current_x = next_x;
                                current_y = next_y;
                                found_next = true;
                                break;
                            }
                        }
                    }
                }
                if found_next {
                    break;
                }
            }

            if !found_next {
                break;
            }
        }

        // Prüfe ob Kontur geschlossen werden kann
        if contour.len() >= 3 {
            let first = contour.vertices[0];
            let last = contour.vertices[contour.len() - 1];
            if first.distance(last) < self.field.cell_size * 1.5 {
                contour.close();
            }
        }

        contour
    }
}

impl<'a> Iterator for MarchingSquaresIterator<'a> {
    type Item = Contour;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((start_x, start_y)) = self.find_next_contour_start() {
            let contour = self.trace_contour(start_x, start_y);
            if !contour.is_empty() {
                return Some(contour);
            }
        }
        None
    }
}

/// Utility-Funktionen für Marching Squares
pub struct MarchingSquares;

impl MarchingSquares {
    /// Extrahiert alle Konturen aus einem Skalarfeld
    pub fn extract_contours(field: &MetaballField, threshold: f32) -> Vec<Contour> {
        let iterator = MarchingSquaresIterator::new(field, threshold);
        iterator.collect()
    }

    /// Extrahiert die erste/größte Kontur
    pub fn extract_main_contour(field: &MetaballField, threshold: f32) -> Option<Contour> {
        let contours = Self::extract_contours(field, threshold);
        contours.into_iter().max_by_key(|c| c.len())
    }

    /// Vereinfacht eine Kontur durch Entfernung naher Punkte
    pub fn simplify_contour(contour: &Contour, tolerance: f32) -> Contour {
        if contour.vertices.len() <= 2 {
            return contour.clone();
        }

        let mut simplified = Contour::new();
        simplified.add_vertex(contour.vertices[0]);

        for i in 1..contour.vertices.len() {
            let current = contour.vertices[i];
            let last = simplified.vertices[simplified.vertices.len() - 1];

            if current.distance(last) >= tolerance {
                simplified.add_vertex(current);
            }
        }

        if contour.is_closed {
            simplified.close();
        }

        simplified
    }

    /// Glättet eine Kontur mit einfachem Moving Average
    pub fn smooth_contour(contour: &Contour, window_size: usize) -> Contour {
        if contour.vertices.len() <= window_size || window_size < 2 {
            return contour.clone();
        }

        let mut smoothed = Contour::new();
        let half_window = window_size / 2;

        for i in 0..contour.vertices.len() {
            let mut sum = Vec2::ZERO;
            let mut count = 0;

            for j in 0..window_size {
                let idx = if contour.is_closed {
                    (i + contour.vertices.len() - half_window + j) % contour.vertices.len()
                } else {
                    (i + j)
                        .saturating_sub(half_window)
                        .min(contour.vertices.len() - 1)
                };

                sum += contour.vertices[idx];
                count += 1;
            }

            smoothed.add_vertex(sum / count as f32);
        }

        if contour.is_closed {
            smoothed.close();
        }

        smoothed
    }

    /// Berechnet die Fläche einer geschlossenen Kontur (Shoelace-Formel)
    pub fn calculate_contour_area(contour: &Contour) -> f32 {
        if !contour.is_closed || contour.vertices.len() < 3 {
            return 0.0;
        }

        let mut area = 0.0;
        let n = contour.vertices.len();

        for i in 0..n {
            let j = (i + 1) % n;
            area += contour.vertices[i].x * contour.vertices[j].y;
            area -= contour.vertices[j].x * contour.vertices[i].y;
        }

        (area * 0.5).abs()
    }

    /// Berechnet den Umfang einer Kontur
    pub fn calculate_contour_perimeter(contour: &Contour) -> f32 {
        if contour.vertices.len() < 2 {
            return 0.0;
        }

        let mut perimeter = 0.0;

        for i in 0..contour.vertices.len() - 1 {
            perimeter += contour.vertices[i].distance(contour.vertices[i + 1]);
        }

        if contour.is_closed && contour.vertices.len() >= 2 {
            perimeter += contour.vertices[contour.vertices.len() - 1].distance(contour.vertices[0]);
        }

        perimeter
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::geometry::metalballs::MetaballField;

    #[test]
    fn test_marching_squares_basic() {
        let mut field = MetaballField::new(10, 10, 1.0);

        // Erstelle ein einfaches Pattern
        for x in 2..8 {
            for y in 2..8 {
                field.set(x, y, 1.5);
            }
        }

        let contours = MarchingSquares::extract_contours(&field, 1.0);
        assert!(!contours.is_empty());
    }

    #[test]
    fn test_contour_operations() {
        let mut contour = Contour::new();
        contour.add_vertex(Vec2::new(0.0, 0.0));
        contour.add_vertex(Vec2::new(1.0, 0.0));
        contour.add_vertex(Vec2::new(1.0, 1.0));
        contour.add_vertex(Vec2::new(0.0, 1.0));
        contour.close();

        let area = MarchingSquares::calculate_contour_area(&contour);
        assert!((area - 1.0).abs() < 1e-6);

        let perimeter = MarchingSquares::calculate_contour_perimeter(&contour);
        assert!((perimeter - 4.0).abs() < 1e-6);
    }
}
