// src/math/geometry/metalballs/marching_squares.rs

use crate::math::geometry::metalballs::field::MetaballField;
use crate::math::utils::constants; // Für EPSILON
use bevy::math::Vec2;

/// Repräsentiert eine extrahierte Konturlinie (Iso-Linie) aus einem Skalarfeld.
#[derive(Debug, Clone, Default)] // Default hinzugefügt
pub struct Contour {
    pub vertices: Vec<Vec2>,
    pub is_closed: bool,
}

impl Contour {
    pub fn new() -> Self {
        Self::default()
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
        if self.vertices.len() >= 3 {
            // Eine geschlossene Kontur braucht mind. 3 Punkte
            self.is_closed = true;
        }
    }

    pub fn len(&self) -> usize {
        self.vertices.len()
    }

    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }
}

/// Implementiert den Marching Squares Algorithmus als Iterator,
/// der Konturen aus einem `MetaballField` extrahiert.
#[derive(PartialEq)]
pub struct MarchingSquaresIterator<'a> {
    field: &'a MetaballField,
    threshold: f32,
    current_x: usize, // Nächste zu prüfende Spalte für einen neuen Konturstart
    current_y: usize, // Nächste zu prüfende Zeile für einen neuen Konturstart
    visited_cells: Vec<Vec<bool>>, // Markiert Zellen (definiert durch linke obere Ecke), die bereits Teil einer verfolgten Kontur waren
    // Die Edge-Table definiert für jede der 16 Zellkonfigurationen,
    // welche Kanten geschnitten werden. Jedes Paar (f32, f32) ist ein (rel_x, rel_y)
    // Offset innerhalb der Zelle, um den Start- und Endpunkt eines Liniensegments zu definieren.
    // 0.0 ist linke/obere Kante, 0.5 ist Mitte der Kante, 1.0 ist rechte/untere Kante.
    // Reihenfolge: (Top-Start, Top-End), (Right-Start, Right-End), (Bottom-Start, Bottom-End), (Left-Start, Left-End)
    // oder Paare von Liniensegmenten.
    edge_table: [[Option<((f32, f32), (f32, f32))>; 2]; 16], // Max 2 Liniensegmente pro Zelle
}

impl<'a> MarchingSquaresIterator<'a> {
    pub fn new(field: &'a MetaballField, threshold: f32) -> Self {
        let visited_cells = if field.width > 0 && field.height > 0 {
            vec![vec![false; field.height - 1]; field.width - 1] // Für Zellen, nicht Punkte
        } else {
            vec![] // Leeres Feld, keine Zellen
        };

        let mut iterator = Self {
            field,
            threshold,
            current_x: 0,
            current_y: 0,
            visited_cells,
            edge_table: Default::default(), // Wird in initialize_edge_table gefüllt
        };
        iterator.initialize_edge_table();
        iterator
    }

    /// Initialisiert die Lookup-Tabelle für Marching Squares.
    /// Definiert Liniensegmente für jede der 16 möglichen Zellkonfigurationen.
    /// Punkte sind relativ zur Zelle (0,0 oben links, 1,1 unten rechts).
    /// Kantenmitten: Top(0.5,0), Right(1,0.5), Bottom(0.5,1), Left(0,0.5)
    fn initialize_edge_table(&mut self) {
        // Konvention: Bit 8: oben-links, Bit 4: oben-rechts, Bit 2: unten-rechts, Bit 1: unten-links
        // Liniensegmente: ((x1,y1), (x2,y2))
        self.edge_table[0] = [None, None]; // 0000
        self.edge_table[1] = [Some(((0.0, 0.5), (0.5, 1.0))), None]; // 0001 (L->B)
        self.edge_table[2] = [Some(((0.5, 1.0), (1.0, 0.5))), None]; // 0010 (B->R)
        self.edge_table[3] = [Some(((0.0, 0.5), (1.0, 0.5))), None]; // 0011 (L->R)
        self.edge_table[4] = [Some(((0.5, 0.0), (1.0, 0.5))), None]; // 0100 (T->R)
        self.edge_table[5] = [
            Some(((0.0, 0.5), (0.5, 0.0))),
            Some(((0.5, 1.0), (1.0, 0.5))),
        ]; // 0101 (L->T, B->R) Ambiguous, Standardauflösung
        self.edge_table[6] = [Some(((0.5, 1.0), (0.5, 0.0))), None]; // 0110 (B->T via R) (R Kante)
        self.edge_table[7] = [Some(((0.0, 0.5), (0.5, 0.0))), None]; // 0111 (L->T)
        self.edge_table[8] = [Some(((0.5, 0.0), (0.0, 0.5))), None]; // 1000 (T->L)
        self.edge_table[9] = [Some(((0.5, 1.0), (0.5, 0.0))), None]; // 1001 (B->T via L) (L Kante)
        self.edge_table[10] = [
            Some(((0.5, 0.0), (1.0, 0.5))),
            Some(((0.0, 0.5), (0.5, 1.0))),
        ]; // 1010 (T->R, L->B) Ambiguous, Standardauflösung
        self.edge_table[11] = [Some(((0.5, 0.0), (1.0, 0.5))), None]; // 1011 (T->R)
        self.edge_table[12] = [Some(((1.0, 0.5), (0.0, 0.5))), None]; // 1100 (R->L) (T Kante)
        self.edge_table[13] = [Some(((0.5, 1.0), (1.0, 0.5))), None]; // 1101 (B->R)
        self.edge_table[14] = [Some(((0.5, 1.0), (0.0, 0.5))), None]; // 1110 (B->L)
        self.edge_table[15] = [None, None]; // 1111
    }

    /// Berechnet den Konfigurations-Index (0-15) für eine Zelle (x, y - obere linke Ecke).
    fn get_cell_configuration(&self, x: usize, y: usize) -> usize {
        let mut config_idx = 0;
        if self.field.get(x, y) >= self.threshold {
            config_idx |= 8;
        } // Oben-Links
        if self.field.get(x + 1, y) >= self.threshold {
            config_idx |= 4;
        } // Oben-Rechts
        if self.field.get(x + 1, y + 1) >= self.threshold {
            config_idx |= 2;
        } // Unten-Rechts
        if self.field.get(x, y + 1) >= self.threshold {
            config_idx |= 1;
        } // Unten-Links
        config_idx
    }

    /// Interpoliert die genaue Position eines Schnittpunkts auf einer Zellkante.
    /// `p1_val`, `p2_val` sind die Feldwerte an den Enden der Kante.
    /// `p1_coord`, `p2_coord` sind die Weltkoordinaten der Enden der Kante.
    fn interpolate_intersection(
        &self,
        p1_val: f32,
        p2_val: f32,
        p1_coord: Vec2,
        p2_coord: Vec2,
    ) -> Vec2 {
        if (p2_val - p1_val).abs() < constants::EPSILON {
            // Werte sind (fast) gleich
            return (p1_coord + p2_coord) * 0.5; // Nimm die Mitte
        }
        // t = (threshold - val1) / (val2 - val1)
        let t = ((self.threshold - p1_val) / (p2_val - p1_val)).clamp(0.0, 1.0);
        p1_coord.lerp(p2_coord, t)
    }

    /// Berechnet die tatsächlichen Schnittpunkte für eine Zelle basierend auf ihrer Konfiguration.
    fn get_cell_edge_intersections(
        &self,
        x_cell: usize,
        y_cell: usize,
        config_idx: usize,
    ) -> Vec<Vec2> {
        if config_idx == 0 || config_idx == 15 {
            return Vec::new();
        }

        let mut segments = Vec::new();
        let cell_size = self.field.cell_size;

        // Eckpunkte der Zelle in Weltkoordinaten
        let p_tl_coord = Vec2::new(x_cell as f32 * cell_size, y_cell as f32 * cell_size);
        let p_tr_coord = Vec2::new((x_cell + 1) as f32 * cell_size, y_cell as f32 * cell_size);
        let p_br_coord = Vec2::new(
            (x_cell + 1) as f32 * cell_size,
            (y_cell + 1) as f32 * cell_size,
        );
        let p_bl_coord = Vec2::new(x_cell as f32 * cell_size, (y_cell + 1) as f32 * cell_size);

        // Feldwerte an den Eckpunkten
        let val_tl = self.field.get(x_cell, y_cell);
        let val_tr = self.field.get(x_cell + 1, y_cell);
        let val_br = self.field.get(x_cell + 1, y_cell + 1);
        let val_bl = self.field.get(x_cell, y_cell + 1);

        for edge_opt in self.edge_table[config_idx].iter() {
            if let Some(((rel_x1, rel_y1), (rel_x2, rel_y2))) = edge_opt {
                let mut p1 = Vec2::ZERO;
                let mut p2 = Vec2::ZERO;

                // Punkt 1 des Segments
                if *rel_y1 == 0.0_f32 {
                    // Top edge
                    p1 = self.interpolate_intersection(val_tl, val_tr, p_tl_coord, p_tr_coord);
                } else if *rel_x1 == 1.0_f32 {
                    // Right edge
                    p1 = self.interpolate_intersection(val_tr, val_br, p_tr_coord, p_br_coord);
                } else if *rel_y1 == 1.0_f32 {
                    // Bottom edge
                    p1 = self.interpolate_intersection(val_bl, val_br, p_bl_coord, p_br_coord);
                } else if *rel_x1 == 0.0_f32 {
                    // Left edge
                    p1 = self.interpolate_intersection(val_tl, val_bl, p_tl_coord, p_bl_coord);
                }

                // Punkt 2 des Segments
                if *rel_y2 == 0.0_f32 {
                    // Top edge
                    p2 = self.interpolate_intersection(val_tl, val_tr, p_tl_coord, p_tr_coord);
                } else if *rel_x2 == 1.0_f32 {
                    // Right edge
                    p2 = self.interpolate_intersection(val_tr, val_br, p_tr_coord, p_br_coord);
                } else if *rel_y2 == 1.0_f32 {
                    // Bottom edge
                    p2 = self.interpolate_intersection(val_bl, val_br, p_bl_coord, p_br_coord);
                } else if *rel_x2 == 0.0_f32 {
                    // Left edge
                    p2 = self.interpolate_intersection(val_tl, val_bl, p_tl_coord, p_bl_coord);
                }
                segments.push(p1);
                segments.push(p2);
            }
        }
        segments
    }

    /// Findet die nächste unbesuchte Zelle, die eine Kontur enthält.
    fn find_next_unvisited_contour_cell(&mut self) -> Option<(usize, usize)> {
        if self.field.width <= 1 || self.field.height <= 1 {
            return None;
        }

        for y in self.current_y..(self.field.height - 1) {
            let start_x = if y == self.current_y {
                self.current_x
            } else {
                0
            };
            for x in start_x..(self.field.width - 1) {
                if !self.visited_cells[x][y] {
                    let config = self.get_cell_configuration(x, y);
                    if config != 0 && config != 15 {
                        // Zelle enthält eine Kontur
                        self.current_x = x; // Merke für den nächsten Aufruf
                        self.current_y = y;
                        return Some((x, y));
                    }
                    self.visited_cells[x][y] = true; // Markiere als besucht, auch wenn leer
                }
            }
            self.current_x = 0; // Für die nächste Zeile von vorne beginnen
        }
        None // Keine weiteren Konturen gefunden
    }

    /// Verfolgt eine einzelne Kontur von einem Startpunkt aus.
    /// Diese Funktion ist komplex und hier nur stark vereinfacht für den Iterator.
    /// Eine vollständige Implementierung würde Kantenverfolgung und das Schließen von Schleifen benötigen.
    fn trace_single_contour(&mut self, start_x: usize, start_y: usize) -> Contour {
        let mut contour = Contour::new();
        // Für eine einfache Implementierung im Iterator geben wir nur die Segmente der aktuellen Zelle zurück.
        // Eine echte Konturverfolgung ist aufwändiger.

        if start_x < self.field.width - 1 && start_y < self.field.height - 1 {
            let config = self.get_cell_configuration(start_x, start_y);
            self.visited_cells[start_x][start_y] = true; // Markiere Zelle als verarbeitet

            if config != 0 && config != 15 {
                let intersections = self.get_cell_edge_intersections(start_x, start_y, config);
                for i in (0..intersections.len()).step_by(2) {
                    if i + 1 < intersections.len() {
                        // Für den Iterator geben wir jedes Segment als separate "Kontur" zurück
                        // oder bauen eine längere Kontur auf, wenn der Algorithmus erweitert wird.
                        // Hier: Nur die Punkte der aktuellen Zelle.
                        if contour.is_empty() {
                            contour.add_vertex(intersections[i]);
                        }
                        contour.add_vertex(intersections[i + 1]);
                    }
                }
                // Vereinfachung: Wir versuchen nicht, die Kontur zu schließen oder zu verfolgen.
                // Die Contour wird als offen betrachtet, es sei denn, sie schließt sich zufällig.
            }
        }
        contour
    }
}

impl<'a> Iterator for MarchingSquaresIterator<'a> {
    type Item = Contour;

    fn next(&mut self) -> Option<Self::Item> {
        if self.field.width <= 1 || self.field.height <= 1 {
            return None;
        }

        while let Some((start_x, start_y)) = self.find_next_unvisited_contour_cell() {
            // Die Logik hier müsste eine vollständige Konturverfolgung sein.
            // Für den Iterator geben wir pro gefundener Zelle mit Kontur deren Segmente zurück.
            // Das ist nicht ideal, aber ein Start. Eine echte Konturverfolgung ist komplex.
            // Hier wird `trace_single_contour` sehr rudimentär sein.
            let contour_segments = self.trace_single_contour(start_x, start_y);

            // Setze den Iterator für den nächsten Suchlauf korrekt:
            if self.current_x < self.field.width - 2 {
                // -2 weil wir Zellen (x,y) bis (width-2, height-2) prüfen
                self.current_x += 1;
            } else {
                self.current_x = 0;
                if self.current_y < self.field.height - 2 {
                    self.current_y += 1;
                } else {
                    // Ende erreicht, aber find_next_unvisited_contour_cell hätte None zurückgeben sollen.
                    // Dieser Fall sollte durch die Schleifenbedingung oben abgedeckt sein.
                    return None;
                }
            }

            if !contour_segments.is_empty() {
                return Some(contour_segments);
            }
        }
        None
    }
}

/// Utility-Funktionen für Marching Squares.
pub struct MarchingSquares;

impl MarchingSquares {
    /// Extrahiert alle Kontursegmente aus einem Skalarfeld.
    /// Beachte: Dies gibt keine zusammenhängenden Konturen zurück, sondern die Segmente pro Zelle.
    /// Für zusammenhängende Konturen wäre ein komplexerer Algorithmus nötig.
    pub fn extract_contour_segments(field: &MetaballField, threshold: f32) -> Vec<Contour> {
        if field.width <= 1 || field.height <= 1 {
            return Vec::new();
        }
        let iterator = MarchingSquaresIterator::new(field, threshold);
        iterator.collect()
    }

    // Die Methoden simplify_contour, smooth_contour, calculate_contour_area, calculate_contour_perimeter
    // können beibehalten werden, aber sie arbeiten am besten mit *geschlossenen, zusammenhängenden* Konturen.
    // Die aktuelle Implementierung des Iterators liefert das nicht unbedingt.

    /// Vereinfacht eine Kontur durch den Ramer-Douglas-Peucker-Algorithmus.
    pub fn simplify_contour_rdp(contour: &Contour, epsilon: f32) -> Contour {
        if contour.vertices.len() <= 2 {
            return contour.clone();
        }
        let simplified_vertices = rdp_simplify(&contour.vertices, epsilon);

        // Stelle sicher, dass geschlossene Konturen geschlossen bleiben, wenn sie noch >2 Punkte haben
        let mut result = Contour::new();
        result.vertices = simplified_vertices;
        if contour.is_closed && result.vertices.len() >= 3 {
            if result.vertices.first() != result.vertices.last() {
                if let Some(first) = result.vertices.first().copied() {
                    result.vertices.push(first);
                }
            }
            result.close();
        } else if result.vertices.len() < 3 {
            result.is_closed = false; // Kann nicht mehr geschlossen sein
        }
        result
    }

    // ... (smooth_contour, calculate_contour_area, calculate_contour_perimeter bleiben ähnlich)
    // aber ihre Nützlichkeit hängt von der Qualität der Eingabekonturen ab.

    /// Berechnet die Fläche einer geschlossenen Kontur (Shoelace-Formel)
    pub fn calculate_contour_area(contour: &Contour) -> f32 {
        if !contour.is_closed || contour.vertices.len() < 3 {
            return 0.0;
        }
        let mut area = 0.0;
        let n = contour.vertices.len(); // Inklusive des duplizierten Endpunkts
        for i in 0..n - 1 {
            // Iteriere bis zum vorletzten Punkt
            let p1 = contour.vertices[i];
            let p2 = contour.vertices[i + 1];
            area += (p1.x * p2.y) - (p2.x * p1.y);
        }
        (area * 0.5).abs()
    }
}

/// Ramer-Douglas-Peucker Algorithmus zur Linienvereinfachung.
fn rdp_simplify(points: &[Vec2], epsilon: f32) -> Vec<Vec2> {
    if points.len() <= 2 {
        return points.to_vec();
    }

    let mut dmax = 0.0;
    let mut index = 0;
    let end = points.len() - 1;

    for i in 1..end {
        let d = perpendicular_distance(points[i], points[0], points[end]);
        if d > dmax {
            index = i;
            dmax = d;
        }
    }

    if dmax > epsilon {
        let mut rec_results1 = rdp_simplify(&points[0..=index], epsilon);
        let rec_results2 = rdp_simplify(&points[index..=end], epsilon);

        rec_results1.pop(); // Entferne Duplikat
        rec_results1.extend(rec_results2);
        rec_results1
    } else {
        vec![points[0], points[end]]
    }
}

/// Senkrechter Abstand eines Punktes zu einem Liniensegment.
fn perpendicular_distance(pt: Vec2, line_start: Vec2, line_end: Vec2) -> f32 {
    let dx = line_end.x - line_start.x;
    let dy = line_end.y - line_start.y;

    if dx == 0.0 && dy == 0.0 {
        // Start und Endpunkt sind gleich
        return pt.distance(line_start);
    }

    let num = (dy * pt.x - dx * pt.y + line_end.x * line_start.y - line_end.y * line_start.x).abs();
    let den = (dy * dy + dx * dx).sqrt();
    num / den
}
