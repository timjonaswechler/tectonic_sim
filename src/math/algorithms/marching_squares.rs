// src/math/algorithms/marching_squares.rs

// Geändert: ScalarField2D Trait verwenden
use crate::math::scalar_field::ScalarField2D;
use crate::math::utils::constants;
use bevy::math::Vec2;

/// Repräsentiert eine extrahierte Konturlinie (Iso-Linie) aus einem Skalarfeld.
#[derive(Debug, Clone, Default)]
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
/// der Konturen aus einem beliebigen `ScalarField2D` extrahiert.
#[derive(PartialEq)] // PartialEq bleibt, da edge_table nun Option<...> enthält
pub struct MarchingSquaresIterator<'a, F: ScalarField2D + ?Sized> {
    field: &'a F,
    threshold: f32,
    current_x: usize,
    current_y: usize,
    visited_cells: Vec<Vec<bool>>,
    edge_table: [[Option<((f32, f32), (f32, f32))>; 2]; 16],
}

impl<'a, F: ScalarField2D + ?Sized> MarchingSquaresIterator<'a, F> {
    pub fn new(field: &'a F, threshold: f32) -> Self {
        let field_width = field.width(); // Methoden vom Trait verwenden
        let field_height = field.height();

        let visited_cells = if field_width > 0 && field_height > 0 {
            vec![vec![false; field_height.saturating_sub(1)]; field_width.saturating_sub(1)]
        } else {
            vec![]
        };

        let mut iterator = Self {
            field,
            threshold,
            current_x: 0,
            current_y: 0,
            visited_cells,
            edge_table: Default::default(),
        };
        iterator.initialize_edge_table();
        iterator
    }

    fn initialize_edge_table(&mut self) {
        // Konvention: Bit 8: oben-links, Bit 4: oben-rechts, Bit 2: unten-rechts, Bit 1: unten-links
        // Liniensegmente: ((x1,y1), (x2,y2)) - relative Kantenpositionen
        // 0.0: linke/obere Kante der Zelle
        // 0.5: Mitte der Kante
        // 1.0: rechte/untere Kante der Zelle
        // Edge Order in lookup: Left(0), Top(1), Right(2), Bottom(3) -> Bit for configuration
        // Our interpretation: Output line segment from point P on edge E1 to point Q on edge E2.
        // P: (relative_x_start, relative_y_start) referring to the edge it's on
        // Q: (relative_x_end, relative_y_end) referring to the edge it's on
        // Statt komplexer Relativ-zu-Relativ, definieren wir relative Punkte auf Zellgrenzen.
        // Top edge: (x, 0), Left: (0, y), Right: (1, y), Bottom: (x, 1)
        // Wo 0.5 auf einer Kante einen Schnittpunkt in der Mitte der Kante bedeutet.

        self.edge_table[0] = [None, None]; // 0000: --
        self.edge_table[1] = [Some(((0.0, 0.5), (0.5, 1.0))), None]; // 0001: L -> B
        self.edge_table[2] = [Some(((0.5, 1.0), (1.0, 0.5))), None]; // 0010: B -> R
        self.edge_table[3] = [Some(((0.0, 0.5), (1.0, 0.5))), None]; // 0011: L -> R
        self.edge_table[4] = [Some(((0.5, 0.0), (1.0, 0.5))), None]; // 0100: T -> R
        self.edge_table[5] = [
            Some(((0.5, 0.0), (0.0, 0.5))), // T -> L (upper segment)
            Some(((0.5, 1.0), (1.0, 0.5))), // B -> R (lower segment) - Ambiguous Case 5
        ];
        self.edge_table[6] = [Some(((0.5, 0.0), (0.5, 1.0))), None]; // 0110: T -> B (verbindet TR und BR)
        self.edge_table[7] = [Some(((0.5, 0.0), (0.0, 0.5))), None]; // 0111: T -> L
        self.edge_table[8] = [Some(((0.0, 0.5), (0.5, 0.0))), None]; // 1000: L -> T
        self.edge_table[9] = [Some(((0.5, 0.0), (0.5, 1.0))), None]; // 1001: T -> B (verbindet TL und BL)
        self.edge_table[10] = [
            Some(((0.5, 0.0), (1.0, 0.5))), // T -> R (upper segment)
            Some(((0.0, 0.5), (0.5, 1.0))), // L -> B (lower segment) - Ambiguous Case 10
        ];
        self.edge_table[11] = [Some(((0.5, 0.0), (1.0, 0.5))), None]; // 1011: T -> R
        self.edge_table[12] = [Some(((0.0, 0.5), (1.0, 0.5))), None]; // 1100: L -> R
        self.edge_table[13] = [Some(((0.5, 1.0), (1.0, 0.5))), None]; // 1101: B -> R
        self.edge_table[14] = [Some(((0.5, 1.0), (0.0, 0.5))), None]; // 1110: B -> L
        self.edge_table[15] = [None, None]; // 1111: ##
    }

    fn get_cell_configuration(&self, x: usize, y: usize) -> usize {
        let mut config_idx = 0;
        if self.field.get_value(x, y) >= self.threshold {
            // field.get_value
            config_idx |= 8; // Top-left
        }
        if self.field.get_value(x + 1, y) >= self.threshold {
            config_idx |= 4; // Top-right
        }
        if self.field.get_value(x + 1, y + 1) >= self.threshold {
            config_idx |= 2; // Bottom-right
        }
        if self.field.get_value(x, y + 1) >= self.threshold {
            config_idx |= 1; // Bottom-left
        }
        config_idx
    }

    fn interpolate_intersection(
        &self,
        val1: f32,
        val2: f32,
        p1_world: Vec2,
        p2_world: Vec2,
    ) -> Vec2 {
        if (val1 - val2).abs() < constants::EPSILON {
            return (p1_world + p2_world) * 0.5; // Midpoint
        }
        let t = ((self.threshold - val1) / (val2 - val1)).clamp(0.0, 1.0);
        p1_world.lerp(p2_world, t)
    }

    // Konvertiert relative Kantenkoordinaten zu Weltkoordinaten.
    // cell_origin_world ist die Weltkoordinate der oberen linken Ecke der Zelle.
    fn get_world_intersection_point(
        &self,
        rel_x: f32,
        rel_y: f32,
        x_cell: usize, // grid cell index x
        y_cell: usize, // grid cell index y
    ) -> Vec2 {
        let cell_size = self.field.cell_size();
        let cell_world_origin = self.field.cell_to_world(x_cell, y_cell);

        // Grid-Eckpunkte für Interpolation
        let p_tl = self.field.get_value(x_cell, y_cell); // Top-left
        let p_tr = self.field.get_value(x_cell + 1, y_cell); // Top-right
        let p_bl = self.field.get_value(x_cell, y_cell + 1); // Bottom-left
        let p_br = self.field.get_value(x_cell + 1, y_cell + 1); // Bottom-right

        // Weltkoordinaten der Zell-Eckpunkte
        let world_tl = cell_world_origin;
        let world_tr = world_tl + Vec2::new(cell_size, 0.0);
        let world_bl = world_tl + Vec2::new(0.0, cell_size);
        let world_br = world_tl + Vec2::new(cell_size, cell_size);

        // Logik für jede der vier Kanten
        if rel_y == 0.0 && rel_x > 0.0 && rel_x < 1.0 {
            // Top edge (ohne Ecken)
            self.interpolate_intersection(p_tl, p_tr, world_tl, world_tr)
        } else if rel_x == 1.0 && rel_y > 0.0 && rel_y < 1.0 {
            // Right edge
            self.interpolate_intersection(p_tr, p_br, world_tr, world_br)
        } else if rel_y == 1.0 && rel_x > 0.0 && rel_x < 1.0 {
            // Bottom edge
            self.interpolate_intersection(p_bl, p_br, world_bl, world_br)
        } else if rel_x == 0.0 && rel_y > 0.0 && rel_y < 1.0 {
            // Left edge
            self.interpolate_intersection(p_tl, p_bl, world_tl, world_bl)
        } else {
            // Dies sollte nicht für Mittenpunkte der Kanten (0.5) passieren
            // Wenn es passiert, ist rel_x oder rel_y 0.0 oder 1.0, was einer Ecke entspricht
            // Hier sollte der Marching Squares Algo eigentlich nicht Kantenpunkte auf Ecken legen, sondern dazwischen.
            // Als Fallback: Lineare Interpolation der Zell-Ecke auf der relativen Position
            cell_world_origin + Vec2::new(rel_x * cell_size, rel_y * cell_size)
        }
    }

    fn get_cell_edge_intersections(
        &self,
        x_cell: usize,
        y_cell: usize,
        config_idx: usize,
    ) -> Vec<Vec2> {
        let segments = Vec::new();
        if config_idx == 0 || config_idx == 15 {
            return segments;
        }

        let mut intersections = Vec::new(); // Wird Paare von Punkten enthalten

        for edge_def_opt in self.edge_table[config_idx].iter() {
            if let Some(((rel_x1, rel_y1), (rel_x2, rel_y2))) = edge_def_opt {
                let world_p1 = self.get_world_intersection_point(*rel_x1, *rel_y1, x_cell, y_cell);
                let world_p2 = self.get_world_intersection_point(*rel_x2, *rel_y2, x_cell, y_cell);
                intersections.push(world_p1);
                intersections.push(world_p2);
            }
        }
        intersections
    }

    fn find_next_unvisited_contour_cell(&mut self) -> Option<(usize, usize)> {
        let field_width = self.field.width();
        let field_height = self.field.height();

        if field_width <= 1 || field_height <= 1 {
            return None;
        }

        for y in self.current_y..field_height.saturating_sub(1) {
            let start_x = if y == self.current_y {
                self.current_x
            } else {
                0
            };
            for x in start_x..field_width.saturating_sub(1) {
                if !self.visited_cells[x][y] {
                    let config = self.get_cell_configuration(x, y);
                    if config != 0 && config != 15 {
                        self.current_x = x;
                        self.current_y = y;
                        return Some((x, y));
                    }
                    self.visited_cells[x][y] = true;
                }
            }
            self.current_x = 0;
        }
        None
    }

    fn trace_single_contour(&mut self, start_x: usize, start_y: usize) -> Contour {
        let mut contour = Contour::new();
        let field_width = self.field.width();
        let field_height = self.field.height();

        if start_x < field_width.saturating_sub(1) && start_y < field_height.saturating_sub(1) {
            let config = self.get_cell_configuration(start_x, start_y);
            self.visited_cells[start_x][start_y] = true;

            if config != 0 && config != 15 {
                let intersections = self.get_cell_edge_intersections(start_x, start_y, config);
                for i in (0..intersections.len()).step_by(2) {
                    if i + 1 < intersections.len() {
                        if contour.is_empty() {
                            contour.add_vertex(intersections[i]);
                        }
                        contour.add_vertex(intersections[i + 1]);
                    }
                }
            }
        }
        contour
    }
}

impl<'a, F: ScalarField2D + ?Sized> Iterator for MarchingSquaresIterator<'a, F> {
    type Item = Contour;

    fn next(&mut self) -> Option<Self::Item> {
        let field_width = self.field.width();
        let field_height = self.field.height();

        if field_width <= 1 || field_height <= 1 {
            return None;
        }

        // Loop for finding the start of a new contour part
        while let Some((start_x, start_y)) = self.find_next_unvisited_contour_cell() {
            let contour_segments = self.trace_single_contour(start_x, start_y);

            // Advance current_x, current_y for the *next* call to find_next_unvisited_contour_cell
            if self.current_x < field_width.saturating_sub(2) {
                self.current_x += 1;
            } else {
                self.current_x = 0;
                if self.current_y < field_height.saturating_sub(2) {
                    self.current_y += 1;
                } else {
                    // Reached the end of the grid to check for new contours.
                    // This logic means the next find_next_unvisited_contour_cell will likely return None
                    // if all previous cells were visited or contained no contours.
                    // Or, if find_next_unvisited_contour_cell still finds something,
                    // current_x/y will be updated by it.
                }
            }

            if !contour_segments.is_empty() {
                return Some(contour_segments);
            }
            // If contour_segments is empty (e.g. traced cell was marked visited but yielded nothing new),
            // loop continues to find next unvisited cell.
        }
        None // No more unvisited cells with contours found
    }
}

pub struct MarchingSquares;

impl MarchingSquares {
    // Generische Funktion, die den Trait ScalarField2D verwendet
    pub fn extract_contour_segments<F: ScalarField2D + ?Sized>(
        field: &F,
        threshold: f32,
    ) -> Vec<Contour> {
        if field.width() <= 1 || field.height() <= 1 {
            return Vec::new();
        }
        let iterator = MarchingSquaresIterator::new(field, threshold);
        iterator.collect()
    }

    pub fn simplify_contour_rdp(contour: &Contour, epsilon: f32) -> Contour {
        if contour.vertices.len() <= 2 {
            return contour.clone();
        }
        let simplified_vertices = rdp_simplify(&contour.vertices, epsilon);

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
            result.is_closed = false;
        }
        result
    }

    pub fn calculate_contour_area(contour: &Contour) -> f32 {
        if !contour.is_closed || contour.vertices.len() < 3 {
            return 0.0;
        }
        let mut area = 0.0;
        let n = contour.vertices.len();
        for i in 0..n.saturating_sub(1) {
            let p1 = contour.vertices[i];
            let p2 = contour.vertices[i + 1];
            area += (p1.x * p2.y) - (p2.x * p1.y);
        }
        (area * 0.5).abs()
    }
}

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

        rec_results1.pop();
        rec_results1.extend(rec_results2);
        rec_results1
    } else {
        vec![points[0], points[end]]
    }
}

fn perpendicular_distance(pt: Vec2, line_start: Vec2, line_end: Vec2) -> f32 {
    let dx = line_end.x - line_start.x;
    let dy = line_end.y - line_start.y;

    if dx.abs() < constants::EPSILON && dy.abs() < constants::EPSILON {
        return pt.distance(line_start);
    }

    let num = (dy * pt.x - dx * pt.y + line_end.x * line_start.y - line_end.y * line_start.x).abs();
    let den = (dy * dy + dx * dx).sqrt();
    if den < constants::EPSILON {
        // Sollte nicht passieren, wenn dx/dy nicht beide ~0 sind
        return pt.distance(line_start);
    }
    num / den
}
