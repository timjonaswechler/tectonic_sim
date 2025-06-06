// src/math/algorithms/point_relaxation/lloyd.rs

use crate::math::{
    error::{MathError, MathResult},
    types::{Bounds2D, SpadePoint},
    utils::constants,
};
// Die LloydConfig und BoundaryHandling enums kommen jetzt aus dem übergeordneten Modul oder config.rs
// Für eine saubere Trennung könnten sie auch hier definiert oder aus einem spezifischen config-Modul importiert werden.
// Da wir sie in `math/point_distribution/voronoi/config.rs` haben, importieren wir von dort.
// Besser wäre es, wenn LloydConfig und BoundaryHandling auch hier oder in einem allgemeineren config-Modul wären.
// Fürs Erste nehmen wir den Pfad an, der existiert, oder definieren sie hier neu, wenn sie unabhängig sein sollen.
// Annahme: Wir definieren sie hier neu, um das Modul eigenständiger zu machen.

use bevy::math::Vec2;
use spade::{DelaunayTriangulation, Triangulation as SpadeTriangulation}; // SpadeTriangulation für insert etc.

/// Konfiguration für den Lloyd-Relaxationsalgorithmus.
#[derive(Debug, Clone)]
pub struct LloydConfig {
    pub max_iterations: usize,
    pub convergence_tolerance_sq: f32, // Quadratische Toleranz für Effizienz
    pub boundary_handling: BoundaryHandling,
}

/// Strategien für die Behandlung von Punkten an den Grenzen während der Lloyd-Relaxation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryHandling {
    Mirror,
    Clamp,
    Wrap,
    // Fixed (wurde entfernt, da schwer allgemein zu implementieren)
}

impl Default for LloydConfig {
    fn default() -> Self {
        Self {
            max_iterations: 15,
            convergence_tolerance_sq: (1e-4 * 1e-4), // Quadrat der Toleranz
            boundary_handling: BoundaryHandling::Clamp,
        }
    }
}

/// Führt Lloyd-Relaxation auf einer Menge von 2D-Punkten durch,
/// um deren Verteilung gleichmäßiger zu gestalten.
pub struct LloydRelaxation {
    config: LloydConfig,
}

/// Statistiken, die während einer Lloyd-Relaxationsdurchlauf gesammelt werden.
#[derive(Debug, Clone, Default)]
pub struct LloydRelaxationStats {
    pub iterations_performed: usize,
    pub total_movement_last_iteration: f32,
    pub max_movement_last_iteration: f32,
    pub converged: bool,
    // Optional: Detailliertere Statistiken pro Iteration
    // pub iteration_details: Vec<LloydIterationDetailStats>,
}

// #[derive(Debug, Clone)]
// pub struct LloydIterationDetailStats { /* ... */ }

impl LloydRelaxation {
    pub fn new(config: LloydConfig) -> Self {
        Self { config }
    }

    /// Wendet Lloyd-Relaxation auf die gegebenen Punkte innerhalb der spezifizierten Grenzen an.
    /// `initial_points`: Die Startpunkte.
    /// `boundary`: Die Grenzen, innerhalb derer die Relaxation stattfindet.
    /// Gibt die relaxierten Punkte und Statistiken zurück.
    pub fn relax_points(
        &self,
        initial_points: &[Vec2],
        boundary: &Bounds2D,
    ) -> MathResult<(Vec<Vec2>, LloydRelaxationStats)> {
        if initial_points.len() < 3 {
            // Delaunay braucht mind. 3 Punkte für nicht-degenerierte Triangulation
            return Err(MathError::InsufficientPoints {
                expected: 3,
                actual: initial_points.len(),
            });
        }
        if !boundary.is_valid() {
            return Err(MathError::InvalidConfiguration {
                message: "Lloyd relaxation boundary is invalid.".to_string(),
            });
        }

        let mut current_points_bevy = initial_points.to_vec();
        let mut stats = LloydRelaxationStats::default();

        for iteration in 0..self.config.max_iterations {
            // 1. Konvertiere zu SpadePoints (nur die innerhalb der Boundary)
            let mut current_points_spade: Vec<SpadePoint> =
                Vec::with_capacity(current_points_bevy.len());
            for &p_bevy in &current_points_bevy {
                if boundary.contains_point(p_bevy) {
                    // Nur Punkte innerhalb der Boundary berücksichtigen
                    current_points_spade.push(SpadePoint::new(p_bevy.x as f32, p_bevy.y as f32));
                } else {
                    // Was tun mit Punkten außerhalb? Momentan werden sie ignoriert für die Triangulation.
                    // Man könnte sie separat behandeln oder an die Boundary clippen *bevor* sie hier reinkommen.
                    // Für die Lloyd-Iteration ist es oft besser, sie sind schon drin.
                    // Hier fügen wir sie nicht zur Triangulation hinzu.
                }
            }

            if current_points_spade.len() < 3 {
                // Nicht genug Punkte für eine valide Triangulation nach dem Filtern
                stats.iterations_performed = iteration;
                // Wenn keine Punkte mehr da sind, oder zu wenige, was ist das Ergebnis?
                // Gebe die Punkte zurück, die noch da sind (eventuell leer).
                return Ok((
                    current_points_bevy
                        .into_iter()
                        .filter(|p| boundary.contains_point(*p))
                        .collect(),
                    stats,
                ));
            }

            // 2. Erstelle Delaunay-Triangulation
            let mut triangulation = DelaunayTriangulation::<SpadePoint>::new();
            for spade_point in &current_points_spade {
                // `insert` kann fehlschlagen, wenn Punkte exakt aufeinanderliegen oder andere Degenerierungen.
                // Wir ignorieren den Fehler hier, da Duplikate keinen großen Einfluss haben sollten.
                let _ = triangulation.insert(*spade_point);
            }

            if triangulation.num_vertices() == 0 {
                // Keine gültigen Vertices in der Triangulation
                stats.iterations_performed = iteration;
                break;
            }

            // 3. Berechne neue Positionen (Zentroide der Voronoi-Zellen)
            let mut new_points_bevy = Vec::with_capacity(current_points_bevy.len());
            let mut max_movement_sq_this_iter = 0.0_f32;
            let mut total_movement_sq_this_iter = 0.0_f32;

            // Zuordnung von originalen Indizes zu Spade-Vertices kann komplex sein, wenn Punkte gefiltert wurden.
            // Einfacher: Iteriere über die Vertices der Triangulation.
            // Und was ist mit den Punkten, die außerhalb der Boundary waren und nicht trianguliert wurden?
            // Die Logik hier muss alle originalen Punkte `current_points_bevy` durchgehen.

            let mut temp_new_positions: Vec<Option<Vec2>> = vec![None; current_points_bevy.len()];

            for (idx, original_point_bevy) in current_points_bevy.iter().enumerate() {
                // Finde den entsprechenden Vertex in der Triangulation.
                // Wenn der Punkt außerhalb der Boundary war, wurde er nicht trianguliert.
                if !boundary.contains_point(*original_point_bevy) {
                    // Behandle Punkte außerhalb der Boundary (z.B. beibehalten oder clippen)
                    temp_new_positions[idx] =
                        Some(self.apply_boundary_handling_to_point(*original_point_bevy, boundary));
                    continue;
                }

                let spade_original_point =
                    SpadePoint::new(original_point_bevy.x, original_point_bevy.y);
                if let Some(vertex_handle) = triangulation.locate_vertex(spade_original_point) {
                    // Berechne Voronoi-Zelle (Liste von Umkreismittelpunkten)
                    let mut voronoi_cell_vertices_spade: Vec<SpadePoint> = Vec::new();
                    for face_handle in triangulation.inner_faces() {
                        if face_handle
                            .vertices()
                            .into_iter()
                            .any(|v| v == vertex_handle)
                        {
                            if let Some(circumcenter) = face_handle.circumcenter().into() {
                                voronoi_cell_vertices_spade
                                    .push(SpadePoint::new(circumcenter.x, circumcenter.y));
                            }
                        }
                    }

                    if voronoi_cell_vertices_spade.len() >= 3 {
                        // Sortiere Voronoi-Vertices um den Generator (original_point)
                        let gen_spade = vertex_handle.position();
                        voronoi_cell_vertices_spade.sort_unstable_by(|a, b| {
                            let angle_a = (a.y - gen_spade.y).atan2(a.x - gen_spade.x);
                            let angle_b = (b.y - gen_spade.y).atan2(b.x - gen_spade.x);
                            angle_a
                                .partial_cmp(&angle_b)
                                .unwrap_or(std::cmp::Ordering::Equal)
                        });

                        // Clippe die Zelle gegen die Boundary (optional, aber oft gut)
                        // Hier implementieren wir Clipping nicht direkt, sondern verlassen uns auf Boundary-Handling später.
                        // Eine robustere Methode wäre, hier die Zelle zu clippen.

                        // Berechne den Zentroid der (ggf. geclippten) Voronoi-Zelle
                        if let Some(centroid_spade) =
                            calculate_polygon_centroid_spade(&voronoi_cell_vertices_spade)
                        {
                            let mut new_pos_bevy = Vec2::new(centroid_spade.x, centroid_spade.y);
                            new_pos_bevy =
                                self.apply_boundary_handling_to_point(new_pos_bevy, boundary);
                            temp_new_positions[idx] = Some(new_pos_bevy);
                        } else {
                            temp_new_positions[idx] = Some(*original_point_bevy); // Behalte alte Position bei Fehler
                        }
                    } else {
                        // Nicht genug Punkte für eine Zelle, behalte alte Position
                        temp_new_positions[idx] = Some(*original_point_bevy);
                    }
                } else {
                    // Punkt nicht in Triangulation gefunden (sollte nicht passieren, wenn er drin war)
                    temp_new_positions[idx] = Some(*original_point_bevy);
                }
            }

            // Sammle neue Punkte und berechne Bewegung
            new_points_bevy.clear(); // Wiederverwenden
            let mut total_valid_points_moved = 0;
            for (idx, opt_new_pos) in temp_new_positions.into_iter().enumerate() {
                let new_pos = opt_new_pos.unwrap_or(current_points_bevy[idx]); // Fallback
                new_points_bevy.push(new_pos);
                let movement_sq = current_points_bevy[idx].distance_squared(new_pos);
                total_movement_sq_this_iter += movement_sq;
                if movement_sq > max_movement_sq_this_iter {
                    max_movement_sq_this_iter = movement_sq;
                }
                if movement_sq > constants::EPSILON_SQUARED {
                    total_valid_points_moved += 1;
                }
            }

            current_points_bevy = new_points_bevy;
            stats.iterations_performed = iteration + 1;
            stats.total_movement_last_iteration = if total_valid_points_moved > 0 {
                (total_movement_sq_this_iter / total_valid_points_moved as f32).sqrt()
            } else {
                0.0
            };
            stats.max_movement_last_iteration = max_movement_sq_this_iter.sqrt();

            if max_movement_sq_this_iter < self.config.convergence_tolerance_sq {
                stats.converged = true;
                break;
            }
        }
        Ok((current_points_bevy, stats))
    }

    /// Wendet die konfigurierte Grenzbehandlung auf einen Punkt an.
    fn apply_boundary_handling_to_point(&self, point: Vec2, boundary: &Bounds2D) -> Vec2 {
        match self.config.boundary_handling {
            BoundaryHandling::Clamp => {
                boundary.closest_point(point) // Clamp zum Rand der Bounds
            }
            BoundaryHandling::Wrap => {
                let mut wrapped_x = point.x;
                let mut wrapped_y = point.y;
                let width = boundary.width();
                let height = boundary.height();

                if width > constants::EPSILON {
                    while wrapped_x < boundary.min.x {
                        wrapped_x += width;
                    }
                    while wrapped_x >= boundary.max.x {
                        wrapped_x -= width;
                    } // >= für max
                } else {
                    wrapped_x = boundary.min.x; // Kein Wrap bei Breite 0
                }
                if height > constants::EPSILON {
                    while wrapped_y < boundary.min.y {
                        wrapped_y += height;
                    }
                    while wrapped_y >= boundary.max.y {
                        wrapped_y -= height;
                    } // >= für max
                } else {
                    wrapped_y = boundary.min.y;
                }
                Vec2::new(wrapped_x, wrapped_y)
            }
            BoundaryHandling::Mirror => {
                let mut mirrored_x = point.x;
                let mut mirrored_y = point.y;
                // Spiegelung an jeder Kante, wenn außerhalb
                if point.x < boundary.min.x {
                    mirrored_x = boundary.min.x + (boundary.min.x - point.x);
                } else if point.x > boundary.max.x {
                    mirrored_x = boundary.max.x - (point.x - boundary.max.x);
                }
                if point.y < boundary.min.y {
                    mirrored_y = boundary.min.y + (boundary.min.y - point.y);
                } else if point.y > boundary.max.y {
                    mirrored_y = boundary.max.y - (point.y - boundary.max.y);
                }

                // Nach dem Spiegeln erneut clippen, falls es über die gegenüberliegende Seite hinausgeht
                // bei sehr kleinen Boundaries. Einfacher ist oft, nur einmal zu spiegeln.
                // Hier: nur einmal spiegeln.
                Vec2::new(mirrored_x, mirrored_y)
            }
        }
    }
}

/// Berechnet den Zentroid eines Polygons, das durch SpadePoints definiert ist.
fn calculate_polygon_centroid_spade(polygon_vertices: &[SpadePoint]) -> Option<SpadePoint> {
    let n = polygon_vertices.len();
    if n < 3 {
        return None;
    } // Nicht genug Punkte für eine Fläche

    let mut area_sum_doubled = 0.0;
    let mut centroid_x_sum = 0.0;
    let mut centroid_y_sum = 0.0;

    for i in 0..n {
        let p1 = polygon_vertices[i];
        let p2 = polygon_vertices[(i + 1) % n]; // Nächster Punkt mit Umlauf

        let cross_term = p1.x * p2.y - p2.x * p1.y;
        area_sum_doubled += cross_term;
        centroid_x_sum += (p1.x + p2.x) * cross_term;
        centroid_y_sum += (p1.y + p2.y) * cross_term;
    }

    if area_sum_doubled.abs() < constants::EPSILON {
        // Polygon ist degeneriert (Linie oder Punkt)
        // Fallback: Arithmetischer Mittelpunkt
        let mut sum_x = 0.0;
        let mut sum_y = 0.0;
        for p in polygon_vertices {
            sum_x += p.x;
            sum_y += p.y;
        }
        return Some(SpadePoint::new(sum_x / n as f32, sum_y / n as f32));
    }

    let six_times_area = 3.0 * area_sum_doubled; // Eigentlich 6 * (Area), Area = 0.5 * area_sum_doubled
    Some(SpadePoint::new(
        centroid_x_sum / six_times_area,
        centroid_y_sum / six_times_area,
    ))
}
