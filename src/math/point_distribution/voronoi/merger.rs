// src/math/point_distribution/voronoi/merger.rs

use crate::math::{
    algorithms::convex_hull::{ConvexHullAlgorithm, ConvexHullComputer},
    error::MathResult,
    point_distribution::voronoi::voronoi_diagram::VoronoiCell,
    prelude::SeedResource,
    utils::constants,
};
use bevy::math::Vec2;
use std::collections::HashSet;

/// Definiert verschiedene Strategien zum Auswählen und Zusammenfügen von Voronoi-Zellen.
#[derive(Debug, Clone, Copy)]
pub enum MergeStrategy {
    /// Fügt Zellen basierend auf Flächenkriterien zusammen.
    AreaBased {
        min_area: f32,
        max_area: f32,
    },
    /// Fügt Zellen basierend auf der Anzahl ihrer Nachbarn zusammen (vereinfacht durch Vertex-Anzahl).
    NeighborhoodBased {
        max_vertices_for_cell: usize,
    },
    /// Wählt eine bestimmte Anzahl der größten Zellen aus.
    LargestCells {
        count: usize,
    },
    /// Wählt Zellen basierend auf Ähnlichkeit ihrer Eigenschaften (z.B. Fläche, Umfang) aus.
    PropertyBased {
        similarity_threshold: f32,
    }, // [0,1]
    /// Wählt eine zufällige Anzahl von Zellen oder spezifische zufällige Zellen aus.
    Random {
        target_count: usize,
    },
    /// Versucht, alle ausgewählten Zellen zu einer einzigen Kontur zu verschmelzen.
    UnionAllSelected, // Neue Strategie
    /// Wählt zufällig eine Zelle mit mindestens 6 Vertices und erweitert die Auswahl auf verbundene Nachbarzellen
    RandomConnectedCells {
        min_vertices: usize,
        max_cells: usize,
    },
}

impl Default for MergeStrategy {
    fn default() -> Self {
        // Eine einzelne große Zelle ist oft das Ziel für z.B. Inselgenerierung
        MergeStrategy::LargestCells { count: 1 }
    }
}

/// Verantwortlich für das Auswählen und Zusammenfügen von Voronoi-Zellen zu größeren Polygonen.
pub struct PolygonMerger {
    // Die VoronoiConfig wird hier nicht mehr direkt benötigt, da die Parameter
    // über die MergeStrategy oder andere Wege hereinkommen sollten.
    // config: VoronoiConfig, // Entfernt, um Abhängigkeit zu reduzieren
    merge_strategy: MergeStrategy,
    // Toleranz für geometrische Operationen, falls benötigt
    tolerance: f32,
}

impl PolygonMerger {
    pub fn new(strategy: MergeStrategy) -> Self {
        Self {
            merge_strategy: strategy,
            tolerance: constants::EPSILON * 10.0,
        }
    }

    // `with_strategy` ist nicht mehr nötig, wenn sie im Konstruktor gesetzt wird.
    // Falls doch:
    // pub fn with_strategy(mut self, strategy: MergeStrategy) -> Self {
    //     self.merge_strategy = strategy;
    //     self
    // }

    /// Wählt Voronoi-Zellen basierend auf der konfigurierten Strategie aus
    /// und versucht, sie zu einem oder mehreren Polygonen (als `Vec<Vec2>`) zusammenzufügen.
    /// Gibt eine Liste von Punktlisten zurück, wobei jede innere Liste ein resultierendes Polygon darstellt.
    pub fn merge_selected_cells(
        &self,
        available_cells: &[VoronoiCell],
        seed_resource: &mut SeedResource,
    ) -> MathResult<Vec<Vec<Vec2>>> {
        // Gibt Vec von Polygon-Punktlisten zurück
        if available_cells.is_empty() {
            return Ok(Vec::new()); // Keine Zellen zum Mergen
        }

        let selected_cells = self.select_cells_to_merge(available_cells, seed_resource)?;

        if selected_cells.is_empty() {
            return Ok(Vec::new()); // Keine Zellen nach Auswahl übrig
        }

        // Basierend auf der Strategie, wie die Zusammenführung erfolgen soll:
        match self.merge_strategy {
            MergeStrategy::UnionAllSelected => {
                // Versuche, alle ausgewählten Zellen zu einem Polygon zu verschmelzen.
                // Eine einfache Methode ist, die konvexe Hülle aller Vertices zu nehmen.
                // Für konkave Ergebnisse wäre ein echter Polygon-Union-Algorithmus nötig.
                let mut all_vertices_from_selected: Vec<Vec2> = Vec::new();
                for cell in selected_cells {
                    all_vertices_from_selected.extend_from_slice(&cell.vertices);
                }

                if all_vertices_from_selected.len() < 3 {
                    return Ok(Vec::new()); // Nicht genug Punkte für ein Polygon
                }

                let hull_computer = ConvexHullComputer::new(ConvexHullAlgorithm::AndrewMonotone)
                    .with_tolerance(self.tolerance);
                let merged_hull_points =
                    hull_computer.compute_hull_points(&all_vertices_from_selected)?;

                if merged_hull_points.len() >= 3 {
                    Ok(vec![merged_hull_points])
                } else {
                    Ok(Vec::new())
                }
            }
            // Andere Strategien könnten jede ausgewählte Zelle einzeln zurückgeben
            // oder spezifischere Merge-Logiken implementieren.
            // Fürs Erste: Die meisten Selektionsstrategien geben einfach die Konturen der ausgewählten Zellen zurück.
            _ => {
                let result_polygons: Vec<Vec<Vec2>> = selected_cells
                    .into_iter()
                    .map(|cell| cell.vertices.clone()) // Clone the vertices to avoid moving from reference
                    .filter(|vertices| vertices.len() >= 3) // Nur gültige Polygone
                    .collect();
                Ok(result_polygons)
            }
        }
    }

    /// Wählt Zellen basierend auf der konfigurierten Strategie aus.
    fn select_cells_to_merge<'a>(
        &self,
        all_cells: &'a [VoronoiCell],
        seed_resource: &mut SeedResource,
    ) -> MathResult<Vec<&'a VoronoiCell>> {
        // Gibt Referenzen auf die Zellen zurück
        if all_cells.is_empty() {
            return Ok(Vec::new());
        }

        match self.merge_strategy {
            MergeStrategy::AreaBased { min_area, max_area } => {
                Ok(all_cells
                    .iter()
                    .filter(|cell| {
                        let area = cell.area(); // area() von VoronoiCell
                        area >= min_area && area <= max_area
                    })
                    .collect())
            }
            MergeStrategy::NeighborhoodBased {
                max_vertices_for_cell,
            } => Ok(all_cells
                .iter()
                .filter(|cell| cell.vertices.len() <= max_vertices_for_cell)
                .collect()),
            MergeStrategy::RandomConnectedCells { min_vertices, max_cells } => {
                Self::select_random_connected_cells(all_cells, seed_resource, min_vertices, max_cells)
            }
            MergeStrategy::LargestCells { count } => {
                let mut sorted_cells: Vec<&VoronoiCell> = all_cells.iter().collect();
                sorted_cells.sort_unstable_by(|a, b| {
                    b.area()
                        .partial_cmp(&a.area())
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                Ok(sorted_cells.into_iter().take(count).collect())
            }
            MergeStrategy::PropertyBased {
                similarity_threshold,
            } => {
                // Beispiel: Wähle Zellen, deren Fläche nahe am Durchschnitt liegt.
                if all_cells.is_empty() {
                    return Ok(Vec::new());
                }
                let avg_area: f32 =
                    all_cells.iter().map(|c| c.area()).sum::<f32>() / all_cells.len() as f32;
                let threshold_val = avg_area * similarity_threshold;

                Ok(all_cells
                    .iter()
                    .filter(|cell| (cell.area() - avg_area).abs() <= threshold_val)
                    .collect())
            }
            MergeStrategy::Random { target_count } => {
                if target_count == 0 {
                    return Ok(Vec::new());
                }
                // Wenn target_count größer/gleich der Anzahl der Zellen ist, gib alle zurück.
                // Erstelle eine Kopie der Referenzen, damit sie gemischt werden können,
                // ohne die ursprüngliche Reihenfolge von all_cells zu ändern, falls das wichtig ist.
                if target_count >= all_cells.len() {
                    return Ok(all_cells.iter().collect());
                }

                // Erstelle einen Vektor von Indizes [0, 1, ..., all_cells.len() - 1]
                let mut indices: Vec<usize> = (0..all_cells.len()).collect();

                // Mische die Indizes mit deiner SeedResource.
                // Du kannst optional eine Kategorie für das Mischen festlegen,
                // um die Auswahl deterministisch für diesen spezifischen "Random"-Merge zu machen.
                seed_resource.shuffle(&mut indices, None);

                // Nimm die ersten target_count Indizes aus dem gemischten Vektor
                // und mappe sie zurück auf die Zellen-Referenzen.
                let selected_cells: Vec<&'a VoronoiCell> = indices
                    .into_iter()
                    .take(target_count)
                    .map(|idx| &all_cells[idx]) // Hole die Referenz auf die Zelle am ausgewählten Index
                    .collect();

                Ok(selected_cells)
            }
            MergeStrategy::UnionAllSelected => {
                // Bei UnionAllSelected werden *alle* initialen Zellen als Kandidaten für die Vereinigung betrachtet,
                // es sei denn, man möchte hier noch eine Vorfilterung.
                // Für diese Implementierung nehmen wir an, dass `UnionAllSelected` bedeutet,
                // dass die Auswahl bereits getroffen wurde und wir sie nur noch vereinigen.
                // Wenn `select_cells_to_merge` für `UnionAllSelected` aufgerufen wird,
                // könnte es einfach alle Zellen zurückgeben.
                Ok(all_cells.iter().collect())
            }
        }
    }

    /// Implementiert die RandomConnectedCells Strategie:
    /// 1. Wählt zufällig eine Zelle mit mindestens min_vertices Vertices
    /// 2. Erweitert die Auswahl schrittweise auf benachbarte Zellen bis max_cells erreicht ist
    fn select_random_connected_cells<'a>(
        all_cells: &'a [VoronoiCell],
        seed_resource: &mut SeedResource,
        min_vertices: usize,
        max_cells: usize,
    ) -> MathResult<Vec<&'a VoronoiCell>> {
        if all_cells.is_empty() || max_cells == 0 {
            return Ok(Vec::new());
        }

        // Finde alle Zellen mit mindestens min_vertices Vertices
        let candidate_cells: Vec<&VoronoiCell> = all_cells
            .iter()
            .filter(|cell| cell.vertices.len() >= min_vertices)
            .collect();

        if candidate_cells.is_empty() {
            return Ok(Vec::new());
        }

        // Wähle zufällig eine Startzelle
        let start_index = seed_resource.next_i32_in_range(0, candidate_cells.len() as i32, None) as usize;
        let start_cell = candidate_cells[start_index];

        let mut selected_cells = Vec::new();
        let mut visited = HashSet::new();
        let mut expansion_candidates = Vec::new();
        
        // Füge die Startzelle hinzu
        selected_cells.push(start_cell);
        visited.insert(start_cell.id);
        expansion_candidates.extend(start_cell.neighbor_ids.iter());

        // Erweitere die Auswahl auf benachbarte Zellen
        while selected_cells.len() < max_cells && !expansion_candidates.is_empty() {
            // Entferne bereits besuchte IDs aus den Kandidaten
            expansion_candidates.retain(|&id| !visited.contains(&id));
            
            if expansion_candidates.is_empty() {
                break;
            }

            // Wähle zufällig eine der verfügbaren Nachbarzellen
            let random_index = seed_resource.next_i32_in_range(0, expansion_candidates.len() as i32, None) as usize;
            let selected_neighbor_id = expansion_candidates.swap_remove(random_index);

            // Finde die entsprechende Zelle
            if let Some(neighbor_cell) = all_cells.iter().find(|c| c.id == selected_neighbor_id) {
                selected_cells.push(neighbor_cell);
                visited.insert(neighbor_cell.id);
                
                // Füge die Nachbarn dieser Zelle als neue Kandidaten hinzu
                expansion_candidates.extend(neighbor_cell.neighbor_ids.iter());
            }
        }

        Ok(selected_cells)
    }
}

// Die alten `MergeUtils` und `MergeQuality` Strukturen könnten in ein
// separates Analyse- oder Utility-Modul für Polygone oder Voronoi-Diagramme passen,
// da sie nicht direkt Teil des Merging-Prozesses sind, sondern diesen bewerten.
