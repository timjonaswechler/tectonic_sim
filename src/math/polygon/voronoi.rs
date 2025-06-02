// In src/math/polygon/voronoi_merged_shape.rs

use bevy::log::warn;
// Notwendige Imports, die wir später im Modul brauchen werden
use rand::Rng; // Für Zufallsauswahl
use rand::SeedableRng; // Für seed_from_u64
use spade::handles::VertexHandle;
use spade::{DelaunayTriangulation, Point2, Triangulation};
use std::collections::{HashMap, HashSet, VecDeque};

// Wir werden Point2<f64> von spade verwenden. Wenn dein Projekt einen eigenen Point2-Typ hat,
// müssten wir hier Konvertierungen vornehmen oder den Typ generisch machen.
// Fürs Erste bleiben wir bei spade::Point2.

// Type alias für interne Verwendung
type InternalPointKey = (i64, i64); // Um Namenskonflikte zu vermeiden, falls PointKey woanders definiert ist

// Hilfsfunktion: vertex_to_key (intern)
fn vertex_to_internal_key(vertex: Point2<f64>) -> InternalPointKey {
    (
        (vertex.x * 1_000_000.0).round() as i64,
        (vertex.y * 1_000_000.0).round() as i64,
    )
}

// Hilfsfunktion: calculate_polygon_centroid (intern)
// (Implementierung von calculate_polygon_centroid hier einfügen)
fn calculate_polygon_centroid(polygon_vertices: &[Point2<f64>]) -> Option<Point2<f64>> {
    if polygon_vertices.len() < 1 {
        return None;
    }
    if polygon_vertices.len() == 1 {
        return Some(polygon_vertices[0]);
    }
    if polygon_vertices.len() == 2 {
        return Some(Point2::new(
            (polygon_vertices[0].x + polygon_vertices[1].x) / 2.0,
            (polygon_vertices[0].y + polygon_vertices[1].y) / 2.0,
        ));
    }
    let mut area = 0.0;
    let mut centroid_x = 0.0;
    let mut centroid_y = 0.0;
    let n = polygon_vertices.len();
    for i in 0..n {
        let p1 = polygon_vertices[i];
        let p2 = polygon_vertices[(i + 1) % n];
        let cross_product = (p1.x * p2.y) - (p2.x * p1.y);
        area += cross_product;
        centroid_x += (p1.x + p2.x) * cross_product;
        centroid_y += (p1.y + p2.y) * cross_product;
    }
    if area.abs() < 1e-9 {
        let sum_x: f64 = polygon_vertices.iter().map(|p| p.x).sum();
        let sum_y: f64 = polygon_vertices.iter().map(|p| p.y).sum();
        return Some(Point2::new(sum_x / n as f64, sum_y / n as f64));
    }
    area *= 0.5;
    centroid_x /= 6.0 * area;
    centroid_y /= 6.0 * area;
    Some(Point2::new(centroid_x, centroid_y))
}

// Hilfsfunktion: center_and_normalize_polygon (intern)
// (Implementierung von center_and_normalize_polygon hier einfügen,
//  stellt sicher, dass sie calculate_polygon_centroid verwendet)
fn core_center_and_normalize_polygon(polygon: &[Point2<f64>]) -> Vec<Point2<f64>> {
    if polygon.len() < 2 {
        if polygon.len() == 1 {
            return vec![Point2::new(0.0, 0.0)];
        } // Einzelner Punkt wird zu (0,0)
        return Vec::new(); // Leeres Polygon oder ungültig
    }

    let centroid = match calculate_polygon_centroid(polygon) {
        Some(c) => c,
        None => {
            // Sollte nicht passieren, wenn polygon.len() >= 2
            let sum_x: f64 = polygon.iter().map(|p| p.x).sum();
            let sum_y: f64 = polygon.iter().map(|p| p.y).sum();
            Point2::new(sum_x / polygon.len() as f64, sum_y / polygon.len() as f64)
        }
    };

    let centered_polygon: Vec<Point2<f64>> = polygon
        .iter()
        .map(|p| Point2::new(p.x - centroid.x, p.y - centroid.y))
        .collect();

    let mut max_abs_coord = 0.0f64;
    for p in &centered_polygon {
        max_abs_coord = max_abs_coord.max(p.x.abs()).max(p.y.abs());
    }

    if max_abs_coord.abs() < 1e-9 {
        // Alle Punkte sind (0,0)
        return centered_polygon
            .iter()
            .map(|_| Point2::new(0.0, 0.0))
            .collect(); // Gib (0,0)-Punkte zurück
    }

    centered_polygon
        .iter()
        .map(|p| Point2::new(p.x / max_abs_coord, p.y / max_abs_coord))
        .collect()
}

// Hilfsfunktion: apply_chaikin_smoothing (intern)
// (Implementierung von apply_chaikin_smoothing hier einfügen)
fn core_apply_chaikin_smoothing(polygon: &[Point2<f64>], iterations: usize) -> Vec<Point2<f64>> {
    if polygon.len() < 3 || iterations == 0 {
        return polygon.to_vec();
    }
    let mut current_polygon_points = polygon.to_vec();
    for _iter in 0..iterations {
        if current_polygon_points.len() < 3 {
            break;
        }
        let mut smoothed_polygon_segment_endpoints = Vec::new();
        let n = current_polygon_points.len();
        if n == 0 {
            break;
        }
        let t = 0.25;
        for i in 0..n {
            let vertex_i_minus_1 = current_polygon_points[(i + n - 1) % n];
            let vertex_i = current_polygon_points[i];
            let vertex_i_plus_1 = current_polygon_points[(i + 1) % n];
            let new_point_on_edge_to_prev_x = (1.0 - t) * vertex_i.x + t * vertex_i_minus_1.x;
            let new_point_on_edge_to_prev_y = (1.0 - t) * vertex_i.y + t * vertex_i_minus_1.y;
            smoothed_polygon_segment_endpoints.push(Point2::new(
                new_point_on_edge_to_prev_x,
                new_point_on_edge_to_prev_y,
            ));
            let new_point_on_edge_to_next_x = (1.0 - t) * vertex_i.x + t * vertex_i_plus_1.x;
            let new_point_on_edge_to_next_y = (1.0 - t) * vertex_i.y + t * vertex_i_plus_1.y;
            smoothed_polygon_segment_endpoints.push(Point2::new(
                new_point_on_edge_to_next_x,
                new_point_on_edge_to_next_y,
            ));
        }
        current_polygon_points = smoothed_polygon_segment_endpoints;
    }
    current_polygon_points
}

// ... (weitere interne Hilfsfunktionen wie get_voronoi_cell_vertices, is_voronoi_cell_complete etc.
//      werden hier folgen. Ihre Namen könnten `core_` als Präfix bekommen, um sie als interne
//      Implementierungsdetails zu kennzeichnen oder einfach klein geschrieben bleiben, wenn sie
//      nicht exportiert werden.)

// --- Kernlogikfunktionen ---
// Diese Funktionen arbeiten auf einem generischen Koordinatensystem und werden von der Hauptfunktion aufgerufen.

fn core_run_lloyds_algorithm(
    initial_points: &[Point2<f64>],
    num_iterations: usize,
    // Die Canvas-Grenzen sind relativ zum Koordinatensystem der `initial_points`
    boundary_min: Point2<f64>,
    boundary_max: Point2<f64>,
) -> DelaunayTriangulation<Point2<f64>> {
    let mut current_triangulation = DelaunayTriangulation::<Point2<f64>>::new();
    for p in initial_points {
        current_triangulation.insert(*p).ok();
    }

    for _iter in 0..num_iterations {
        if current_triangulation.num_vertices() == 0 {
            break;
        }

        let mut next_generation_points = Vec::with_capacity(current_triangulation.num_vertices());

        for vertex_handle in current_triangulation.vertices() {
            let old_pos = vertex_handle.position();
            // core_get_voronoi_cell_vertices wird `boundary_min/max` nicht direkt verwenden,
            // aber `core_is_voronoi_cell_complete` wird es.
            let voronoi_cell_polygon =
                core_get_voronoi_cell_vertices(vertex_handle, boundary_min, boundary_max);

            if !core_is_voronoi_cell_complete(&voronoi_cell_polygon, boundary_min, boundary_max)
                || voronoi_cell_polygon.len() < 3
            {
                next_generation_points.push(old_pos);
                continue;
            }

            if let Some(centroid) = calculate_polygon_centroid(&voronoi_cell_polygon) {
                let new_x = centroid
                    .x
                    .max(boundary_min.x + 1e-3)
                    .min(boundary_max.x - 1e-3);
                let new_y = centroid
                    .y
                    .max(boundary_min.y + 1e-3)
                    .min(boundary_max.y - 1e-3);
                next_generation_points.push(Point2::new(new_x, new_y));
            } else {
                next_generation_points.push(old_pos);
            }
        }

        if next_generation_points.is_empty() {
            break;
        }

        current_triangulation = DelaunayTriangulation::<Point2<f64>>::new();
        let mut unique_inserted_keys = HashSet::new();
        for p_new in next_generation_points {
            if unique_inserted_keys.insert(vertex_to_internal_key(p_new)) {
                current_triangulation.insert(p_new).ok();
            }
        }
    }
    current_triangulation
}

// Diese Funktion benötigt jetzt die Boundaries, um die Voronoi-Zellen ggf. zu clippen oder
// um Circumcenter, die weit außerhalb liegen (Zeichen für Randzellen), zu behandeln.
// Für eine Version, die auf [-1,1] normalisierten Daten arbeitet, könnten die Boundaries implizit sein.
// Aber für die *korrekte* Voronoi-Zell-Erzeugung, besonders für Randzellen,
// ist das Clipping eigentlich notwendig. `spade` gibt uns Circumcenter, die sehr weit weg sein können.
fn core_get_voronoi_cell_vertices(
    vertex_handle: spade::handles::VertexHandle<Point2<f64>>,
    _boundary_min: Point2<f64>, // Vorerst nicht direkt für Clipping verwendet, aber für Kontext
    _boundary_max: Point2<f64>,
) -> Vec<Point2<f64>> {
    let mut voronoi_vertices_map = HashMap::new();
    let generator_pos = vertex_handle.position();

    for edge in vertex_handle.out_edges() {
        let face = edge.face();
        if face.is_outer() {
            continue;
        }
        // Hole die drei Eckpunkte des Dreiecks (face)
        let mut triangle_points = Vec::new();
        if let Some(starting_edge) = face.adjacent_edge() {
            let mut current_edge = starting_edge;
            loop {
                triangle_points.push(current_edge.from().position());
                current_edge = current_edge.next();
                if current_edge == starting_edge || triangle_points.len() > 3 {
                    break;
                }
            }
        }
        if triangle_points.len() == 3 {
            let p1 = triangle_points[0];
            let p2 = triangle_points[1];
            let p3 = triangle_points[2];
            if let Some(circumcenter) = calculate_circumcenter(p1, p2, p3) {
                let key = (
                    (circumcenter.x * 1e7).round() as i64,
                    (circumcenter.y * 1e7).round() as i64,
                );
                voronoi_vertices_map.insert(key, circumcenter);
            }
        }
    }

    let mut voronoi_vertices: Vec<Point2<f64>> = voronoi_vertices_map.values().cloned().collect();

    voronoi_vertices.sort_by(|a, b| {
        let angle_a = (a.y - generator_pos.y).atan2(a.x - generator_pos.x);
        let angle_b = (b.y - generator_pos.y).atan2(b.x - generator_pos.x);
        angle_a
            .partial_cmp(&angle_b)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    voronoi_vertices
}

fn core_is_voronoi_cell_complete(
    voronoi_vertices: &[Point2<f64>],
    boundary_min: Point2<f64>,
    boundary_max: Point2<f64>,
) -> bool {
    if voronoi_vertices.len() < 3 {
        return false;
    }
    // Eine Zelle gilt als "nicht komplett" (im Sinne von unendlich), wenn einer ihrer
    // Vertices sehr weit von den definierten Grenzen entfernt ist.
    let large_factor = 5.0; // Wie weit außerhalb ist "zu weit"
    let margin_x = (boundary_max.x - boundary_min.x) * large_factor;
    let margin_y = (boundary_max.y - boundary_min.y) * large_factor;

    for vertex in voronoi_vertices {
        if vertex.x < boundary_min.x - margin_x
            || vertex.x > boundary_max.x + margin_x
            || vertex.y < boundary_min.y - margin_y
            || vertex.y > boundary_max.y + margin_y
        {
            return false; // Vertex ist extrem weit außerhalb
        }
    }
    true
}

fn core_find_inner_generator_points(
    triangulation: &DelaunayTriangulation<Point2<f64>>,
    boundary_min: Point2<f64>,
    boundary_max: Point2<f64>,
    min_dist_to_boundary_factor: f64, // z.B. 0.1 für 10% der kürzeren Seite
) -> Vec<Point2<f64>> {
    let mut inner_points = Vec::new();
    let width = boundary_max.x - boundary_min.x;
    let height = boundary_max.y - boundary_min.y;
    let min_abs_distance = width.min(height) * min_dist_to_boundary_factor;

    for vertex_handle in triangulation.vertices() {
        let pos = vertex_handle.position();

        let dist_ok = pos.x > boundary_min.x + min_abs_distance
            && pos.x < boundary_max.x - min_abs_distance
            && pos.y > boundary_min.y + min_abs_distance
            && pos.y < boundary_max.y - min_abs_distance;

        let voronoi_cell_polygon =
            core_get_voronoi_cell_vertices(vertex_handle, boundary_min, boundary_max);
        let cell_complete =
            core_is_voronoi_cell_complete(&voronoi_cell_polygon, boundary_min, boundary_max);
        let neighbor_count = vertex_handle.out_edges().count();

        if dist_ok && cell_complete && neighbor_count >= 5 {
            inner_points.push(pos);
        }
    }
    inner_points
}

fn core_find_delaunay_neighbors(
    triangulation: &DelaunayTriangulation<Point2<f64>>,
    target_generator_pos: Point2<f64>,
) -> Vec<Point2<f64>> {
    let mut neighbors = Vec::new();
    let tolerance = 1e-9;
    // Finde den VertexHandle. In einer größeren Anwendung wäre eine Map von Pos zu Handle effizienter.
    if let Some(vh) = triangulation.vertices().find(|v| {
        let vp = v.position();
        (vp.x - target_generator_pos.x).abs() < tolerance
            && (vp.y - target_generator_pos.y).abs() < tolerance
    }) {
        for edge in vh.out_edges() {
            neighbors.push(edge.to().position());
        }
    }
    neighbors
}

fn core_select_connected_inner_cells(
    triangulation: &DelaunayTriangulation<Point2<f64>>,
    inner_generator_points: &[Point2<f64>],
    target_size: usize,
    rng: &mut impl Rng, // Generischer Rng
) -> Vec<Point2<f64>> {
    if inner_generator_points.is_empty() || target_size == 0 {
        return Vec::new();
    }

    let start_index = rng.random_range(0..inner_generator_points.len());
    let start_point = inner_generator_points[start_index];

    let mut selected_points = vec![start_point];
    let mut candidates_queue = VecDeque::new();
    let mut visited_keys = HashSet::new();

    visited_keys.insert(vertex_to_internal_key(start_point));

    for neighbor_pos in core_find_delaunay_neighbors(triangulation, start_point) {
        // Prüfe, ob der Nachbar auch ein "inner_generator_point" ist
        if inner_generator_points
            .iter()
            .any(|p| (p.x - neighbor_pos.x).abs() < 1e-9 && (p.y - neighbor_pos.y).abs() < 1e-9)
        {
            if visited_keys.insert(vertex_to_internal_key(neighbor_pos)) {
                candidates_queue.push_back(neighbor_pos);
            }
        }
    }

    while selected_points.len() < target_size && !candidates_queue.is_empty() {
        let candidate_idx_in_queue = rng.random_range(0..candidates_queue.len());
        let chosen_point = candidates_queue.remove(candidate_idx_in_queue).unwrap();
        selected_points.push(chosen_point);

        for neighbor_pos in core_find_delaunay_neighbors(triangulation, chosen_point) {
            if inner_generator_points
                .iter()
                .any(|p| (p.x - neighbor_pos.x).abs() < 1e-9 && (p.y - neighbor_pos.y).abs() < 1e-9)
            {
                if visited_keys.insert(vertex_to_internal_key(neighbor_pos)) {
                    candidates_queue.push_back(neighbor_pos);
                }
            }
        }
    }
    selected_points
}

fn core_merge_selected_voronoi_polygons(
    triangulation: &DelaunayTriangulation<Point2<f64>>,
    selected_generator_points: &[Point2<f64>],
    boundary_min: Point2<f64>,
    boundary_max: Point2<f64>,
) -> Vec<Point2<f64>> {
    if selected_generator_points.len() < 1 {
        return Vec::new();
    }

    let mut generator_point_to_handle_map: HashMap<InternalPointKey, VertexHandle<Point2<f64>>> =
        HashMap::new();
    for vh in triangulation.vertices() {
        generator_point_to_handle_map.insert(vertex_to_internal_key(vh.position()), vh);
    }

    let valid_selected_handles: Vec<VertexHandle<Point2<f64>>> = selected_generator_points
        .iter()
        .filter_map(|p| generator_point_to_handle_map.get(&vertex_to_internal_key(*p)))
        .cloned()
        .collect();

    if valid_selected_handles.is_empty() {
        return Vec::new();
    }

    let mut edge_counts: HashMap<(InternalPointKey, InternalPointKey), usize> = HashMap::new();
    let mut point_key_to_original_point: HashMap<InternalPointKey, Point2<f64>> = HashMap::new();

    for vertex_handle in &valid_selected_handles {
        let cell_polygon_points =
            core_get_voronoi_cell_vertices(*vertex_handle, boundary_min, boundary_max);
        if cell_polygon_points.len() < 3 {
            continue;
        }

        for i in 0..cell_polygon_points.len() {
            let p1 = cell_polygon_points[i];
            let p2 = cell_polygon_points[(i + 1) % cell_polygon_points.len()];
            let key1 = vertex_to_internal_key(p1);
            let key2 = vertex_to_internal_key(p2);
            point_key_to_original_point.entry(key1).or_insert(p1);
            point_key_to_original_point.entry(key2).or_insert(p2);
            let edge = if key1 < key2 {
                (key1, key2)
            } else {
                (key2, key1)
            };
            *edge_counts.entry(edge).or_insert(0) += 1;
        }
    }

    let mut boundary_graph: HashMap<InternalPointKey, Vec<InternalPointKey>> = HashMap::new();
    let mut all_boundary_keys = HashSet::new();

    for ((k1, k2), count) in edge_counts {
        if count == 1 {
            boundary_graph.entry(k1).or_default().push(k2);
            boundary_graph.entry(k2).or_default().push(k1);
            all_boundary_keys.insert(k1);
            all_boundary_keys.insert(k2);
        }
    }

    if boundary_graph.is_empty() || all_boundary_keys.is_empty() {
        return Vec::new();
    }

    let mut ordered_polygon_keys = Vec::new();
    let mut visited_points_in_path = HashSet::new();
    let mut current_key = *all_boundary_keys.iter().next().unwrap();
    let path_start_key = current_key;

    // Pfadfindungslogik (vereinfacht für Übersicht, die komplexere Version mit Winkel-Sortierung ist besser)
    for _ in 0..all_boundary_keys.len() + 2 {
        ordered_polygon_keys.push(current_key);
        visited_points_in_path.insert(current_key);

        let neighbors = match boundary_graph.get(&current_key) {
            Some(n) => n,
            None => break,
        };
        let mut next_key_opt = None;
        for &neighbor_key in neighbors {
            if neighbor_key == path_start_key && ordered_polygon_keys.len() > 2 {
                next_key_opt = Some(neighbor_key);
                break;
            }
            if !visited_points_in_path.contains(&neighbor_key) {
                next_key_opt = Some(neighbor_key);
                break;
            }
        }
        if let Some(next_key) = next_key_opt {
            current_key = next_key;
        } else {
            break;
        }
        if current_key == path_start_key && ordered_polygon_keys.len() > 2 {
            break;
        }
        if ordered_polygon_keys.len() >= all_boundary_keys.len() + 1
            && current_key != path_start_key
        {
            return Vec::new();
        } // Fehler: Pfad zu lang, nicht geschlossen
    }
    if ordered_polygon_keys.last() != Some(&path_start_key) && current_key == path_start_key {
        // Polygon geschlossen, aber der Startpunkt wurde nicht als letztes Element hinzugefügt, das ist ok.
    } else if ordered_polygon_keys.last() != Some(&path_start_key)
        || ordered_polygon_keys.len() <= 2
    {
        return Vec::new(); // Nicht korrekt geschlossen oder zu klein
    }

    ordered_polygon_keys
        .into_iter()
        .map(|key| *point_key_to_original_point.get(&key).unwrap())
        .collect()
}

/// Configuration for generating the Voronoi merged shape.
#[derive(Debug, Clone)]
pub struct VoronoiMergedShapeConfig {
    /// Number of initial random points to generate in the [-1, 1] normalized square.
    pub num_initial_points: usize,
    /// Number of Lloyd's relaxation iterations.
    pub num_lloyd_iterations: usize,
    /// Target number of connected inner cells to select for merging.
    pub target_connected_cell_group_size: usize,
    /// Number of Chaikin smoothing iterations to apply to the final merged polygon.
    pub chaikin_smoothing_iterations: usize,
    /// Factor to determine "inner" points (e.g., 0.1 means 10% from boundary).
    pub inner_point_boundary_factor: f64,
    // Seed für den RNG, um reproduzierbare Ergebnisse zu ermöglichen (optional)
    pub rng_seed: Option<u64>,
}

impl Default for VoronoiMergedShapeConfig {
    fn default() -> Self {
        let mut thread_rng = rand::rng();
        let mut rng = rand::rngs::StdRng::from_rng(&mut thread_rng);

        Self {
            num_initial_points: rng.random_range(100..=1000), // Default zwischen 100 und 1000
            num_lloyd_iterations: 10,
            target_connected_cell_group_size: rng.random_range(3..=7),
            chaikin_smoothing_iterations: 1,
            inner_point_boundary_factor: 0.15, // 15% vom Rand nach "inner"
            rng_seed: None,
        }
    }
}

/// Generates a smoothed, merged Voronoi shape, normalized to the [-1, 1] range.
///
/// The process involves:
/// 1. Generating initial random points in a [-1, 1] square.
/// 2. Applying Lloyd's relaxation to these points.
/// 3. Identifying "inner" Voronoi cells based on their distance from the [-1, 1] boundary.
/// 4. Selecting a connected group of these inner cells.
/// 5. Merging the Voronoi polygons of these selected cells into a single polygon.
/// 6. Applying Chaikin smoothing to this merged polygon.
/// 7. Centering and normalizing the final smoothed polygon to the [-1, 1] range.
///
/// Returns a `Vec<Point2<f64>>` representing the vertices of the final polygon,
/// or an empty Vec if a polygon could not be successfully generated.
pub fn generate_normalized_voronoi_merged_shape(
    config: &VoronoiMergedShapeConfig,
) -> Vec<Point2<f64>> {
    // Initialisiere RNG
    let mut rng = if let Some(seed) = config.rng_seed {
        rand::rngs::StdRng::seed_from_u64(seed)
    } else {
        let mut thread_rng = rand::rng();
        rand::rngs::StdRng::from_rng(&mut thread_rng)
    };
    if config.num_initial_points < 100 {
        warn!(
            "VoronoiMergedShapeConfig: num_initial_points should be at least 100 for triangulation."
        );
        return Vec::new();
    }
    if config.target_connected_cell_group_size < 1 {
        warn!("VoronoiMergedShapeConfig: target_connected_cell_group_size should be at least 1.");
        return Vec::new();
    }
    if config.chaikin_smoothing_iterations < 1 {
        warn!("VoronoiMergedShapeConfig: chaikin_smoothing_iterations should be at least 1.");
        return Vec::new();
    }
    if config.inner_point_boundary_factor <= 0.0 || config.inner_point_boundary_factor >= 1.0 {
        warn!(
            "VoronoiMergedShapeConfig: inner_point_boundary_factor should be between 0.0 and 1.0."
        );
        return Vec::new();
    }
    if config.num_lloyd_iterations < 1 {
        warn!("VoronoiMergedShapeConfig: num_lloyd_iterations should be at least 1.");
        return Vec::new();
    }
    if config.target_connected_cell_group_size < 1 {
        warn!("VoronoiMergedShapeConfig: target_connected_cell_group_size should be at least 1.");
        return Vec::new();
    }
    // 1. Erzeuge initiale Zufallspunkte im Bereich [-1, 1] x [-1, 1]
    let mut initial_points = Vec::new();
    let boundary_min = Point2::new(-1.0, -1.0);
    let boundary_max = Point2::new(1.0, 1.0);

    for _ in 0..config.num_initial_points {
        let x = rng.random_range(boundary_min.x..boundary_max.x);
        let y = rng.random_range(boundary_min.y..boundary_max.y);
        initial_points.push(Point2::new(x, y));
    }
    if initial_points.len() < 3 {
        return Vec::new();
    }

    // 2. Lloyd's Algorithmus anwenden (arbeitet auf den [-1,1] koordinierten Punkten)
    let triangulation = core_run_lloyds_algorithm(
        &initial_points,
        config.num_lloyd_iterations,
        boundary_min,
        boundary_max,
    );
    if triangulation.num_vertices() < 3 {
        return Vec::new();
    }

    // 3. Innere Generatorpunkte finden
    let inner_generator_points = core_find_inner_generator_points(
        &triangulation,
        boundary_min,
        boundary_max,
        config.inner_point_boundary_factor,
    );
    if inner_generator_points.is_empty() {
        return Vec::new();
    }

    // 4. Verbundene innere Zellen auswählen
    let selected_generator_points = core_select_connected_inner_cells(
        &triangulation,
        &inner_generator_points,
        config.target_connected_cell_group_size,
        &mut rng,
    );
    if selected_generator_points.is_empty() {
        return Vec::new();
    }

    // 5. Ausgewählte Voronoi-Polygone verschmelzen
    let mut merged_polygon_vertices = core_merge_selected_voronoi_polygons(
        &triangulation,
        &selected_generator_points,
        boundary_min,
        boundary_max,
    );
    if merged_polygon_vertices.len() < 3 {
        return Vec::new();
    }

    // 6. Chaikin Smoothing
    merged_polygon_vertices = core_apply_chaikin_smoothing(
        &merged_polygon_vertices,
        config.chaikin_smoothing_iterations,
    );
    if merged_polygon_vertices.len() < 3 {
        return Vec::new();
    }

    // 7. Zentrieren und Normalisieren des finalen Polygons
    // Da alle vorherigen Schritte bereits im [-1,1]-Raum (oder relativ dazu) stattfanden,
    // sollte das `merged_polygon_vertices` bereits ungefähr in diesem Raum liegen.
    // Die Normalisierung stellt sicher, dass es exakt zentriert und auf [-1,1] skaliert ist.
    core_center_and_normalize_polygon(&merged_polygon_vertices)
}

fn calculate_circumcenter(
    p1: Point2<f64>,
    p2: Point2<f64>,
    p3: Point2<f64>,
) -> Option<Point2<f64>> {
    let ax = p1.x;
    let ay = p1.y;
    let bx = p2.x;
    let by = p2.y;
    let cx = p3.x;
    let cy = p3.y;

    let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));

    if d.abs() < 1e-10 {
        return None; // Punkte sind kollinear
    }

    let ux = ((ax * ax + ay * ay) * (by - cy)
        + (bx * bx + by * by) * (cy - ay)
        + (cx * cx + cy * cy) * (ay - by))
        / d;
    let uy = ((ax * ax + ay * ay) * (cx - bx)
        + (bx * bx + by * by) * (ax - cx)
        + (cx * cx + cy * cy) * (bx - ax))
        / d;

    Some(Point2::new(ux, uy))
}
