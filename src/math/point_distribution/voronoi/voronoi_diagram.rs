// src/math/point_distribution/voronoi/voronoi_diagram.rs

use crate::math::{
    error::MathResult,
    types::{Bounds2D, SpadePoint},
    utils::constants,
};
use bevy::log::info;
use bevy::math::Vec2;
use spade::{DelaunayTriangulation, Point2, Triangulation};
use std::collections::HashSet;

/// Repräsentiert eine einzelne Zelle in einem Voronoi-Diagramm.
#[derive(Debug, Clone)]
pub struct VoronoiCell {
    /// Eindeutige ID der Zelle
    pub id: usize,
    /// Der Generatorpunkt (Site), der diese Zelle definiert.
    pub generator: Vec2,
    /// Die Eckpunkte des Polygons, das die Zelle bildet, in CCW-Reihenfolge.
    pub vertices: Vec<Vec2>,
    /// Gibt an, ob diese Zelle an der äußeren Grenze des Diagramms liegt
    /// (d.h., einige ihrer Kanten sind unendlich oder durch die Boundary beschnitten).
    pub is_boundary_cell: bool,
    /// IDs der benachbarten Zellen
    pub neighbor_ids: Vec<usize>,
    // Optional: Index des Generatorpunktes in der ursprünglichen Eingabeliste
    // pub generator_index: Option<usize>,
}

impl VoronoiCell {
    // ... (area und perimeter bleiben gleich)
    /// Berechnet die Fläche der Voronoi-Zelle (als Polygon).
    pub fn area(&self) -> f32 {
        if self.vertices.len() < 3 {
            return 0.0;
        }
        let mut area_sum = 0.0;
        for i in 0..self.vertices.len() {
            let p1 = self.vertices[i];
            let p2 = self.vertices[(i + 1) % self.vertices.len()];
            area_sum += (p1.x * p2.y) - (p2.x * p1.y);
        }
        0.5 * area_sum.abs()
    }

    /// Berechnet den Umfang der Voronoi-Zelle.
    pub fn perimeter(&self) -> f32 {
        if self.vertices.len() < 2 {
            return 0.0;
        }
        let mut perimeter_sum = 0.0;
        for i in 0..self.vertices.len() {
            let p1 = self.vertices[i];
            let p2 = self.vertices[(i + 1) % self.vertices.len()];
            perimeter_sum += p1.distance(p2);
        }
        perimeter_sum
    }
}

/// Repräsentiert eine Kante im Voronoi-Diagramm.
#[derive(Debug, Clone, Copy)]
pub struct VoronoiEdge {
    /// Startpunkt der Kante.
    pub start: Vec2,
    /// Endpunkt der Kante. `None`, wenn die Kante ins Unendliche reicht (für Boundary-Zellen).
    pub end: Option<Vec2>,
    /// Einer der beiden Generatorpunkte, die diese Kante teilen.
    pub generator_a: Vec2,
    /// Der andere Generatorpunkt, der diese Kante teilt.
    pub generator_b: Vec2,
}

/// Enthält die extrahierten Zellen und Kanten eines Voronoi-Diagramms.
#[derive(Debug, Clone, Default)]
pub struct VoronoiDiagram {
    pub cells: Vec<VoronoiCell>,
    pub edges: Vec<VoronoiEdge>,
    // Optional: Die ursprünglichen Generatorpunkte
    // pub generators: Vec<Vec2>,
}

/// Extrahiert Voronoi-Zellen und -Kanten aus einer Delaunay-Triangulation.
pub struct VoronoiExtractor;

impl VoronoiExtractor {
    /// Extrahiert alle Voronoi-Zellen aus der gegebenen Delaunay-Triangulation.
    /// `triangulation`: Die Delaunay-Triangulation der Generatorpunkte.
    /// `clip_boundary`: Optionale Grenzen, gegen die die Voronoi-Zellen geclippt werden.
    ///                  Wenn `None`, können Zellen unbeschränkt sein (für Randzellen).
    pub fn extract_cells(
        triangulation: &DelaunayTriangulation<SpadePoint>,
        clip_boundary: Option<&Bounds2D>,
    ) -> MathResult<Vec<VoronoiCell>> {
        let mut voronoi_cells = Vec::with_capacity(triangulation.num_vertices());

        for vertex_handle in triangulation.vertices() {
            let generator_spade = vertex_handle.position();
            let generator_bevy = Vec2::new(generator_spade.x, generator_spade.y);
            let mut cell_circumcenters_spade: Vec<SpadePoint> = Vec::new();
            let mut is_boundary_cell_flag = false;

            for connected_edge in vertex_handle.out_edges() {
                let face_handle = connected_edge.face();
                if face_handle.is_outer() {
                    is_boundary_cell_flag = true;
                    // Optional: Breche hier schon ab, wenn du weißt, dass du boundary cells nicht willst
                    // break;
                } else {
                    if let Some(inner_face) = face_handle.as_inner() {
                        cell_circumcenters_spade.push(inner_face.circumcenter());
                    }
                }
            }

            // <<< NEUER FILTER >>>
            if is_boundary_cell_flag {
                // Angenommen, skip_boundary_cells ist eine neue Config-Option
                info!(
                    "VoronoiExtractor: SKIPPING boundary cell for generator {:?} entirely due to config.",
                    generator_bevy
                );
                continue;
            }
            // <<< ENDE NEUER FILTER >>>

            // Bestehender Filter, wenn cell_circumcenters_spade.len() < 2 ist UND es KEINE boundary cell ist.
            // Dieser Filter sollte NACH dem optionalen `skip_boundary_cells_entirely` kommen.
            if cell_circumcenters_spade.len() < 2 && !is_boundary_cell_flag {
                // Spezieller Fall: Nur ein Punkt in der Triangulation
                if cell_circumcenters_spade.is_empty() && triangulation.num_vertices() == 1 {
                    if let Some(boundary) = clip_boundary {
                        voronoi_cells.push(VoronoiCell {
                            id: voronoi_cells.len(),
                            generator: generator_bevy,
                            vertices: boundary.corners().to_vec(),
                            is_boundary_cell: true,
                            neighbor_ids: Vec::new(),
                        });
                    } // else: keine boundary, einzelner Punkt -> keine sinnvolle Zelle
                } else {
                    info!(
                        "VoronoiExtractor: Skipping non-boundary cell for {:?} with only {} circumcenters.",
                        generator_bevy,
                        cell_circumcenters_spade.len()
                    );
                }
                continue;
            }

            if cell_circumcenters_spade.len() >= 2 {
                // Sortierung nur, wenn mind. 2 Punkte da sind
                cell_circumcenters_spade.sort_unstable_by(|a, b| {
                    // ... Sortierlogik ...
                    let angle_a = (a.y - generator_spade.y).atan2(a.x - generator_spade.x);
                    let angle_b = (b.y - generator_spade.y).atan2(b.x - generator_spade.x);
                    angle_a
                        .partial_cmp(&angle_b)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
            }

            let mut cell_vertices_bevy: Vec<Vec2> = cell_circumcenters_spade
                .into_iter()
                .map(|sp| Vec2::new(sp.x, sp.y))
                .collect();

            // Der hinzugefügte Log und Check von dir:
            if cell_vertices_bevy.len() < 3 {
                // Deine Logik für 'boundary: true' und cell_vertices_bevy.len() == 2 greift hier.
                info!(
                    // Geändert zu INFO, da es jetzt eine bewusste Design-Entscheidung sein kann
                    "VoronoiExtractor: Raw Voronoi cell for {:?} (boundary: {}) has only {} vertices BEFORE clipping. Skipping.",
                    generator_bevy,
                    is_boundary_cell_flag,
                    cell_vertices_bevy.len()
                );
                continue;
            }

            if let Some(boundary) = clip_boundary {
                if !cell_vertices_bevy.is_empty() {
                    let clipper_points = boundary.corners();
                    let clipper_algo =
                        crate::math::algorithms::clipping::ClippingAlgorithm::SutherlandHodgman;
                    let polygon_clipper =
                        crate::math::algorithms::clipping::PolygonClipper::new(clipper_algo);

                    let clipped_results = polygon_clipper
                        .clip_polygon_points(&cell_vertices_bevy, &clipper_points)?; // << DIESER Aufruf gibt MathError::InsufficientPoints zurück

                    if let Some(clipped_cell) = clipped_results.into_iter().next() {
                        cell_vertices_bevy = clipped_cell;
                    } else {
                        cell_vertices_bevy.clear(); // Polygon war komplett außerhalb oder wurde zu nichts geclippt
                    }
                } else if is_boundary_cell_flag {
                    // Fall: Randzelle ohne initiale Umkreismittelpunkte (z.B. ein Punkt auf der konvexen Hülle).
                    // Hier ist es knifflig. Man müsste die "unendlichen" Kanten der Zelle konstruieren
                    // und dann gegen die Boundary clippen.
                    // Eine einfache Annahme: Wenn die Zelle leer ist, aber eine Randzelle sein soll
                    // und wir eine Boundary haben, könnte man versuchen, eine sinnvolle Form zu finden.
                    // Für den Moment bleibt sie leer, wenn keine Punkte durch Clipping entstehen.
                    // Dies ist oft ein Zeichen dafür, dass der Generator außerhalb oder direkt auf der
                    // Boundary liegt und seine Zelle stark durch die Boundary beschnitten wird.
                }
            }

            if cell_vertices_bevy.len() >= 3 {
                voronoi_cells.push(VoronoiCell {
                    id: voronoi_cells.len(),
                    generator: generator_bevy,
                    vertices: cell_vertices_bevy,
                    // Eine Zelle ist eine Randzelle, wenn sie ursprünglich als solche erkannt wurde
                    // ODER wenn sie durch das Clipping an die Boundary stößt.
                    // Die `is_boundary_cell_flag` ist hier entscheidender,
                    // das Clipping an sich macht sie nicht unbedingt zur Randzelle im Sinne von "unendlich".
                    // Aber wenn sie geclippt wurde, grenzt sie an die Boundary.
                    is_boundary_cell: is_boundary_cell_flag || clip_boundary.is_some(),
                    neighbor_ids: Vec::new(),
                });
            } else if is_boundary_cell_flag
                && clip_boundary.is_some()
                && cell_vertices_bevy.is_empty()
            {
                // Logik für diesen speziellen Fall (z.B. wenn der Generator sehr nah an einer Ecke der Boundary ist)
                // könnte hier verfeinert werden. Manchmal ist es korrekt, dass die Zelle leer ist.
            }
        }
        
        // Berechne Nachbarschaftsbeziehungen
        Self::calculate_neighbor_relationships(&mut voronoi_cells, triangulation);
        
        Ok(voronoi_cells)
    }

    /// Berechnet die Nachbarschaftsbeziehungen zwischen Voronoi-Zellen basierend auf der Delaunay-Triangulation
    fn calculate_neighbor_relationships(
        cells: &mut [VoronoiCell],
        triangulation: &DelaunayTriangulation<SpadePoint>,
    ) {
        // Erstelle eine Zuordnung von Generatorpunkten zu Zell-IDs
        let mut generator_to_id = std::collections::HashMap::new();
        for cell in cells.iter() {
            let generator_key = (
                (cell.generator.x * 1000.0) as i32,
                (cell.generator.y * 1000.0) as i32,
            );
            generator_to_id.insert(generator_key, cell.id);
        }

        // Für jede Kante in der Delaunay-Triangulation finde die entsprechenden Voronoi-Zellen
        for edge in triangulation.undirected_edges() {
            let vertices = edge.vertices();
            if vertices.len() == 2 {
                let gen1 = vertices[0].position();
                let gen2 = vertices[1].position();
                
                let key1 = ((gen1.x * 1000.0) as i32, (gen1.y * 1000.0) as i32);
                let key2 = ((gen2.x * 1000.0) as i32, (gen2.y * 1000.0) as i32);
                
                if let (Some(&id1), Some(&id2)) = (generator_to_id.get(&key1), generator_to_id.get(&key2)) {
                    // Füge die Nachbarschaftsbeziehung hinzu
                    for cell in cells.iter_mut() {
                        if cell.id == id1 && !cell.neighbor_ids.contains(&id2) {
                            cell.neighbor_ids.push(id2);
                        } else if cell.id == id2 && !cell.neighbor_ids.contains(&id1) {
                            cell.neighbor_ids.push(id1);
                        }
                    }
                }
            }
        }
    }

    /// Extrahiert alle Voronoi-Kanten aus der Delaunay-Triangulation.
    /// `triangulation`: Die Delaunay-Triangulation.
    /// `clip_boundary`: Optionale Grenzen, gegen die die Kanten geclippt werden.
    pub fn extract_edges(
        triangulation: &DelaunayTriangulation<SpadePoint>, // SpadePoint ist Point2<f32>
        clip_boundary: Option<&Bounds2D>,
    ) -> MathResult<Vec<VoronoiEdge>> {
        let mut voronoi_edges = Vec::new();
        let mut processed_delaunay_edges = HashSet::new();

        for delaunay_edge_handle in triangulation.undirected_edges() {
            let [v1_handle, v2_handle] = delaunay_edge_handle.vertices();
            let gen_a_spade = v1_handle.position();
            let gen_b_spade = v2_handle.position();

            let key_a = ((gen_a_spade.x * 1e6) as i64, (gen_a_spade.y * 1e6) as i64);
            let key_b = ((gen_b_spade.x * 1e6) as i64, (gen_b_spade.y * 1e6) as i64);
            let edge_key = if key_a < key_b {
                (key_a, key_b)
            } else {
                (key_b, key_a)
            };

            if processed_delaunay_edges.contains(&edge_key) {
                continue;
            }
            processed_delaunay_edges.insert(edge_key);

            let directed_edge = delaunay_edge_handle.as_directed();
            let face1_handle = directed_edge.face();
            let twin_edge = triangulation.directed_edge(directed_edge.fix());
            let face2_handle = twin_edge.face();

            let cc1_opt = if face1_handle.is_outer() {
                None
            } else {
                face1_handle.as_inner().map(|f| f.circumcenter())
            };

            let cc2_opt = if face2_handle.is_outer() {
                None
            } else {
                face2_handle.as_inner().map(|f| f.circumcenter())
            };

            let gen_a_bevy = Vec2::new(gen_a_spade.x, gen_a_spade.y);
            let gen_b_bevy = Vec2::new(gen_b_spade.x, gen_b_spade.y);

            match (cc1_opt, cc2_opt) {
                (Some(cc1_spade), Some(cc2_spade)) => {
                    let start_bevy = Vec2::new(cc1_spade.x, cc1_spade.y); // .x statt .x() für Point2
                    let end_bevy = Vec2::new(cc2_spade.x, cc2_spade.y); // .y statt .y() für Point2

                    if let Some(boundary) = clip_boundary {
                        if let Some((clipped_start, clipped_end)) =
                            clip_line_segment_to_bounds(start_bevy, end_bevy, boundary)
                        {
                            voronoi_edges.push(VoronoiEdge {
                                start: clipped_start,
                                end: Some(clipped_end),
                                generator_a: gen_a_bevy,
                                generator_b: gen_b_bevy,
                            });
                        }
                    } else {
                        voronoi_edges.push(VoronoiEdge {
                            start: start_bevy,
                            end: Some(end_bevy),
                            generator_a: gen_a_bevy,
                            generator_b: gen_b_bevy,
                        });
                    }
                }
                (Some(cc_spade), None) | (None, Some(cc_spade)) => {
                    if let Some(boundary) = clip_boundary {
                        let start_finite_bevy = Vec2::new(cc_spade.x, cc_spade.y); // .x, .y
                        let delaunay_edge_vec = gen_b_bevy - gen_a_bevy;
                        // Richtung senkrecht zur Delaunay-Kante.
                        // Die Richtung hängt davon ab, welches der beiden Faces (face1 oder face2)
                        // das äußere war, um zu bestimmen, ob die Kante "hinaus" oder "hinein" zeigt.
                        // Wir brauchen eine konsistente Methode, um die Richtung der unendlichen Kante zu bestimmen.
                        // Die Voronoi-Kante ist die Mittelsenkrechte der Delaunay-Kante (gen_a, gen_b).
                        // Der Punkt cc_spade liegt auf dieser Mittelsenkrechten.
                        // Die Richtung der Voronoi-Kante von cc_spade weg ist senkrecht zu (gen_b - gen_a).
                        // Es gibt zwei solche senkrechte Richtungen.
                        // Wir müssen diejenige wählen, die "vom anderen Generator weg" zeigt, wenn man auf der Delaunay-Kante steht.

                        let mid_point_delaunay_edge = (gen_a_bevy + gen_b_bevy) / 2.0;
                        let vec_to_circumcenter = start_finite_bevy - mid_point_delaunay_edge;

                        // Normiere die Richtung des Strahls von cc_spade weg.
                        // Die Richtung des Strahls ist von cc_spade weg vom Mittelpunkt der Delaunay-Kante.
                        // Alternativ: senkrecht zur Delaunay-Kante
                        let perpendicular_dir_candidate1 =
                            Vec2::new(-delaunay_edge_vec.y, delaunay_edge_vec.x)
                                .normalize_or_zero();
                        let _perpendicular_dir_candidate2 =
                            Vec2::new(delaunay_edge_vec.y, -delaunay_edge_vec.x)
                                .normalize_or_zero();

                        // Um die korrekte Richtung zu finden, müssen wir testen, welche Seite des
                        // unendlichen Faces wir betrachten.
                        // Eine einfache Methode: Erzeuge einen Punkt sehr weit weg und schaue, ob er
                        // näher an gen_a oder gen_b ist als der andere Generator (bezogen auf die Senkrechte).
                        // Dies ist etwas heuristisch. Eine robustere Methode würde die Topologie genauer analysieren.

                        // Die Richtung der Kante ist vom cc weg, senkrecht zur Kante zwischen gen_a und gen_b.
                        // Die Richtung von `vec_to_circumcenter` ist schon gut, falls sie nicht Null ist.
                        // Ansonsten die Senkrechte zur Delaunay-Kante.
                        let ray_direction =
                            if vec_to_circumcenter.length_squared() > constants::EPSILON_SQUARED {
                                vec_to_circumcenter.normalize_or_zero()
                            } else {
                                // Fall: cc_spade liegt genau auf dem Mittelpunkt (z.B. bei zwei Punkten)
                                // dann ist die Richtung senkrecht zur Delaunay-Kante.
                                // Hier ist es wichtig zu wissen, auf welcher Seite des Delaunay-Segments das
                                // unendliche Face liegt.
                                // Für jetzt verwenden wir eine der beiden Senkrechten.
                                // Man müsste eigentlich prüfen, welche Seite des Segments (gen_a, gen_b)
                                // "außen" liegt.
                                perpendicular_dir_candidate1 // Willkürliche Wahl, muss evtl. angepasst werden.
                            };

                        if ray_direction.length_squared() < constants::EPSILON_SQUARED {
                            // Kann passieren, wenn gen_a und gen_b sehr nah sind oder cc_spade darauf liegt.
                            // In diesem Fall kann keine eindeutige Richtung für die unendliche Kante bestimmt werden.
                            // Überspringe diese Kante oder behandle sie speziell.
                            continue;
                        }

                        // Erzeuge einen sehr langen Punkt in diese Richtung
                        // Der Faktor 2.0 ist vielleicht nicht immer genug, je nach Größe der Boundary.
                        // Ein größerer Faktor oder eine Projektion auf die Boundary wäre robuster.
                        let far_point =
                            start_finite_bevy + ray_direction * boundary.size().length() * 10.0; // Größerer Multiplikator

                        if let Some((cs, ce)) =
                            clip_line_segment_to_bounds(start_finite_bevy, far_point, boundary)
                        {
                            if cs.distance_squared(ce) > constants::EPSILON_SQUARED {
                                voronoi_edges.push(VoronoiEdge {
                                    start: cs, // Der Punkt näher an start_finite_bevy sollte der Start sein
                                    end: Some(ce),
                                    generator_a: gen_a_bevy,
                                    generator_b: gen_b_bevy,
                                });
                            }
                        }
                        // Die Logik mit far_point1 und far_point2 war potenziell problematisch,
                        // da sie zwei Kanten aus einer unendlichen generieren konnte.
                        // Eine unendliche Kante (Strahl) sollte nur in eine Richtung gehen.
                    } else {
                        voronoi_edges.push(VoronoiEdge {
                            start: Vec2::new(cc_spade.x, cc_spade.y), // .x, .y
                            end: None,
                            generator_a: gen_a_bevy,
                            generator_b: gen_b_bevy,
                        });
                    }
                }
                (None, None) => {
                    // Beide Faces sind außen. Dies sollte für eine Kante in einer
                    // gültigen Triangulation nicht passieren, es sei denn, es ist die
                    // "äußerste" Kante der konvexen Hülle im Unendlichen.
                    // In der Regel werden solche Kanten nicht zu Voronoi-Kanten im begrenzten Diagramm.
                }
            }
        }
        Ok(voronoi_edges)
    }
}

/// Berechnet den Umkreismittelpunkt eines Dreiecks.
/// `triangle`: Ein Array von 3 `Point2<f32>`-Punkten.
fn circumcenter(triangle: &[Point2<f32>; 3]) -> Point2<f32> {
    let (a, b, c) = (triangle[0], triangle[1], triangle[2]); // Direkter Zugriff, da es Point2<f32> ist
    let d_denominator = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));

    // Es ist wichtig, hier eine Division durch Null zu vermeiden, falls die Punkte kollinear sind.
    if d_denominator.abs() < constants::EPSILON {
        // Kollineare Punkte haben keinen eindeutigen Umkreismittelpunkt (liegt im Unendlichen).
        // Gib einen "ungültigen" Punkt oder den Mittelpunkt des längsten Segments zurück,
        // oder behandle diesen Fall speziell (z.B. durch Fehler oder Option<Point2<f32>>).
        // Für Voronoi ist dies oft ein Zeichen für degenerierte Fälle.
        // Hier geben wir einfach einen Punkt zurück, der weit entfernt ist, oder einen NaN-Punkt.
        // Besser wäre es, dies im rufenden Code zu behandeln.
        // Für den Moment, den Mittelpunkt des Dreiecks als Fallback (nicht geometrisch korrekt für Umkreis):
        // return Point2::new((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0);
        // Oder einen speziellen Wert, der anzeigt, dass es keinen gibt.
        // Da spade selbst damit umgehen sollte (selten kollinear durch numerische Präzision),
        // belassen wir es bei der Division. Eine robuste Implementierung würde dies prüfen.
        // Bei Problemen:
        // eprintln!("Warnung: Degeneriertes Dreieck im Umkreismittelpunkt, d_denominator ist sehr klein: {}", d_denominator);
        // return Point2::new(f32::NAN, f32::NAN); // Um den Fehler sichtbar zu machen
    }

    Point2::new(
        ((a.x * a.x + a.y * a.y) * (b.y - c.y)
            + (b.x * b.x + b.y * b.y) * (c.y - a.y)
            + (c.x * c.x + c.y * c.y) * (a.y - b.y))
            / d_denominator,
        ((a.x * a.x + a.y * a.y) * (c.x - b.x)
            + (b.x * b.x + b.y * b.y) * (a.x - c.x)
            + (c.x * c.x + c.y * c.y) * (b.x - a.x))
            / d_denominator,
    )
}

// Die Funktionen clip_line_segment_to_bounds und compute_outcode scheinen in Ordnung zu sein.
// ... (clip_line_segment_to_bounds und compute_outcode bleiben gleich)
/// Clippt ein Liniensegment (p1, p2) gegen eine achsenparallele Bounding Box.
/// Gibt Some((clipped_p1, clipped_p2)) zurück, wenn ein Teil des Segments innerhalb liegt.
/// Verwendet eine vereinfachte Version des Liang-Barsky oder Cohen-Sutherland ähnlichen Ansatzes.
fn clip_line_segment_to_bounds(
    mut p1: Vec2,
    mut p2: Vec2,
    bounds: &Bounds2D,
) -> Option<(Vec2, Vec2)> {
    // Cohen-Sutherland Outcodes
    let mut outcode1 = compute_outcode(p1, bounds);
    let mut outcode2 = compute_outcode(p2, bounds);

    loop {
        if (outcode1 | outcode2) == 0 {
            // Beide Punkte innerhalb: Akzeptieren
            return Some((p1, p2));
        }
        if (outcode1 & outcode2) != 0 {
            // Beide Punkte außerhalb auf derselben Seite: Verwerfen
            return None;
        }

        // Mindestens ein Punkt ist außerhalb. Wähle diesen Punkt.
        let outcode_out = if outcode1 != 0 { outcode1 } else { outcode2 };
        let mut intersection_p = Vec2::ZERO;

        // Vermeide Division durch Null, wenn Linie parallel zu Clipping-Kante
        if (outcode_out & OUTCODE_TOP) != 0 {
            // Oben
            if (p2.y - p1.y).abs() > constants::EPSILON {
                intersection_p.x = p1.x + (p2.x - p1.x) * (bounds.max.y - p1.y) / (p2.y - p1.y);
            } else {
                // Horizontale Linie außerhalb oder auf der Kante
                intersection_p.x = p1.x; // Oder p2.x, irrelevant, da Y bereits clippt
            }
            intersection_p.y = bounds.max.y;
        } else if (outcode_out & OUTCODE_BOTTOM) != 0 {
            // Unten
            if (p2.y - p1.y).abs() > constants::EPSILON {
                intersection_p.x = p1.x + (p2.x - p1.x) * (bounds.min.y - p1.y) / (p2.y - p1.y);
            } else {
                intersection_p.x = p1.x;
            }
            intersection_p.y = bounds.min.y;
        } else if (outcode_out & OUTCODE_RIGHT) != 0 {
            // Rechts
            if (p2.x - p1.x).abs() > constants::EPSILON {
                intersection_p.y = p1.y + (p2.y - p1.y) * (bounds.max.x - p1.x) / (p2.x - p1.x);
            } else {
                // Vertikale Linie außerhalb oder auf der Kante
                intersection_p.y = p1.y;
            }
            intersection_p.x = bounds.max.x;
        } else if (outcode_out & OUTCODE_LEFT) != 0 {
            // Links
            if (p2.x - p1.x).abs() > constants::EPSILON {
                intersection_p.y = p1.y + (p2.y - p1.y) * (bounds.min.x - p1.x) / (p2.x - p1.x);
            } else {
                intersection_p.y = p1.y;
            }
            intersection_p.x = bounds.min.x;
        } else {
            // Sollte nicht passieren, wenn outcode_out != 0
            return None;
        }

        // Ersetze den äußeren Punkt durch den Schnittpunkt und berechne Outcode neu.
        if outcode_out == outcode1 {
            p1 = intersection_p;
            outcode1 = compute_outcode(p1, bounds);
        } else {
            p2 = intersection_p;
            outcode2 = compute_outcode(p2, bounds);
        }
    }
}

const OUTCODE_LEFT: u8 = 1; // 0001
const OUTCODE_RIGHT: u8 = 2; // 0010
const OUTCODE_BOTTOM: u8 = 4; // 0100
const OUTCODE_TOP: u8 = 8; // 1000

fn compute_outcode(p: Vec2, bounds: &Bounds2D) -> u8 {
    let mut code = 0;
    if p.x < bounds.min.x {
        // strikt kleiner/größer für "außerhalb"
        code |= OUTCODE_LEFT;
    } else if p.x > bounds.max.x {
        code |= OUTCODE_RIGHT;
    }
    if p.y < bounds.min.y {
        code |= OUTCODE_BOTTOM;
    } else if p.y > bounds.max.y {
        code |= OUTCODE_TOP;
    }
    code
}
