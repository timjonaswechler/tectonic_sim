// src/debug/visualization/svg.rs
use crate::math::{
    point_distribution::voronoi::VoronoiCell,
    types::{Bounds2D, SpadePoint as MathSpadePoint},
};
use bevy::log::info;
use bevy::math::Vec2;
use std::borrow::Cow;
use std::io::Write;

// ===================================================================================
// 1. HILFS-STRUCT für die SVG-Erstellung
// ===================================================================================
/// Ein Helfer zum Erstellen einer SVG-Datei.
struct SvgBuilder {
    content: String,
    // Relative Größen, die vom Builder berechnet werden
    stroke_w_normal: f64,
    stroke_w_thin: f64,
    point_radius: f64,
    font_size_coord: f64,
}

impl SvgBuilder {
    /// Erstellt ein neues SVG-Grundgerüst mit Header, Stil und Hintergrund.
    fn new(display_bounds: &Bounds2D, svg_pixel_size: f64) -> Self {
        let viewbox_min_x = display_bounds.min.x as f64;
        let viewbox_min_y = display_bounds.min.y as f64;
        let viewbox_width = display_bounds.width() as f64;
        let viewbox_height = display_bounds.height() as f64;

        let stroke_w_normal = (viewbox_width + viewbox_height) / 2.0 * 0.005;
        let stroke_w_thin = (viewbox_width + viewbox_height) / 2.0 * 0.002;
        let point_radius = (viewbox_width + viewbox_height) / 2.0 * 0.004;
        let font_size_coord = point_radius * 1.1;

        let content = format!(
            r#"<?xml version="1.0" encoding="UTF-8"?>
<svg width="{svg_pixel_size}" height="{svg_pixel_size}" viewBox="{viewbox_min_x} {viewbox_min_y} {viewbox_width} {viewbox_height}" xmlns="http://www.w3.org/2000/svg">
  <style>
    .background {{ fill: #f0f0f0; fill-opacity: 1.0; }}
    .display-bounds {{ fill: none; stroke: #cccccc; stroke-width: {stroke_w_thin}; stroke-dasharray: 5,5; }}
    .clip-bounds {{ fill: none; stroke: #888888; stroke-width: {stroke_w_thin}; stroke-dasharray: 2,2; }}
    .initial-point {{ fill: #ffaaaa; stroke: #cc0000; stroke-width: {stroke_w_thin}; }}
    .relaxed-point {{ fill: #aaccff; stroke: #0000cc; stroke-width: {stroke_w_thin}; }}
    .raw-voronoi-cell {{ fill: rgba(255, 220, 150, 0.3); stroke: #ffaa00; stroke-width: {stroke_w_thin}; }}
    .clipped-voronoi-cell {{ fill: rgba(150, 255, 150, 0.5); stroke: #00aa00; stroke-width: {stroke_w_normal}; }}
    .merged-polygon {{ fill: rgba(200, 150, 255, 0.7); stroke: #5500aa; stroke-width: {stroke_w_normal}; }}
    .merged-polygon-coord {{ 
        font-family: monospace; 
        font-size: {font_size_coord:.3}px; 
        fill: #000000; 
        stroke: white; 
        stroke-width: {stroke_w_thin:.3}; 
        paint-order: stroke fill;
        text-anchor: middle;
        dominant-baseline: middle;
    }}
  </style>
  <rect x="{viewbox_min_x}" y="{viewbox_min_y}" width="{viewbox_width}" height="{viewbox_height}" class="background" />
"#,
        );

        Self {
            content,
            stroke_w_normal,
            stroke_w_thin,
            point_radius,
            font_size_coord,
        }
    }

    /// Zeichnet ein Polygon.
    fn draw_polygon(&mut self, vertices: &[Vec2], class: &str) {
        if vertices.len() < 2 {
            return;
        }
        let points_str: String = vertices
            .iter()
            .map(|p| format!("{:.3},{:.3}", p.x, p.y))
            .collect::<Vec<_>>()
            .join(" ");
        self.content.push_str(&format!(
            r#"  <polygon points="{}" class="{}" />
"#,
            points_str, class
        ));
    }

    /// Zeichnet einen Kreis.
    fn draw_circle(&mut self, center: &Vec2, radius: f64, class: &str) {
        self.content.push_str(&format!(
            r#"  <circle cx="{:.3}" cy="{:.3}" r="{:.3}" class="{}" />
"#,
            center.x, center.y, radius, class
        ));
    }

    /// Zeichnet Text.
    fn draw_text(&mut self, pos: &Vec2, text: &str, class: &str) {
        self.content.push_str(&format!(
            r#"  <text x="{:.3}" y="{:.3}" class="{}">{}</text>
"#,
            pos.x, pos.y, class, text
        ));
    }

    /// Zeichnet ein Rechteck.
    fn draw_rect(&mut self, bounds: &Bounds2D, class: &str) {
        self.content.push_str(&format!(
            r#"  <rect x="{}" y="{}" width="{}" height="{}" class="{}" />
"#,
            bounds.min.x,
            bounds.min.y,
            bounds.width(),
            bounds.height(),
            class
        ));
    }

    /// Speichert die SVG-Datei und schließt die Tags.
    fn save(mut self, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
        self.content.push_str("</svg>");
        let mut file = std::fs::File::create(filename)?;
        file.write_all(self.content.as_bytes())?;
        info!("Debug SVG '{}' wurde erstellt.", filename);
        Ok(())
    }
}

// ===================================================================================
// NEUE HILFSFUNKTION: Bounding Box für Polygone berechnen
// ===================================================================================
/// Berechnet die kleinste Bounding Box, die alle übergebenen Polygone umschließt.
fn calculate_bounds_for_polygons(polygons: &[Vec<Vec2>]) -> Option<Bounds2D> {
    let mut min_p = Vec2::new(f32::MAX, f32::MAX);
    let mut max_p = Vec2::new(f32::MIN, f32::MIN);
    let mut has_points = false;

    for poly in polygons {
        for point in poly {
            min_p.x = min_p.x.min(point.x);
            min_p.y = min_p.y.min(point.y);
            max_p.x = max_p.x.max(point.x);
            max_p.y = max_p.y.max(point.y);
            has_points = true;
        }
    }

    if has_points {
        Some(Bounds2D {
            min: min_p,
            max: max_p,
        })
    } else {
        None
    }
}

// ===================================================================================
// GEÄNDERTE FOKUSSIERTE FUNKTION
// ===================================================================================
/// Erstellt eine SVG-Datei, die die finalen, gemergten Polygone anzeigt.
///
/// # Arguments
/// * `filename` - Der Dateipfad für die zu erstellende SVG.
/// * `display_bounds` - Die Grenzen der SVG-Leinwand ("ViewBox").
/// * `merged_polygons` - Die Polygone, die gezeichnet werden sollen.
/// * `svg_pixel_size` - Die Größe der SVG in Pixeln (Breite und Höhe).
/// * `scale_to_fit` - Wenn `true`, werden die Polygone skaliert und zentriert,
///   um die `display_bounds` maximal auszufüllen, unter Beibehaltung des Seitenverhältnisses.
#[cfg(debug_assertions)]
pub fn create_merged_polygon_svg(
    filename: &str,
    display_bounds: &Bounds2D,
    merged_polygons: &[Vec<Vec2>],
    svg_pixel_size: f64,
    scale_to_fit: bool, // NEUER Parameter
) -> Result<(), Box<dyn std::error::Error>> {
    let mut svg = SvgBuilder::new(display_bounds, svg_pixel_size);
    svg.draw_rect(display_bounds, "display-bounds");

    // Cow (Clone-on-Write) ist hier ideal: Wir vermeiden das Klonen der Polygon-Daten,
    // wenn keine Skalierung nötig ist (Cow::Borrowed). Nur wenn wir skalieren,
    // erzeugen wir neue Daten und speichern sie in Cow::Owned.
    let polygons_to_draw: Cow<[Vec<Vec2>]> = if scale_to_fit {
        // 1. Finde die Bounding Box der aktuellen Polygone
        if let Some(poly_bounds) = calculate_bounds_for_polygons(merged_polygons) {
            // Verhindere Division durch Null, wenn das Polygon keine Fläche hat
            if poly_bounds.width() < 1e-6 || poly_bounds.height() < 1e-6 {
                Cow::Borrowed(merged_polygons) // Keine Skalierung möglich, Original verwenden
            } else {
                // 2. Berechne Skalierungsfaktoren
                let scale_x = display_bounds.width() / poly_bounds.width();
                let scale_y = display_bounds.height() / poly_bounds.height();
                // Nimm den kleineren Faktor, um das Seitenverhältnis zu wahren
                let scale = scale_x.min(scale_y) * 0.95; // 0.95 für einen kleinen Rand

                // 3. Berechne Zentren für die Verschiebung
                let poly_center = poly_bounds.center();
                let display_center = display_bounds.center();

                // 4. Transformiere alle Punkte
                let scaled_polygons: Vec<Vec<Vec2>> = merged_polygons
                    .iter()
                    .map(|poly| {
                        poly.iter()
                            .map(|&p| {
                                // Formel: v_final = (v_original - poly_center) * scale + display_center
                                (p - poly_center) * scale + display_center
                            })
                            .collect()
                    })
                    .collect();

                Cow::Owned(scaled_polygons)
            }
        } else {
            Cow::Borrowed(merged_polygons) // Keine Polygone zum Skalieren
        }
    } else {
        Cow::Borrowed(merged_polygons) // Skalierung nicht gewünscht
    };

    // Zeichne die (ggf. skalierten) Polygone und ihre Koordinaten
    for poly_vertices in polygons_to_draw.iter() {
        svg.draw_polygon(poly_vertices, "merged-polygon");

        for vertex in poly_vertices {
            // Die angezeigten Koordinaten sind nun die transformierten Koordinaten
            let coord_text = format!("({:.1},{:.1})", vertex.x, vertex.y);
            svg.draw_text(vertex, &coord_text, "merged-polygon-coord");
        }
    }

    svg.save(filename)
}

// ===================================================================================
// 3. UMGESCHRIEBENE "ALTE" FUNKTION, DIE JETZT DEN BUILDER NUTZT
// ===================================================================================
/// Erstellt eine umfassende Debug-SVG mit allen Schritten des Voronoi-Prozesses.
#[cfg(debug_assertions)]
pub fn create_voronoi_debug_svg(
    filename: &str,
    display_bounds: &Bounds2D,
    clip_bounds: Option<&Bounds2D>,
    initial_points: Option<&[Vec2]>,
    relaxed_points: &[Vec2],
    voronoi_cells_raw: Option<&[Vec<Vec2>]>,
    voronoi_cells_clipped: Option<&[VoronoiCell]>,
    merged_polygons: Option<&[Vec<Vec2>]>,
    svg_pixel_size: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut svg = SvgBuilder::new(display_bounds, svg_pixel_size);

    // Zeichne display_bounds
    svg.draw_rect(display_bounds, "display-bounds");

    // Zeichne clip_bounds
    if let Some(cb) = clip_bounds {
        svg.draw_rect(cb, "clip-bounds");
    }

    // Zeichne Roh-Voronoi-Zellen
    if let Some(cells) = voronoi_cells_raw {
        for cell_vertices in cells {
            svg.draw_polygon(cell_vertices, "raw-voronoi-cell");
        }
    }

    // Zeichne geclippte Voronoi-Zellen
    if let Some(cells) = voronoi_cells_clipped {
        for cell in cells {
            svg.draw_polygon(&cell.vertices, "clipped-voronoi-cell");
        }
    }

    // Zeichne gemergte Polygone und ihre Koordinaten
    if let Some(polys) = merged_polygons {
        for poly_vertices in polys {
            svg.draw_polygon(poly_vertices, "merged-polygon");
            for vertex in poly_vertices {
                let coord_text = format!("({:.1},{:.1})", vertex.x, vertex.y);
                svg.draw_text(vertex, &coord_text, "merged-polygon-coord");
            }
        }
    }

    // Zeichne initiale Punkte
    if let Some(points) = initial_points {
        for p in points {
            svg.draw_circle(p, svg.point_radius * 0.7, "initial-point");
        }
    }

    // Zeichne relaxierte Punkte
    for p in relaxed_points {
        svg.draw_circle(p, svg.point_radius, "relaxed-point");
    }

    svg.save(filename)
}
