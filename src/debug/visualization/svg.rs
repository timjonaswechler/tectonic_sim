use spade::Point2;
use std::io::Write;
// In deiner main.rs oder einem dedizierten Debug-Modul:

// Diese Funktion wird nur kompiliert, wenn das "debug_voronoi_svg" Feature aktiv ist.
#[cfg(debug_assertions)]
pub fn create_debug_svg_for_normalized_polygon(
    polygon_points: &[Point2<f64>], // Erwartet normalisierte Punkte [-1,1]
    filename: &str,
    svg_pixel_size: f64,
    padding_factor: f64, // z.B. 0.05 für 5% Padding im ViewBox
) -> Result<(), Box<dyn std::error::Error>> {
    if polygon_points.is_empty() {
        println!("Debug SVG: Polygon ist leer, erstelle kein SVG.");
        return Ok(());
    }

    let mut svg_content = String::new();

    // ViewBox so einstellen, dass der Bereich [-1,1] plus Padding abgedeckt wird
    let viewbox_min = -1.0 - padding_factor;
    let viewbox_size = 2.0 + 2.0 * padding_factor;

    svg_content.push_str(&format!(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<svg width="{}" height="{}" viewBox="{} {} {} {}" xmlns="http://www.w3.org/2000/svg">
  <defs>
    <style>
      .normalized-polygon {{ fill: rgba(100, 100, 255, 0.7); stroke: #333399; stroke-width: {}; /* Strichstärke relativ zur ViewBox */ }}
      .background {{ fill: #f8f8f8; }}
    </style>
  </defs>
  <rect x="{}" y="{}" width="{}" height="{}" class="background" />
"#,
        svg_pixel_size, svg_pixel_size, // SVG-Größe in Pixeln
        viewbox_min, viewbox_min, viewbox_size, viewbox_size, // ViewBox-Einstellungen
        0.02, // Beispiel für stroke-width
        viewbox_min, viewbox_min, viewbox_size, viewbox_size // Hintergrund-Rechteck
    ));

    let points_str: String = polygon_points
        .iter()
        .map(|p| format!("{:.6},{:.6}", p.x, p.y)) // Koordinaten sind schon normalisiert
        .collect::<Vec<_>>()
        .join(" ");

    svg_content.push_str(&format!(
        r#"  <polygon points="{}" class="normalized-polygon" />
"#,
        points_str
    ));

    // Optional: Achsen oder Ursprung für Debugging
    svg_content
        .push_str(r#"  <line x1="-1" y1="0" x2="1" y2="0" stroke="grey" stroke-width="0.01" />"#);
    svg_content
        .push_str(r#"  <line x1="0" y1="-1" x2="0" y2="1" stroke="grey" stroke-width="0.01" />"#);
    svg_content.push_str(r#"  <circle cx="0" cy="0" r="0.03" fill="red" />"#);

    svg_content.push_str("</svg>");

    let mut file = std::fs::File::create(filename)?;
    file.write_all(svg_content.as_bytes())?;
    println!("Debug SVG '{}' wurde erstellt.", filename);
    Ok(())
}
