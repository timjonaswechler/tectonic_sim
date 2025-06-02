// In src/math/polygon/voronoi.rs (oder ein neues Modul wie src/math/polygon/projection.rs)
use crate::math::sphere::algorithms::{
    create_orthogonal_tangent_vectors, rotate_vector_around_axis,
};
use bevy::prelude::{Vec2, Vec3}; // Verwende Bevy's Vec2 für Konsistenz, wenn möglich
use spade::Point2;

/// Projiziert ein normalisiertes 2D-Polygon (Koordinaten im Bereich [-1, 1])
/// auf die Oberfläche einer Kugel.
///
/// # Arguments
/// * `normalized_polygon_2d` - Die Eckpunkte des 2D-Polygons (x,y in [-1,1]).
/// * `sphere_center_point` - Der 3D-Punkt auf der Kugeloberfläche, der das Zentrum der Projektion sein soll (normalisiert).
/// * `sphere_radius` - Der Radius der Kugel.
/// * `max_angular_extent_degrees` - Die maximale Winkelausdehnung des Polygons auf der Kugel (Durchmesser).
///                                  Ein Punkt des 2D-Polygons bei (1,0) wird so skaliert, dass er
///                                  `max_angular_extent_degrees / 2.0` vom `sphere_center_point` entfernt ist.
///
/// # Returns
/// Eine Liste von 3D-Punkten, die das Polygon auf der Kugeloberfläche darstellen.
pub fn project_normalized_2d_polygon_to_sphere(
    normalized_polygon_2d: &[Point2<f64>], // Input von Voronoi
    sphere_center_point_normalized: Vec3,  // Zentrum auf Einheitskugel
    sphere_radius: f32,
    max_angular_extent_degrees: f32,
) -> Vec<Vec3> {
    if normalized_polygon_2d.is_empty() {
        return Vec::new();
    }

    // 1. Skalierungsfaktor bestimmen
    // Der weiteste Punkt im normalisierten 2D-Polygon ist bei sqrt(1^2+1^2) = sqrt(2) vom Ursprung entfernt.
    // Aber wir wollen, dass ein Punkt bei (1,0) oder (0,1) einem bestimmten Winkel entspricht.
    // Die maximale Ausdehnung ist der "Durchmesser" des Polygons auf der Kugel.
    // Der Radius der Ausdehnung ist also `max_angular_extent_degrees / 2.0`.
    let max_angular_radius_rad = (max_angular_extent_degrees / 2.0).to_radians();

    // Ein 2D-Input-Punkt (1.0, 0.0) im normalisierten Polygon soll auf der Kugel
    // einen Winkelabstand von `max_angular_radius_rad` zum `sphere_center_point` haben.
    // Die Bogenlänge auf der Einheitskugel ist gleich dem Winkel in Radiant.
    // Also entspricht die Koordinate 1.0 im 2D-Polygon einem Winkel von `max_angular_radius_rad`.
    // Der Skalierungsfaktor `s` ist also `max_angular_radius_rad`.
    // Ein 2D-Punkt (x_2d, y_2d) wird zu einem Winkel-Offset (x_2d * s, y_2d * s)
    let angle_scale_factor = max_angular_radius_rad;

    // 2. Lokales Koordinatensystem am Kugelmittelpunkt
    // `sphere_center_point_normalized` ist die Normale der Tangentialebene.
    let (tangent_u, tangent_v) = create_orthogonal_tangent_vectors(sphere_center_point_normalized);
    // tangent_u könnte als "lokale X-Achse / Ost" und tangent_v als "lokale Y-Achse / Nord" in der Ebene betrachtet werden.

    let mut projected_vertices_3d = Vec::with_capacity(normalized_polygon_2d.len());

    for p2d_f64 in normalized_polygon_2d {
        let p2d = Vec2::new(p2d_f64.x as f32, p2d_f64.y as f32); // Konvertiere zu Vec2<f32>

        // 3. Umwandlung der 2D-Koordinaten in kleine Winkel-Offsets
        // x_angle_offset entspricht einer Rotation um tangent_v (lokale Y-Achse der Tangentialebene)
        // y_angle_offset entspricht einer Rotation um tangent_u (lokale X-Achse der Tangentialebene)
        let x_angle_offset_rad = p2d.x * angle_scale_factor;
        let y_angle_offset_rad = p2d.y * angle_scale_factor;

        // 4. Punkt auf Kugel durch zwei Rotationen erzeugen:
        // Starte am Kugelmittelpunkt.
        // Rotiere zuerst um `tangent_v` (lokale Y-Achse der Ebene) um `x_angle_offset_rad`.
        // Dies bewegt den Punkt entlang der `tangent_u`-Richtung (lokale X-Achse der Ebene).
        let point_after_first_rotation = rotate_vector_around_axis(
            sphere_center_point_normalized, // Der zu rotierende Punkt ist der Zentrumspunkt selbst
            tangent_v,                      // Rotationsachse
            x_angle_offset_rad,             // Winkel
        );

        // Rotiere nun diesen neuen Punkt um die *ursprüngliche* `tangent_u` Achse.
        // ACHTUNG: Die Rotationsachse muss relativ zum aktuellen System bleiben.
        // Besser: Konstruiere den Punkt direkt über sphärische Addition oder eine andere Projektionsmethode.

        // Alternative: Tangentialebenen-Approximation für kleine Winkel (weniger genau, aber einfacher)
        // P_on_tangent_plane = sphere_center_point_normalized * sphere_radius
        //                      + tangent_u * (p2d.x * max_angular_radius_rad * sphere_radius) -- falsch, Winkel zu Distanz
        //                      + tangent_v * (p2d.y * max_angular_radius_rad * sphere_radius);
        // P_on_sphere = P_on_tangent_plane.normalize() * sphere_radius;

        // Robusterer Ansatz mit Rotationen oder slerp:
        // Wir möchten einen Punkt, der vom `sphere_center_point_normalized` aus
        // um `y_angle_offset_rad` entlang `tangent_v` ("lokal Nord") und
        // um `x_angle_offset_rad` entlang `tangent_u` ("lokal Ost") verschoben ist.

        // Methode: Punkt auf Tangentialebene konstruieren und dann auf Kugel normalisieren.
        // Die `angle_scale_factor` gibt den Winkel an, nicht die Distanz auf der Tangentialebene.
        // Für kleine Winkel: Distanz auf Tangentialebene ≈ Winkel * Radius.
        // Hier arbeiten wir auf der Einheitskugel (Radius 1) für die Offsets.
        let dist_on_tangent_plane_u = x_angle_offset_rad; // tan(x_angle) * 1, für kleine Winkel approx x_angle
        let dist_on_tangent_plane_v = y_angle_offset_rad; // tan(y_angle) * 1, für kleine Winkel approx y_angle

        // Punkt in der Tangentialebene relativ zum Kugelmittelpunkt (angenommen Kugelzentrum ist 0,0,0)
        let point_on_tangent_offset =
            tangent_u * dist_on_tangent_plane_u + tangent_v * dist_on_tangent_plane_v;

        // Der Punkt, der auf die Kugel projiziert werden soll:
        let point_to_project = sphere_center_point_normalized + point_on_tangent_offset;

        // Auf Kugeloberfläche normalisieren und skalieren
        let final_3d_point_on_sphere = point_to_project.normalize() * sphere_radius;

        projected_vertices_3d.push(final_3d_point_on_sphere);
    }

    projected_vertices_3d
}
