use bevy::prelude::*;

/// Komponente, um ein Polygon auf der Kugeloberfläche zu definieren.
#[derive(Component, Debug)]
pub struct PolygonVisual {
    /// Die Punkte des Polygons (sollten normalisierte Vektoren sein).
    pub points: Vec<Vec3>,
    /// Die Farbe des Polygons.
    pub color: Color,
}

impl Default for PolygonVisual {
    fn default() -> Self {
        Self {
            points: Vec::new(),
            color: Color::YELLOW,
        }
    }
}

impl PolygonVisual {
    pub fn new(points: Vec<Vec3>, color: Color) -> Self {
        Self { points, color }
    }
}

/// System, das alle `PolygonVisual`-Entitäten nimmt und sie als Gizmos zeichnet.
pub fn draw_polygons_system(query: Query<&PolygonVisual>, mut gizmos: Gizmos) {
    for polygon in query.iter() {
        // Zeichne das Polygon Linie für Linie
        if polygon.points.len() < 3 {
            warn!("PolygonVisual requires at least 3 points to form a polygon.");
            continue;
        }
        // Schließe das Polygon, indem der erste Punkt am Ende hinzugefügt wird
        let mut closed_points = polygon.points.clone();
        closed_points.push(closed_points[0]);
        for i in 0..polygon.points.len() {
            // Zeichne eine Linie zwischen jedem Punkt und dem nächsten
            let start = closed_points[i];
            let end = closed_points[i + 1];
            gizmos.line(start, end, polygon.color);
        }
    }
}
