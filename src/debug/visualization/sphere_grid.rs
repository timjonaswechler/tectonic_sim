// src/debug/visualization/sphere_grid.rs
use crate::physics::sim; // Stelle sicher, dass dieser Pfad zu deinen sim::resources passt
use bevy::prelude::*;
use std::f32::consts::PI;

pub fn draw_sphere_grid_gizmos(
    mut gizmos: Gizmos,
    sim_params: Res<sim::resources::SimulationParameters>, // Zugriff auf den Planetenradius
) {
    let radius = sim_params.planet_radius_km; // Verwende den Radius aus den Simulationsparametern

    if radius <= 0.0 {
        // Sicherheitscheck, falls Radius noch nicht initialisiert oder ungültig
        return;
    }

    let num_latitude_lines = 17; // Anzahl der Breitenkreise (ohne Pole)
    let num_longitude_lines = 35; // Anzahl der Längengrade
    let segments_per_longitude_line = 32; // Auflösung der Längengrad-Linien

    let lat_color = Color::rgba(0.0, 1.0, 0.0, 0.3); // Grün mit etwas Transparenz
    let lon_color = Color::rgba(0.0, 0.5, 1.0, 0.3); // Blau mit etwas Transparenz
    let axis_color_x = Color::RED;
    let axis_color_y = Color::LIME_GREEN;
    let axis_color_z = Color::CYAN;

    // --- Breitenkreise (Latitude) ---
    // Diese sind Kreise parallel zur XZ-Ebene.
    for i in 1..=num_latitude_lines {
        let fraction = i as f32 / (num_latitude_lines + 1) as f32;
        let lat_angle_rad = fraction * PI - PI / 2.0; // von -PI/2 bis +PI/2 (ohne Pole)

        let y_coord = radius * lat_angle_rad.sin();
        let circle_radius_in_xz_plane = radius * lat_angle_rad.cos();

        if circle_radius_in_xz_plane > 0.001 * radius {
            // Skaliere die Toleranz mit dem Radius
            gizmos.circle(
                Vec3::new(0.0, y_coord, 0.0),
                Direction3d::Y, // Normale zeigt entlang der Y-Achse
                circle_radius_in_xz_plane,
                lat_color,
            );
        }
    }

    // --- Längengrade (Longitude) ---
    // Diese sind Halbkreise von Pol zu Pol.
    for i in 0..num_longitude_lines {
        let lon_angle_rad = (i as f32 / num_longitude_lines as f32) * 2.0 * PI;

        let mut points = Vec::with_capacity(segments_per_longitude_line + 1);
        for j in 0..=segments_per_longitude_line {
            let segment_fraction = j as f32 / segments_per_longitude_line as f32;
            let polar_angle_rad = segment_fraction * PI - PI / 2.0; // von -PI/2 (Südpol) bis +PI/2 (Nordpol)

            // Kugelkoordinaten zu Kartesisch (Y ist Höhe):
            // x = r * cos(polar) * cos(longitude_around_y)
            // y = r * sin(polar)
            // z = r * cos(polar) * sin(longitude_around_y)
            let x = radius * polar_angle_rad.cos() * lon_angle_rad.cos();
            let y = radius * polar_angle_rad.sin();
            let z = radius * polar_angle_rad.cos() * lon_angle_rad.sin();
            points.push(Vec3::new(x, y, z));
        }
        gizmos.linestrip(points, lon_color);
    }

    // --- Optional: Koordinatenachsen zeichnen (skaliert mit Radius) ---
    let axis_length = radius;
    gizmos.line(Vec3::ZERO, Vec3::X * axis_length, axis_color_x);
    gizmos.line(Vec3::ZERO, Vec3::Y * axis_length, axis_color_y);
    gizmos.line(Vec3::ZERO, Vec3::Z * axis_length, axis_color_z);
}
