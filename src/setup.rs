// ./src/setup.rs
use crate::{debug::visualization::normal_vector::NormalArrowVisual, physics::sim};
use bevy::prelude::*;
use bevy_panorbit_camera::PanOrbitCamera;

pub fn setup_scene(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut sim_params: ResMut<sim::resources::SimulationParameters>,
) {
    // Kugel als Platzhalter für den Planeten
    commands.spawn(PbrBundle {
        mesh: meshes.add(Sphere::new(sim_params.planet_radius_km).mesh().build()),
        material: materials.add(StandardMaterial {
            base_color: Color::rgb(0.3, 0.5, 0.3),
            metallic: 0.2,
            perceptual_roughness: 0.7,
            ..default()
        }),
        transform: Transform::from_xyz(0.0, 0.0, 0.0),
        ..default()
    });

    // Licht
    commands.spawn(PointLightBundle {
        point_light: PointLight {
            shadows_enabled: true,
            intensity: 10_000_000.,
            range: 100.0,
            ..default()
        },
        transform: Transform::from_xyz(4.0, 8.0, 4.0),
        ..default()
    });

    // Kamera
    commands.spawn((
        Camera3dBundle {
            transform: Transform::from_xyz(-2.0, 2.5, 5.0).looking_at(Vec3::ZERO, Vec3::Y),
            ..default()
        },
        PanOrbitCamera {
            button_orbit: MouseButton::Right,
            button_pan: MouseButton::Middle,
            radius: Some(5.0), // Start Entfernung
            ..default()
        },
    ));
    // --- Test-Pfeile spawnen ---
    // Pfeil am Nordpol
    let north_pole_origin = Vec3::new(0.0, sim_params.planet_radius_km, 0.0);
    commands.spawn(NormalArrowVisual {
        origin: north_pole_origin,
        direction: north_pole_origin.normalize(), // Normale ist der normalisierte Vektor vom Kugelzentrum zum Punkt
        length: 0.3,
        color: Color::RED,
    });

    // Pfeil am Äquator (X-Achse)
    let equator_x_origin = Vec3::new(sim_params.planet_radius_km, 0.0, 0.0);
    commands.spawn(NormalArrowVisual {
        origin: equator_x_origin,
        direction: equator_x_origin.normalize(),
        length: 0.3,
        color: Color::GREEN,
    });

    // Pfeil am Äquator (Z-Achse)
    let equator_z_origin = Vec3::new(0.0, 0.0, sim_params.planet_radius_km);
    commands.spawn(NormalArrowVisual {
        origin: equator_z_origin,
        direction: equator_z_origin.normalize(),
        length: 0.3,
        color: Color::BLUE,
    });
}
