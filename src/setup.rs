// ./src/setup.rs
use crate::debug::visualization::normal_vector::NormalArrowVisual;
use bevy::prelude::*;
use bevy_panorbit_camera::PanOrbitCamera;

const SPHERE_RADIUS: f32 = 1.0;

pub fn setup_scene(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // Kugel als Platzhalter für den Planeten
    commands.spawn(PbrBundle {
        mesh: meshes.add(Sphere::new(SPHERE_RADIUS).mesh().build()), // Standard Kugelradius 1.0
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
    let north_pole_origin = Vec3::new(0.0, SPHERE_RADIUS, 0.0);
    commands.spawn(NormalArrowVisual {
        origin: north_pole_origin,
        direction: north_pole_origin.normalize(), // Normale ist der normalisierte Vektor vom Kugelzentrum zum Punkt
        length: 0.3,
        color: Color::RED,
    });

    // Pfeil am Äquator (X-Achse)
    let equator_x_origin = Vec3::new(SPHERE_RADIUS, 0.0, 0.0);
    commands.spawn(NormalArrowVisual {
        origin: equator_x_origin,
        direction: equator_x_origin.normalize(),
        length: 0.3,
        color: Color::GREEN,
    });

    // Pfeil am Äquator (Z-Achse)
    let equator_z_origin = Vec3::new(0.0, 0.0, SPHERE_RADIUS);
    commands.spawn(NormalArrowVisual {
        origin: equator_z_origin,
        direction: equator_z_origin.normalize(),
        length: 0.3,
        color: Color::BLUE,
    });

    // Ein schräger Pfeil
    let angled_origin = Vec3::new(1.0, 1.0, 0.0).normalize() * SPHERE_RADIUS;
    commands.spawn(NormalArrowVisual {
        origin: angled_origin,
        direction: angled_origin.normalize(),
        length: 0.4,
        color: Color::ORANGE,
    });
}
