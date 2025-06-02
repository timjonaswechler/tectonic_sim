use bevy::prelude::*;

/// Komponente, um einen Normalenpfeil auf der Kugeloberfläche zu definieren.
#[derive(Component, Debug)]
pub struct NormalArrowVisual {
    /// Der Ursprung des Pfeils auf der Kugeloberfläche (Weltkoordinaten).
    pub origin: Vec3,
    /// Die Richtung des Pfeils (sollte ein normalisierter Vektor sein).
    pub direction: Vec3,
    /// Die Länge des Pfeils.
    pub length: f32,
    /// Die Farbe des Pfeils.
    pub color: Color,
}

impl Default for NormalArrowVisual {
    fn default() -> Self {
        Self {
            origin: Vec3::ZERO,
            direction: Vec3::Y, // Standardmäßig nach oben
            length: 0.5,
            color: Color::YELLOW,
        }
    }
}

impl NormalArrowVisual {
    pub fn new(origin: Vec3, length: f32, color: Color) -> Self {
        Self {
            origin,
            direction: origin.normalize(), // Stelle sicher, dass die Richtung normalisiert ist
            length,
            color,
        }
    }
}

/// System, das alle `NormalArrowVisual`-Entitäten nimmt und sie als Gizmos zeichnet.
pub fn draw_normal_arrows_system(query: Query<&NormalArrowVisual>, mut gizmos: Gizmos) {
    for arrow_data in query.iter() {
        let tip = arrow_data.origin + arrow_data.direction * arrow_data.length;
        // `gizmos.arrow` zeichnet einen Pfeil vom ersten zum zweiten Punkt.
        // Beachte: `arrow` hat standardmäßig eine feste Kopfgröße.
        // Für mehr Kontrolle über den Kopf müsstest du ihn manuell aus Linien zeichnen.
        gizmos.arrow(arrow_data.origin, tip, arrow_data.color);

        // Alternative: Nur eine Linie, wenn der Pfeilkopf von gizmos.arrow stört
        // gizmos.line(arrow_data.origin, tip, arrow_data.color);
    }
}
