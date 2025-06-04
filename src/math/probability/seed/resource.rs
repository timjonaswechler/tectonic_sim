// /Users/tim-jonaswechler/GitHub-Projekte/tectonic_sim/src/math/probability/seed/resource.rs
use bevy::prelude::*;
use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
// Annahme: Float ist f64, wie in noise.rs definiert
// use crate::math::probability::noise::noise::Float as NoiseFloat; // Pfad ggf. anpassen

// Struktur, um den Zustand des Spiral-Cursors zu speichern
#[derive(Debug, Clone, Copy)]
struct SpiralCursor {
    x: i64,            // Aktuelle diskrete X-Position auf dem Gitter
    y: i64,            // Aktuelle diskrete Y-Position auf dem Gitter
    dx: i64,           // Bewegungsrichtung X (1, 0, -1, 0)
    dy: i64,           // Bewegungsrichtung Y (0, 1, 0, -1)
    segment_len: i64,  // Aktuelle Länge des geraden Segments der Spirale
    segment_pass: i64, // Wie viele Schritte im aktuellen Segment schon gemacht wurden
    turn_counter: i64, // Zählt, wie oft schon abgebogen wurde in der aktuellen "Schicht" der Spirale
}

impl Default for SpiralCursor {
    fn default() -> Self {
        Self {
            x: 0,
            y: 0,
            dx: 1, // Startet mit Bewegung nach rechts
            dy: 0,
            segment_len: 1,
            segment_pass: 0,
            turn_counter: 0,
        }
    }
}

impl SpiralCursor {
    /// Bewegt den Cursor zum nächsten Punkt auf der Spirale.
    fn next_point(&mut self) -> (i64, i64) {
        let current_x = self.x;
        let current_y = self.y;

        self.x += self.dx;
        self.y += self.dy;
        self.segment_pass += 1;

        if self.segment_pass == self.segment_len {
            self.segment_pass = 0;
            // Drehe gegen den Uhrzeigersinn
            let temp_dx = self.dx;
            self.dx = -self.dy;
            self.dy = temp_dx;
            self.turn_counter += 1;
            if self.turn_counter == 2 {
                // Nach zwei Drehungen (z.B. rechts -> hoch -> links)
                self.turn_counter = 0;
                self.segment_len += 1; // Verlängere das Segment
            }
        }
        (current_x, current_y)
    }
}

#[derive(Resource, Debug, Clone)]
pub struct SeedResource {
    pub seed: u64,
    pub probabilities: super::noise::Noise, // Deine Noise-Instanz

    // Für kategorisierte sequentielle Werte mit Spirale und Z-Ebene
    // Key: Kategorie-Name (String), Value: SpiralCursor
    categorical_spiral_cursors: HashMap<String, SpiralCursor>,

    // Standard Z-Offset Multiplikator für Kategorien
    category_z_multiplier: f64,
}

impl SeedResource {
    pub fn from_seed(seed_val: u64) -> Self {
        let mut noise_gen = super::noise::Noise::default();
        noise_gen.set_seed(Some(seed_val as i64));
        // Da du sagst, Frequenz ist 1 und step_size ist 1, setzen wir das hier so.
        // Das bedeutet, wir verwenden die Noise-Koordinaten direkt.
        noise_gen.set_frequency(Some(1.0));
        noise_gen.set_fractal_octaves(Some(1)); // Oder was auch immer hier passend ist

        Self {
            seed: seed_val,
            probabilities: noise_gen,
            categorical_spiral_cursors: HashMap::new(),
            category_z_multiplier: 10.0, // Sorgt für genügend Abstand zwischen Z-Ebenen der Kategorien
        }
    }

    pub fn from_text<S: AsRef<str>>(text: S) -> Self {
        let mut hasher = DefaultHasher::new();
        text.as_ref().hash(&mut hasher);
        let seed_val = hasher.finish();
        Self::from_seed(seed_val)
    }

    /// Liefert den "nächsten" Noise-Wert für eine bestimmte Kategorie
    /// unter Verwendung einer eigenen Z-Ebene und eines spiralförmigen Cursors in X/Y.
    pub fn next_value_for_category(&mut self, category: &str) -> f64 {
        let spiral_cursor = self
            .categorical_spiral_cursors
            .entry(category.to_string())
            .or_insert_with(SpiralCursor::default);

        let (cx, cy) = spiral_cursor.next_point();

        // Z-Koordinate basierend auf dem Hash der Kategorie bestimmen
        let mut category_hasher = DefaultHasher::new();
        category.hash(&mut category_hasher);
        // Der Hash wird verwendet, um eine "einzigartige" Z-Ebene zu finden.
        // Modulo und Multiplikator, um die Z-Werte etwas zu streuen und Kollisionen unwahrscheinlicher zu machen.
        // Die `probabilities` Instanz wird hier ihren eigenen `seed` (vom Hauptseed abgeleitet) verwenden.
        // Die Z-Koordinate selbst ist relativ zum internen Seed der Noise-Instanz.
        let cz = (category_hasher.finish() % 1000) as f64 * self.category_z_multiplier;

        // Verwende get_noise_3d, da jede Kategorie ihre eigene Z-Ebene hat.
        // Die Typkonvertierung von i64 (Cursor) zu f64 (Float in Noise) ist wichtig.
        let value = self.probabilities.get_noise_3d(
            cx as super::noise::noise::Float, // Float ist f64 in noise.rs
            cy as super::noise::noise::Float,
            cz as super::noise::noise::Float,
        );

        value
    }

    /// Liefert den "nächsten" unspezifischen Noise-Wert.
    /// Nutzt die Kategorie "default" implizit.
    pub fn next_value(&mut self) -> f64 {
        self.next_value_for_category("__default__") // Eine Standard-Kategorie
    }

    /// Ändert den globalen Seed und setzt alle Cursor zurück.
    pub fn reset_with_new_seed(&mut self, new_seed_val: u64) {
        self.seed = new_seed_val;
        self.probabilities.set_seed(Some(new_seed_val as i64));
        // Ggf. Frequenz/Oktaven etc. hier anpassen, wenn sie vom Seed abhängen sollen.
        // Da Frequenz=1 fest ist, brauchen wir das hier nicht.

        // Cursor zurücksetzen
        self.categorical_spiral_cursors.clear();
        // info!("SeedResource resettet mit neuem Seed: {}", new_seed_val);
    }
}

impl Default for SeedResource {
    fn default() -> Self {
        let seed_number = rand::random::<u64>();
        Self::from_seed(seed_number)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_text_seed_consistency() {
        let s1 = SeedResource::from_text("abc");
        let s2 = SeedResource::from_text("abc");
        assert_eq!(s1.seed, s2.seed);
    }

    #[test]
    fn test_numeric_seed() {
        let num = 1337u64;
        let s = SeedResource::from_seed(num);
        assert_eq!(s.seed, num);
    }

    #[test]
    fn test_spiral_cursor_movement() {
        let mut cursor = SpiralCursor::default();
        let points = [
            (0, 0),
            (1, 0),
            (1, 1),
            (0, 1),
            (-1, 1),
            (-1, 0),
            (-1, -1),
            (0, -1),
            (1, -1),
            (2, -1),
            (2, 0),
            (2, 1),
            (2, 2),
            (1, 2),
            (0, 2),
            (-1, 2),
            (-2, 2),
            (-2, 1),
            (-2, 0),
            (-2, -1),
            (-2, -2),
            (-1, -2),
            (0, -2),
            (1, -2),
            (2, -2),
        ];
        for &expected_point in points.iter() {
            assert_eq!(cursor.next_point(), expected_point);
        }
    }

    #[test]
    fn test_categorized_spiral_values() {
        let mut sr = SeedResource::default();
        let cat1_val1 = sr.next_value_for_category("terrain");
        let cat2_val1 = sr.next_value_for_category("weather");
        let cat1_val2 = sr.next_value_for_category("terrain");

        // Es ist sehr wahrscheinlich, aber nicht 100% garantiert, dass sie verschieden sind,
        // besonders bei Noise-Funktionen. Aber für unterschiedliche Cursor/Kategorien ist es das Ziel.
        if sr.probabilities.noise_type != super::super::noise::types::NoiseType::Value
            && sr.probabilities.noise_type != super::super::noise::types::NoiseType::ValueCubic
        {
            // Value Noise kann bei (0,0,z) für verschiedene z anfangs sehr ähnlich sein, wenn z klein ist
            // oder der Z-Offset nicht groß genug ist.
            assert_ne!(
                cat1_val1, cat2_val1,
                "Werte aus verschiedenen Kategorien sollten (meistens) unterschiedlich sein."
            );
        }
        assert_ne!(
            cat1_val1, cat1_val2,
            "Aufeinanderfolgende Werte derselben Kategorie sollten (meistens) unterschiedlich sein."
        );

        // Generiere ein paar Werte, um zu sehen, ob der Cursor sich bewegt
        let mut previous_val = cat1_val2;
        for _ in 0..5 {
            let next_val = sr.next_value_for_category("terrain");
            assert_ne!(
                previous_val, next_val,
                "Weitere aufeinanderfolgende Werte sollten (meistens) unterschiedlich sein."
            );
            previous_val = next_val;
        }
    }

    #[test]
    fn test_default_next_value_uses_spiral() {
        let mut sr = SeedResource::default();
        let val1 = sr.next_value();
        let val2 = sr.next_value();
        assert_ne!(
            val1, val2,
            "Aufeinanderfolgende next_value() mit Spiralcursor sollten (meistens) unterschiedlich sein."
        );
    }
}
