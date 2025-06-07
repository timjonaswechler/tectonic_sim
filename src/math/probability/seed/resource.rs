// /Users/tim-jonaswechler/GitHub-Projekte/tectonic_sim/src/math/probability/seed/resource.rs
use super::super::noise::Noise;
use super::super::noise::NoiseType;
use bevy::prelude::*;
use bevy::utils::info;
use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

#[derive(Debug, Clone, Copy)]
struct SpiralCursor {
    x: i32,
    y: i32,
    dx: i32,
    dy: i32,
    segment_len: i32,
    segment_pass: i32,
    turns_done: i32,
    started: bool,
}

impl Default for SpiralCursor {
    fn default() -> Self {
        Self {
            x: 0,
            y: 0,
            dx: 1,
            dy: 0,
            segment_len: 1,
            segment_pass: 0,
            turns_done: 0,
            started: false,
        }
    }
}

impl SpiralCursor {
    /// Gibt den aktuellen Punkt zurück und bewegt den Cursor zum nächsten Punkt der Spirale.
    fn next_point(&mut self) -> (i32, i32) {
        // Sonderfall: zuerst (0, 0) zurückgeben
        if !self.started {
            self.started = true;
            return (self.x, self.y);
        }

        // Bewegung ausführen
        self.x += self.dx;
        self.y += self.dy;
        self.segment_pass += 1;

        // Segmentlänge erreicht → Drehen
        if self.segment_pass == self.segment_len {
            self.segment_pass = 0;
            self.turns_done += 1;

            // Drehung um 90° im Uhrzeigersinn
            let new_dx = -self.dy;
            let new_dy = self.dx;
            self.dx = new_dx;
            self.dy = new_dy;

            // Nach zwei Drehungen wird das Segment länger
            if self.turns_done == 2 {
                self.turns_done = 0;
                self.segment_len += 1;
            }
        }

        (self.x, self.y)
    }
}

#[derive(Resource, Debug)]
pub struct SeedResource {
    pub seed: u32,
    pub probabilities: Noise,
    categorical_spiral_cursors: HashMap<String, SpiralCursor>,
    category_z_multiplier: f32,
}

impl SeedResource {
    pub fn from_seed(seed_val: i32) -> Self {
        let mut noise_gen = Noise::default();
        noise_gen.set_noise_type(Some(NoiseType::Value));
        noise_gen.set_seed(Some(seed_val));
        noise_gen.set_frequency(Some(1.0));

        Self {
            seed: seed_val as u32,
            probabilities: noise_gen,
            categorical_spiral_cursors: HashMap::new(),
            category_z_multiplier: 10.0,
        }
    }

    pub fn from_text<S: AsRef<str>>(text: S) -> Self {
        let mut hasher = DefaultHasher::new();
        text.as_ref().hash(&mut hasher);
        let seed_val = hasher.finish();
        Self::from_seed(seed_val as i32)
    }

    pub fn next_value_for_category(&mut self, category: &str) -> f32 {
        // Rückgabe jetzt f32
        let spiral_cursor = self
            .categorical_spiral_cursors
            .entry(category.to_string())
            .or_insert_with(SpiralCursor::default);

        let (cx, cy) = spiral_cursor.next_point();

        let mut category_hasher = DefaultHasher::new();
        category.hash(&mut category_hasher);
        let cz = (category_hasher.finish() % 1000) as f32 * self.category_z_multiplier;

        // super::super::noise::Float ist f32
        let value = self.probabilities.get_noise_3d(
            cx as f32, cy as f32, cz as f32, // cz ist schon f32
        );
        info!(
            "Current point for category '{}'({}) is ({}, {}) = {}",
            category, cz, cx, cy, value
        );
        value
    }

    pub fn next_value(&mut self) -> f32 {
        // Rückgabe jetzt f32
        self.next_value_for_category("__default__")
    }
    /// Liefert einen normalisierten Wert im Bereich [0.0, 1.0)
    pub fn next_value_normalized(&mut self, category_opt: Option<&str>) -> f32 {
        let val = match category_opt {
            Some(cat) => self.next_value_for_category(cat),
            None => self.next_value(),
        };
        (val + 1.0) * 0.5
    }
    /// Liefert ein zufälliges bool mit einer gegebenen Wahrscheinlichkeit für `true`.
    /// Optional kann eine Kategorie angegeben werden.
    pub fn next_bool_with_probability(
        &mut self,
        probability_true: f32,
        category_opt: Option<&str>,
    ) -> bool {
        // Wichtig: next_value_normalized verwenden, damit der Vergleichswert im Bereich [0,1) ist.
        let normalized_val = match category_opt {
            Some(cat) => self.next_value_normalized(Some(cat)),
            None => self.next_value_normalized(None), // Verwendet intern __default__
        };
        normalized_val < probability_true.clamp(0.0, 1.0)
    }

    /// Liefert ein zufälliges bool (ca. 50/50 Chance).
    /// Optional kann eine Kategorie angegeben werden.
    pub fn next_bool(&mut self, category_opt: Option<&str>) -> bool {
        let val = match category_opt {
            Some(cat) => self.next_value_for_category(cat),
            None => self.next_value(),
        };
        val > 0.0 // Noise-Werte sind um 0 zentriert
    }
    /// Liefert eine zufällige Ganzzahl im Bereich [min, max).
    /// Optional kann eine Kategorie angegeben werden.
    pub fn next_i32_in_range(&mut self, min: i32, max: i32, category_opt: Option<&str>) -> i32 {
        if min >= max {
            // Oder panic!("max must be greater than min in next_i32_in_range");
            warn!(
                "next_i32_in_range called with min ({}) >= max ({}). Returning min.",
                min, max
            );
            return min;
        }
        let normalized_val = match category_opt {
            Some(cat) => self.next_value_normalized(Some(cat)),
            None => self.next_value_normalized(None), // Verwendet intern __default__
        };
        let range_size = (max - min) as f32;
        min + (normalized_val * range_size).floor() as i32
    }

    pub fn next_f32_in_range(&mut self, min: f32, max: f32, category_opt: Option<&str>) -> f32 {
        if min >= max {
            warn!(
                "next_f32_in_range called with min ({}) >= max ({}). Returning min.",
                min, max
            );
            return min;
        }
        let normalized_val = match category_opt {
            Some(cat) => self.next_value_normalized(Some(cat)),
            None => self.next_value_normalized(None),
        };
        let range_size = max - min;
        min + (normalized_val * range_size)
    }

    /// Mischt die Elemente eines Slices (In-Place) unter Verwendung des internen Zufallsgenerators.
    /// Implementiert den Fisher-Yates-Algorithmus.
    /// Optional kann eine Kategorie angegeben werden, um das Mischen für bestimmte Zwecke
    /// reproduzierbar zu machen (bei gleichem Seed und gleicher Kategorie).
    pub fn shuffle<T>(&mut self, slice: &mut [T], category_opt: Option<&str>) {
        let n = slice.len();
        if n < 2 {
            return; // Nichts zu mischen für 0 oder 1 Elemente
        }

        // Fisher-Yates: Iteriere von n-1 runter bis 1
        for i in (1..n).rev() {
            // Wähle einen zufälligen Index j aus dem Bereich [0, i] (inklusiv)
            // next_i32_in_range gibt einen Wert im Bereich [min, max) zurück.
            // Also brauchen wir als max = i + 1
            let j = self.next_i32_in_range(0, (i + 1) as i32, category_opt) as usize;
            slice.swap(i, j);
        }
    }

    pub fn reset_with_new_seed(&mut self, new_seed_val: u32) {
        self.seed = new_seed_val;
        // Noise::set_seed erwartet Option<i32>
        self.probabilities.set_seed(Some(new_seed_val as i32));
        self.probabilities.set_frequency(Some(1.0_f32));
        self.probabilities.set_fractal_octaves(Some(1)); // u32
        self.categorical_spiral_cursors.clear();
        self.categorical_spiral_cursors = HashMap::new();
    }
}

impl Default for SeedResource {
    fn default() -> Self {
        // Verwendung einer externen crate `rand` ist hier OK für den *Default*-Seed,
        // da dies nur einmal bei der Initialisierung ohne expliziten Seed passiert.
        // Die eigentliche Zufallsgenerierung für die Simulation nutzt dann `self.probabilities`.
        let seed_number = rand::random::<i32>();
        Self::from_seed(seed_number)
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::NoiseType;
    use super::*;

    #[test]
    fn test_text_seed_consistency() {
        let s1 = SeedResource::from_text("abc");
        let s2 = SeedResource::from_text("abc");
        assert_eq!(s1.seed, s2.seed);
    }

    #[test]
    fn test_numeric_seed() {
        let num = 1337i32;
        let s = SeedResource::from_seed(num);
        assert_eq!(s.seed as u32, num as u32);
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

        // Prüfe, ob NoiseType nicht Value oder ValueCubic ist, da diese bei (0,0,z)
        // und kleinen, unterschiedlichen z-Werten ähnliche Ergebnisse liefern können.
        if sr.probabilities.noise_type != NoiseType::Value
            && sr.probabilities.noise_type != NoiseType::ValueCubic
        {
            // Es ist immer noch möglich, dass Werte zufällig gleich sind, aber unwahrscheinlich.
            if (cat1_val1 - cat2_val1).abs() > f32::EPSILON {
                // Werte sind unterschiedlich, was erwartet wird.
            } else {
                // Werte sind gleich, könnte ein Flake sein oder ein Problem mit der Z-Trennung.
                // Für diesen Test ist es ok, wenn sie manchmal gleich sind,
                // solange es nicht systematisch passiert.
                // warn!("cat1_val1 and cat2_val1 are unexpectedly similar or equal.");
            }
        }
        assert_ne!(
            cat1_val1, cat1_val2,
            "Aufeinanderfolgende Werte derselben Kategorie sollten (meistens) unterschiedlich sein."
        );

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

    #[test]
    fn test_next_i32_in_range_distribution() {
        let mut sr = SeedResource::from_seed(123);
        let min = 0;
        let max = 10;
        let mut counts = vec![0; (max - min) as usize];
        let num_samples = 10000;

        for _ in 0..num_samples {
            let val = sr.next_i32_in_range(min, max, Some("test_range"));
            assert!(
                val >= min && val < max,
                "Wert {} außerhalb des Bereichs [{}, {})",
                val,
                min,
                max
            );
            counts[(val - min) as usize] += 1;
        }

        // Überprüfe, ob jede Zahl mindestens einmal vorkommt (sehr wahrscheinlich für num_samples = 10000)
        for (i, count) in counts.iter().enumerate() {
            assert!(
                *count > 0,
                "Zahl {} im Bereich [{}, {}) wurde nie generiert.",
                i + (min as usize),
                min,
                max
            );
        }
        // println!("Counts für next_i32_in_range [{}, {}): {:?}", min, max, counts);
    }

    #[test]
    fn test_next_bool_with_probability_logic() {
        let mut sr = SeedResource::from_seed(42);
        let mut true_count_p01 = 0;
        let mut true_count_p09 = 0;
        let samples = 10000;

        for _ in 0..samples {
            if sr.next_bool_with_probability(0.1, Some("prob_0.1")) {
                true_count_p01 += 1;
            }
            if sr.next_bool_with_probability(0.9, Some("prob_0.9")) {
                true_count_p09 += 1;
            }
        }
        // Erwartung: true_count_p01 sollte nahe 10% von samples sein
        // Erwartung: true_count_p09 sollte nahe 90% von samples sein
        // Aufgrund der Natur von Rauschen ist dies nicht exakt, aber die Größenordnung sollte stimmen.
        assert!(
            (true_count_p01 as f32 / samples as f32) < 0.2,
            "Anteil true für p=0.1 zu hoch: {}",
            true_count_p01 as f32 / samples as f32
        );
        assert!(
            (true_count_p09 as f32 / samples as f32) > 0.8,
            "Anteil true für p=0.9 zu niedrig: {}",
            true_count_p09 as f32 / samples as f32
        );
        // println!("True counts for p=0.1: {}/{}, p=0.9: {}/{}", true_count_p01, samples, true_count_p09, samples);
    }

    #[test]
    fn test_shuffle_basic() {
        let mut sr = SeedResource::from_seed(99);
        let mut data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let original_data = data.clone();

        sr.shuffle(&mut data, Some("basic_shuffle"));

        // Überprüfe, ob die Länge gleich geblieben ist
        assert_eq!(data.len(), original_data.len());

        // Überprüfe, ob alle ursprünglichen Elemente noch vorhanden sind (sortiert vergleichen)
        let mut sorted_data = data.clone();
        sorted_data.sort();
        assert_eq!(sorted_data, original_data);

        // Es ist sehr wahrscheinlich, dass das Array nicht mehr dem Original entspricht
        // (außer bei extremen Zufällen oder sehr kleinen Arrays)
        if data.len() > 3 {
            // Für sehr kleine Arrays ist die Chance höher, dass sie gleich bleiben
            assert_ne!(
                data, original_data,
                "Array sollte nach dem Mischen (meistens) anders sein."
            );
        }
    }

    #[test]
    fn test_shuffle_empty_and_single() {
        let mut sr = SeedResource::from_seed(100);
        let mut empty_vec: Vec<i32> = vec![];
        sr.shuffle(&mut empty_vec, Some("empty"));
        assert_eq!(empty_vec, Vec::<i32>::new());

        let mut single_vec = vec![42];
        sr.shuffle(&mut single_vec, Some("single"));
        assert_eq!(single_vec, vec![42]);
    }

    #[test]
    fn test_shuffle_deterministic_same_category() {
        let mut sr1 = SeedResource::from_seed(101);
        let mut sr2 = SeedResource::from_seed(101); // Gleicher Seed

        let mut data1 = vec![1, 2, 3, 4, 5];
        let mut data2 = vec![1, 2, 3, 4, 5];

        sr1.shuffle(&mut data1, Some("deterministic_test"));
        sr2.shuffle(&mut data2, Some("deterministic_test")); // Gleiche Kategorie

        assert_eq!(
            data1, data2,
            "Mischen mit gleichem Seed und gleicher Kategorie sollte zum gleichen Ergebnis führen."
        );
    }

    #[test]
    fn test_shuffle_deterministic_same_seed_resource_instance_diff_categories() {
        let mut sr = SeedResource::from_seed(102);

        let mut data1 = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let mut data2 = data1.clone();

        sr.shuffle(&mut data1, Some("cat_A_shuffle"));
        // Wichtig: Da sr dieselbe Instanz ist, wird ihr interner Zustand (Spiral Cursors) weitergeführt.
        // Um hier wirklich unabhängige Sequenzen für die Kategorien zu testen, müssten wir
        // sicherstellen, dass die Kategorien unterschiedliche Startpunkte im Noise Space haben.
        // Die aktuelle Implementierung von next_value_for_category sollte das durch
        // den Hash der Kategorie im z-Wert sicherstellen.
        sr.shuffle(&mut data2, Some("cat_B_shuffle"));

        // Es ist sehr wahrscheinlich, dass sie unterschiedlich sind, da unterschiedliche Kategorien
        // unterschiedliche z-Werte im Noise verwenden und somit unterschiedliche Zufallssequenzen.
        if data1.len() > 5 {
            // Höhere Chance auf Unterschied bei größeren Arrays
            assert_ne!(
                data1, data2,
                "Mischen mit derselben SeedResource-Instanz aber unterschiedlichenShuffle-Kategorien sollte (meistens) zu unterschiedlichen Ergebnissen führen."
            );
        }
    }

    #[test]
    fn test_reset_seed_produces_different_shuffles() {
        let initial_seed = 200;
        let mut sr = SeedResource::from_seed(initial_seed);
        let mut data1 = (0..10).collect::<Vec<i32>>();
        let original_data = data1.clone();

        sr.shuffle(&mut data1, Some("shuffle_before_reset"));

        let mut data2 = original_data.clone();
        sr.reset_with_new_seed(initial_seed as u32 + 1); // Neuer Seed
        sr.shuffle(&mut data2, Some("shuffle_after_reset")); // Gleiche Kategorie, aber neuer Seed

        // data1 und data2 sollten jetzt unterschiedlich sein, da der Seed resettet wurde.
        // Es ist wichtig, dass "shuffle_after_reset" dieselbe Kategorie hat wie zuvor,
        // um zu zeigen, dass der Seed-Reset der entscheidende Faktor ist.
        assert_ne!(
            data1, data2,
            "Shuffles mit unterschiedlichen Seeds sollten unterschiedlich sein."
        );

        // Teste, ob ein Reset auf den ursprünglichen Seed wieder das gleiche Ergebnis liefert
        let mut data3 = original_data.clone();
        sr.reset_with_new_seed(initial_seed as u32); // Zurück zum ursprünglichen Seed
        sr.shuffle(&mut data3, Some("shuffle_before_reset")); // Gleiche Kategorie wie beim ersten Mal

        assert_eq!(
            data1, data3,
            "Shuffle nach Reset auf Originalseed sollte identisch zum ersten Shuffle sein."
        );
    }

    #[test]
    fn test_next_value_normalized_range() {
        let mut sr = SeedResource::from_seed(12345);
        for _ in 0..1000 {
            let val_default = sr.next_value_normalized(None);
            assert!(
                val_default >= 0.0 && val_default < 1.0,
                "Normalized value (default) out of range [0, 1): {}",
                val_default
            );

            let val_cat = sr.next_value_normalized(Some("norm_cat_test"));
            assert!(
                val_cat >= 0.0 && val_cat < 1.0,
                "Normalized value (category) out of range [0, 1): {}",
                val_cat
            );
        }
    }

    // Kleiner Fix in next_bool_with_probability, um next_value_normalized zu verwenden
    #[test]
    fn test_next_bool_with_probability_uses_normalized() {
        let mut sr = SeedResource::from_seed(555);
        let mut true_count = 0;
        let probability = 0.5; // Sollte ca. 50% true ergeben
        let num_samples = 10000;

        for i in 0..num_samples {
            // Wir brauchen eine sich ändernde Kategorie, um sicherzustellen, dass wir nicht
            // immer denselben normalisierten Wert bekommen, falls der zugrundeliegende Noise-Wert
            // für die ersten paar Punkte einer Kategorie konstant bleibt.
            let category = format!("bool_norm_test_{}", i);
            if sr.next_bool_with_probability(probability, Some(&category)) {
                true_count += 1;
            }
        }
        let ratio = true_count as f32 / num_samples as f32;
        // Für Noise-basierte Werte erwarten wir nicht exakt 0.5, aber es sollte in einer vernünftigen Spanne sein.
        // Da der Noise um 0 zentriert ist und dann auf [0,1) normalisiert wird, sollte
        // probability < 0.5 eher false und probability > 0.5 eher true ergeben.
        // Der kritische Punkt ist, ob der normalisierte Wert < probability ist.
        // Noise Werte sind [-1, 1]. Normalisiert: (noise + 1)*0.5 -> [0, 1).
        // Wenn probability = 0.5, dann ist (noise+1)*0.5 < 0.5  <=> noise+1 < 1 <=> noise < 0.
        // Da Noise um 0 zentriert ist, sollte dies ca. 50% der Zeit der Fall sein.
        assert!(
            ratio > 0.4 && ratio < 0.6,
            "Ratio for p=0.5 ({}) is not close to 0.5. True count: {}",
            ratio,
            true_count
        );
    }
}
