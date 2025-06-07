//! Defines the `SeedResource` and related structures for seed-based random number generation.
//!
//! The `SeedResource` uses a noise function and a spiral pattern to generate
//! deterministic sequences of pseudo-random numbers based on an initial seed.
//! It supports categorized random number generation, allowing for different, independent
//! sequences of random numbers for various purposes within an application.

use super::super::noise::Noise;
use super::super::noise::NoiseType;
use bevy::prelude::*;
use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// Represents a cursor that moves in a spiral pattern on a 2D grid.
///
/// The spiral starts at (0,0) and expands outwards. This is used by `SeedResource`
/// to pick unique 2D coordinates for generating noise values for different categories
/// or subsequent calls.
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
    /// Creates a new `SpiralCursor` initialized at the starting state (0,0), ready to move.
    fn default() -> Self {
        Self {
            x: 0,
            y: 0,
            dx: 1, // Start by moving right
            dy: 0,
            segment_len: 1,
            segment_pass: 0,
            turns_done: 0,
            started: false,
        }
    }
}

impl SpiralCursor {
    /// Returns the current point of the cursor and then advances it to the next point in the spiral.
    ///
    /// The sequence of points generated is:
    /// (0,0), (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1), (1,-1), (2,-1), ...
    fn next_point(&mut self) -> (i32, i32) {
        // Special case: return (0,0) first
        if !self.started {
            self.started = true;
            return (self.x, self.y);
        }

        // Execute movement
        self.x += self.dx;
        self.y += self.dy;
        self.segment_pass += 1;

        // Reached end of current segment -> turn
        if self.segment_pass == self.segment_len {
            self.segment_pass = 0;
            self.turns_done += 1;

            // Rotate 90 degrees clockwise
            let new_dx = -self.dy;
            let new_dy = self.dx;
            self.dx = new_dx;
            self.dy = new_dy;

            // After two turns, the segment length increases
            if self.turns_done == 2 {
                self.turns_done = 0;
                self.segment_len += 1;
            }
        }

        (self.x, self.y)
    }
}

/// A Bevy `Resource` that manages seed-based random number generation.
///
/// It uses a noise function (`Noise`) and `SpiralCursor`s to generate deterministic
/// pseudo-random values. The generation can be categorized, allowing for independent
/// sequences of random numbers for different uses, all derived from the same initial seed.
///
/// The resource can be initialized with a specific `u32` seed or a string, which will be hashed
/// to a seed. It provides methods to get raw noise values, normalized values, booleans,
/// integers within a range, and floats within a range. It also includes a Fisher-Yates shuffle method.
#[derive(Resource, Debug)]
pub struct SeedResource {
    /// The current master seed value.
    pub seed: u32,
    /// The noise generator instance used to produce random values.
    pub probabilities: Noise,
    /// Stores a `SpiralCursor` for each category of random number generation.
    /// This allows different categories to have their own independent sequence of noise coordinates.
    categorical_spiral_cursors: HashMap<String, SpiralCursor>,
    /// A multiplier applied to the hashed category value to determine the z-coordinate for 3D noise.
    /// This helps in separating noise values for different categories.
    category_z_multiplier: f32,
}

impl SeedResource {
    /// Creates a new `SeedResource` initialized with a specific integer seed value.
    ///
    /// # Arguments
    ///
    /// * `seed_val`: The `i32` seed value. It will be cast to `u32` for storage.
    pub fn from_seed(seed_val: i32) -> Self {
        let mut noise_gen = Noise::default();
        noise_gen.set_noise_type(Some(NoiseType::Value)); // Consider making this configurable
        noise_gen.set_seed(Some(seed_val));
        noise_gen.set_frequency(Some(1.0)); // Consider making this configurable

        Self {
            seed: seed_val as u32,
            probabilities: noise_gen,
            categorical_spiral_cursors: HashMap::new(),
            category_z_multiplier: 10.0, // Factor to separate categories in z-dimension for noise
        }
    }

    /// Creates a new `SeedResource` initialized with a seed derived from a text string.
    ///
    /// The text string is hashed to generate a `u64` value, which is then cast to `i32`
    /// to create the seed.
    ///
    /// # Arguments
    ///
    /// * `text`: A string slice or any type that can be referenced as a string slice (`AsRef<str>`).
    pub fn from_text<S: AsRef<str>>(text: S) -> Self {
        let mut hasher = DefaultHasher::new();
        text.as_ref().hash(&mut hasher);
        let seed_val = hasher.finish();
        Self::from_seed(seed_val as i32) // Note: u64 to i32 cast might lose information
    }

    /// Generates the next pseudo-random `f32` value for a specific category.
    ///
    /// Values are typically in the range `[-1.0, 1.0]`, depending on the noise function.
    /// Each category uses its own `SpiralCursor` to ensure a unique sequence of values.
    ///
    /// # Arguments
    ///
    /// * `category`: A string slice representing the category for which to generate the value.
    pub fn next_value_for_category(&mut self, category: &str) -> f32 {
        let spiral_cursor = self
            .categorical_spiral_cursors
            .entry(category.to_string())
            .or_insert_with(SpiralCursor::default);

        let (cx, cy) = spiral_cursor.next_point();

        // Generate a z-coordinate based on the category hash to differentiate categories in 3D noise space
        let mut category_hasher = DefaultHasher::new();
        category.hash(&mut category_hasher);
        // Modulo and multiplier to distribute categories along the z-axis
        let cz = (category_hasher.finish() % 1000) as f32 * self.category_z_multiplier;

        let value = self.probabilities.get_noise_3d(cx as f32, cy as f32, cz);
        value
    }

    /// Generates the next pseudo-random `f32` value using the default category ("__default__").
    ///
    /// Values are typically in the range `[-1.0, 1.0]`.
    pub fn next_value(&mut self) -> f32 {
        self.next_value_for_category("__default__")
    }

    /// Generates the next pseudo-random `f32` value, normalized to the range `[0.0, 1.0)`.
    ///
    /// Optionally, a category can be specified. If `None`, the default category is used.
    /// The normalization assumes the raw noise value is in `[-1.0, 1.0]`.
    ///
    /// # Arguments
    ///
    /// * `category_opt`: An `Option<&str>` specifying the category. `None` uses the default.
    pub fn next_value_normalized(&mut self, category_opt: Option<&str>) -> f32 {
        let val = match category_opt {
            Some(cat) => self.next_value_for_category(cat),
            None => self.next_value(),
        };
        (val + 1.0) * 0.5 // Normalize from [-1,1] to [0,1)
    }

    /// Generates a random `bool` with a specified probability for `true`.
    ///
    /// It uses a normalized random value and compares it against `probability_true`.
    /// Optionally, a category can be specified.
    ///
    /// # Arguments
    ///
    /// * `probability_true`: The probability (0.0 to 1.0) that the function returns `true`. Clamped to `[0.0, 1.0]`.
    /// * `category_opt`: An `Option<&str>` specifying the category. `None` uses the default.
    pub fn next_bool_with_probability(
        &mut self,
        probability_true: f32,
        category_opt: Option<&str>,
    ) -> bool {
        let normalized_val = self.next_value_normalized(category_opt);
        normalized_val < probability_true.clamp(0.0, 1.0)
    }

    /// Generates a random `bool` with an approximately 50/50 chance of being `true` or `false`.
    ///
    /// This is determined by checking if the raw noise value (typically `[-1.0, 1.0]`) is greater than `0.0`.
    /// Optionally, a category can be specified.
    ///
    /// # Arguments
    ///
    /// * `category_opt`: An `Option<&str>` specifying the category. `None` uses the default.
    pub fn next_bool(&mut self, category_opt: Option<&str>) -> bool {
        let val = match category_opt {
            Some(cat) => self.next_value_for_category(cat),
            None => self.next_value(),
        };
        val > 0.0 // Noise values are centered around 0
    }

    /// Generates a random `i32` integer within the specified range `[min, max)`.
    ///
    /// The upper bound `max` is exclusive. If `min >= max`, `min` is returned and a warning is logged.
    /// Optionally, a category can be specified.
    ///
    /// # Arguments
    ///
    /// * `min`: The minimum inclusive value of the range.
    /// * `max`: The maximum exclusive value of the range.
    /// * `category_opt`: An `Option<&str>` specifying the category. `None` uses the default.
    pub fn next_i32_in_range(&mut self, min: i32, max: i32, category_opt: Option<&str>) -> i32 {
        if min >= max {
            warn!(
                "next_i32_in_range called with min ({}) >= max ({}). Returning min.",
                min, max
            );
            return min;
        }
        let normalized_val = self.next_value_normalized(category_opt);
        let range_size = (max - min) as f32;
        min + (normalized_val * range_size).floor() as i32
    }

    /// Generates a random `f32` floating-point number within the specified range `[min, max)`.
    ///
    /// If `min >= max`, `min` is returned and a warning is logged.
    /// Note: Due to floating point precision, `max` might occasionally be inclusive.
    /// Optionally, a category can be specified.
    ///
    /// # Arguments
    ///
    /// * `min`: The minimum inclusive value of the range.
    /// * `max`: The maximum exclusive value of the range.
    /// * `category_opt`: An `Option<&str>` specifying the category. `None` uses the default.
    pub fn next_f32_in_range(&mut self, min: f32, max: f32, category_opt: Option<&str>) -> f32 {
        if min >= max {
            warn!(
                "next_f32_in_range called with min ({}) >= max ({}). Returning min.",
                min, max
            );
            return min;
        }
        let normalized_val = self.next_value_normalized(category_opt);
        let range_size = max - min;
        min + (normalized_val * range_size)
    }

    /// Shuffles the elements of a slice in-place using the Fisher-Yates algorithm.
    ///
    /// The shuffling is deterministic based on the current seed and optional category.
    /// If the slice has fewer than 2 elements, it is not modified.
    ///
    /// # Arguments
    ///
    /// * `slice`: A mutable slice of elements to be shuffled.
    /// * `category_opt`: An `Option<&str>` specifying the category for random number generation
    ///   during shuffling. `None` uses the default category.
    pub fn shuffle<T>(&mut self, slice: &mut [T], category_opt: Option<&str>) {
        let n = slice.len();
        if n < 2 {
            return; // Nothing to shuffle for 0 or 1 elements
        }

        // Fisher-Yates: Iterate from n-1 down to 1
        for i in (1..n).rev() {
            // Choose a random index j from the range [0, i] (inclusive)
            // next_i32_in_range returns a value in [min, max), so max needs to be i + 1
            let j = self.next_i32_in_range(0, (i + 1) as i32, category_opt) as usize;
            slice.swap(i, j);
        }
    }

    /// Resets the `SeedResource` with a new seed value.
    ///
    /// This changes the master seed, re-initializes the underlying noise generator with this new seed,
    /// and clears all categorical spiral cursors. Subsequent random number generations will
    /// produce different sequences based on this new seed.
    ///
    /// # Arguments
    ///
    /// * `new_seed_val`: The `u32` value for the new seed.
    pub fn reset_with_new_seed(&mut self, new_seed_val: u32) {
        self.seed = new_seed_val;
        self.probabilities.set_seed(Some(new_seed_val as i32));
        // Reset other noise parameters to a default state if desired
        self.probabilities.set_frequency(Some(1.0_f32));
        self.probabilities.set_fractal_octaves(Some(1)); // Example: reset octaves
        self.categorical_spiral_cursors.clear(); // Reset cursors for categories
    }
}

impl Default for SeedResource {
    /// Creates a `SeedResource` with a default seed.
    ///
    /// The default seed is generated using `rand::random::<i32>()`. This means that
    /// each default-initialized `SeedResource` will likely have a different random sequence
    /// unless an explicit seed is set later.
    fn default() -> Self {
        // Using an external crate `rand` here is acceptable for the *default* seed,
        // as this typically happens once at initialization if no explicit seed is provided.
        // The actual simulation's random generation then uses `self.probabilities` (the noise generator).
        let seed_number = rand::random::<i32>();
        Self::from_seed(seed_number)
    }
}
