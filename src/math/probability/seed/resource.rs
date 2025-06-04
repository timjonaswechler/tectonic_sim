use bevy::prelude::*;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

#[derive(Resource, Debug, Clone)]
pub struct SeedResource {
    pub seed: u64,
    pub probabilities: super::noise::Noise,
}

impl SeedResource {
    pub fn from_seed(seed: u64) -> Self {
        Self {
            seed,
            probabilities: super::noise::Noise::new()
                .set_seed(seed)
                .set_frequency(1.0)
                .set_fractal_octaves(1),
        }
    }

    pub fn from_text<S: AsRef<str>>(text: S) -> Self {
        let mut hasher = DefaultHasher::new();
        text.as_ref().hash(&mut hasher);
        let seed = hasher.finish();
        Self::from_seed(seed)
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
}
