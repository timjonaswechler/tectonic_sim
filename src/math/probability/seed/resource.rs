use bevy::prelude::*;

#[derive(Resource, Debug, Clone)]
pub struct SeedResource {
    pub seed: u64,
    pub probabilities: super::noise::Noise,
}

impl Default for SeedResource {
    fn default() -> Self {
        let seed_number = rand::random::<u64>();
        Self {
            seed: seed_number,
            probabilities: super::noise::Noise::new()
                .set_seed(seed_number)
                .set_frequency(1.0)
                .set_fractal_octaves(1),
        }
    }
}
