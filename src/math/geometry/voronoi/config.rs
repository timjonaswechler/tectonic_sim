// src/math/geometry/voronoi/config.rs
use crate::math::error::MathResult;

#[derive(Debug, Clone)]
pub struct VoronoiConfig {
    pub initial_points: usize,
    pub lloyd_iterations: usize,
    pub target_cells: usize,
    pub smoothing_iterations: usize,
    pub boundary_factor: f64,
    pub seed: Option<u64>,
}

impl VoronoiConfig {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_points(mut self, count: usize) -> Self {
        self.initial_points = count;
        self
    }

    pub fn with_lloyd_iterations(mut self, iterations: usize) -> Self {
        self.lloyd_iterations = iterations;
        self
    }

    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    pub fn validate(&self) -> MathResult<()> {
        if self.initial_points < 3 {
            return Err(MathError::InvalidConfiguration {
                message: "Need at least 3 initial points".to_string(),
            });
        }

        if self.target_cells == 0 {
            return Err(MathError::InvalidConfiguration {
                message: "Target cells must be greater than 0".to_string(),
            });
        }

        if !(0.0..1.0).contains(&self.boundary_factor) {
            return Err(MathError::InvalidConfiguration {
                message: "Boundary factor must be between 0.0 and 1.0".to_string(),
            });
        }

        Ok(())
    }
}

impl Default for VoronoiConfig {
    fn default() -> Self {
        Self {
            initial_points: 500,
            lloyd_iterations: 15,
            target_cells: 5,
            smoothing_iterations: 1,
            boundary_factor: 0.15,
            seed: None,
        }
    }
}
