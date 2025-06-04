// src/math/geometry/voronoi/generator.rs
use super::config::VoronoiConfig;
use super::lloyd::LloydRelaxation;
use super::merger::PolygonMerger;
use crate::math::{
    error::{MathError, MathResult},
    geometry::polygon::smoothing::ChaikinSmoother,
    types::*,
};
use spade::Point2;

pub struct VoronoiGenerator {
    config: VoronoiConfig,
}

impl VoronoiGenerator {
    pub fn new(config: VoronoiConfig) -> MathResult<Self> {
        config.validate()?;
        Ok(Self { config })
    }

    pub fn generate(&self) -> MathResult<Vec<Point2<f64>>> {
        // 1. Initiale Punkte generieren
        let initial_points = self.generate_initial_points()?;

        // 2. Lloyd-Relaxation anwenden
        let relaxation = LloydRelaxation::new(&self.config);
        let triangulation = relaxation.apply(&initial_points)?;

        // 3. Innere Zellen finden und verschmelzen
        let merger = PolygonMerger::new(&self.config);
        let merged_polygon = merger.merge_cells(&triangulation)?;

        // 4. Glätten
        let smoother = ChaikinSmoother::new(self.config.smoothing_iterations);
        let smoothed = smoother.apply(&merged_polygon)?;

        // 5. Normalisieren
        Ok(self.normalize_polygon(&smoothed))
    }

    fn generate_initial_points(&self) -> MathResult<Vec<Point2<f64>>> {
        use rand::Rng;

        let mut rng = self.create_rng();
        let mut points = Vec::with_capacity(self.config.initial_points);

        for _ in 0..self.config.initial_points {
            let x = rng.gen_range(-1.0..=1.0);
            let y = rng.gen_range(-1.0..=1.0);
            points.push(Point2::new(x, y));
        }

        Ok(points)
    }

    fn normalize_polygon(&self, polygon: &[Point2<f64>]) -> Vec<Point2<f64>> {
        // Implementierung der Normalisierung
        // (vereinfacht hier)
        polygon.to_vec()
    }

    fn create_rng(&self) -> impl rand::Rng {
        use rand::SeedableRng;

        match self.config.seed {
            Some(seed) => rand::rngs::StdRng::seed_from_u64(seed),
            None => rand::rngs::StdRng::from_entropy(),
        }
    }
}
