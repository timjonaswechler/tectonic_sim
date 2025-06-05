// src/math/geometry/voronoi/generator.rs
use super::config::VoronoiConfig;
use super::lloyd::LloydRelaxation;
use super::merger::PolygonMerger;
use crate::math::probability::SeedResource;

use crate::math::{error::MathResult, geometry::polygon::smoothing::ChaikinSmoother};
use crate::next_random_value;
use spade::Point2;

pub struct VoronoiGenerator {
    config: VoronoiConfig,
}

impl VoronoiGenerator {
    pub fn new(config: VoronoiConfig) -> MathResult<Self> {
        config.validate()?;
        Ok(Self { config })
    }

    pub fn generate(&self, seed_resource: Res<SeedResource>) -> MathResult<Vec<Point2<f64>>> {
        // 1. Initiale Punkte generieren
        let initial_points = self.generate_initial_points(seed_resource)?;

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

    fn generate_initial_points(
        &self,
        seed_resource: Res<SeedResource>,
    ) -> MathResult<Vec<Point2<f64>>> {
        let mut points = Vec::with_capacity(self.config.initial_points);

        for _ in 0..self.config.initial_points {
            let x = next_random_value!(seed_resource);
            let y = next_random_value!(seed_resource);
            points.push(Point2::new(x, y));
        }

        Ok(points)
    }

    fn normalize_polygon(&self, polygon: &[Point2<f64>]) -> Vec<Point2<f64>> {
        // Implementierung der Normalisierung
        // (vereinfacht hier)
        polygon.to_vec()
    }
}
