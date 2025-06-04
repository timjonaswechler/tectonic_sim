
use crate::math::{error::*, types::*, utils::*, geometry::sphere::coordinates::*};
use rand::{Rng, SeedableRng};
use std::f32::consts::{PI, TAU};

/// Verschiedene Sampling-Methoden für Kugel-Oberflächen
#[derive(Debug, Clone, Copy)]
pub enum SamplingMethod {
    /// Uniform verteilte Punkte
    Uniform,
    /// Fibonacci-Spirale (gleichmäßige Verteilung)
    Fibonacci,
    /// Stratified Sampling (Grid-basiert)
    Stratified,
    /// Poisson Disk Sampling (minimaler Abstand)
    PoissonDisk { min_distance: f32 },
}

/// Sphere-Sampler für verschiedene Sampling-Strategien
pub struct SphereSampler {
    method: SamplingMethod,
    radius: f32,
    rng: rand::rngs::StdRng,
}

impl SphereSampler {
    /// Erstellt einen neuen Sampler
    pub fn new(method: SamplingMethod, radius: f32) -> Self {
        Self {
            method,
            radius,
            rng: rand::rngs::StdRng::from_entropy(),
        }
    }
    
    /// Erstellt einen Sampler mit festem Seed
    pub fn with_seed(method: SamplingMethod, radius: f32, seed: u64) -> Self {
        Self {
            method,
            radius,
            rng: rand::rngs::StdRng::seed_from_u64(seed),
        }
    }
    
    /// Generiert eine bestimmte Anzahl von Punkten
    pub fn sample_points(&mut self, count: usize) -> Vec<Point3D> {
        match self.method {
            SamplingMethod::Uniform => self.uniform_sampling(count),
            SamplingMethod::Fibonacci => self.fibonacci_sampling(count),
            SamplingMethod::Stratified => self.stratified_sampling(count),
            SamplingMethod::PoissonDisk { min_distance } => {
                self.poisson_disk_sampling(count, min_distance)
            }
        }
    }
    
    /// Generiert einen einzelnen zufälligen Punkt
    pub fn sample_single_point(&mut self) -> Point3D {
        match self.method {
            SamplingMethod::Uniform => self.uniform_point(),
            _ => self.sample_points(1)[0],
        }
    }
    
    /// Uniform verteilte Punkte (Marsaglia-Methode)
    fn uniform_sampling(&mut self, count: usize) -> Vec<Point3D> {
        (0..count).map(|_| self.uniform_point()).collect()
    }
    
    fn uniform_point(&mut self) -> Point3D {
        loop {
            let x1 = self.rng.gen_range(-1.0..=1.0);
            let x2 = self.rng.gen_range(-1.0..=1.0);
            let length_sq = x1 * x1 + x2 * x2;
            
            if length_sq < 1.0 && length_sq > constants::EPSILON {
                let factor = 2.0 * (1.0 - length_sq).sqrt();
                return Point3D::new(
                    2.0 * x1 * (1.0 - length_sq).sqrt() * self.radius,
                    2.0 * x2 * (1.0 - length_sq).sqrt() * self.radius,
                    (1.0 - 2.0 * length_sq) * self.radius,
                );
            }
        }
    }
    
    /// Fibonacci-Spirale Sampling
    fn fibonacci_sampling(&mut self, count: usize) -> Vec<Point3D> {
        let mut points = Vec::with_capacity(count);
        let golden_ratio = (1.0 + 5.0_f32.sqrt()) * 0.5;
        
        for i in 0..count {
            let theta = TAU * (i as f32 / golden_ratio);
            let phi_cos = 1.0 - 2.0 * (i as f32 / (count - 1) as f32);
            let phi_sin = (1.0 - phi_cos * phi_cos).sqrt();
            
            let x = phi_sin * theta.cos();
            let y = phi_sin * theta.sin();
            let z = phi_cos;
            
            points.push(Point3D::new(x * self.radius, y * self.radius, z * self.radius));
        }
        
        points
    }
    
    /// Stratified Sampling (Grid-basiert)
    fn stratified_sampling(&mut self, count: usize) -> Vec<Point3D> {
        let side_count = (count as f32).sqrt().ceil() as usize;
        let mut points = Vec::with_capacity(side_count * side_count);
        
        for i in 0..side_count {
            for j in 0..side_count {
                // Jitter für bessere Verteilung
                let u = (i as f32 + self.rng.gen::<f32>()) / side_count as f32;
                let v = (j as f32 + self.rng.gen::<f32>()) / side_count as f32;
                
                // Equal-area mapping auf Kugel
                let theta = TAU * u;
                let phi = (2.0 * v - 1.0).acos();
                
                let sin_phi = phi.sin();
                let x = sin_phi * theta.cos();
                let y = sin_phi * theta.sin();
                let z = phi.cos();
                
                points.push(Point3D::new(x * self.radius, y * self.radius, z * self.radius));
                
                if points.len() >= count {
                    break;
                }
            }
            if points.len() >= count {
                break;
            }
        }
        
        points.truncate(count);
        points
    }
    
    /// Poisson Disk Sampling
    fn poisson_disk_sampling(&mut self, target_count: usize, min_distance: f32) -> Vec<Point3D> {
        let mut points = Vec::new();
        let mut active_list = Vec::new();
        
        // Startpunkt
        let first_point = self.uniform_point();
        points.push(first_point);
        active_list.push(0);
        
        let max_attempts = 30;
        
        while !active_list.is_empty() && points.len() < target_count {
            let active_index = self.rng.gen_range(0..active_list.len());
            let point_index = active_list[active_index];
            let current_point = points[point_index];
            
            let mut found_valid = false;
            
            for _ in 0..max_attempts {
                let candidate = self.generate_candidate_around(current_point, min_distance);
                
                if self.is_valid_candidate(&candidate, &points, min_distance) {
                    points.push(candidate);
                    active_list.push(points.len() - 1);
                    found_valid = true;
                    break;
                }
            }
            
            if !found_valid {
                active_list.swap_remove(active_index);
            }
        }
        
        points
    }
    
    fn generate_candidate_around(&mut self, point: Point3D, min_distance: f32) -> Point3D {
        // Generiere Punkt in Kugel-Kappe um den gegebenen Punkt
        let normal = point.normalize();
        let (tangent1, tangent2) = CoordinateConverter::sphere_tangents(point);
        
        // Zufällige Richtung in Tangentialebene
        let angle = self.rng.gen_range(0.0..TAU);
        let radius_factor = self.rng.gen_range(min_distance..min_distance * 2.0);
        
        let offset = tangent1 * (angle.cos() * radius_factor) + 
                     tangent2 * (angle.sin() * radius_factor);
        
        let candidate = point + offset;
        CoordinateConverter::project_to_sphere(candidate, self.radius)
    }
    
    fn is_valid_candidate(&self, candidate: &Point3D, existing_points: &[Point3D], min_distance: f32) -> bool {
        for point in existing_points {
            let distance = candidate.distance(*point);
            if distance < min_distance {
                return false;
            }
        }
        true
    }
}

/// Spezielle Sampling-Patterns
pub struct SamplingPatterns;

impl SamplingPatterns {
    /// Erzeugt Punkte auf den Eckpunkten eines Ikosaeders
    pub fn icosahedron_vertices(radius: f32) -> Vec<Point3D> {
        let phi = constants::GOLDEN_RATIO;
        let a = 1.0 / (phi * phi + 1.0).sqrt();
        let b = phi * a;
        
        let vertices = vec![
            Point3D::new(-a, b, 0.0), Point3D::new(a, b, 0.0),
            Point3D::new(-a, -b, 0.0), Point3D::new(a, -b, 0.0),
            Point3D::new(0.0, -a, b), Point3D::new(0.0, a, b),
            Point3D::new(0.0, -a, -b), Point3D::new(0.0, a, -b),
            Point3D::new(b, 0.0, -a), Point3D::new(b, 0.0, a),
            Point3D::new(-b, 0.0, -a), Point3D::new(-b, 0.0, a),
        ];
        
        vertices.into_iter()
            .map(|v| v.normalize() * radius)
            .collect()
    }
    
    /// Erzeugt Punkte durch geodätische Subdivision eines Ikosaeders
    pub fn geodesic_sphere(radius: f32, subdivisions: usize) -> Vec<Point3D> {
        let mut vertices = Self::icosahedron_vertices(radius);
        
        for _ in 0..subdivisions {
            vertices = Self::subdivide_geodesic(vertices, radius);
        }
        
        vertices
    }
    
    fn subdivide_geodesic(vertices: Vec<Point3D>, radius: f32) -> Vec<Point3D> {
        // Vereinfachte geodätische Subdivision
        // In einer vollständigen Implementierung würde man die Dreiecks-Topologie verwalten
        let mut new_vertices = vertices.clone();
        
        for i in 0..vertices.len() {
            for j in (i + 1)..vertices.len() {
                let midpoint = (vertices[i] + vertices[j]) * 0.5;
                let projected = CoordinateConverter::project_to_sphere(midpoint, radius);
                new_vertices.push(projected);
            }
        }
        
        new_vertices
    }
    
    /// Erzeugt Punkte entlang Breitengrad-Kreisen
    pub fn latitude_circles(radius: f32, num_circles: usize, points_per_circle: usize) -> Vec<Point3D> {
        let mut points = Vec::new();
        
        for i in 0..num_circles {
            let latitude = PI * (i as f32 / (num_circles - 1) as f32) - PI * 0.5;
            let circle_radius = radius * latitude.cos();
            let height = radius * latitude.sin();
            
            for j in 0..points_per_circle {
                let longitude = TAU * (j as f32 / points_per_circle as f32);
                
                let x = circle_radius * longitude.cos();
                let z = circle_radius * longitude.sin();
                let y = height;
                
                points.push(Point3D::new(x, y, z));
            }
        }
        
        points
    }
    
    /// Generiert zufällige Punkte in einem Kugel-Segment (Calotte)
    pub fn spherical_cap(radius: f32, height: f32, count: usize, seed: Option<u64>) -> Vec<Point3D> {
        let mut rng = match seed {
            Some(s) => rand::rngs::StdRng::seed_from_u64(s),
            None => rand::rngs::StdRng::from_entropy(),
        };
        
        let mut points = Vec::with_capacity(count);
        let cos_max_polar = (radius - height) / radius;
        
        for _ in 0..count {
            // Uniform sampling in spherical cap
            let cos_polar = rng.gen_range(cos_max_polar..=1.0);
            let polar = cos_polar.acos();
            let azimuth = rng.gen_range(0.0..TAU);
            
            let spherical = SphericalCoordinates::new(radius, polar, azimuth).unwrap();
            points.push(spherical.to_cartesian());
        }
        
        points
    }
}

/// Sampling-Utilities für spezielle Anwendungen
pub struct SamplingUtils;

impl SamplingUtils {
    /// Berechnet die lokale Dichte von Punkten
    pub fn calculate_density(points: &[Point3D], query_point: Point3D, radius: f32) -> f32 {
        let neighbors = points.iter()
            .filter(|&&p| query_point.distance(p) <= radius)
            .count();
        
        neighbors as f32 / (4.0 * PI * radius * radius)
    }
    
    /// Findet die k nächsten Nachbarn
    pub fn k_nearest_neighbors(points: &[Point3D], query_point: Point3D, k: usize) -> Vec<(Point3D, f32)> {
        let mut distances: Vec<(Point3D, f32)> = points.iter()
            .map(|&p| (p, query_point.distance(p)))
            .collect();
        
        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        distances.truncate(k);
        distances
    }
    
    /// Berechnet Voronoi-Zell-Volumen näherungsweise
    pub fn approximate_voronoi_volume(points: &[Point3D], point_index: usize, radius: f32) -> f32 {
        if point_index >= points.len() {
            return 0.0;
        }
        
        let point = points[point_index];
        let neighbors = Self::k_nearest_neighbors(points, point, 10);
        
        if neighbors.is_empty() {
            return 4.0 * PI * radius * radius * radius / 3.0 / points.len() as f32;
        }
        
        let avg_neighbor_distance = neighbors.iter()
            .map(|(_, dist)| *dist)
            .sum::<f32>() / neighbors.len() as f32;
        
        4.0 * PI * (avg_neighbor_distance * 0.5).powi(3) / 3.0
    }
    
    /// Prüft ob Punkte-Verteilung uniform ist (Lloyd-Index)
    pub fn uniformity_measure(points: &[Point3D], radius: f32) -> f32 {
        let expected_neighbors = points.len() as f32 / (4.0 * PI);
        let mut variance_sum = 0.0;
        
        for &point in points {
            let actual_neighbors = Self::k_nearest_neighbors(points, point, 20).len() as f32;
            let diff = actual_neighbors - expected_neighbors;
            variance_sum += diff * diff;
        }
        
        (variance_sum / points.len() as f32).sqrt() / expected_neighbors
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_uniform_sampling() {
        let mut sampler = SphereSampler::with_seed(SamplingMethod::Uniform, 1.0, 42);
        let points = sampler.sample_points(100);
        
        assert_eq!(points.len(), 100);
        
        // Alle Punkte sollten auf der Einheitskugel liegen
        for point in &points {
            assert!(comparison::nearly_equal(point.length(), 1.0));
        }
    }
    
    #[test]
    fn test_fibonacci_sampling() {
        let mut sampler = SphereSampler::new(SamplingMethod::Fibonacci, 2.0);
        let points = sampler.sample_points(50);
        
        assert_eq!(points.len(), 50);
        
        for point in &points {
            assert!(comparison::nearly_equal(point.length(), 2.0));
        }
    }
    
    #[test]
    fn test_icosahedron_vertices() {
        let vertices = SamplingPatterns::icosahedron_vertices(1.0);
        assert_eq!(vertices.len(), 12);
        
        for vertex in &vertices {
            assert!(comparison::nearly_equal(vertex.length(), 1.0));
        }
    }
    
    #[test]
    fn test_sampling_utils() {
        let points = vec![
            Point3D::new(1.0, 0.0, 0.0),
            Point3D::new(0.0, 1.0, 0.0),
            Point3D::new(0.0, 0.0, 1.0),
        ];
        
        let query = Point3D::new(0.5, 0.5, 0.5);
        let neighbors = SamplingUtils::k_nearest_neighbors(&points, query, 2);
        
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors[0].1 <= neighbors[1].1); // Sortiert nach Distanz
    }
}