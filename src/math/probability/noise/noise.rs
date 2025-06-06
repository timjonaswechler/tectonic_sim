// Importiere Typen und Konstanten aus den anderen Moduldateien
use super::constants::*;
use super::types::*;
use bevy::prelude::Resource;

// Float und FloatOps müssen hier bekannt sein.
// Sie werden im übergeordneten fastnoise_mod/mod.rs oder src/lib.rs definiert/importiert
// und dann hier mit `use super::{Float, FloatOps};` oder `use crate::{Float, FloatOps};` importiert.
// Für dieses Beispiel nehmen wir an, sie sind im crate root oder `fastnoise_mod/mod.rs`
// und wir importieren sie relativ.
// pub(crate) type Float = f32;

// #[cfg(feature = "libm")]
// use num_traits::float::FloatCore as FloatOps; // FloatCore für no_std

// #[cfg(all(feature = "std", not(feature = "libm")))]
// #[cfg(feature = "f32")]
// use f32 as FloatOps;

// #[cfg(all(not(feature = "std"), not(feature = "libm")))]
// compile_error!("`fastnoise-lite` crate: either the \"std\" or \"libm\" feature must be enabled");
// #[cfg(all(not(feature = "std"), not(feature = "libm")))]
// use Float as FloatOps; // Dummy
// #[cfg(all(not(feature = "std"), not(feature = "libm")))]
// impl DummyFloatExt for Float {}
// #[cfg(all(not(feature = "std"), not(feature = "libm")))]
// trait DummyFloatExt: Sized {
//     fn sqrt(self) -> Self {
//         unimplemented!()
//     }
//     fn trunc(self) -> Self {
//         unimplemented!()
//     }
//     fn abs(self) -> Self {
//         unimplemented!()
//     }
// }
use f32 as FloatOps;
use f32 as Float;

#[derive(Clone, Debug, Resource)] // Resource Trait für Bevy
pub struct Noise {
    pub seed: i32,
    pub frequency: f32,
    pub noise_type: NoiseType,
    pub rotation_type_3d: RotationType3D,
    transform_type_3d: TransformType3D, // ist private, bleibt so

    pub fractal_type: FractalType,
    pub octaves: u32,
    pub lacunarity: f32,
    pub gain: f32,
    pub weighted_strength: f32,
    pub ping_pong_strength: f32,

    fractal_bounding: f32, // ist private, bleibt so

    pub cellular_distance_function: CellularDistanceFunction,
    pub cellular_return_type: CellularReturnType,
    pub cellular_jitter_modifier: f32,

    pub domain_warp_type: DomainWarpType,
    warp_transform_type_3d: TransformType3D, // ist private, bleibt so
    pub domain_warp_amp: f32,
}

impl Default for Noise {
    fn default() -> Self {
        Self {
            seed: 1337,
            frequency: 1.0,
            noise_type: NoiseType::OpenSimplex2,
            rotation_type_3d: RotationType3D::None,
            transform_type_3d: TransformType3D::DefaultOpenSimplex2,

            fractal_type: FractalType::None,
            octaves: 3,
            lacunarity: 2.,
            gain: 0.5,
            weighted_strength: 0.,
            ping_pong_strength: 2.,

            fractal_bounding: 1. / 1.75,

            cellular_distance_function: CellularDistanceFunction::EuclideanSq,
            cellular_return_type: CellularReturnType::Distance,
            cellular_jitter_modifier: 1.,

            domain_warp_type: DomainWarpType::OpenSimplex2,
            warp_transform_type_3d: TransformType3D::DefaultOpenSimplex2,
            domain_warp_amp: 1.,
        }
    }
}

impl Noise {
    // =====================
    // Constructor functions
    // =====================

    /// # Constructor
    ///
    /// Create new FastNoise object with the default seed of `1337`.
    pub fn new() -> Self {
        Self::default()
    }

    // ============================
    // Setter convenience functions
    // ============================

    /// Sets seed used for all noise types.
    ///
    /// If set to [`None`], it is reset to its default: `1337`.
    pub fn set_seed(&mut self, seed: i32) {
        self.seed = seed;
    }

    /// Sets frequency used for all noise types.
    ///
    /// If set to [`None`], it is reset to its default: `0.01`.
    pub fn set_frequency(&mut self, frequency: f32) {
        self.frequency = frequency;
    }

    /// Sets noise algorithm used for [`get_noise_2d`](Self::get_noise_2d)/[`get_noise_3d`](Self::get_noise_3d).
    ///
    /// If set to [`None`], it is reset to its default: [`NoiseType::OpenSimplex2`].
    pub fn set_noise_type(&mut self, noise_type: Option<NoiseType>) {
        self.noise_type = noise_type.unwrap_or(Self::default().noise_type);
        self.update_transform_type_3d();
    }

    /// Sets domain rotation type for 3D Noise and 3D DomainWarp.
    /// Can aid in reducing directional artifacts when sampling a 2D plane in 3D.
    ///
    /// If set to [`None`], it is reset to its default: [`RotationType3D::None`].
    pub fn set_rotation_type_3d(&mut self, rotation_type_3d: Option<RotationType3D>) {
        self.rotation_type_3d = rotation_type_3d.unwrap_or(Self::default().rotation_type_3d);
        self.update_transform_type_3d();
        self.update_warp_transform_type_3d();
    }

    /// Sets method for combining octaves in all fractal noise types.
    ///
    /// If set to [`None`], it is reset to its default: [`FractalType::None`].
    ///
    /// Note: [`FractalType::DomainWarpProgressive`]/[`FractalType::DomainWarpIndependent`] only affects [`domain_warp_2d`](Self::domain_warp_2d).
    pub fn set_fractal_type(&mut self, fractal_type: Option<FractalType>) {
        self.fractal_type = fractal_type.unwrap_or(Self::default().fractal_type);
    }

    /// Sets octave count for all fractal noise types.
    ///
    /// If set to [`None`], it is reset to its default: `3`.
    pub fn set_fractal_octaves(&mut self, octaves: Option<u32>) {
        self.octaves = octaves.unwrap_or(Self::default().octaves);
        self.calculate_fractal_bounding();
    }

    /// Sets octave lacunarity for all fractal noise types.
    ///
    /// If set to [`None`], it is reset to its default: `2.0`.
    pub fn set_fractal_lacunarity(&mut self, lacunarity: Option<f32>) {
        self.lacunarity = lacunarity.unwrap_or(Self::default().lacunarity);
    }

    /// Sets octave gain for all fractal noise types.
    ///
    /// If set to [`None`], it is reset to its default: `0.5`.
    pub fn set_fractal_gain(&mut self, gain: Option<f32>) {
        self.gain = gain.unwrap_or(Self::default().gain);
        self.calculate_fractal_bounding();
    }

    /// Sets octave weighting for all none DomainWarp fractal types.
    ///
    /// If set to [`None`], it is reset to its default: `0.0`.
    ///
    /// Note: Keep between 0..1 to maintain -1..1 output bounding.
    pub fn set_fractal_weighted_strength(&mut self, weighted_strength: Option<f32>) {
        self.weighted_strength = weighted_strength.unwrap_or(Self::default().weighted_strength);
    }

    /// Sets strength of the fractal ping pong effect.
    ///
    /// If set to [`None`], it is reset to its default: `2.0`.
    pub fn set_fractal_ping_pong_strength(&mut self, ping_pong_strength: Option<f32>) {
        self.ping_pong_strength = ping_pong_strength.unwrap_or(Self::default().ping_pong_strength);
    }

    /// Sets distance function used in cellular noise calculations.
    ///
    /// If set to [`None`], it is reset to its default: [`CellularDistanceFunction::EuclideanSq`].
    pub fn set_cellular_distance_function(
        &mut self,
        cellular_distance_function: Option<CellularDistanceFunction>,
    ) {
        self.cellular_distance_function =
            cellular_distance_function.unwrap_or(Self::default().cellular_distance_function);
    }

    /// Sets return type from cellular noise calculations.
    ///
    /// If set to [`None`], it is reset to its default: [`CellularReturnType::Distance`].
    pub fn set_cellular_return_type(&mut self, cellular_return_type: Option<CellularReturnType>) {
        self.cellular_return_type =
            cellular_return_type.unwrap_or(Self::default().cellular_return_type);
    }

    /// Sets the maximum distance a cellular point can move from its grid position.
    ///
    /// If set to [`None`], it is reset to its default: `1.0`.
    ///
    /// Note: Setting this higher than 1 will cause artifacts.
    pub fn set_cellular_jitter(&mut self, cellular_jitter: Option<f32>) {
        self.cellular_jitter_modifier =
            cellular_jitter.unwrap_or(Self::default().cellular_jitter_modifier);
    }

    /// Sets the warp algorithm when using [`domain_warp_2d`](Self::domain_warp_2d).
    ///
    /// If set to [`None`], it is reset to its default: [`DomainWarpType::OpenSimplex2`].
    pub fn set_domain_warp_type(&mut self, domain_warp_type: Option<DomainWarpType>) {
        self.domain_warp_type = domain_warp_type.unwrap_or(Self::default().domain_warp_type);
        self.update_warp_transform_type_3d();
    }

    /// Sets the maximum warp distance from original position when using [`domain_warp_2d`](Self::domain_warp_2d).
    ///
    /// If set to [`None`], it is reset to its default: `1.0`.
    pub fn set_domain_warp_amp(&mut self, domain_warp_amp: Option<f32>) {
        self.domain_warp_amp = domain_warp_amp.unwrap_or(Self::default().domain_warp_amp);
    }

    // =========================
    // Noise generator functions
    // =========================

    /// 2D noise at given position using current settings.
    ///
    /// Noise output bounded between -1..1.
    ///
    /// Example usage:
    /// ```rs
    /// let noise = get_noise_2d(x, y); // Value in the -1..1 range
    /// let noise = (noise + 1.) / 2.; // Consider remapping it to the 0..1 range
    /// ```
    pub fn get_noise_2d(&self, x: Float, y: Float) -> f32 {
        let (x, y) = self.transform_noise_coordinate_2d(x, y);

        match self.fractal_type {
            FractalType::FBm => self.gen_fractal_fbm_2d(x, y),
            FractalType::Ridged => self.gen_fractal_ridged_2d(x, y),
            FractalType::PingPong => self.gen_fractal_ping_pong_2d(x, y),
            _ => self.gen_noise_single_2d(self.seed, x, y),
        }
    }

    /// 3D noise at given position using current settings.
    ///
    /// Noise output is bounded between -1..1.
    ///
    /// Example usage:
    /// ```rs
    /// let noise = get_noise_3d(x, y, z); // Value in the -1..1 range
    /// let noise = (noise + 1.) / 2.; // Consider remapping it to the 0..1 range
    /// ```
    pub fn get_noise_3d(&self, x: Float, y: Float, z: Float) -> f32 {
        let (x, y, z) = self.transform_noise_coordinate_3d(x, y, z);

        match self.fractal_type {
            FractalType::FBm => self.gen_fractal_fbm_3d(x, y, z),
            FractalType::Ridged => self.gen_fractal_ridged_3d(x, y, z),
            FractalType::PingPong => self.gen_fractal_ping_pong_3d(x, y, z),
            _ => self.gen_noise_single_3d(self.seed, x, y, z),
        }
    }

    // ===============================================
    // Domain warp coordinate transformation functions
    // ===============================================

    /// 2D warps the input position using current domain warp settings.
    ///
    /// Example usage:
    /// ```rs
    /// let (x, y) = domain_warp_2d(x, y);
    /// let noise = get_noise_2d(x, y); // Value in the -1..1 range
    /// let noise = (noise + 1.) / 2.; // Consider remapping it to the 0..1 range
    /// ```
    pub fn domain_warp_2d(&self, x: Float, y: Float) -> (Float, Float) {
        match self.fractal_type {
            FractalType::DomainWarpProgressive => self.domain_warp_fractal_progressive_2d(x, y),
            FractalType::DomainWarpIndependent => self.domain_warp_fractal_independent_2d(x, y),
            _ => self.domain_warp_single_2d(x, y),
        }
    }

    /// 3D warps the input position using current domain warp settings.
    ///
    /// Example usage:
    /// ```rs
    /// let (x, y, z) = domain_warp_3d(x, y, z);
    /// let noise = get_noise_3d(x, y, z); // Value in the -1..1 range
    /// let noise = (noise + 1.) / 2.; // Consider remapping it to the 0..1 range
    /// ```
    pub fn domain_warp_3d(&self, x: Float, y: Float, z: Float) -> (Float, Float, Float) {
        match self.fractal_type {
            FractalType::DomainWarpProgressive => self.domain_warp_fractal_progressive_3d(x, y, z),
            FractalType::DomainWarpIndependent => self.domain_warp_fractal_independent_3d(x, y, z),
            _ => self.domain_warp_single_3d(x, y, z),
        }
    }

    // ==============
    // Math functions
    // ==============

    #[inline(always)]
    fn fast_floor(f: Float) -> i32 {
        if f >= 0. { f as i32 } else { f as i32 - 1 }
    }

    #[inline(always)]
    fn fast_round(f: Float) -> i32 {
        if f >= 0. {
            (f + 0.5) as i32
        } else {
            (f - 0.5) as i32
        }
    }

    #[inline(always)]
    fn lerp(a: f32, b: f32, t: f32) -> f32 {
        a + t * (b - a)
    }

    #[inline(always)]
    fn interp_hermite(t: f32) -> f32 {
        t * t * (t * -2. + 3.)
    }

    #[inline(always)]
    fn interp_quintic(t: f32) -> f32 {
        t * t * t * (t * (t * 6. - 15.) + 10.)
    }

    #[inline(always)]
    fn cubic_lerp(a: f32, b: f32, c: f32, d: f32, t: f32) -> f32 {
        let p = (d - c) - (a - b);
        t * t * t * p + t * t * ((a - b) - p) + t * (c - a) + b
    }

    #[inline(always)]
    fn ping_pong(t: f32) -> f32 {
        let t = t - FloatOps::trunc(t * 0.5) * 2.;

        if t < 1. { t } else { 2. - t }
    }

    fn calculate_fractal_bounding(&mut self) {
        let gain = FloatOps::abs(self.gain);
        let mut amp = gain;
        let mut amp_fractal = 1.;
        for _ in 1..self.octaves {
            amp_fractal += amp;
            amp *= gain;
        }
        self.fractal_bounding = 1. / amp_fractal;
    }

    // ==============
    // Hash functions
    // ==============

    #[inline(always)]
    fn hash_2d(seed: i32, x_primed: i32, y_primed: i32) -> i32 {
        let hash = seed ^ x_primed ^ y_primed;
        hash.wrapping_mul(0x27d4eb2d)
    }

    #[inline(always)]
    fn hash_3d(seed: i32, x_primed: i32, y_primed: i32, z_primed: i32) -> i32 {
        let hash = seed ^ x_primed ^ y_primed ^ z_primed;
        hash.wrapping_mul(0x27d4eb2d)
    }

    // ===================================
    // Internal noise generator algorithms
    // ===================================

    #[inline(always)]
    fn val_coord_2d(seed: i32, x_primed: i32, y_primed: i32) -> f32 {
        let hash = Self::hash_2d(seed, x_primed, y_primed);
        let hash = hash.wrapping_mul(hash);
        let hash = hash ^ (hash << 19);
        hash as f32 * (1. / 2147483648.)
    }

    #[inline(always)]
    fn val_coord_3d(seed: i32, x_primed: i32, y_primed: i32, z_primed: i32) -> f32 {
        let hash = Self::hash_3d(seed, x_primed, y_primed, z_primed);

        let hash = hash.wrapping_mul(hash);
        let hash = hash ^ (hash << 19);
        hash as f32 * (1. / 2147483648.)
    }

    #[inline(always)]
    fn grad_coord_2d(seed: i32, x_primed: i32, y_primed: i32, xd: f32, yd: f32) -> f32 {
        let hash = Self::hash_2d(seed, x_primed, y_primed);
        let hash = hash ^ (hash >> 15);
        let hash = hash & (127 << 1);

        let xg = GRADIENTS_2D[hash as usize];
        let yg = GRADIENTS_2D[(hash | 1) as usize];

        xd * xg + yd * yg
    }

    #[inline(always)]
    fn grad_coord_3d(
        seed: i32,
        x_primed: i32,
        y_primed: i32,
        z_primed: i32,
        xd: f32,
        yd: f32,
        zd: f32,
    ) -> f32 {
        let hash = Self::hash_3d(seed, x_primed, y_primed, z_primed);
        let hash = hash ^ (hash >> 15);
        let hash = hash & (63 << 2);

        let xg = GRADIENTS_3D[hash as usize];
        let yg = GRADIENTS_3D[(hash | 1) as usize];
        let zg = GRADIENTS_3D[(hash | 2) as usize];

        xd * xg + yd * yg + zd * zg
    }

    #[inline(always)]
    fn grad_coord_out_2d(seed: i32, x_primed: i32, y_primed: i32) -> (f32, f32) {
        let hash = Self::hash_2d(seed, x_primed, y_primed) & (255 << 1);

        let xo = RAND_VECS_2D[hash as usize];
        let yo = RAND_VECS_2D[(hash | 1) as usize];

        (xo, yo)
    }

    #[inline(always)]
    fn grad_coord_out_3d(
        seed: i32,
        x_primed: i32,
        y_primed: i32,
        z_primed: i32,
    ) -> (f32, f32, f32) {
        let hash = Self::hash_3d(seed, x_primed, y_primed, z_primed) & (255 << 2);

        let xo = RAND_VECS_3D[hash as usize];
        let yo = RAND_VECS_3D[(hash | 1) as usize];
        let zo = RAND_VECS_3D[(hash | 2) as usize];

        (xo, yo, zo)
    }

    #[inline(always)]
    fn grad_coord_dual_2d(seed: i32, x_primed: i32, y_primed: i32, xd: f32, yd: f32) -> (f32, f32) {
        let hash = Self::hash_2d(seed, x_primed, y_primed);
        let index1 = hash & (127 << 1);
        let index2 = (hash >> 7) & (255 << 1);

        let xg = GRADIENTS_2D[index1 as usize];
        let yg = GRADIENTS_2D[(index1 | 1) as usize];
        let value = xd * xg + yd * yg;

        let xgo = RAND_VECS_2D[index2 as usize];
        let ygo = RAND_VECS_2D[(index2 | 1) as usize];

        let xo = value * xgo;
        let yo = value * ygo;

        (xo, yo)
    }

    #[inline(always)]
    fn grad_coord_dual_3d(
        seed: i32,
        x_primed: i32,
        y_primed: i32,
        z_primed: i32,
        xd: f32,
        yd: f32,
        zd: f32,
    ) -> (f32, f32, f32) {
        let hash = Self::hash_3d(seed, x_primed, y_primed, z_primed);
        let index1 = hash & (63 << 2);
        let index2 = (hash >> 6) & (255 << 2);

        let xg = GRADIENTS_3D[index1 as usize];
        let yg = GRADIENTS_3D[(index1 | 1) as usize];
        let zg = GRADIENTS_3D[(index1 | 2) as usize];
        let value = xd * xg + yd * yg + zd * zg;

        let xgo = RAND_VECS_3D[index2 as usize];
        let ygo = RAND_VECS_3D[(index2 | 1) as usize];
        let zgo = RAND_VECS_3D[(index2 | 2) as usize];

        let xo = value * xgo;
        let yo = value * ygo;
        let zo = value * zgo;

        (xo, yo, zo)
    }

    // Generic noise gen

    fn gen_noise_single_2d(&self, seed: i32, x: Float, y: Float) -> f32 {
        match self.noise_type {
            NoiseType::OpenSimplex2 => self.single_simplex_2d(seed, x, y),
            NoiseType::OpenSimplex2S => self.single_open_simplex_2s_2d(seed, x, y),
            NoiseType::Cellular => self.single_cellular_2d(seed, x, y),
            NoiseType::Perlin => self.single_perlin_2d(seed, x, y),
            NoiseType::ValueCubic => self.single_value_cubic_2d(seed, x, y),
            NoiseType::Value => self.single_value_2d(seed, x, y),
        }
    }

    fn gen_noise_single_3d(&self, seed: i32, x: Float, y: Float, z: Float) -> f32 {
        match self.noise_type {
            NoiseType::OpenSimplex2 => self.single_open_simplex_2(seed, x, y, z),
            NoiseType::OpenSimplex2S => self.single_open_simplex_2s_3d(seed, x, y, z),
            NoiseType::Cellular => self.single_cellular_3d(seed, x, y, z),
            NoiseType::Perlin => self.single_perlin_3d(seed, x, y, z),
            NoiseType::ValueCubic => self.single_value_cubic_3d(seed, x, y, z),
            NoiseType::Value => self.single_value_3d(seed, x, y, z),
        }
    }

    // Noise Coordinate Transforms (frequency, and possible skew or rotation)

    #[inline(always)]
    fn transform_noise_coordinate_2d(&self, x: Float, y: Float) -> (Float, Float) {
        let mut x = x * self.frequency as Float;
        let mut y = y * self.frequency as Float;

        match self.noise_type {
            NoiseType::OpenSimplex2 | NoiseType::OpenSimplex2S => {
                let sqrt3 = 1.7320508075688772935274463415059;
                let f2 = 0.5 * (sqrt3 - 1.);
                let t = (x + y) * f2;

                x += t;
                y += t;
            }
            _ => {}
        }

        (x, y)
    }

    #[inline(always)]
    fn transform_noise_coordinate_3d(&self, x: Float, y: Float, z: Float) -> (Float, Float, Float) {
        let mut x = x * self.frequency as Float;
        let mut y = y * self.frequency as Float;
        let mut z = z * self.frequency as Float;

        match self.transform_type_3d {
            TransformType3D::ImproveXYPlanes => {
                let xy = x + y;
                let s2 = xy * -0.211324865405187;
                z *= 0.577350269189626;
                x += s2 - z;
                y = y + s2 - z;
                z += xy * 0.577350269189626;
            }
            TransformType3D::ImproveXZPlanes => {
                let xz = x + z;
                let s2 = xz * -0.211324865405187;
                y *= 0.577350269189626;
                x += s2 - y;
                z += s2 - y;
                y += xz * 0.577350269189626;
            }
            TransformType3D::DefaultOpenSimplex2 => {
                let r3 = 2. / 3.;
                let r = (x + y + z) * r3; // Rotation, not skew
                x = r - x;
                y = r - y;
                z = r - z;
            }
            _ => {}
        }

        (x, y, z)
    }

    fn update_transform_type_3d(&mut self) {
        match self.rotation_type_3d {
            RotationType3D::ImproveXYPlanes => {
                self.transform_type_3d = TransformType3D::ImproveXYPlanes;
            }
            RotationType3D::ImproveXZPlanes => {
                self.transform_type_3d = TransformType3D::ImproveXZPlanes;
            }
            _ => match self.noise_type {
                NoiseType::OpenSimplex2 | NoiseType::OpenSimplex2S => {
                    self.transform_type_3d = TransformType3D::DefaultOpenSimplex2;
                }
                _ => {
                    self.transform_type_3d = TransformType3D::None;
                }
            },
        }
    }

    // Domain Warp Coordinate Transforms

    #[inline(always)]
    fn transform_domain_warp_coordinate_2d(&self, x: Float, y: Float) -> (Float, Float) {
        let mut x = x;
        let mut y = y;

        match self.domain_warp_type {
            DomainWarpType::OpenSimplex2 | DomainWarpType::OpenSimplex2Reduced => {
                let sqrt3 = 1.7320508075688772935274463415059;
                let f2 = 0.5 * (sqrt3 - 1.);
                let t = (x + y) * f2;

                x += t;
                y += t;
            }
            _ => {}
        }

        (x, y)
    }

    #[inline(always)]
    fn transform_domain_warp_coordinate_3d(
        &self,
        x: Float,
        y: Float,
        z: Float,
    ) -> (Float, Float, Float) {
        let mut x = x;
        let mut y = y;
        let mut z = z;

        match self.warp_transform_type_3d {
            TransformType3D::ImproveXYPlanes => {
                let xy = x + y;
                let s2 = xy * -0.211324865405187;
                z *= 0.577350269189626;
                x += s2 - z;
                y = y + s2 - z;
                z += xy * 0.577350269189626;
            }
            TransformType3D::ImproveXZPlanes => {
                let xz = x + z;
                let s2 = xz * -0.211324865405187;
                y *= 0.577350269189626;
                x += s2 - y;
                z += s2 - y;
                y += xz * 0.577350269189626;
            }
            TransformType3D::DefaultOpenSimplex2 => {
                let r3 = 2. / 3.;
                let r = (x + y + z) * r3; // Rotation, not skew
                x = r - x;
                y = r - y;
                z = r - z;
            }
            _ => {}
        }

        (x, y, z)
    }

    fn update_warp_transform_type_3d(&mut self) {
        match self.rotation_type_3d {
            RotationType3D::ImproveXYPlanes => {
                self.warp_transform_type_3d = TransformType3D::ImproveXYPlanes;
            }
            RotationType3D::ImproveXZPlanes => {
                self.warp_transform_type_3d = TransformType3D::ImproveXZPlanes;
            }
            _ => match self.domain_warp_type {
                DomainWarpType::OpenSimplex2 | DomainWarpType::OpenSimplex2Reduced => {
                    self.warp_transform_type_3d = TransformType3D::DefaultOpenSimplex2;
                }
                _ => {
                    self.warp_transform_type_3d = TransformType3D::None;
                }
            },
        }
    }

    // Fractal FBm

    fn gen_fractal_fbm_2d(&self, x: Float, y: Float) -> f32 {
        let mut x = x;
        let mut y = y;

        let mut seed = self.seed;
        let mut sum = 0.;
        let mut amp = self.fractal_bounding;

        for _ in 0..self.octaves {
            let noise = self.gen_noise_single_2d(seed, x, y);

            seed += 1;
            sum += noise * amp;
            amp *= Self::lerp(1., (noise + 1.).min(2.) * 0.5, self.weighted_strength);

            x *= self.lacunarity as Float;
            y *= self.lacunarity as Float;
            amp *= self.gain;
        }

        sum
    }

    fn gen_fractal_fbm_3d(&self, x: Float, y: Float, z: Float) -> f32 {
        let mut x = x;
        let mut y = y;
        let mut z = z;

        let mut seed = self.seed;
        let mut sum = 0.;
        let mut amp = self.fractal_bounding;

        for _ in 0..self.octaves {
            let noise = self.gen_noise_single_3d(seed, x, y, z);

            seed += 1;
            sum += noise * amp;
            amp *= Self::lerp(1., (noise + 1.) * 0.5, self.weighted_strength);

            x *= self.lacunarity as Float;
            y *= self.lacunarity as Float;
            z *= self.lacunarity as Float;
            amp *= self.gain;
        }

        sum
    }

    // Fractal Ridged

    fn gen_fractal_ridged_2d(&self, x: Float, y: Float) -> f32 {
        let mut x = x;
        let mut y = y;

        let mut seed = self.seed;
        let mut sum = 0.;
        let mut amp = self.fractal_bounding;

        for _ in 0..self.octaves {
            let noise = FloatOps::abs(self.gen_noise_single_2d(seed, x, y));
            seed += 1;

            sum += (noise * -2. + 1.) * amp;
            amp *= Self::lerp(1., 1. - noise, self.weighted_strength);

            x *= self.lacunarity as Float;
            y *= self.lacunarity as Float;
            amp *= self.gain;
        }

        sum
    }

    fn gen_fractal_ridged_3d(&self, x: Float, y: Float, z: Float) -> f32 {
        let mut x = x;
        let mut y = y;
        let mut z = z;

        let mut seed = self.seed;
        let mut sum = 0.;
        let mut amp = self.fractal_bounding;

        for _ in 0..self.octaves {
            let noise = FloatOps::abs(self.gen_noise_single_3d(seed, x, y, z));
            seed += 1;

            sum += (noise * -2. + 1.) * amp;
            amp *= Self::lerp(1., 1. - noise, self.weighted_strength);

            x *= self.lacunarity as Float;
            y *= self.lacunarity as Float;
            z *= self.lacunarity as Float;
            amp *= self.gain;
        }

        sum
    }

    // Fractal PingPong

    fn gen_fractal_ping_pong_2d(&self, x: Float, y: Float) -> f32 {
        let mut x = x;
        let mut y = y;

        let mut seed = self.seed;
        let mut sum = 0.;
        let mut amp = self.fractal_bounding; // DEBUG: This is 0.400000006, should be 0.4

        for _ in 0..self.octaves {
            let noise = Self::ping_pong(
                (self.gen_noise_single_2d(seed, x, y) + 1.) * self.ping_pong_strength,
            );
            seed += 1;

            sum += (noise - 0.5) * 2. * amp;
            amp *= Self::lerp(1., noise, self.weighted_strength);

            x *= self.lacunarity as Float;
            y *= self.lacunarity as Float;
            amp *= self.gain;
        }

        sum
    }

    fn gen_fractal_ping_pong_3d(&self, x: Float, y: Float, z: Float) -> f32 {
        let mut x = x;
        let mut y = y;
        let mut z = z;

        let mut seed = self.seed;
        let mut sum = 0.;
        let mut amp = self.fractal_bounding;

        for _ in 0..self.octaves {
            let noise = Self::ping_pong(
                (self.gen_noise_single_3d(seed, x, y, z) + 1.) * self.ping_pong_strength,
            );
            seed += 1;

            sum += (noise - 0.5) * 2. * amp;
            amp *= Self::lerp(1., noise, self.weighted_strength);

            x *= self.lacunarity as Float;
            y *= self.lacunarity as Float;
            z *= self.lacunarity as Float;
            amp *= self.gain;
        }

        sum
    }

    // Simplex/OpenSimplex2 Noise

    fn single_simplex_2d(&self, seed: i32, x: Float, y: Float) -> f32 {
        // 2D OpenSimplex2 case uses the same algorithm as ordinary Simplex.

        let sqrt3 = 1.7320508075688772935274463415059;
        let g2 = (3. - sqrt3) / 6.;

        /*
         * --- Skew moved to TransformNoiseCoordinateXY method ---
         * let f2 = 0.5 * (sqrt3 - 1.);
         * let s = (x + y) * f2;
         * x += s; y += s;
         */

        let i = Self::fast_floor(x);
        let j = Self::fast_floor(y);
        #[allow(clippy::unnecessary_cast)]
        let xi = (x - i as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let yi = (y - j as Float) as f32;

        let t = (xi + yi) * g2;
        let x0 = xi - t;
        let y0 = yi - t;

        let i = i.wrapping_mul(PRIME_X);
        let j = j.wrapping_mul(PRIME_Y);

        let a = 0.5 - x0 * x0 - y0 * y0;
        let n0 = if a <= 0. {
            0.
        } else {
            (a * a) * (a * a) * Self::grad_coord_2d(seed, i, j, x0, y0)
        };

        let c = (2. * (1. - 2. * g2) * (1. / g2 - 2.)) * t
            + ((-2. * (1. - 2. * g2) * (1. - 2. * g2)) + a);
        let n2 = if c <= 0. {
            0.
        } else {
            let x2 = x0 + (2. * g2 - 1.);
            let y2 = y0 + (2. * g2 - 1.);

            (c * c)
                * (c * c)
                * Self::grad_coord_2d(
                    seed,
                    i.wrapping_add(PRIME_X),
                    j.wrapping_add(PRIME_Y),
                    x2,
                    y2,
                )
        };

        let n1 = if y0 > x0 {
            let x1 = x0 + g2;
            let y1 = y0 + (g2 - 1.);
            let b = 0.5 - x1 * x1 - y1 * y1;

            if b <= 0. {
                0.
            } else {
                (b * b) * (b * b) * Self::grad_coord_2d(seed, i, j.wrapping_add(PRIME_Y), x1, y1)
            }
        } else {
            let x1 = x0 + (g2 - 1.);
            let y1 = y0 + g2;
            let b = 0.5 - x1 * x1 - y1 * y1;

            if b <= 0. {
                0.
            } else {
                (b * b) * (b * b) * Self::grad_coord_2d(seed, i.wrapping_add(PRIME_X), j, x1, y1)
            }
        };

        (n0 + n1 + n2) * 99.83685446303647
    }

    fn single_open_simplex_2(&self, seed: i32, x: Float, y: Float, z: Float) -> f32 {
        // 3D OpenSimplex2 case uses two offset rotated cube grids.

        /*
         * --- Rotation moved to TransformNoiseCoordinateXYZ method ---
         * let r3 = 2. / 3.;
         * let r = (x + y + z) * r3; // Rotation, not skew
         * x = r - x; y = r - y; z = r - z;
         */

        let mut seed = seed;

        let i = Self::fast_round(x);
        let j = Self::fast_round(y);
        let k = Self::fast_round(z);
        #[allow(clippy::unnecessary_cast)]
        let mut x0 = (x - i as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let mut y0 = (y - j as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let mut z0 = (z - k as Float) as f32;

        let mut x_n_sign = (-1. - x0) as i32 | 1;
        let mut y_n_sign = (-1. - y0) as i32 | 1;
        let mut z_n_sign = (-1. - z0) as i32 | 1;

        let mut ax0 = x_n_sign as f32 * -x0;
        let mut ay0 = y_n_sign as f32 * -y0;
        let mut az0 = z_n_sign as f32 * -z0;

        let mut i = i.wrapping_mul(PRIME_X);
        let mut j = j.wrapping_mul(PRIME_Y);
        let mut k = k.wrapping_mul(PRIME_Z);

        let mut value = 0.;
        let mut a = (0.6 - x0 * x0) - (y0 * y0 + z0 * z0);

        let mut l = 0;
        loop {
            if a > 0. {
                value += (a * a) * (a * a) * Self::grad_coord_3d(seed, i, j, k, x0, y0, z0);
            }

            if ax0 >= ay0 && ax0 >= az0 {
                let b = a + ax0 + ax0;
                if b > 1. {
                    let b = b - 1.;
                    value += (b * b)
                        * (b * b)
                        * Self::grad_coord_3d(
                            seed,
                            i.wrapping_sub(x_n_sign.wrapping_mul(PRIME_X)),
                            j,
                            k,
                            x0 + x_n_sign as f32,
                            y0,
                            z0,
                        );
                }
            } else if ay0 > ax0 && ay0 >= az0 {
                let b = a + ay0 + ay0;
                if b > 1. {
                    let b = b - 1.;
                    value += (b * b)
                        * (b * b)
                        * Self::grad_coord_3d(
                            seed,
                            i,
                            j.wrapping_sub(y_n_sign.wrapping_mul(PRIME_Y)),
                            k,
                            x0,
                            y0 + y_n_sign as f32,
                            z0,
                        );
                }
            } else {
                let b = a + az0 + az0;
                if b > 1. {
                    let b = b - 1.;
                    value += (b * b)
                        * (b * b)
                        * Self::grad_coord_3d(
                            seed,
                            i,
                            j,
                            k.wrapping_sub(z_n_sign.wrapping_mul(PRIME_Z)),
                            x0,
                            y0,
                            z0 + z_n_sign as f32,
                        );
                }
            }

            if l == 1 {
                break;
            }

            ax0 = 0.5 - ax0;
            ay0 = 0.5 - ay0;
            az0 = 0.5 - az0;

            x0 = x_n_sign as f32 * ax0;
            y0 = y_n_sign as f32 * ay0;
            z0 = z_n_sign as f32 * az0;

            a = a + (0.75 - ax0) - (ay0 + az0);

            i = i.wrapping_add((x_n_sign >> 1) & PRIME_X);
            j = j.wrapping_add((y_n_sign >> 1) & PRIME_Y);
            k = k.wrapping_add((z_n_sign >> 1) & PRIME_Z);

            x_n_sign = -x_n_sign;
            y_n_sign = -y_n_sign;
            z_n_sign = -z_n_sign;

            seed = !seed;

            l += 1;
        }

        value * 32.69428253173828125
    }

    // OpenSimplex2S Noise

    fn single_open_simplex_2s_2d(&self, seed: i32, x: Float, y: Float) -> f32 {
        // 2D OpenSimplex2S case is a modified 2D simplex noise.

        let sqrt3 = 1.7320508075688772935274463415059;
        let g2 = (3. - sqrt3) / 6.;

        /*
         * --- Skew moved to TransformNoiseCoordinateXY method ---
         * let f2 = 0.5 * (sqrt3 - 1);
         * let s = (x + y) * f2;
         * x += s; y += s;
         */

        let i = Self::fast_floor(x);
        let j = Self::fast_floor(y);
        #[allow(clippy::unnecessary_cast)]
        let xi = (x - i as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let yi = (y - j as Float) as f32;

        let i = i.wrapping_mul(PRIME_X);
        let j = j.wrapping_mul(PRIME_Y);
        let i1 = i.wrapping_add(PRIME_X);
        let j1 = j.wrapping_add(PRIME_Y);

        let t = (xi + yi) * g2;
        let x0 = xi - t;
        let y0 = yi - t;

        let a0 = (2. / 3.) - x0 * x0 - y0 * y0;
        let mut value = (a0 * a0) * (a0 * a0) * Self::grad_coord_2d(seed, i, j, x0, y0);

        let a1 = (2. * (1. - 2. * g2) * (1. / g2 - 2.)) * t
            + ((-2. * (1. - 2. * g2) * (1. - 2. * g2)) + a0);
        let x1 = x0 - (1. - 2. * g2);
        let y1 = y0 - (1. - 2. * g2);
        value += (a1 * a1) * (a1 * a1) * Self::grad_coord_2d(seed, i1, j1, x1, y1);

        // Nested conditionals were faster than compact bit logic/arithmetic.
        let xmyi = xi - yi;
        if t > g2 {
            if xi + xmyi > 1. {
                let x2 = x0 + (3. * g2 - 2.);
                let y2 = y0 + (3. * g2 - 1.);
                let a2 = (2. / 3.) - x2 * x2 - y2 * y2;
                if a2 > 0. {
                    value += (a2 * a2)
                        * (a2 * a2)
                        * Self::grad_coord_2d(
                            seed,
                            i.wrapping_add(PRIME_X << 1),
                            j.wrapping_add(PRIME_Y),
                            x2,
                            y2,
                        )
                }
            } else {
                let x2 = x0 + g2;
                let y2 = y0 + (g2 - 1.);
                let a2 = (2. / 3.) - x2 * x2 - y2 * y2;
                if a2 > 0. {
                    value += (a2 * a2)
                        * (a2 * a2)
                        * Self::grad_coord_2d(seed, i, j.wrapping_add(PRIME_Y), x2, y2)
                }
            }

            if yi - xmyi > 1. {
                let x3 = x0 + (3. * g2 - 1.);
                let y3 = y0 + (3. * g2 - 2.);
                let a3 = (2. / 3.) - x3 * x3 - y3 * y3;
                if a3 > 0. {
                    value += (a3 * a3)
                        * (a3 * a3)
                        * Self::grad_coord_2d(
                            seed,
                            i.wrapping_add(PRIME_X),
                            j.wrapping_add(PRIME_Y << 1),
                            x3,
                            y3,
                        )
                }
            } else {
                let x3 = x0 + (g2 - 1.);
                let y3 = y0 + g2;
                let a3 = (2. / 3.) - x3 * x3 - y3 * y3;
                if a3 > 0. {
                    value += (a3 * a3)
                        * (a3 * a3)
                        * Self::grad_coord_2d(seed, i.wrapping_add(PRIME_X), j, x3, y3)
                }
            }
        } else {
            if xi + xmyi < 0. {
                let x2 = x0 + (1. - g2);
                let y2 = y0 - g2;
                let a2 = (2. / 3.) - x2 * x2 - y2 * y2;
                if a2 > 0. {
                    value += (a2 * a2)
                        * (a2 * a2)
                        * Self::grad_coord_2d(seed, i.wrapping_sub(PRIME_X), j, x2, y2)
                }
            } else {
                let x2 = x0 + (g2 - 1.);
                let y2 = y0 + g2;
                let a2 = (2. / 3.) - x2 * x2 - y2 * y2;
                if a2 > 0. {
                    value += (a2 * a2)
                        * (a2 * a2)
                        * Self::grad_coord_2d(seed, i.wrapping_add(PRIME_X), j, x2, y2)
                }
            }

            if yi < xmyi {
                let x2 = x0 - g2;
                let y2 = y0 - (g2 - 1.);
                let a2 = (2. / 3.) - x2 * x2 - y2 * y2;
                if a2 > 0. {
                    value += (a2 * a2)
                        * (a2 * a2)
                        * Self::grad_coord_2d(seed, i, j.wrapping_sub(PRIME_Y), x2, y2)
                }
            } else {
                let x2 = x0 + g2;
                let y2 = y0 + (g2 - 1.);
                let a2 = (2. / 3.) - x2 * x2 - y2 * y2;
                if a2 > 0. {
                    value += (a2 * a2)
                        * (a2 * a2)
                        * Self::grad_coord_2d(seed, i, j.wrapping_add(PRIME_Y), x2, y2)
                }
            }
        }

        value * 18.24196194486065
    }

    fn single_open_simplex_2s_3d(&self, seed: i32, x: Float, y: Float, z: Float) -> f32 {
        // 3D OpenSimplex2S case uses two offset rotated cube grids.

        /*
         * --- Rotation moved to TransformNoiseCoordinateXYZ method ---
         * let R3 = 2. / 3.;
         * let r = (x + y + z) * R3; // Rotation, not skew
         * x = r - x; y = r - y; z = r - z;
        	*/

        let i = Self::fast_floor(x);
        let j = Self::fast_floor(y);
        let k = Self::fast_floor(z);
        #[allow(clippy::unnecessary_cast)]
        let xi = (x - i as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let yi = (y - j as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let zi = (z - k as Float) as f32;

        let i = i.wrapping_mul(PRIME_X);
        let j = j.wrapping_mul(PRIME_Y);
        let k = k.wrapping_mul(PRIME_Z);
        let seed2 = seed + 1293373;

        let x_n_mask = (-0.5 - xi) as i32;
        let y_n_mask = (-0.5 - yi) as i32;
        let z_n_mask = (-0.5 - zi) as i32;

        let x0 = xi + x_n_mask as f32;
        let y0 = yi + y_n_mask as f32;
        let z0 = zi + z_n_mask as f32;
        let a0 = 0.75 - x0 * x0 - y0 * y0 - z0 * z0;
        let mut value = (a0 * a0)
            * (a0 * a0)
            * Self::grad_coord_3d(
                seed,
                i.wrapping_add(x_n_mask & PRIME_X),
                j.wrapping_add(y_n_mask & PRIME_Y),
                k.wrapping_add(z_n_mask & PRIME_Z),
                x0,
                y0,
                z0,
            );

        let x1 = xi - 0.5;
        let y1 = yi - 0.5;
        let z1 = zi - 0.5;
        let a1 = 0.75 - x1 * x1 - y1 * y1 - z1 * z1;
        value += (a1 * a1)
            * (a1 * a1)
            * Self::grad_coord_3d(
                seed2,
                i.wrapping_add(PRIME_X),
                j.wrapping_add(PRIME_Y),
                k.wrapping_add(PRIME_Z),
                x1,
                y1,
                z1,
            );

        let x_a_flip_mask_0 = ((x_n_mask | 1) << 1) as f32 * x1;
        let y_a_flip_mask_0 = ((y_n_mask | 1) << 1) as f32 * y1;
        let z_a_flip_mask_0 = ((z_n_mask | 1) << 1) as f32 * z1;
        let x_a_flip_mask_1 = (-2 - (x_n_mask << 2)) as f32 * x1 - 1.;
        let y_a_flip_mask_1 = (-2 - (y_n_mask << 2)) as f32 * y1 - 1.;
        let z_a_flip_mask_1 = (-2 - (z_n_mask << 2)) as f32 * z1 - 1.;

        let mut skip_5 = false;
        let a2 = x_a_flip_mask_0 + a0;
        if a2 > 0. {
            let x2 = x0 - (x_n_mask | 1) as f32;
            let y2 = y0;
            let z2 = z0;
            value += (a2 * a2)
                * (a2 * a2)
                * Self::grad_coord_3d(
                    seed,
                    i.wrapping_add(!x_n_mask & PRIME_X),
                    j.wrapping_add(y_n_mask & PRIME_Y),
                    k.wrapping_add(z_n_mask & PRIME_Z),
                    x2,
                    y2,
                    z2,
                );
        } else {
            let a3 = y_a_flip_mask_0 + z_a_flip_mask_0 + a0;
            if a3 > 0. {
                let x3 = x0;
                let y3 = y0 - (y_n_mask | 1) as f32;
                let z3 = z0 - (z_n_mask | 1) as f32;
                value += (a3 * a3)
                    * (a3 * a3)
                    * Self::grad_coord_3d(
                        seed,
                        i.wrapping_add(x_n_mask & PRIME_X),
                        j.wrapping_add(!y_n_mask & PRIME_Y),
                        k.wrapping_add(!z_n_mask & PRIME_Z),
                        x3,
                        y3,
                        z3,
                    );
            }

            let a4 = x_a_flip_mask_1 + a1;
            if a4 > 0. {
                let x4 = (x_n_mask | 1) as f32 + x1;
                let y4 = y1;
                let z4 = z1;
                value += (a4 * a4)
                    * (a4 * a4)
                    * Self::grad_coord_3d(
                        seed2,
                        i.wrapping_add(x_n_mask & PRIME_X_2),
                        j.wrapping_add(PRIME_Y),
                        k.wrapping_add(PRIME_Z),
                        x4,
                        y4,
                        z4,
                    );
                skip_5 = true;
            }
        }

        let mut skip_9 = false;
        let a6 = y_a_flip_mask_0 + a0;
        if a6 > 0. {
            let x6 = x0;
            let y6 = y0 - (y_n_mask | 1) as f32;
            let z6 = z0;
            value += (a6 * a6)
                * (a6 * a6)
                * Self::grad_coord_3d(
                    seed,
                    i.wrapping_add(x_n_mask & PRIME_X),
                    j.wrapping_add(!y_n_mask & PRIME_Y),
                    k.wrapping_add(z_n_mask & PRIME_Z),
                    x6,
                    y6,
                    z6,
                );
        } else {
            let a7 = x_a_flip_mask_0 + z_a_flip_mask_0 + a0;
            if a7 > 0. {
                let x7 = x0 - (x_n_mask | 1) as f32;
                let y7 = y0;
                let z7 = z0 - (z_n_mask | 1) as f32;
                value += (a7 * a7)
                    * (a7 * a7)
                    * Self::grad_coord_3d(
                        seed,
                        i.wrapping_add(!x_n_mask & PRIME_X),
                        j.wrapping_add(y_n_mask & PRIME_Y),
                        k.wrapping_add(!z_n_mask & PRIME_Z),
                        x7,
                        y7,
                        z7,
                    );
            }

            let a8 = y_a_flip_mask_1 + a1;
            if a8 > 0. {
                let x8 = x1;
                let y8 = (y_n_mask | 1) as f32 + y1;
                let z8 = z1;
                value += (a8 * a8)
                    * (a8 * a8)
                    * Self::grad_coord_3d(
                        seed2,
                        i.wrapping_add(PRIME_X),
                        j.wrapping_add(y_n_mask & (PRIME_Y << 1)),
                        k.wrapping_add(PRIME_Z),
                        x8,
                        y8,
                        z8,
                    );
                skip_9 = true;
            }
        }

        let mut skip_d = false;
        let a_a = z_a_flip_mask_0 + a0;
        if a_a > 0. {
            let x_a = x0;
            let y_a = y0;
            let z_a = z0 - (z_n_mask | 1) as f32;
            value += (a_a * a_a)
                * (a_a * a_a)
                * Self::grad_coord_3d(
                    seed,
                    i.wrapping_add(x_n_mask & PRIME_X),
                    j.wrapping_add(y_n_mask & PRIME_Y),
                    k.wrapping_add(!z_n_mask & PRIME_Z),
                    x_a,
                    y_a,
                    z_a,
                );
        } else {
            let a_b = x_a_flip_mask_0 + y_a_flip_mask_0 + a0;
            if a_b > 0. {
                let x_b = x0 - (x_n_mask | 1) as f32;
                let y_b = y0 - (y_n_mask | 1) as f32;
                let z_b = z0;
                value += (a_b * a_b)
                    * (a_b * a_b)
                    * Self::grad_coord_3d(
                        seed,
                        i.wrapping_add(!x_n_mask & PRIME_X),
                        j.wrapping_add(!y_n_mask & PRIME_Y),
                        k.wrapping_add(z_n_mask & PRIME_Z),
                        x_b,
                        y_b,
                        z_b,
                    );
            }

            let a_c = z_a_flip_mask_1 + a1;
            if a_c > 0. {
                let x_c = x1;
                let y_c = y1;
                let z_c = (z_n_mask | 1) as f32 + z1;
                value += (a_c * a_c)
                    * (a_c * a_c)
                    * Self::grad_coord_3d(
                        seed2,
                        i.wrapping_add(PRIME_X),
                        j.wrapping_add(PRIME_Y),
                        k.wrapping_add(z_n_mask & (PRIME_Z << 1)),
                        x_c,
                        y_c,
                        z_c,
                    );
                skip_d = true;
            }
        }

        if !skip_5 {
            let a5 = y_a_flip_mask_1 + z_a_flip_mask_1 + a1;
            if a5 > 0. {
                let x5 = x1;
                let y5 = (y_n_mask | 1) as f32 + y1;
                let z5 = (z_n_mask | 1) as f32 + z1;
                value += (a5 * a5)
                    * (a5 * a5)
                    * Self::grad_coord_3d(
                        seed2,
                        i.wrapping_add(PRIME_X),
                        j.wrapping_add(y_n_mask & (PRIME_Y << 1)),
                        k.wrapping_add(z_n_mask & (PRIME_Z << 1)),
                        x5,
                        y5,
                        z5,
                    );
            }
        }

        if !skip_9 {
            let a9 = x_a_flip_mask_1 + z_a_flip_mask_1 + a1;
            if a9 > 0. {
                let x9 = (x_n_mask | 1) as f32 + x1;
                let y9 = y1;
                let z9 = (z_n_mask | 1) as f32 + z1;
                value += (a9 * a9)
                    * (a9 * a9)
                    * Self::grad_coord_3d(
                        seed2,
                        i.wrapping_add(x_n_mask & PRIME_X_2),
                        j.wrapping_add(PRIME_Y),
                        k.wrapping_add(z_n_mask & (PRIME_Z << 1)),
                        x9,
                        y9,
                        z9,
                    );
            }
        }

        if !skip_d {
            let a_d = x_a_flip_mask_1 + y_a_flip_mask_1 + a1;
            if a_d > 0. {
                let x_d = (x_n_mask | 1) as f32 + x1;
                let y_d = (y_n_mask | 1) as f32 + y1;
                let z_d = z1;
                value += (a_d * a_d)
                    * (a_d * a_d)
                    * Self::grad_coord_3d(
                        seed2,
                        i.wrapping_add(x_n_mask & (PRIME_X << 1)),
                        j.wrapping_add(y_n_mask & (PRIME_Y << 1)),
                        k.wrapping_add(PRIME_Z),
                        x_d,
                        y_d,
                        z_d,
                    );
            }
        }

        value * 9.046026385208288
    }

    // Cellular Noise

    fn single_cellular_2d(&self, seed: i32, x: Float, y: Float) -> f32 {
        let xr = Self::fast_round(x);
        let yr = Self::fast_round(y);

        let mut distance0 = f32::MAX;
        let mut distance1 = f32::MAX;
        let mut closest_hash = 0;

        let cellular_jitter = 0.43701595 * self.cellular_jitter_modifier;

        let mut x_primed = (xr - 1).wrapping_mul(PRIME_X);
        let y_primed_base = (yr - 1).wrapping_mul(PRIME_Y);

        match self.cellular_distance_function {
            CellularDistanceFunction::Euclidean | CellularDistanceFunction::EuclideanSq => {
                for xi in xr - 1..=xr + 1 {
                    let mut y_primed = y_primed_base;

                    for yi in yr - 1..=yr + 1 {
                        let hash = Self::hash_2d(seed, x_primed, y_primed);
                        let idx = hash & (255 << 1);

                        #[allow(clippy::unnecessary_cast)]
                        let vec_x =
                            (xi as Float - x) as f32 + RAND_VECS_2D[idx as usize] * cellular_jitter;
                        #[allow(clippy::unnecessary_cast)]
                        let vec_y = (yi as Float - y) as f32
                            + RAND_VECS_2D[(idx | 1) as usize] * cellular_jitter;

                        let new_distance = vec_x * vec_x + vec_y * vec_y;

                        distance1 = distance1.min(new_distance).max(distance0);
                        if new_distance < distance0 {
                            distance0 = new_distance;
                            closest_hash = hash;
                        }
                        y_primed = y_primed.wrapping_add(PRIME_Y);
                    }
                    x_primed = x_primed.wrapping_add(PRIME_X);
                }
            }
            CellularDistanceFunction::Manhattan => {
                for xi in xr - 1..=xr + 1 {
                    let mut y_primed = y_primed_base;

                    for yi in yr - 1..=yr + 1 {
                        let hash = Self::hash_2d(seed, x_primed, y_primed);
                        let idx = hash & (255 << 1);

                        #[allow(clippy::unnecessary_cast)]
                        let vec_x =
                            (xi as Float - x) as f32 + RAND_VECS_2D[idx as usize] * cellular_jitter;
                        #[allow(clippy::unnecessary_cast)]
                        let vec_y = (yi as Float - y) as f32
                            + RAND_VECS_2D[(idx | 1) as usize] * cellular_jitter;

                        let new_distance = FloatOps::abs(vec_x) + FloatOps::abs(vec_y);

                        distance1 = distance1.min(new_distance).max(distance0);
                        if new_distance < distance0 {
                            distance0 = new_distance;
                            closest_hash = hash;
                        }
                        y_primed = y_primed.wrapping_add(PRIME_Y);
                    }
                    x_primed = x_primed.wrapping_add(PRIME_X);
                }
            }
            CellularDistanceFunction::Hybrid => {
                for xi in xr - 1..=xr + 1 {
                    let mut y_primed = y_primed_base;

                    for yi in yr - 1..=yr + 1 {
                        let hash = Self::hash_2d(seed, x_primed, y_primed);
                        let idx = hash & (255 << 1);

                        #[allow(clippy::unnecessary_cast)]
                        let vec_x =
                            (xi as Float - x) as f32 + RAND_VECS_2D[idx as usize] * cellular_jitter;
                        #[allow(clippy::unnecessary_cast)]
                        let vec_y = (yi as Float - y) as f32
                            + RAND_VECS_2D[(idx | 1) as usize] * cellular_jitter;

                        let new_distance = (FloatOps::abs(vec_x) + FloatOps::abs(vec_y))
                            + (vec_x * vec_x + vec_y * vec_y);

                        distance1 = distance1.min(new_distance).max(distance0);
                        if new_distance < distance0 {
                            distance0 = new_distance;
                            closest_hash = hash;
                        }
                        y_primed = y_primed.wrapping_add(PRIME_Y);
                    }
                    x_primed = x_primed.wrapping_add(PRIME_X);
                }
            }
        }

        if self.cellular_distance_function == CellularDistanceFunction::Euclidean
            && self.cellular_return_type >= CellularReturnType::Distance
        {
            distance0 = FloatOps::sqrt(distance0);

            if self.cellular_return_type >= CellularReturnType::Distance2 {
                distance1 = FloatOps::sqrt(distance1);
            }
        }

        match self.cellular_return_type {
            CellularReturnType::CellValue => closest_hash as f32 * (1. / 2147483648.),
            CellularReturnType::Distance => distance0 - 1.,
            CellularReturnType::Distance2 => distance1 - 1.,
            CellularReturnType::Distance2Add => (distance1 + distance0) * 0.5 - 1.,
            CellularReturnType::Distance2Sub => distance1 - distance0 - 1.,
            CellularReturnType::Distance2Mul => distance1 * distance0 * 0.5 - 1.,
            CellularReturnType::Distance2Div => distance0 / distance1 - 1.,
        }
    }

    fn single_cellular_3d(&self, seed: i32, x: Float, y: Float, z: Float) -> f32 {
        let xr = Self::fast_round(x);
        let yr = Self::fast_round(y);
        let zr = Self::fast_round(z);

        let mut distance0 = f32::MAX;
        let mut distance1 = f32::MAX;
        let mut closest_hash = 0;

        let cellular_jitter = 0.39614353 * self.cellular_jitter_modifier;

        let mut x_primed = (xr - 1).wrapping_mul(PRIME_X);
        let y_primed_base = (yr - 1).wrapping_mul(PRIME_Y);
        let z_primed_base = (zr - 1).wrapping_mul(PRIME_Z);

        match self.cellular_distance_function {
            CellularDistanceFunction::Euclidean | CellularDistanceFunction::EuclideanSq => {
                for xi in xr - 1..=xr + 1 {
                    let mut y_primed = y_primed_base;

                    for yi in yr - 1..=yr + 1 {
                        let mut z_primed = z_primed_base;

                        for zi in zr - 1..=zr + 1 {
                            let hash = Self::hash_3d(seed, x_primed, y_primed, z_primed);
                            let idx = hash & (255 << 2);

                            #[allow(clippy::unnecessary_cast)]
                            let vec_x = (xi as Float - x) as f32
                                + RAND_VECS_3D[idx as usize] * cellular_jitter;
                            #[allow(clippy::unnecessary_cast)]
                            let vec_y = (yi as Float - y) as f32
                                + RAND_VECS_3D[(idx | 1) as usize] * cellular_jitter;
                            #[allow(clippy::unnecessary_cast)]
                            let vec_z = (zi as Float - z) as f32
                                + RAND_VECS_3D[(idx | 2) as usize] * cellular_jitter;

                            let new_distance = vec_x * vec_x + vec_y * vec_y + vec_z * vec_z;

                            distance1 = distance1.min(new_distance).max(distance0);
                            if new_distance < distance0 {
                                distance0 = new_distance;
                                closest_hash = hash;
                            }
                            z_primed = z_primed.wrapping_add(PRIME_Z);
                        }
                        y_primed = y_primed.wrapping_add(PRIME_Y);
                    }
                    x_primed = x_primed.wrapping_add(PRIME_X);
                }
            }
            CellularDistanceFunction::Manhattan => {
                for xi in xr - 1..=xr + 1 {
                    let mut y_primed = y_primed_base;

                    for yi in yr - 1..=yr + 1 {
                        let mut z_primed = z_primed_base;

                        for zi in zr - 1..=zr + 1 {
                            let hash = Self::hash_3d(seed, x_primed, y_primed, z_primed);
                            let idx = hash & (255 << 2);

                            #[allow(clippy::unnecessary_cast)]
                            let vec_x = (xi as Float - x) as f32
                                + RAND_VECS_3D[idx as usize] * cellular_jitter;
                            #[allow(clippy::unnecessary_cast)]
                            let vec_y = (yi as Float - y) as f32
                                + RAND_VECS_3D[(idx | 1) as usize] * cellular_jitter;
                            #[allow(clippy::unnecessary_cast)]
                            let vec_z = (zi as Float - z) as f32
                                + RAND_VECS_3D[(idx | 2) as usize] * cellular_jitter;

                            let new_distance =
                                FloatOps::abs(vec_x) + FloatOps::abs(vec_y) + FloatOps::abs(vec_z);

                            distance1 = distance1.min(new_distance).max(distance0);
                            if new_distance < distance0 {
                                distance0 = new_distance;
                                closest_hash = hash;
                            }
                            z_primed = z_primed.wrapping_add(PRIME_Z);
                        }
                        y_primed = y_primed.wrapping_add(PRIME_Y);
                    }
                    x_primed = x_primed.wrapping_add(PRIME_X);
                }
            }
            CellularDistanceFunction::Hybrid => {
                for xi in xr - 1..=xr + 1 {
                    let mut y_primed = y_primed_base;

                    for yi in yr - 1..=yr + 1 {
                        let mut z_primed = z_primed_base;

                        for zi in zr - 1..=zr + 1 {
                            let hash = Self::hash_3d(seed, x_primed, y_primed, z_primed);
                            let idx = hash & (255 << 2);

                            #[allow(clippy::unnecessary_cast)]
                            let vec_x = (xi as Float - x) as f32
                                + RAND_VECS_3D[idx as usize] * cellular_jitter;
                            #[allow(clippy::unnecessary_cast)]
                            let vec_y = (yi as Float - y) as f32
                                + RAND_VECS_3D[(idx | 1) as usize] * cellular_jitter;
                            #[allow(clippy::unnecessary_cast)]
                            let vec_z = (zi as Float - z) as f32
                                + RAND_VECS_3D[(idx | 2) as usize] * cellular_jitter;

                            let new_distance = (FloatOps::abs(vec_x)
                                + FloatOps::abs(vec_y)
                                + FloatOps::abs(vec_z))
                                + (vec_x * vec_x + vec_y * vec_y + vec_z * vec_z);

                            distance1 = distance1.min(new_distance).max(distance0);
                            if new_distance < distance0 {
                                distance0 = new_distance;
                                closest_hash = hash;
                            }
                            z_primed = z_primed.wrapping_add(PRIME_Z);
                        }
                        y_primed = y_primed.wrapping_add(PRIME_Y);
                    }
                    x_primed = x_primed.wrapping_add(PRIME_X);
                }
            }
        }

        if self.cellular_distance_function == CellularDistanceFunction::Euclidean
            && self.cellular_return_type >= CellularReturnType::Distance
        {
            distance0 = FloatOps::sqrt(distance0);

            if self.cellular_return_type >= CellularReturnType::Distance2 {
                distance1 = FloatOps::sqrt(distance1);
            }
        }

        match self.cellular_return_type {
            CellularReturnType::CellValue => closest_hash as f32 * (1. / 2147483648.),
            CellularReturnType::Distance => distance0 - 1.,
            CellularReturnType::Distance2 => distance1 - 1.,
            CellularReturnType::Distance2Add => (distance1 + distance0) * 0.5 - 1.,
            CellularReturnType::Distance2Sub => distance1 - distance0 - 1.,
            CellularReturnType::Distance2Mul => distance1 * distance0 * 0.5 - 1.,
            CellularReturnType::Distance2Div => distance0 / distance1 - 1.,
        }
    }

    // Perlin Noise

    fn single_perlin_2d(&self, seed: i32, x: Float, y: Float) -> f32 {
        let x0 = Self::fast_floor(x);
        let y0 = Self::fast_floor(y);

        #[allow(clippy::unnecessary_cast)]
        let xd0 = (x - x0 as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let yd0 = (y - y0 as Float) as f32;
        let xd1 = xd0 - 1.;
        let yd1 = yd0 - 1.;

        let xs = Self::interp_quintic(xd0);
        let ys = Self::interp_quintic(yd0);

        let x0 = x0.wrapping_mul(PRIME_X);
        let y0 = y0.wrapping_mul(PRIME_Y);
        let x1 = x0.wrapping_add(PRIME_X);
        let y1 = y0.wrapping_add(PRIME_Y);

        let xf0 = Self::lerp(
            Self::grad_coord_2d(seed, x0, y0, xd0, yd0),
            Self::grad_coord_2d(seed, x1, y0, xd1, yd0),
            xs,
        );
        let xf1 = Self::lerp(
            Self::grad_coord_2d(seed, x0, y1, xd0, yd1),
            Self::grad_coord_2d(seed, x1, y1, xd1, yd1),
            xs,
        );

        Self::lerp(xf0, xf1, ys) * 1.4247691104677813
    }

    fn single_perlin_3d(&self, seed: i32, x: Float, y: Float, z: Float) -> f32 {
        let x0 = Self::fast_floor(x);
        let y0 = Self::fast_floor(y);
        let z0 = Self::fast_floor(z);

        #[allow(clippy::unnecessary_cast)]
        let xd0 = (x - x0 as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let yd0 = (y - y0 as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let zd0 = (z - z0 as Float) as f32;
        let xd1 = xd0 - 1.;
        let yd1 = yd0 - 1.;
        let zd1 = zd0 - 1.;

        let xs = Self::interp_quintic(xd0);
        let ys = Self::interp_quintic(yd0);
        let zs = Self::interp_quintic(zd0);

        let x0 = x0.wrapping_mul(PRIME_X);
        let y0 = y0.wrapping_mul(PRIME_Y);
        let z0 = z0.wrapping_mul(PRIME_Z);
        let x1 = x0.wrapping_add(PRIME_X);
        let y1 = y0.wrapping_add(PRIME_Y);
        let z1 = z0.wrapping_add(PRIME_Z);

        let xf00 = Self::lerp(
            Self::grad_coord_3d(seed, x0, y0, z0, xd0, yd0, zd0),
            Self::grad_coord_3d(seed, x1, y0, z0, xd1, yd0, zd0),
            xs,
        );
        let xf10 = Self::lerp(
            Self::grad_coord_3d(seed, x0, y1, z0, xd0, yd1, zd0),
            Self::grad_coord_3d(seed, x1, y1, z0, xd1, yd1, zd0),
            xs,
        );
        let xf01 = Self::lerp(
            Self::grad_coord_3d(seed, x0, y0, z1, xd0, yd0, zd1),
            Self::grad_coord_3d(seed, x1, y0, z1, xd1, yd0, zd1),
            xs,
        );
        let xf11 = Self::lerp(
            Self::grad_coord_3d(seed, x0, y1, z1, xd0, yd1, zd1),
            Self::grad_coord_3d(seed, x1, y1, z1, xd1, yd1, zd1),
            xs,
        );

        let yf0 = Self::lerp(xf00, xf10, ys);
        let yf1 = Self::lerp(xf01, xf11, ys);

        Self::lerp(yf0, yf1, zs) * 0.964921414852142333984375
    }

    // Value Cubic Noise

    fn single_value_cubic_2d(&self, seed: i32, x: Float, y: Float) -> f32 {
        let x1 = Self::fast_floor(x);
        let y1 = Self::fast_floor(y);

        #[allow(clippy::unnecessary_cast)]
        let xs = (x - x1 as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let ys = (y - y1 as Float) as f32;

        let x1 = x1.wrapping_mul(PRIME_X);
        let y1 = y1.wrapping_mul(PRIME_Y);
        let x0 = x1.wrapping_sub(PRIME_X);
        let y0 = y1.wrapping_sub(PRIME_Y);
        let x2 = x1.wrapping_add(PRIME_X);
        let y2 = y1.wrapping_add(PRIME_Y);
        let x3 = x1.wrapping_add(PRIME_X_2);
        let y3 = y1.wrapping_add(PRIME_Y_2);

        Self::cubic_lerp(
            Self::cubic_lerp(
                Self::val_coord_2d(seed, x0, y0),
                Self::val_coord_2d(seed, x1, y0),
                Self::val_coord_2d(seed, x2, y0),
                Self::val_coord_2d(seed, x3, y0),
                xs,
            ),
            Self::cubic_lerp(
                Self::val_coord_2d(seed, x0, y1),
                Self::val_coord_2d(seed, x1, y1),
                Self::val_coord_2d(seed, x2, y1),
                Self::val_coord_2d(seed, x3, y1),
                xs,
            ),
            Self::cubic_lerp(
                Self::val_coord_2d(seed, x0, y2),
                Self::val_coord_2d(seed, x1, y2),
                Self::val_coord_2d(seed, x2, y2),
                Self::val_coord_2d(seed, x3, y2),
                xs,
            ),
            Self::cubic_lerp(
                Self::val_coord_2d(seed, x0, y3),
                Self::val_coord_2d(seed, x1, y3),
                Self::val_coord_2d(seed, x2, y3),
                Self::val_coord_2d(seed, x3, y3),
                xs,
            ),
            ys,
        ) * (1. / (1.5 * 1.5))
    }

    fn single_value_cubic_3d(&self, seed: i32, x: Float, y: Float, z: Float) -> f32 {
        let x1 = Self::fast_floor(x);
        let y1 = Self::fast_floor(y);
        let z1 = Self::fast_floor(z);

        #[allow(clippy::unnecessary_cast)]
        let xs = (x - x1 as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let ys = (y - y1 as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let zs = (z - z1 as Float) as f32;

        let x1 = x1.wrapping_mul(PRIME_X);
        let y1 = y1.wrapping_mul(PRIME_Y);
        let z1 = z1.wrapping_mul(PRIME_Z);

        let x0 = x1.wrapping_sub(PRIME_X);
        let y0 = y1.wrapping_sub(PRIME_Y);
        let z0 = z1.wrapping_sub(PRIME_Z);
        let x2 = x1.wrapping_add(PRIME_X);
        let y2 = y1.wrapping_add(PRIME_Y);
        let z2 = z1.wrapping_add(PRIME_Z);
        let x3 = x1.wrapping_add(PRIME_X_2);
        let y3 = y1.wrapping_add(PRIME_Y_2);
        let z3 = z1.wrapping_add(PRIME_Z_2);

        Self::cubic_lerp(
            Self::cubic_lerp(
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y0, z0),
                    Self::val_coord_3d(seed, x1, y0, z0),
                    Self::val_coord_3d(seed, x2, y0, z0),
                    Self::val_coord_3d(seed, x3, y0, z0),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y1, z0),
                    Self::val_coord_3d(seed, x1, y1, z0),
                    Self::val_coord_3d(seed, x2, y1, z0),
                    Self::val_coord_3d(seed, x3, y1, z0),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y2, z0),
                    Self::val_coord_3d(seed, x1, y2, z0),
                    Self::val_coord_3d(seed, x2, y2, z0),
                    Self::val_coord_3d(seed, x3, y2, z0),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y3, z0),
                    Self::val_coord_3d(seed, x1, y3, z0),
                    Self::val_coord_3d(seed, x2, y3, z0),
                    Self::val_coord_3d(seed, x3, y3, z0),
                    xs,
                ),
                ys,
            ),
            Self::cubic_lerp(
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y0, z1),
                    Self::val_coord_3d(seed, x1, y0, z1),
                    Self::val_coord_3d(seed, x2, y0, z1),
                    Self::val_coord_3d(seed, x3, y0, z1),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y1, z1),
                    Self::val_coord_3d(seed, x1, y1, z1),
                    Self::val_coord_3d(seed, x2, y1, z1),
                    Self::val_coord_3d(seed, x3, y1, z1),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y2, z1),
                    Self::val_coord_3d(seed, x1, y2, z1),
                    Self::val_coord_3d(seed, x2, y2, z1),
                    Self::val_coord_3d(seed, x3, y2, z1),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y3, z1),
                    Self::val_coord_3d(seed, x1, y3, z1),
                    Self::val_coord_3d(seed, x2, y3, z1),
                    Self::val_coord_3d(seed, x3, y3, z1),
                    xs,
                ),
                ys,
            ),
            Self::cubic_lerp(
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y0, z2),
                    Self::val_coord_3d(seed, x1, y0, z2),
                    Self::val_coord_3d(seed, x2, y0, z2),
                    Self::val_coord_3d(seed, x3, y0, z2),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y1, z2),
                    Self::val_coord_3d(seed, x1, y1, z2),
                    Self::val_coord_3d(seed, x2, y1, z2),
                    Self::val_coord_3d(seed, x3, y1, z2),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y2, z2),
                    Self::val_coord_3d(seed, x1, y2, z2),
                    Self::val_coord_3d(seed, x2, y2, z2),
                    Self::val_coord_3d(seed, x3, y2, z2),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y3, z2),
                    Self::val_coord_3d(seed, x1, y3, z2),
                    Self::val_coord_3d(seed, x2, y3, z2),
                    Self::val_coord_3d(seed, x3, y3, z2),
                    xs,
                ),
                ys,
            ),
            Self::cubic_lerp(
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y0, z3),
                    Self::val_coord_3d(seed, x1, y0, z3),
                    Self::val_coord_3d(seed, x2, y0, z3),
                    Self::val_coord_3d(seed, x3, y0, z3),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y1, z3),
                    Self::val_coord_3d(seed, x1, y1, z3),
                    Self::val_coord_3d(seed, x2, y1, z3),
                    Self::val_coord_3d(seed, x3, y1, z3),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y2, z3),
                    Self::val_coord_3d(seed, x1, y2, z3),
                    Self::val_coord_3d(seed, x2, y2, z3),
                    Self::val_coord_3d(seed, x3, y2, z3),
                    xs,
                ),
                Self::cubic_lerp(
                    Self::val_coord_3d(seed, x0, y3, z3),
                    Self::val_coord_3d(seed, x1, y3, z3),
                    Self::val_coord_3d(seed, x2, y3, z3),
                    Self::val_coord_3d(seed, x3, y3, z3),
                    xs,
                ),
                ys,
            ),
            zs,
        ) * (1. / (1.5 * 1.5 * 1.5))
    }

    // Value Noise

    fn single_value_2d(&self, seed: i32, x: Float, y: Float) -> f32 {
        let x0 = Self::fast_floor(x);
        let y0 = Self::fast_floor(y);

        #[allow(clippy::unnecessary_cast)]
        let xs = Self::interp_hermite((x - x0 as Float) as f32);
        #[allow(clippy::unnecessary_cast)]
        let ys = Self::interp_hermite((y - y0 as Float) as f32);

        let x0 = x0.wrapping_mul(PRIME_X);
        let y0 = y0.wrapping_mul(PRIME_Y);
        let x1 = x0.wrapping_add(PRIME_X);
        let y1 = y0.wrapping_add(PRIME_Y);

        let xf0 = Self::lerp(
            Self::val_coord_2d(seed, x0, y0),
            Self::val_coord_2d(seed, x1, y0),
            xs,
        );
        let xf1 = Self::lerp(
            Self::val_coord_2d(seed, x0, y1),
            Self::val_coord_2d(seed, x1, y1),
            xs,
        );

        Self::lerp(xf0, xf1, ys)
    }

    fn single_value_3d(&self, seed: i32, x: Float, y: Float, z: Float) -> f32 {
        let x0 = Self::fast_floor(x);
        let y0 = Self::fast_floor(y);
        let z0 = Self::fast_floor(z);

        #[allow(clippy::unnecessary_cast)]
        let xs = Self::interp_hermite((x - x0 as Float) as f32);
        #[allow(clippy::unnecessary_cast)]
        let ys = Self::interp_hermite((y - y0 as Float) as f32);
        #[allow(clippy::unnecessary_cast)]
        let zs = Self::interp_hermite((z - z0 as Float) as f32);

        let x0 = x0.wrapping_mul(PRIME_X);
        let y0 = y0.wrapping_mul(PRIME_Y);
        let z0 = z0.wrapping_mul(PRIME_Z);
        let x1 = x0.wrapping_add(PRIME_X);
        let y1 = y0.wrapping_add(PRIME_Y);
        let z1 = z0.wrapping_add(PRIME_Z);

        let xf00 = Self::lerp(
            Self::val_coord_3d(seed, x0, y0, z0),
            Self::val_coord_3d(seed, x1, y0, z0),
            xs,
        );
        let xf10 = Self::lerp(
            Self::val_coord_3d(seed, x0, y1, z0),
            Self::val_coord_3d(seed, x1, y1, z0),
            xs,
        );
        let xf01 = Self::lerp(
            Self::val_coord_3d(seed, x0, y0, z1),
            Self::val_coord_3d(seed, x1, y0, z1),
            xs,
        );
        let xf11 = Self::lerp(
            Self::val_coord_3d(seed, x0, y1, z1),
            Self::val_coord_3d(seed, x1, y1, z1),
            xs,
        );

        let yf0 = Self::lerp(xf00, xf10, ys);
        let yf1 = Self::lerp(xf01, xf11, ys);

        Self::lerp(yf0, yf1, zs)
    }

    // Domain Warp

    #[allow(clippy::too_many_arguments)]
    fn do_single_domain_warp_2d(
        &self,
        seed: i32,
        amp: f32,
        freq: f32,
        x: Float,
        y: Float,
        xr: Float,
        yr: Float,
    ) -> (Float, Float) {
        match self.domain_warp_type {
            DomainWarpType::OpenSimplex2 => self.single_domain_warp_simplex_gradient_2d(
                seed,
                amp * 38.283687591552734375,
                freq,
                x,
                y,
                xr,
                yr,
                false,
            ),
            DomainWarpType::OpenSimplex2Reduced => self.single_domain_warp_simplex_gradient_2d(
                seed,
                amp * 16.,
                freq,
                x,
                y,
                xr,
                yr,
                true,
            ),
            DomainWarpType::BasicGrid => {
                self.single_domain_warp_basic_grid_2d(seed, amp, freq, x, y, xr, yr)
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn do_single_domain_warp_3d(
        &self,
        seed: i32,
        amp: f32,
        freq: f32,
        x: Float,
        y: Float,
        z: Float,
        xr: Float,
        yr: Float,
        zr: Float,
    ) -> (Float, Float, Float) {
        match self.domain_warp_type {
            DomainWarpType::OpenSimplex2 => self.single_domain_warp_open_simplex_2_gradient(
                seed,
                amp * 32.69428253173828125,
                freq,
                x,
                y,
                z,
                xr,
                yr,
                zr,
                false,
            ),
            DomainWarpType::OpenSimplex2Reduced => self.single_domain_warp_open_simplex_2_gradient(
                seed,
                amp * 7.71604938271605,
                freq,
                x,
                y,
                z,
                xr,
                yr,
                zr,
                true,
            ),
            DomainWarpType::BasicGrid => {
                self.single_domain_warp_basic_grid_3d(seed, amp, freq, x, y, z, xr, yr, zr)
            }
        }
    }

    // Domain Warp Single Wrapper

    fn domain_warp_single_2d(&self, x: Float, y: Float) -> (Float, Float) {
        let seed = self.seed;
        let amp = self.domain_warp_amp * self.fractal_bounding;
        let freq = self.frequency;

        let (xs, ys) = self.transform_domain_warp_coordinate_2d(x, y);

        self.do_single_domain_warp_2d(seed, amp, freq, xs, ys, x, y)
    }

    fn domain_warp_single_3d(&self, x: Float, y: Float, z: Float) -> (Float, Float, Float) {
        let seed = self.seed;
        let amp = self.domain_warp_amp * self.fractal_bounding;
        let freq = self.frequency;

        let (xs, ys, zs) = self.transform_domain_warp_coordinate_3d(x, y, z);

        self.do_single_domain_warp_3d(seed, amp, freq, xs, ys, zs, x, y, z)
    }

    // Domain Warp Fractal Progressive

    fn domain_warp_fractal_progressive_2d(&self, x: Float, y: Float) -> (Float, Float) {
        let mut x = x;
        let mut y = y;

        let mut seed = self.seed;
        let mut amp = self.domain_warp_amp * self.fractal_bounding;
        let mut freq = self.frequency;

        for _ in 0..self.octaves {
            let (xs, ys) = self.transform_domain_warp_coordinate_2d(x, y);

            (x, y) = self.do_single_domain_warp_2d(seed, amp, freq, xs, ys, x, y);

            seed += 1;
            amp *= self.gain;
            freq *= self.lacunarity;
        }

        (x, y)
    }

    fn domain_warp_fractal_progressive_3d(
        &self,
        x: Float,
        y: Float,
        z: Float,
    ) -> (Float, Float, Float) {
        let mut x = x;
        let mut y = y;
        let mut z = z;

        let mut seed = self.seed;
        let mut amp = self.domain_warp_amp * self.fractal_bounding;
        let mut freq = self.frequency;

        for _ in 0..self.octaves {
            let (xs, ys, zs) = self.transform_domain_warp_coordinate_3d(x, y, z);

            (x, y, z) = self.do_single_domain_warp_3d(seed, amp, freq, xs, ys, zs, x, y, z);

            seed += 1;
            amp *= self.gain;
            freq *= self.lacunarity;
        }

        (x, y, z)
    }

    // Domain Warp Fractal Independant

    fn domain_warp_fractal_independent_2d(&self, x: Float, y: Float) -> (Float, Float) {
        let mut x = x;
        let mut y = y;

        let (xs, ys) = self.transform_domain_warp_coordinate_2d(x, y);

        let mut seed = self.seed;
        let mut amp = self.domain_warp_amp * self.fractal_bounding;
        let mut freq = self.frequency;

        for _ in 0..self.octaves {
            (x, y) = self.do_single_domain_warp_2d(seed, amp, freq, xs, ys, x, y);

            seed += 1;
            amp *= self.gain;
            freq *= self.lacunarity;
        }

        (x, y)
    }

    fn domain_warp_fractal_independent_3d(
        &self,
        x: Float,
        y: Float,
        z: Float,
    ) -> (Float, Float, Float) {
        let mut x = x;
        let mut y = y;
        let mut z = z;

        let (xs, ys, zs) = self.transform_domain_warp_coordinate_3d(x, y, z);

        let mut seed = self.seed;
        let mut amp = self.domain_warp_amp * self.fractal_bounding;
        let mut freq = self.frequency;

        for _ in 0..self.octaves {
            (x, y, z) = self.do_single_domain_warp_3d(seed, amp, freq, xs, ys, zs, x, y, z);

            seed += 1;
            amp *= self.gain;
            freq *= self.lacunarity;
        }

        (x, y, z)
    }

    // Domain Warp Basic Grid

    #[allow(clippy::too_many_arguments)]
    fn single_domain_warp_basic_grid_2d(
        &self,
        seed: i32,
        warp_amp: f32,
        frequency: f32,
        x: Float,
        y: Float,
        xr: Float,
        yr: Float,
    ) -> (Float, Float) {
        let xf = x * frequency as Float;
        let yf = y * frequency as Float;

        let x0 = Self::fast_floor(xf);
        let y0 = Self::fast_floor(yf);

        #[allow(clippy::unnecessary_cast)]
        let xs = Self::interp_hermite((xf - x0 as Float) as f32);
        #[allow(clippy::unnecessary_cast)]
        let ys = Self::interp_hermite((yf - y0 as Float) as f32);

        let x0 = x0.wrapping_mul(PRIME_X);
        let y0 = y0.wrapping_mul(PRIME_Y);
        let x1 = x0.wrapping_add(PRIME_X);
        let y1 = y0.wrapping_add(PRIME_Y);

        let hash0 = Self::hash_2d(seed, x0, y0) & (255 << 1);
        let hash1 = Self::hash_2d(seed, x1, y0) & (255 << 1);

        let lx0x = Self::lerp(
            RAND_VECS_2D[hash0 as usize],
            RAND_VECS_2D[hash1 as usize],
            xs,
        );
        let ly0x = Self::lerp(
            RAND_VECS_2D[(hash0 | 1) as usize],
            RAND_VECS_2D[(hash1 | 1) as usize],
            xs,
        );

        let hash0 = Self::hash_2d(seed, x0, y1) & (255 << 1);
        let hash1 = Self::hash_2d(seed, x1, y1) & (255 << 1);

        let lx1x = Self::lerp(
            RAND_VECS_2D[hash0 as usize],
            RAND_VECS_2D[hash1 as usize],
            xs,
        );
        let ly1x = Self::lerp(
            RAND_VECS_2D[(hash0 | 1) as usize],
            RAND_VECS_2D[(hash1 | 1) as usize],
            xs,
        );

        let xr = xr + (Self::lerp(lx0x, lx1x, ys) * warp_amp) as Float;
        let yr = yr + (Self::lerp(ly0x, ly1x, ys) * warp_amp) as Float;

        (xr, yr)
    }

    #[allow(clippy::too_many_arguments)]
    fn single_domain_warp_basic_grid_3d(
        &self,
        seed: i32,
        warp_amp: f32,
        frequency: f32,
        x: Float,
        y: Float,
        z: Float,
        xr: Float,
        yr: Float,
        zr: Float,
    ) -> (Float, Float, Float) {
        let xf = x * frequency as Float;
        let yf = y * frequency as Float;
        let zf = z * frequency as Float;

        let x0 = Self::fast_floor(xf);
        let y0 = Self::fast_floor(yf);
        let z0 = Self::fast_floor(zf);

        #[allow(clippy::unnecessary_cast)]
        let xs = Self::interp_hermite((xf - x0 as Float) as f32);
        #[allow(clippy::unnecessary_cast)]
        let ys = Self::interp_hermite((yf - y0 as Float) as f32);
        #[allow(clippy::unnecessary_cast)]
        let zs = Self::interp_hermite((zf - z0 as Float) as f32);

        let x0 = x0.wrapping_mul(PRIME_X);
        let y0 = y0.wrapping_mul(PRIME_Y);
        let z0 = z0.wrapping_mul(PRIME_Z);
        let x1 = x0.wrapping_add(PRIME_X);
        let y1 = y0.wrapping_add(PRIME_Y);
        let z1 = z0.wrapping_add(PRIME_Z);

        let hash0 = Self::hash_3d(seed, x0, y0, z0) & (255 << 2);
        let hash1 = Self::hash_3d(seed, x1, y0, z0) & (255 << 2);

        let lx0x = Self::lerp(
            RAND_VECS_3D[hash0 as usize],
            RAND_VECS_3D[hash1 as usize],
            xs,
        );
        let ly0x = Self::lerp(
            RAND_VECS_3D[(hash0 | 1) as usize],
            RAND_VECS_3D[(hash1 | 1) as usize],
            xs,
        );
        let lz0x = Self::lerp(
            RAND_VECS_3D[(hash0 | 2) as usize],
            RAND_VECS_3D[(hash1 | 2) as usize],
            xs,
        );

        let hash0 = Self::hash_3d(seed, x0, y1, z0) & (255 << 2);
        let hash1 = Self::hash_3d(seed, x1, y1, z0) & (255 << 2);

        let lx1x = Self::lerp(
            RAND_VECS_3D[hash0 as usize],
            RAND_VECS_3D[hash1 as usize],
            xs,
        );
        let ly1x = Self::lerp(
            RAND_VECS_3D[(hash0 | 1) as usize],
            RAND_VECS_3D[(hash1 | 1) as usize],
            xs,
        );
        let lz1x = Self::lerp(
            RAND_VECS_3D[(hash0 | 2) as usize],
            RAND_VECS_3D[(hash1 | 2) as usize],
            xs,
        );

        let lx0y = Self::lerp(lx0x, lx1x, ys);
        let ly0y = Self::lerp(ly0x, ly1x, ys);
        let lz0y = Self::lerp(lz0x, lz1x, ys);

        let hash0 = Self::hash_3d(seed, x0, y0, z1) & (255 << 2);
        let hash1 = Self::hash_3d(seed, x1, y0, z1) & (255 << 2);

        let lx0x = Self::lerp(
            RAND_VECS_3D[hash0 as usize],
            RAND_VECS_3D[hash1 as usize],
            xs,
        );
        let ly0x = Self::lerp(
            RAND_VECS_3D[(hash0 | 1) as usize],
            RAND_VECS_3D[(hash1 | 1) as usize],
            xs,
        );
        let lz0x = Self::lerp(
            RAND_VECS_3D[(hash0 | 2) as usize],
            RAND_VECS_3D[(hash1 | 2) as usize],
            xs,
        );

        let hash0 = Self::hash_3d(seed, x0, y1, z1) & (255 << 2);
        let hash1 = Self::hash_3d(seed, x1, y1, z1) & (255 << 2);

        let lx1x = Self::lerp(
            RAND_VECS_3D[hash0 as usize],
            RAND_VECS_3D[hash1 as usize],
            xs,
        );
        let ly1x = Self::lerp(
            RAND_VECS_3D[(hash0 | 1) as usize],
            RAND_VECS_3D[(hash1 | 1) as usize],
            xs,
        );
        let lz1x = Self::lerp(
            RAND_VECS_3D[(hash0 | 2) as usize],
            RAND_VECS_3D[(hash1 | 2) as usize],
            xs,
        );

        let xr = xr + (Self::lerp(lx0y, Self::lerp(lx0x, lx1x, ys), zs) * warp_amp) as Float;
        let yr = yr + (Self::lerp(ly0y, Self::lerp(ly0x, ly1x, ys), zs) * warp_amp) as Float;
        let zr = zr + (Self::lerp(lz0y, Self::lerp(lz0x, lz1x, ys), zs) * warp_amp) as Float;

        (xr, yr, zr)
    }

    // Domain Warp Simplex/OpenSimplex2
    #[allow(clippy::too_many_arguments)]
    fn single_domain_warp_simplex_gradient_2d(
        &self,
        seed: i32,
        warp_amp: f32,
        frequency: f32,
        x: Float,
        y: Float,
        xr: Float,
        yr: Float,
        out_frad_only: bool,
    ) -> (Float, Float) {
        const SQRT3: f32 = 1.7320508075688772935274463415059;
        const G2: f32 = (3. - SQRT3) / 6.;

        let x = x * frequency as Float;
        let y = y * frequency as Float;

        /*
         * --- Skew moved to TransformNoiseCoordinateXY method ---
         * let f2 = 0.5 * (sqrt3 - 1);
         * let s = (x + y) * f2;
         * x += s; y += s;
         */

        let i = Self::fast_floor(x);
        let j = Self::fast_floor(y);
        #[allow(clippy::unnecessary_cast)]
        let xi = (x - i as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let yi = (y - j as Float) as f32;

        let t = (xi + yi) * G2;
        let x0 = xi - t;
        let y0 = yi - t;

        let i = i.wrapping_mul(PRIME_X);
        let j = j.wrapping_mul(PRIME_Y);

        let mut vx = 0.;
        let mut vy = 0.;

        let a = 0.5 - x0 * x0 - y0 * y0;
        if a > 0. {
            let aaaa = (a * a) * (a * a);
            let (xo, yo) = if out_frad_only {
                Self::grad_coord_out_2d(seed, i, j)
            } else {
                Self::grad_coord_dual_2d(seed, i, j, x0, y0)
            };
            vx += aaaa * xo;
            vy += aaaa * yo;
        }

        let c = (2. * (1. - 2. * G2) * (1. / G2 - 2.)) * t
            + ((-2. * (1. - 2. * G2) * (1. - 2. * G2)) + a);
        if c > 0. {
            let x2 = x0 + (2. * G2 - 1.);
            let y2 = y0 + (2. * G2 - 1.);
            let cccc = (c * c) * (c * c);
            let (xo, yo) = if out_frad_only {
                Self::grad_coord_out_2d(seed, i.wrapping_add(PRIME_X), j.wrapping_add(PRIME_Y))
            } else {
                Self::grad_coord_dual_2d(
                    seed,
                    i.wrapping_add(PRIME_X),
                    j.wrapping_add(PRIME_Y),
                    x2,
                    y2,
                )
            };
            vx += cccc * xo;
            vy += cccc * yo;
        }

        if y0 > x0 {
            let x1 = x0 + G2;
            let y1 = y0 + (G2 - 1.);
            let b = 0.5 - x1 * x1 - y1 * y1;
            if b > 0. {
                let bbbb = (b * b) * (b * b);
                let (xo, yo) = if out_frad_only {
                    Self::grad_coord_out_2d(seed, i, j.wrapping_add(PRIME_Y))
                } else {
                    Self::grad_coord_dual_2d(seed, i, j.wrapping_add(PRIME_Y), x1, y1)
                };
                vx += bbbb * xo;
                vy += bbbb * yo;
            }
        } else {
            let x1 = x0 + (G2 - 1.);
            let y1 = y0 + G2;
            let b = 0.5 - x1 * x1 - y1 * y1;
            if b > 0. {
                let bbbb = (b * b) * (b * b);
                let (xo, yo) = if out_frad_only {
                    Self::grad_coord_out_2d(seed, i.wrapping_add(PRIME_X), j)
                } else {
                    Self::grad_coord_dual_2d(seed, i.wrapping_add(PRIME_X), j, x1, y1)
                };
                vx += bbbb * xo;
                vy += bbbb * yo;
            }
        }

        let xr = xr + (vx * warp_amp) as Float;
        let yr = yr + (vy * warp_amp) as Float;

        (xr, yr)
    }

    #[allow(clippy::too_many_arguments)]
    fn single_domain_warp_open_simplex_2_gradient(
        &self,
        seed: i32,
        warp_amp: f32,
        frequency: f32,
        x: Float,
        y: Float,
        z: Float,
        xr: Float,
        yr: Float,
        zr: Float,
        out_grad_only: bool,
    ) -> (Float, Float, Float) {
        let mut seed = seed;

        let x = x * frequency as Float;
        let y = y * frequency as Float;
        let z = z * frequency as Float;

        /*
         * --- Rotation moved to TransformDomainWarpCoordinate method ---
         * let r3 = 2. / 3.;
         * let r = (x + y + z) * r3; // Rotation, not skew
         * x = r - x; y = r - y; z = r - z;
         */

        let i = Self::fast_round(x);
        let j = Self::fast_round(y);
        let k = Self::fast_round(z);
        #[allow(clippy::unnecessary_cast)]
        let mut x0 = (x - i as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let mut y0 = (y - j as Float) as f32;
        #[allow(clippy::unnecessary_cast)]
        let mut z0 = (z - k as Float) as f32;

        let mut x_n_sign = (-x0 - 1.) as i32 | 1;
        let mut y_n_sign = (-y0 - 1.) as i32 | 1;
        let mut z_n_sign = (-z0 - 1.) as i32 | 1;

        let mut ax0 = x_n_sign as f32 * -x0;
        let mut ay0 = y_n_sign as f32 * -y0;
        let mut az0 = z_n_sign as f32 * -z0;

        let mut i = i.wrapping_mul(PRIME_X);
        let mut j = j.wrapping_mul(PRIME_Y);
        let mut k = k.wrapping_mul(PRIME_Z);

        let mut vx = 0.;
        let mut vy = 0.;
        let mut vz = 0.;

        let mut a = (0.6 - x0 * x0) - (y0 * y0 + z0 * z0);
        let mut l = 0;
        loop {
            if a > 0. {
                let aaaa = (a * a) * (a * a);
                let (xo, yo, zo) = if out_grad_only {
                    Self::grad_coord_out_3d(seed, i, j, k)
                } else {
                    Self::grad_coord_dual_3d(seed, i, j, k, x0, y0, z0)
                };
                vx += aaaa * xo;
                vy += aaaa * yo;
                vz += aaaa * zo;
            }

            let mut b = a;
            let mut i1 = i;
            let mut j1 = j;
            let mut k1 = k;
            let mut x1 = x0;
            let mut y1 = y0;
            let mut z1 = z0;

            if ax0 >= ay0 && ax0 >= az0 {
                x1 += x_n_sign as f32;
                b = b + ax0 + ax0;
                i1 = i1.wrapping_sub(x_n_sign.wrapping_mul(PRIME_X));
            } else if ay0 > ax0 && ay0 >= az0 {
                y1 += y_n_sign as f32;
                b = b + ay0 + ay0;
                j1 = j1.wrapping_sub(y_n_sign.wrapping_mul(PRIME_Y));
            } else {
                z1 += z_n_sign as f32;
                b = b + az0 + az0;
                k1 = k1.wrapping_sub(z_n_sign.wrapping_mul(PRIME_Z));
            }

            if b > 1. {
                b -= 1.;
                let bbbb = (b * b) * (b * b);
                let (xo, yo, zo) = if out_grad_only {
                    Self::grad_coord_out_3d(seed, i1, j1, k1)
                } else {
                    Self::grad_coord_dual_3d(seed, i1, j1, k1, x1, y1, z1)
                };
                vx += bbbb * xo;
                vy += bbbb * yo;
                vz += bbbb * zo;
            }

            if l == 1 {
                break;
            }

            ax0 = 0.5 - ax0;
            ay0 = 0.5 - ay0;
            az0 = 0.5 - az0;

            x0 = x_n_sign as f32 * ax0;
            y0 = y_n_sign as f32 * ay0;
            z0 = z_n_sign as f32 * az0;

            a += (0.75 - ax0) - (ay0 + az0);

            i = i.wrapping_add((x_n_sign >> 1) & PRIME_X);
            j = j.wrapping_add((y_n_sign >> 1) & PRIME_Y);
            k = k.wrapping_add((z_n_sign >> 1) & PRIME_Z);

            x_n_sign = -x_n_sign;
            y_n_sign = -y_n_sign;
            z_n_sign = -z_n_sign;

            seed += 1293373;

            l += 1;
        }

        let xr = xr + (vx * warp_amp) as Float;
        let yr = yr + (vy * warp_amp) as Float;
        let zr = zr + (vz * warp_amp) as Float;

        (xr, yr, zr)
    }
}
