// src/math/mod.rs

// ... (bestehende Deklarationen für error, utils, types, geometry) ...
pub mod algorithms; // Neu strukturiert
pub mod error;
pub mod geometry;
pub mod point_distribution; // Neu strukturiert
pub mod probability;
pub mod types;
pub mod utils; // Bleibt vorerst unberührt

pub use bevy::math::{IVec2, IVec3, Mat2, Mat3, Mat4, Quat, UVec2, UVec3, Vec2, Vec3};
pub use error::{MathError, MathResult};
pub use types::{
    Bool2, Bool3, Bounds2D, Bounds3D, SpadePoint, Vec2Builder, Vector2DExt, Vector3DExt,
    bevy_to_spade_points, spade_to_bevy_points,
};

pub mod prelude {
    // Fehler und Ergebnis
    pub use super::error::{MathError, MathResult};

    // Utilities
    pub use super::utils::{
        angles, comparison, constants, easing, numerical, random, simple_geometry,
    };

    // Kern-Typen und Vektor-Erweiterungen
    pub use super::types::{
        Bool2, Bool3, Bounds2D, Bounds3D, SpadePoint, Vec2Builder, Vector2DExt, Vector3DExt,
        bevy_to_spade_points, spade_to_bevy_points,
    };
    pub use bevy::math::{IVec2, IVec3, Mat2, Mat3, Mat4, Quat, UVec2, UVec3, Vec2, Vec3};

    // Geometrie-Modul Exporte
    pub use super::geometry::{
        // Metalballs
        metalballs::{
            CombinationMode,
            CompositeInfluence,
            Contour,
            FieldInfluence,
            GridConfig,
            InfluenceUtils,
            LineInfluence,
            MarchingSquares,
            MarchingSquaresIterator,
            MetaballFalloff,
            MetaballField,
            MetaballSource,
            Metaballs,
            MetaballsBuilder,
            PointInfluence,
            PolygonFalloff as MetalballPolygonFalloff, // Alias um Kollision zu vermeiden
            PolygonInfluence,
            PolygonInfluenceBuilder,
        },
        // Polygon
        polygon::{
            AffineTransform,
            AffineTransformBuilder,
            AffineTransformable,
            ConformalType,
            ConvexHullComputer, // Behalten, da nicht explizit als verschoben markiert durch Fehler
            Polygon,
            Polygon2DProjector,
            Polygon3DProjector,
            PolygonBuilder,
            PolygonDistortion,
            PolygonDistortionEffects,
            PolygonDistortionFn,
            PolygonDistortionType,
            PolygonOrientation,
            PolygonProjection2DType,
            PolygonProjectionEffects,
            PolygonProperties,
            PolygonTriangle, // Bleibt hier laut Compiler-Vorschlag (crate::math::geometry::PolygonTriangle)
            PolygonValidator,
            ProjectionBatch,
            ProjectionComplex,
            QuickValidation,
            ShapeGenerators, // ShapeGenerators hier, wenn es Polygon-spezifisch ist
            SineAxis,
            ValidationError,
            ValidationLevel,
            ValidationReport,
            ValidationSuggestion,
            ValidationWarning,
        },
        // Sphere
        sphere::{
            CoordinateConverter, GeographicCoordinates, SamplingPatterns, SphereProjectionType,
            SphereProjectionUtils, SphereProjector, SphereRotation, SphereSampler,
            SphereSamplingMethod, SphereUVMapper, SphericalCoordinates,
        },
    };

    // Algorithmen-Modul Exporte
    pub use super::algorithms::{
        // Von geometry::polygon hierher verschobene Elemente
        BooleanOpType,
        ClippingAlgorithm,
        ConvexHullAlgorithm,
        PolygonBooleanOps,
        PolygonClipper,
        PolygonTriangulator,
        // Smoothing - korrigierte Pfade
        SmoothingAlgorithmParameters, // Direkt unter algorithms laut Compiler-Vorschlag
        TriangulationAlgorithm,

        // Point Relaxation - korrigierte Pfade
        point_relaxation::lloyd::LloydRelaxationStats,
        point_relaxation::{
            BoundaryHandling as LloydBoundaryHandling, LloydConfig, LloydRelaxation,
        },
        smoothing::{
            BezierSmoother, CatmullRomSmoother, ChaikinSmoother, Smoothing, SmoothingParameters,
        },
    };

    // Point Distribution Exporte
    pub use super::point_distribution::voronoi::{
        MergeStrategy as VoronoiMergeStrategy, // Alias falls generische Namen existieren
        PolygonMerger as VoronoiPolygonMerger,
        VoronoiBuilder,
        VoronoiCell,
        VoronoiConfig, // LloydConfig ist in VoronoiConfig enthalten oder wird intern erstellt
        VoronoiDiagram,
        VoronoiEdge,
        VoronoiExtractor,
    };

    pub use super::probability::{SeedPlugin, SeedResource};
}
