// src/math/geometry/mod.rs

// Deklaration der Haupt-Geometriemodule
pub mod metalballs;
pub mod polygon;
pub mod sphere; // Jetzt deklariert

// Re-Exporte für einen schnellen Zugriff auf die Kern-Geometrietypen,
// falls man nicht das gesamte `math::prelude` importieren möchte.

// Polygon-Exporte
pub use self::polygon::{
    Polygon,
    PolygonBuilder,
    PolygonOrientation,
    PolygonProperties,
    PolygonValidator,
    QuickValidation,
    ShapeGenerators,
    ValidationError,
    ValidationLevel,
    ValidationReport,
    ValidationSuggestion,
    ValidationWarning,
    operations::{
        boolean::{
            BooleanAnalysis, BooleanOperation, BooleanOperations, BooleanStats, PolygonBoolean,
        },
        clipping::{
            ClippingAlgorithm, ClippingAnalysis, ClippingOperations, ClippingStats, PolygonClipper,
        },
        convex_hull::{ConvexHullAlgorithm, ConvexHullComputer, ConvexityDefect, ConvexityUtils},
        triangulation::{
            PolygonTriangle, PolygonTriangulator, TriangulationAlgorithm, TriangulationUtils,
        },
    },
    // Transformationen und Operationen für Polygon
    transformations::{
        affine::{
            AffineTransform, AffineTransformable, TransformBuilder as AffineTransformBuilder,
        },
        distortion::{
            PolygonDistortion, PolygonDistortionEffects, PolygonDistortionFn,
            PolygonDistortionType, SineAxis,
        },
        projection::{
            Complex as ProjectionComplex, ConformalType, Polygon2DProjector, Polygon3DProjector,
            PolygonProjection2DType, ProjectionBatch,
            ProjectionEffects as PolygonProjectionEffects,
        },
    },
};

// Sphere-Exporte
pub use self::sphere::{
    CoordinateConverter, GeographicCoordinates, SamplingPatterns, SphereProjectionType,
    SphereProjectionUtils, SphereProjector, SphereRotation, SphereSampler, SphereSamplingMethod,
    SphereUVMapper, SphericalCoordinates,
};

// Metalballs-Exporte
pub use self::metalballs::{
    GridConfig, MetaballField, Metaballs, MetaballsBuilder,
    influence::{
        CombinationMode, CompositeInfluence, FieldInfluence, InfluenceUtils, LineInfluence,
        PointInfluence,
    },
    marching_squares::{Contour, MarchingSquares, MarchingSquaresIterator},
    objects::{
        // Konkrete Influence-Objekte
        Metaball,
        MetaballFalloff,
        PolygonFalloff,
        PolygonInfluence,
        PolygonInfluenceBuilder,
    },
};
