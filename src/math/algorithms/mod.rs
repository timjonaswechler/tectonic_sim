// src/math/algorithms/mod.rs

//! # Mathematical Algorithms Module
//!
//! This module consolidates various computational geometry algorithms
//! used within the simulation. It is organized into sub-modules,
//! each focusing on a specific category of algorithms.
//!
//! ## Available Algorithm Categories:
//!
//! - `boolean`: For performing boolean operations (union, intersection, difference, XOR) on polygons.
//! - `clipping`: For clipping polygons against other polygons or rectangular bounds.
//! - `convex_hull`: For computing the convex hull of a set of points.
//! - `smoothing`: For smoothing polylines and polygon boundaries using various techniques.
//! - `triangulation`: For triangulating polygons, a common prerequisite for other geometric operations.
//! - `point_relaxation`: For algorithms like Lloyd's relaxation to distribute points more evenly.
//!
//! Key components from these sub-modules are re-exported for easier access.

// Declaration of the different algorithm categories
pub mod boolean;
pub mod clipping;
pub mod convex_hull;
pub mod metaballs;
pub mod smoothing;
pub mod triangulation; // For Boolean operations on polygon point lists

// Placeholder for future algorithm modules
pub mod point_relaxation;
// pub mod pathfinding;
// pub mod optimization;

// Re-exports for direct access to algorithm modules
// or specific algorithms and their configurations.

// Smoothing Algorithms and Traits
pub use self::smoothing::{
    AlgorithmParameters as SmoothingAlgorithmParameters, // Specific parameters for smoothers
    BezierSmoother,
    CatmullRomSmoother,
    ChaikinSmoother,
    Smoothing,           // Main trait
    SmoothingParameters, // General parameters (if still used)
};

// Convex Hull
pub use self::convex_hull::{
    ConvexHullAlgorithm, // Enum for selecting the hull algorithm
    // ConvexityUtils and ConvexityDefect were removed as they are more related to polygon analysis
    ConvexHullComputer, // Struct to perform convex hull computation
};

// Triangulation
pub use self::triangulation::{
    PolygonTriangulator,    // Struct to perform triangulation
    Triangle,               // The triangle structure used by triangulators
    TriangulationAlgorithm, // Enum for selecting the triangulation algorithm
    TriangulationUtils,     // Utility functions related to triangulation
};

// Clipping
pub use self::clipping::{
    ClippingAlgorithm,  // Enum for selecting the clipping algorithm
    ClippingOperations, // Utility functions for common clipping tasks
    // ClippingStats, ClippingAnalysis, // If needed and implemented
    PolygonClipper, // Struct to perform polygon clipping
};

// Boolean Operations
pub use self::boolean::{
    BooleanOpType, // Enum for the type of boolean operation
    // BooleanOperations (utility struct) was removed as PolygonBooleanOps handles this
    // BooleanAnalysis, BooleanStats, // If needed and implemented
    PolygonBooleanOps, // Renamed from PolygonBoolean to distinguish from Polygon structure
};

// Point Relaxation (e.g., Lloyd's Algorithm)
pub use self::point_relaxation::{
    BoundaryHandling, // Enum for how to handle points near boundaries during relaxation
    // LloydStatistics, // Is defined in lloyd.rs and re-exported here
    LloydConfig,     // Configuration for Lloyd's relaxation
    LloydRelaxation, // Struct to perform Lloyd's relaxation
};
