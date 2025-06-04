// src/math/types/mod.rs
pub mod bounds;
pub mod point;
pub mod vector;

pub use bounds::*;
pub use point::*;
pub use vector::*;

// Re-export häufig verwendete externe Typen
pub use bevy::math::{Vec2, Vec3};
pub use spade::Point2;

// Einheitliche Typen für das gesamte Modul
pub type Point2D = Vec2;
pub type Point3D = Vec3;
pub type SpadePoint = Point2<f64>;
