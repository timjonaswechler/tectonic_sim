// src/math/geometry/sphere/mod.rs

// Deklaration der Untermodule für Kugel-spezifische Funktionalität
pub mod coordinates;
pub mod projection;
pub mod rotation;
pub mod sampling;

// Re-Exporte für den einfachen Zugriff auf die wichtigsten Kugel-Elemente
pub use self::coordinates::{CoordinateConverter, GeographicCoordinates, SphericalCoordinates};
pub use self::projection::{
    SphereProjectionType,
    SphereProjectionUtils, // Hilfsfunktionen für Projektionen
    SphereProjector,
    SphereUVMapper, // UV-Mapping Utilities
};
pub use self::rotation::{
    // BevyQuat wird direkt aus math::prelude oder bevy::math importiert
    SphereRotation, // Hilfsfunktionen für Rotationen mit BevyQuat
};
pub use self::sampling::{
    SamplingPatterns, // Spezielle Punktmuster
    // SamplingUtils, // Wurde entfernt, Funktionalität ggf. woanders
    SphereSampler,
    SphereSamplingMethod,
};
