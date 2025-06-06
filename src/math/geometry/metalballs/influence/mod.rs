// src/math/geometry/metalballs/influence/mod.rs

// Diese Datei existierte vorher nicht explizit im Plan, aber sie ist nötig, um
// die influence-Typen zu bündeln.
// Dein alter `math/geometry/metalballs/influence/mod.rs` war umfassender.
// Wir bauen ihn hier schrittweise neu auf.

pub mod composite;
pub mod line;
pub mod metalball; // Metaball ist jetzt unter objects/metalball.rs
pub mod point;
pub mod polygon; // PolygonInfluence ist jetzt unter objects/polygon.rs
pub mod traits; // Der FieldInfluence Trait

pub use self::composite::*;
pub use self::line::*;
pub use self::metalball::*;
pub use self::point::*;
pub use self::polygon::*;
pub use self::traits::*; // FieldInfluence, InfluenceUtils
