use thiserror::Error;

#[derive(Error, Debug)]
pub enum VoronoiBuildError {
    #[error(
        "Insufficient points at step '{step}': expected at least {expected}, got {actual}. Context: {context}"
    )]
    InsufficientPointsContext {
        step: String,
        expected: usize,
        actual: usize,
        context: String, // Zusätzlicher Text, z.B. Bound-Größen
    },
    #[error(
        "Triangulation failed during step '{step}': {reason}. Points involved: {point_count}, Triangulation vertices: {tri_vertices}"
    )]
    TriangulationFailedContext {
        step: String,
        reason: String,
        point_count: usize,
        tri_vertices: usize,
    },
    #[error("Geometric failure during Voronoi step '{step}': {operation}. Additional info: {info}")]
    GeometricFailureContext {
        step: String,
        operation: String,
        info: String,
    },
}
