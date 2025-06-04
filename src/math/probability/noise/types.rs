#[derive(Copy, Clone, Debug, PartialEq)]
pub enum NoiseType {
    OpenSimplex2,
    OpenSimplex2S,
    Cellular,
    Perlin,
    ValueCubic,
    Value,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum RotationType3D {
    None,
    ImproveXYPlanes,
    ImproveXZPlanes,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum FractalType {
    None,
    FBm,
    Ridged,
    PingPong,
    DomainWarpProgressive,
    DomainWarpIndependent,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum CellularDistanceFunction {
    Euclidean,
    EuclideanSq,
    Manhattan,
    Hybrid,
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub enum CellularReturnType {
    CellValue = 0,
    Distance = 1,
    Distance2 = 2,
    Distance2Add = 3,
    Distance2Sub = 4,
    Distance2Mul = 5,
    Distance2Div = 6,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum DomainWarpType {
    OpenSimplex2,
    OpenSimplex2Reduced,
    BasicGrid,
}

// Interne Enum, die vielleicht nicht public sein muss,
// aber wenn sie in generator.rs gebraucht wird, dann `pub(crate)` oder `pub`
#[derive(Copy, Clone, Debug, PartialEq)]
pub(crate) enum TransformType3D {
    None,
    ImproveXYPlanes,
    ImproveXZPlanes,
    DefaultOpenSimplex2,
}
