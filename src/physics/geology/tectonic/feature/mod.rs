pub mod generic;
// // Topological features.
// TopologicalClosedPlateBoundary
// TopologicalSlabBoundary
// TopologicalNetwork
// UnclassifiedTopologicalFeature

// // Reconstruction features.
// TotalReconstructionSequence
// AbsoluteReferenceFrame

// // Artificial features.
// ClosedPlateBoundary
// ClosedContinentalBoundary
// InferredPaleoBoundary
// OldPlatesGridMark
// MeshNode
pub mod flowline;
pub use flowline::*;
// MotionPath
// PolygonCentroidPoint
// DisplacementPoint
// PoliticalBoundary
// SmallCircle

// // Rock units.
// BasicRockUnit
// RockUnit_carbonate
// RockUnit_siliciclastic
// RockUnit_evaporite
// RockUnit_organic
// RockUnit_chemical
// RockUnit_plutonic
// RockUnit_volcanic
// RockUnit_metamorphic
// RockUnit_indeterminate_igneous
// FossilCollection_small
// FossilCollection_medium
// FossilCollection_large

// // Abstract Geological Plane & Contact features.
// GeologicalPlane
// FoldPlane
// Fault
// TerraneBoundary
// Unconformity
// UnknownContact

// // Tectonic sections.
pub mod continental_rift;
pub use continental_rift::*;
pub mod mid_ocean_ridge;
pub use mid_ocean_ridge::*;
// SubductionZone
// OrogenicBelt
// Transform
// FractureZone
// PassiveContinentalBoundary

// // Fields.
// Bathymetry
// Topography
// Gravimetry
// Magnetics
// GlobalElevation
// OceanicAge
// CrustalThickness
// DynamicTopography
// MantleDensity
// HeatFlow
// SedimentThickness
// Roughness
// SpreadingRate
// SpreadingAsymmetry
// Stress

// // Tangible features.
// Isochron
// MagneticAnomalyIdentification
// MagneticAnomalyShipTrack
// FractureZoneIdentification
// Suture
pub mod island_arc;
pub use island_arc::*;
// HotSpot
// HotSpotTrail
// Seamount
// SlabEdge
// Volcano
// Pluton
// Ophiolite
// NavdatSampleMafic
// NavdatSampleIntermediate
// NavdatSampleFelsicLow
// NavdatSampleFelsicHigh
// AseismicRidge
// Coastline
pub mod craton;
pub use craton::*;
// LargeIgneousProvince
// Basin
pub mod continental_crust;
pub use continental_crust::*;
pub mod oceanic_crust;
pub use oceanic_crust::*;
// TransitionalCrust
// ContinentalFragment
// GeologicalLineation
// PseudoFault
// VirtualGeomagneticPole
// UnclassifiedFeature

// // Rasters.
// Raster

// // 3D scalar fields.
// ScalarField3D

// // In einem Bevy System:
// use bevy::prelude::*;
// use crate::physics::geology::tectonic::feature::{
//     common_components::{PlateId, AppearsAt, Vertices, Center, Area},
//     craton::Craton,
//     continental_rift::{ContinentalRift, ActiveContinentalRift},
// };

// fn setup_features(mut commands: Commands) {
//     // Erstelle einen Craton
//     commands.spawn((
//         Craton, // Marker
//         Name::new("Craton_Alpha"), // Bevy's Name Komponente
//         PlateId(1),
//         AppearsAt(1500.0), // Vor 1500 Millionen Jahren
//         Vertices(vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0), Vec3::new(5.0, 10.0, 0.0)]),
//         Center(Vec3::new(5.0, 3.33, 0.0)), // Kann auch später berechnet werden
//         Area(50.0), // Kann auch später berechnet werden
//     ));

//     // Erstelle einen aktiven kontinentalen Rift
//     commands.spawn((
//         ContinentalRift, // Basis-Marker
//         ActiveContinentalRift, // Zustands-Marker
//         Name::new("Great_Rift_Valley"),
//         PlateId(2), // Gehört zu Platte 2 (oder spaltet Platte X)
//         AppearsAt(30.0),
//         Vertices(vec![Vec3::new(100.0, 50.0, 0.0), Vec3::new(120.0, 55.0, 0.0), Vec3::new(140.0, 45.0, 0.0)]),
//         // Center und Area könnten für einen Rift weniger relevant sein oder anders interpretiert werden
//     ));

//     // Erstelle eine Flowline
//     // Hier könntest du common_components::FlowlinePath oder einfach nur Vertices verwenden
//     use crate::physics::geology::tectonic::feature::flowline::Flowline;
//     use crate::physics::geology::tectonic::feature::common_components::FlowlinePath; // Beispiel

//     commands.spawn((
//         Flowline,
//         Name::new("Flowline_R1_to_CratonA_Vtx3"),
//         PlateId(1), // Gehört zu der Platte, auf der es "fließt"
//         AppearsAt(25.0),
//         FlowlinePath { // Verwendung der spezifischeren Komponente
//             start_rift_entity: None, // Könnte Entity ID eines Rift-Features sein
//             end_crust_vertex_index: None, // Könnte Index eines Vertex sein
//             path_vertices: vec![Vec3::new(110.0, 52.0, 0.0), Vec3::new(80.0, 20.0, 0.0)],
//         }
//         // Alternativ nur:
//         // Vertices(vec![Vec3::new(110.0, 52.0, 0.0), Vec3::new(80.0, 20.0, 0.0)])
//     ));
// }
