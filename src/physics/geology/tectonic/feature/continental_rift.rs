use bevy::prelude::*;

// Marker component
#[derive(Component, Debug, Clone)]
pub struct ContinentalRift;

// Zusätzliche Marker für den Zustand des Rifts
#[derive(Component, Debug, Clone)]
pub struct ActiveContinentalRift;

#[derive(Component, Debug, Clone)]
pub struct FailedContinentalRift;
