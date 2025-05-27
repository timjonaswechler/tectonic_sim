use bevy::prelude::*;

// Hier kommt die Definition von SimulationParameters hin.
// Die Definition ist dieselbe wie im vorherigen Vorschlag,
// nur der Ort der Datei ist jetzt korrigiert.

#[derive(Resource, Debug)]
pub struct SimulationParameters {
    // --- Zeitsteuerung ---
    pub advanced_simulation_time_ma: f64,
    pub display_time_ma: f64,
    pub time_step_ma: f64,
    pub paused: bool,
    pub execute_single_step_request: bool,
    pub execute_step_backward_request: bool,

    // --- Globale Planeten-/Physik-Parameter ---
    pub planet_radius_km: f32,
    pub oceanic_crust_density_kg_m3: f32,
    pub continental_crust_density_kg_m3: f32,
    pub mantle_density_kg_m3: f32,
    pub mantle_viscosity_pas: f64,
    pub mantle_convection_strength_factor: f32,

    // --- Parameter für Plattengenerierung und -verhalten ---
    pub target_plate_count: u32,
    pub min_initial_plate_surface_area_ratio: f32,
    pub max_initial_plate_surface_area_ratio: f32,
    pub initial_continental_crust_probability: f32,
    pub strong_convergence_threshold_m_per_year: f32,
    pub rock_fracture_stress_threshold_pa: f64,
    pub mountain_uplift_factor: f32,
    pub seafloor_spreading_rate_km2_per_year: f32,

    // --- Visualisierung & Debug ---
    pub geology_vel_scale: f32,
    pub vertex_vel_scale: f32,
    pub convergence_threshold_factor: f32,
    pub shear_threshold_factor: f32,
    pub show_plate_boundaries: bool,
    pub show_plate_velocity_arrows: bool,
}

impl Default for SimulationParameters {
    fn default() -> Self {
        Self {
            // Zeit
            advanced_simulation_time_ma: 0.0,
            display_time_ma: 0.0,
            time_step_ma: 0.1,
            paused: true,
            execute_single_step_request: false,
            execute_step_backward_request: false,

            // Planet/Physik
            planet_radius_km: 6371.0,
            oceanic_crust_density_kg_m3: 2900.0,
            continental_crust_density_kg_m3: 2700.0,
            mantle_density_kg_m3: 3300.0,
            mantle_viscosity_pas: 1e21,
            mantle_convection_strength_factor: 1.0,

            // Platten
            target_plate_count: 12,
            min_initial_plate_surface_area_ratio: 0.01,
            max_initial_plate_surface_area_ratio: 0.20,
            initial_continental_crust_probability: 0.3,
            strong_convergence_threshold_m_per_year: 0.02,
            rock_fracture_stress_threshold_pa: 100e6,
            mountain_uplift_factor: 0.1,
            seafloor_spreading_rate_km2_per_year: 2.5,

            // Vis/Debug
            geology_vel_scale: 1.0,
            vertex_vel_scale: 1.0,
            convergence_threshold_factor: 0.7,
            shear_threshold_factor: 0.7,
            show_plate_boundaries: true,
            show_plate_velocity_arrows: true,
        }
    }
}
