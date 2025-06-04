/// Gemeinsame Basiskonfiguration für Grid-basierte Berechnungen.
#[derive(Debug, Clone)]
pub struct GridConfig {
    pub width: usize,
    pub height: usize,
    pub cell_size: f32,
    pub threshold: f32,
}

impl Default for GridConfig {
    fn default() -> Self {
        Self {
            width: 100,
            height: 100,
            cell_size: 8.0,
            threshold: 1.0,
        }
    }
}
