/// Zweidimensionales Raster für Skalarwerte.
/// Speicherung zeilenweise (row-major).
#[derive(Debug, Clone)]
pub struct MetaballField {
    data: Vec<f32>,
    pub width: usize,
    pub height: usize,
    pub cell_size: f32,
}

impl MetaballField {
    pub fn new(width: usize, height: usize, cell: f32) -> Self {
        Self {
            data: vec![0.0; width * height],
            width,
            height,
            cell_size: cell,
        }
    }

    #[inline]
    fn idx(&self, x: usize, y: usize) -> usize {
        y * self.width + x
    }

    pub fn get(&self, x: usize, y: usize) -> f32 {
        self.data[self.idx(x, y)]
    }

    pub fn set(&mut self, x: usize, y: usize, v: f32) {
        let i = self.idx(x, y);
        self.data[i] = v;
    }
}
