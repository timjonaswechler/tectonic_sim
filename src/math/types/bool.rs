// Hilfsstrukturen für komponentenweise boolesche Ergebnisse
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bool2 {
    pub x: bool,
    pub y: bool,
}
impl Bool2 {
    pub fn new(x: bool, y: bool) -> Self {
        Self { x, y }
    }
    pub fn all(self) -> bool {
        self.x && self.y
    }
    pub fn any(self) -> bool {
        self.x || self.y
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bool3 {
    pub x: bool,
    pub y: bool,
    pub z: bool,
}
impl Bool3 {
    pub fn new(x: bool, y: bool, z: bool) -> Self {
        Self { x, y, z }
    }
    pub fn all(self) -> bool {
        self.x && self.y && self.z
    }
    pub fn any(self) -> bool {
        self.x || self.y || self.z
    }
}
