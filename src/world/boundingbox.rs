use std::fmt;

#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct BoundingBox {
    pub x0: usize,
    pub y0: usize,
    pub z0: usize,
    pub x1: usize,
    pub y1: usize,
    pub z1: usize,
}

impl fmt::Display for BoundingBox {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Box: ({}, {}, {}) - ({}, {}, {})",
            self.x0, self.y0, self.z0, self.x1, self.y1, self.z1
        )
    }
}

impl BoundingBox {
    pub fn new(x0: usize, y0: usize, z0: usize, x1: usize, y1: usize, z1: usize) -> BoundingBox {
        BoundingBox {
            x0,
            y0,
            z0,
            x1,
            y1,
            z1,
        }
    }

    /// Returns true if two bounding boxes intersect or touch.
    pub fn intersects(&self, rhs: &BoundingBox) -> bool {
        !(self.x0 > rhs.x1
            || rhs.x0 > self.x1
            || self.y0 > rhs.y1
            || rhs.y0 > self.y1
            || self.z0 > rhs.z1
            || rhs.z0 > self.z0)
    }
}
