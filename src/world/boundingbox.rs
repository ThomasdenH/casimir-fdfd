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

#[derive(Debug, Fail)]
pub enum BoundingBoxError {
    #[fail(display = "expansion outside domain: {} expanded by {}", bbox, distance)]
    ExpansionOutsideDomain {
        bbox: BoundingBox,
        distance: usize
    }
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

    pub fn inside(&self, rhs: &BoundingBox) -> bool {
        self.x0 >= rhs.x0
            && self.y0 >= rhs.y0
            && self.z0 >= rhs.z0
            && self.x1 < rhs.x1
            && self.y1 < rhs.y1
            && self.z1 < rhs.z1
    }

    pub fn expanded(&self, distance: usize) -> Result<BoundingBox, BoundingBoxError> {
        if distance > self.x0 || distance > self.y0 || distance > self.z0 {
            Err(BoundingBoxError::ExpansionOutsideDomain {
                bbox: *self,
                distance
            })
        } else {
            Ok(BoundingBox::new(
                self.x0 - distance,
                self.y0 - distance,
                self.z0 - distance,
                self.x1 + distance,
                self.y1 + distance,
                self.z1 + distance
            ))
        }
    }
}
