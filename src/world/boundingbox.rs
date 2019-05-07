use std::fmt;

/// A box aligned with the grid.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct BoundingBox {
    pub x0: usize,
    pub y0: usize,
    pub z0: usize,
    pub x1: usize,
    pub y1: usize,
    pub z1: usize,
}

/// This error is returned when expanding the `BoundingBox` below the 0 boundary.
#[derive(Debug, Fail)]
#[fail(
    display = "expansion outside domain: {} expanded by {}",
    bbox, distance
)]
pub struct ExpansionOutsideDomainError {
    bbox: BoundingBox,
    distance: usize,
}

impl fmt::Display for BoundingBox {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Box: ({}, {}, {}) - ({}, {}, {})",
            self.x0, self.y0, self.z0, self.x1, self.y1, self.z1
        )
    }
}

impl BoundingBox {
    /// Create a new `BoundingBox`. The smallest coordinates are (`x0`, `y0`, `z0`), and the largest
    /// coordinates are (`x1`, `y1`, `z1`).
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

    /// Returns true if two bounding boxes intersect or touch. The order of the boxes does not
    /// matter.
    pub fn intersects(&self, rhs: &BoundingBox) -> bool {
        !(self.x0 > rhs.x1
            || rhs.x0 > self.x1
            || self.y0 > rhs.y1
            || rhs.y0 > self.y1
            || self.z0 > rhs.z1
            || rhs.z0 > self.z1)
    }

    /// Returns whether this box is inside the `rhs`. The smaller coordinates are allowed to be at
    /// the boundary, while the bigger coordinates must be smaller than the surrounding box to be
    /// inside.
    pub fn inside(&self, rhs: &BoundingBox) -> bool {
        self.x0 >= rhs.x0
            && self.y0 >= rhs.y0
            && self.z0 >= rhs.z0
            && self.x1 < rhs.x1
            && self.y1 < rhs.y1
            && self.z1 < rhs.z1
    }

    /// Return an expanded version of this box. This function returns a
    /// `Result<BoundingBox, ExpansionOutsideDomainError>` because the expansion could cause an
    /// underflow.
    pub fn expanded(&self, distance: usize) -> Result<BoundingBox, ExpansionOutsideDomainError> {
        if distance > self.x0 || distance > self.y0 || distance > self.z0 {
            Err(ExpansionOutsideDomainError {
                bbox: *self,
                distance,
            })
        } else {
            Ok(BoundingBox::new(
                self.x0 - distance,
                self.y0 - distance,
                self.z0 - distance,
                self.x1 + distance,
                self.y1 + distance,
                self.z1 + distance,
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::world::boundingbox::BoundingBox;

    #[test]
    fn test_inside() {
        let outer = BoundingBox::new(0, 0, 0, 10, 10, 10);
        let inner = BoundingBox::new(0, 0, 0, 9, 9, 9);
        assert!(inner.inside(&outer));
        assert!(!outer.inside(&inner));
    }

    #[test]
    fn test_inside_not() {
        let outer = BoundingBox::new(0, 0, 0, 10, 10, 10);
        let inner = BoundingBox::new(0, 0, 0, 10, 10, 10);
        assert!(!inner.inside(&outer));
        assert!(!outer.inside(&inner));
    }

    #[test]
    fn test_intersects() {
        let first = BoundingBox::new(0, 0, 0, 4, 4, 4);
        let second = BoundingBox::new(2, 2, 2, 5, 5, 5);
        assert!(first.intersects(&second));
        assert!(second.intersects(&first));
    }

    #[test]
    fn test_intersects_not() {
        let first = BoundingBox::new(0, 0, 0, 4, 4, 4);
        let second = BoundingBox::new(5, 0, 0, 6, 3, 3);
        assert!(!first.intersects(&second));
        assert!(!second.intersects(&first));
    }

    #[test]
    fn test_expanded() {
        let first = BoundingBox::new(3, 3, 3, 4, 4, 4);
        let expected = BoundingBox::new(0, 0, 0, 7, 7, 7);
        assert_eq!(first.expanded(3).unwrap(), expected);
        assert!(first.expanded(4).is_err());
    }
}
