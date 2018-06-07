use scalarfield::ScalarField;

const OBJECT_PERMITIVITY: f32 = 1000_000.0;

pub struct BoundingBox {
    x0: usize,
    y0: usize,
    z0: usize,
    x1: usize,
    y1: usize,
    z1: usize,
}

impl BoundingBox {
    fn new(x0: usize, y0: usize, z0: usize, x1: usize, y1: usize, z1: usize) -> BoundingBox {
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
    fn intersects(&self, rhs: &BoundingBox) -> bool {
        !(self.x0 > rhs.x1 || rhs.x0 > self.x1 || self.y0 > rhs.y1 || rhs.y0 > self.y1
            || self.z0 > rhs.z1 || rhs.z0 > self.z0)
    }
}

/// A struct representing the geometry
pub struct World {
    nx: usize,
    ny: usize,
    nz: usize,
    bboxes: Vec<BoundingBox>,
    perm: ScalarField,
    mag: ScalarField,
}

impl World {
    fn new(nx: usize, ny: usize, nz: usize) -> World {
        World {
            nx,
            ny,
            nz,
            bboxes: Vec::new(),
            perm: ScalarField::ones(nx, ny, nz),
            mag: ScalarField::ones(nx, ny, nz),
        }
    }

    fn add_box(&mut self, x0: usize, y0: usize, z0: usize, x1: usize, y1: usize, z1: usize) {
        assert!(x0 < x1 && y0 < y1 && z0 < z1 && x1 < self.nx && y1 < self.ny && z1 < self.nz);
        let bbox = BoundingBox::new(x0, y0, z0, x1, y1, z1);
        assert!(!self.bboxes.iter().any(|b| b.intersects(&bbox)));
        self.bboxes.push(bbox);
        for x in x0..x1 {
            for y in y0..y1 {
                for z in z0..z1 {
                    self.perm
                        .set(x as isize, y as isize, z as isize, OBJECT_PERMITIVITY);
                }
            }
        }
    }
}
