use greenfunctions::stress_tensor;
use nalgebra::*;
use scalarfield::ScalarField;
use pbr::ProgressBar;

const OBJECT_PERMITIVITY: f32 = 1000_000.0;

pub struct BoundingBox {
    pub x0: usize,
    pub y0: usize,
    pub z0: usize,
    pub x1: usize,
    pub y1: usize,
    pub z1: usize,
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
    inv_perm: ScalarField,
    inv_mag: ScalarField,
}

impl World {
    pub fn new(nx: usize, ny: usize, nz: usize) -> World {
        World {
            nx,
            ny,
            nz,
            bboxes: Vec::new(),
            perm: ScalarField::ones(nx, ny, nz),
            mag: ScalarField::ones(nx, ny, nz),
            inv_perm: ScalarField::ones(nx, ny, nz),
            inv_mag: ScalarField::ones(nx, ny, nz),
        }
    }

    pub fn add_box(&mut self, x0: usize, y0: usize, z0: usize, x1: usize, y1: usize, z1: usize) {
        assert!(x0 < x1 && y0 < y1 && z0 < z1 && x1 < self.nx && y1 < self.ny && z1 < self.nz);
        let bbox = BoundingBox::new(x0, y0, z0, x1, y1, z1);
        assert!(!self.bboxes.iter().any(|b| b.intersects(&bbox)));
        self.bboxes.push(bbox);
        for x in x0..x1 {
            for y in y0..y1 {
                for z in z0..z1 {
                    self.perm
                        .set(x as isize, y as isize, z as isize, OBJECT_PERMITIVITY);
                    self.inv_perm
                        .set(x as isize, y as isize, z as isize, 1.0 / OBJECT_PERMITIVITY);
                }
            }
        }
    }

    pub fn force_on(&self, i: usize) -> Vector3<f32> {
        // The maximum frequency is given by (speed of light) / (grid element size)

        let start_freq = 0.01;
        let end_freq = 1.0;
        let start_force = self.force_on_for_freq(i, start_freq);
        let end_force = self.force_on_for_freq(i, end_freq);

        let force = self.integrate_force_between_frequencies(i, start_freq, end_freq, start_force, end_force);

        force
    }

    pub fn integrate_force_between_frequencies(
        &self,
        i: usize,
        start_frequency: f32,
        end_frequency: f32,
        start_value: Vector3<f32>,
        end_value: Vector3<f32>
    ) -> Vector3<f32> {
        // Do a recursive integration. The function should be smooth.
        if (start_value - end_value).norm() < 0.001 {
            // The difference is small enough to do a line approximation
            0.5 * (start_value + end_value) * (end_frequency - start_frequency)
        } else {
            let middle_frequency = 0.5 * (start_frequency + end_frequency);
            let middle_value = self.force_on_for_freq(i, middle_frequency);
            self.integrate_force_between_frequencies(
                i,
                start_frequency,
                middle_frequency,
                start_value,
                middle_value
            )
                + self.integrate_force_between_frequencies(
                    i,
                    middle_frequency,
                    end_frequency,
                    middle_value,
                    end_value
                )
        }
    }

    fn force_on_for_freq(&self, i: usize, frequency: f32) -> Vector3<f32> {
        let mut progress = ProgressBar::new(2 * (self.nx * self.ny + self.ny * self.nz + self.nz * self.nx) as u64);

        let size = Vector3::new(self.nx, self.ny, self.nz);
        let mut force = Vector3::new(0.0, 0.0, 0.0);

        // Integrate the force over the faces of the cube
        let bbox = &self.bboxes[i];

        for x in bbox.x0..bbox.x1 {
            for y in bbox.y0..bbox.y1 {
                // Top face
                force += stress_tensor(
                    frequency,
                    Point3::new(x, y, bbox.z1 + 1),
                    &self.perm,
                    &self.mag,
                    &self.inv_perm,
                    &self.inv_mag,
                    size,
                ) * Vector3::new(0.0, 0.0, 1.0);
                progress.inc();

                // Bottom face
                force += stress_tensor(
                    frequency,
                    Point3::new(x, y, bbox.z0 - 1),
                    &self.perm,
                    &self.mag,
                    &self.inv_perm,
                    &self.inv_mag,
                    size,
                ) * Vector3::new(0.0, 0.0, -1.0);
                progress.inc();
            }
        }

        for x in bbox.x0..bbox.x1 {
            for z in bbox.z0..bbox.z1 {
                // Front
                force += stress_tensor(
                    frequency,
                    Point3::new(x, bbox.y1 + 1, z),
                    &self.perm,
                    &self.mag,
                    &self.inv_perm,
                    &self.inv_mag,
                    size,
                ) * Vector3::new(0.0, 1.0, 0.0);
                progress.inc();

                // Back
                force += stress_tensor(
                    frequency,
                    Point3::new(x, bbox.y0 - 1, z),
                    &self.perm,
                    &self.mag,
                    &self.inv_perm,
                    &self.inv_mag,
                    size,
                ) * Vector3::new(0.0, -1.0, 0.0);
                progress.inc();
            }
        }

        for y in bbox.y0..bbox.y1 {
            for z in bbox.z0..bbox.z1 {
                // Front
                force += stress_tensor(
                    frequency,
                    Point3::new(bbox.x1 + 1, y, z),
                    &self.perm,
                    &self.mag,
                    &self.inv_perm,
                    &self.inv_mag,
                    size,
                ) * Vector3::new(1.0, 0.0, 0.0);
                progress.inc();

                // Back
                force += stress_tensor(
                    frequency,
                    Point3::new(bbox.x0 - 1, y, z),
                    &self.perm,
                    &self.mag,
                    &self.inv_perm,
                    &self.inv_mag,
                    size,
                ) * Vector3::new(-1.0, 0.0, 0.0);
                progress.inc();
            }
        }

        progress.finish();

        force
    }
}
