use greenfunctions::stress_tensor;
use nalgebra::*;
use scalarfield::ScalarField;
use rayon::iter::*;
use std::f64::consts::PI;

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
    bboxes: Vec<BoundingBox>
}

impl World {
    pub fn new(nx: usize, ny: usize, nz: usize) -> World {
        World {
            nx,
            ny,
            nz,
            bboxes: Vec::new()
        }
    }

    pub fn add_box(&mut self, x0: usize, y0: usize, z0: usize, x1: usize, y1: usize, z1: usize) {
        debug_assert!(x0 < x1 && y0 < y1 && z0 < z1 && x1 < self.nx && y1 < self.ny && z1 < self.nz);
        let bbox = BoundingBox::new(x0, y0, z0, x1, y1, z1);
        debug_assert!(!self.bboxes.iter().any(|b| b.intersects(&bbox)));
        self.bboxes.push(bbox);
    }

    pub fn force_on(&self, i: usize) -> Vector3<f64> {
        // The maximum frequency is given by (speed of light) / (grid element size)
        let start_freq = 0.001;
        let end_freq = 1.0;
        let start_force = self.force_on_for_freq(i, start_freq);
        let end_force = self.force_on_for_freq(i, end_freq);

        self.integrate_force_between_frequencies(i, start_freq, end_freq, start_force, end_force)
    }

    pub fn integrate_force_between_frequencies(
        &self,
        i: usize,
        start_frequency: f64,
        end_frequency: f64,
        start_value: Vector3<f64>,
        end_value: Vector3<f64>
    ) -> Vector3<f64> {
        // Do a recursive integration. The function should be smooth.
        if (start_value - end_value).norm() < 0.05 {
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

    fn force_on_for_freq(&self, i: usize, frequency: f64) -> Vector3<f64> {
        let mut perm = ScalarField::ones(self.nx, self.ny, self.nz);
        let mag = ScalarField::ones(self.nx, self.ny, self.nz);
        let permitivity = World::permitivity(frequency);
        println!("{}", permitivity);
        for bbox in &self.bboxes {
            for x in bbox.x0..bbox.x1 {
                for y in bbox.y0..bbox.y1 {
                    for z in bbox.z0..bbox.z1 {
                        perm.set(x as isize, y as isize, z as isize, permitivity);
                    }
                }
            }
        }
        let mut inv_perm = perm.clone();
        inv_perm.multiplicative_invert();

        let size = Vector3::new(self.nx, self.ny, self.nz);
        let mut force = Vector3::new(0.0, 0.0, 0.0);

        // Integrate the force over the faces of the cube
        let bbox = &self.bboxes[i];

        println!("Force on for frequency {}.", frequency);

        let perm = &perm;
        let mag = &mag;
        let inv_perm = &inv_perm;

        force += (bbox.x0..bbox.x1).into_par_iter().flat_map(|x|
            (bbox.y0..bbox.y1).into_par_iter().map(move |y| {
                // Top face
                let a = stress_tensor(
                    frequency,
                    Point3::new(x, y, bbox.z1 + 1),
                    &perm,
                    &mag,
                    &inv_perm,
                    &mag,
                    size,
                ) * Vector3::new(0.0, 0.0, 1.0);

                let b = stress_tensor(
                    frequency,
                    Point3::new(x, y, bbox.z0 - 1),
                    &perm,
                    &mag,
                    &inv_perm,
                    &mag,
                    size,
                ) * Vector3::new(0.0, 0.0, -1.0);

                a + b
            })
        ).sum::<Vector3<f64>>();

        force += (bbox.x0..bbox.x1).into_par_iter().flat_map(|x|
            (bbox.z0..bbox.z1).into_par_iter().map(move |z| {
                // Front
                let a = stress_tensor(
                    frequency,
                    Point3::new(x, bbox.y1 + 1, z),
                    &perm,
                    &mag,
                    &inv_perm,
                    &mag,
                    size,
                ) * Vector3::new(0.0, 1.0, 0.0);

                // Back
                let b = stress_tensor(
                    frequency,
                    Point3::new(x, bbox.y0 - 1, z),
                    &perm,
                    &mag,
                    &inv_perm,
                    &mag,
                    size,
                ) * Vector3::new(0.0, -1.0, 0.0);

                a + b
            })
        ).sum::<Vector3<f64>>();

        force += (bbox.y0..bbox.y1).into_par_iter().flat_map(|y|
            (bbox.z0..bbox.z1).into_par_iter().map(move |z| {
                // Right
                let a = stress_tensor(
                    frequency,
                    Point3::new(bbox.x1 + 1, y, z),
                    &perm,
                    &mag,
                    &inv_perm,
                    &mag,
                    size,
                ) * Vector3::new(1.0, 0.0, 0.0);

                // Left
                let b = stress_tensor(
                    frequency,
                    Point3::new(bbox.x0 - 1, y, z),
                    &perm,
                    &mag,
                    &inv_perm,
                    &mag,
                    size,
                ) * Vector3::new(-1.0, 0.0, 0.0);

                a + b
            })
        ).sum::<Vector3<f64>>();

        println!("Force for frequency {}: ({}, {}, {})", frequency, force.x, force.y, force.z);

        force
    }

    fn permitivity(freq: f64) -> f64 {
        let omega_p = 7.79;
        let omega_tau = 48.8;
        let mut total = 0.0;
        for i in 0.. {
            let omega = f64::from(i) * 0.1;
            let added = (omega_p * omega_p * omega_tau)
                / (omega * omega + omega_tau * omega_tau)
                / (omega * omega + freq * freq) * 0.1;
            total += added;
            if added < 0.001 {
                break;
            }
        }
        1.0 + 2.0 / PI * total
    }
}
