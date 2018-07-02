use config::SimulationConfig;
use fields::{ScalarField, VectorField};
use greenfunctions::operator::{Operator, OperatorType};
use nalgebra::*;
use pbr::ProgressBar;
use std::f32::consts::PI;
use std::io::Stdout;
use std::sync::{Arc, Mutex};

#[derive(Eq, PartialEq, Copy, Clone, Hash, Debug)]
pub enum Direction {
    X,
    Y,
    Z,
    NegX,
    NegY,
    NegZ,
}

impl Direction {
    fn vector(self) -> Vector3<f32> {
        match self {
            Direction::X => Vector3::new(1.0, 0.0, 0.0),
            Direction::Y => Vector3::new(0.0, 1.0, 0.0),
            Direction::Z => Vector3::new(0.0, 0.0, 1.0),
            Direction::NegX => Vector3::new(-1.0, 0.0, 0.0),
            Direction::NegY => Vector3::new(0.0, -1.0, 0.0),
            Direction::NegZ => Vector3::new(0.0, 0.0, -1.0),
        }
    }
}

pub struct CosineBasis<'a> {
    p0: Point3<usize>,
    p1: Point3<usize>,
    frequency: f32,
    normal: Direction,
    permitivity: &'a ScalarField,
    simulation_config: &'a SimulationConfig,
    progress_bar: Option<Arc<Mutex<ProgressBar<Stdout>>>>,
}

impl<'a> CosineBasis<'a> {
    /// Will construct a new cosine basis with the following points as start and end points.
    pub fn new(
        p0: Point3<usize>,
        p1: Point3<usize>,
        frequency: f32,
        permitivity: &'a ScalarField,
        simulation_config: &'a SimulationConfig,
        normal: Direction,
    ) -> CosineBasis<'a> {
        CosineBasis {
            p0,
            p1,
            frequency,
            permitivity,
            normal,
            simulation_config,
            progress_bar: None,
        }
    }

    pub fn with_progress_bar(
        mut self,
        progress_bar: Option<Arc<Mutex<ProgressBar<Stdout>>>>,
    ) -> CosineBasis<'a> {
        self.progress_bar = progress_bar;
        self
    }

    pub fn force(&self) -> Vector3<f32> {
        let (amax, bmax) = match self.normal {
            Direction::X | Direction::NegX => (self.p1.y - self.p0.y, self.p1.z - self.p0.z),
            Direction::Y | Direction::NegY => (self.p1.x - self.p0.x, self.p1.z - self.p0.z),
            Direction::Z | Direction::NegZ => (self.p1.x - self.p0.x, self.p1.y - self.p0.y),
        };

        let mut total_force = Vector3::new(0.0, 0.0, 0.0);
        let mut remaining = amax * bmax;
        let mut count = 0;

        for n_total in 0..=(amax + bmax) {
            let b_start = n_total.min(bmax);
            let b_end = 0.max(n_total as i64 - amax as i64) as usize;

            let mut difference = 0.0;

            for nb in b_start..=b_end {
                let na = n_total - nb;
                let force = self.force_for_basis(na, nb);
                total_force += force;
                difference += force.norm();
                count += 1;

                if let Some(ref progress_bar) = self.progress_bar {
                    progress_bar.lock().unwrap().inc();
                    remaining -= 1;
                }
            }

            if difference / ((b_start - b_end + 1) as f32) < total_force.norm() / count as f32 {
                if let Some(ref progress_bar) = self.progress_bar {
                    progress_bar.lock().unwrap().add(remaining as u64);
                }
                return total_force;
            }
        }
        total_force
    }

    fn get_source(&self, na: usize, nb: usize, polarization: Direction) -> VectorField {
        let mut source_field = VectorField::new(self.permitivity.size());
        match self.normal {
            Direction::X | Direction::NegX => {
                let dy = self.p1.y - self.p0.y;
                let dz = self.p1.z - self.p0.z;
                let vector = 2.0 / ((dy * dz) as f32).sqrt() * polarization.vector();
                for y in self.p0.y..self.p1.y {
                    for z in self.p0.z..self.p1.z {
                        source_field[(self.p0.x, y, z)] = vector
                            * (na as f32 * PI * y as f32 / dy as f32).cos()
                            * (nb as f32 * PI * z as f32 / dz as f32).cos();
                    }
                }
            }
            Direction::Y | Direction::NegY => {
                let dx = self.p1.x - self.p0.x;
                let dz = self.p1.z - self.p0.z;
                let vector = 2.0 / ((dx * dz) as f32).sqrt() * polarization.vector();
                for x in self.p0.x..self.p1.x {
                    for z in self.p0.z..self.p1.z {
                        source_field[(x, self.p0.y, z)] = vector
                            * (na as f32 * PI * x as f32 / dx as f32).cos()
                            * (nb as f32 * PI * z as f32 / dz as f32).cos();
                    }
                }
            }
            Direction::Z | Direction::NegZ => {
                let dx = self.p1.x - self.p0.x;
                let dy = self.p1.y - self.p0.y;
                let vector = 2.0 / ((dx * dy) as f32).sqrt() * polarization.vector();
                for x in self.p0.x..self.p1.x {
                    for y in self.p0.y..self.p1.y {
                        source_field[(x, y, self.p0.z)] = vector
                            * (na as f32 * PI * x as f32 / dx as f32).cos()
                            * (nb as f32 * PI * y as f32 / dy as f32).cos();
                    }
                }
            }
        }
        source_field
    }

    pub fn force_for_basis(&self, na: usize, nb: usize) -> Vector3<f32> {
        self.stress_tensor(na, nb) * self.normal.vector()
    }

    pub fn stress_tensor(&self, na: usize, nb: usize) -> Matrix3<f32> {
        let electric_tensor = self.green_tensor(na, nb, OperatorType::Electric);
        let magnetic_tensor = self.green_tensor(na, nb, OperatorType::Magnetic);

        self.frequency * self.frequency / PI
            * ((magnetic_tensor - Matrix3::from_diagonal_element(0.5) * magnetic_tensor.trace())
                * (electric_tensor - Matrix3::from_diagonal_element(0.5) * electric_tensor.trace()))
    }

    pub fn green_tensor(&self, na: usize, nb: usize, operator_type: OperatorType) -> Matrix3<f32> {
        Matrix3::from_columns(&[
            self.green_function(&self.get_source(na, nb, Direction::X), operator_type),
            self.green_function(&self.get_source(na, nb, Direction::Y), operator_type),
            self.green_function(&self.get_source(na, nb, Direction::Z), operator_type),
        ])
    }

    pub fn green_function(
        &self,
        source: &VectorField,
        operator_type: OperatorType,
    ) -> Vector3<f32> {
        // The operator
        let a = Operator::new(self.frequency, self.permitivity, operator_type);
        let size = self.permitivity.size();

        let mut x = VectorField::new(size);
        let mut r = source - &a * x.clone();
        let mut p = r.clone();
        let mut rsold = &r * &r;

        // In theory the conjugate gradient method should converge in N steps. In practice,it converges
        // much quicker.
        let volume = size.x * size.y * size.z;

        let mut temp1 = VectorField::new(size);
        let mut temp2 = VectorField::new(size);
        let mut a_p = VectorField::new(size);

        for _ in 0..volume {
            let (next_a_p, next_temp1, next_temp2) =
                a.mul_with_temps(p.clone_to(a_p), temp1, temp2);
            a_p = next_a_p;
            temp1 = next_temp1;
            temp2 = next_temp2;

            let alpha = rsold / (&p * &a_p);
            x += alpha * &p;
            r -= alpha * &a_p;
            let rsnew = &r * &r;
            if rsnew.sqrt() < self.simulation_config.fdfd_convergence {
                break;
            }
            p = (rsnew / rsold) * p + &r;
            rsold = rsnew;
        }

        // Integrate over the surface
        let mut green = Vector3::new(0.0, 0.0, 0.0);
        match self.normal {
            Direction::X | Direction::NegX => {
                for py in self.p0.y..self.p1.y {
                    for pz in self.p0.z..self.p1.z {
                        green += x[(self.p0.x, py, pz)];
                    }
                }
            }
            Direction::Y | Direction::NegY => {
                for px in self.p0.x..self.p1.x {
                    for pz in self.p0.z..self.p1.z {
                        green += x[(px, self.p0.y, pz)];
                    }
                }
            }
            Direction::Z | Direction::NegZ => {
                for px in self.p0.x..self.p1.x {
                    for py in self.p0.y..self.p1.y {
                        green += x[(px, py, self.p0.z)];
                    }
                }
            }
        }
        green
    }
}
