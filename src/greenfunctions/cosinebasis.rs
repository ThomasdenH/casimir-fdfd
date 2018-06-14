use greenfunctions::operator::{Operator, OperatorType};
use nalgebra::*;
use rayon::iter::*;
use scalarfield::ScalarField;
use std::f64::consts::PI;
use vectorfield::VectorField;

#[derive(Eq, PartialEq, Copy, Clone, Hash, Debug)]
enum Direction {
    X,
    Y,
    Z,
}

impl Direction {
    fn vector(self) -> Vector3<f64> {
        match self {
            Direction::X => Vector3::new(1.0, 0.0, 0.0),
            Direction::Y => Vector3::new(0.0, 1.0, 0.0),
            Direction::Z => Vector3::new(0.0, 0.0, 1.0),
        }
    }
}

pub struct CosineBasis<'a> {
    p0: Point3<usize>,
    p1: Point3<usize>,
    frequency: f64,
    normal: Direction,
    permitivity: &'a ScalarField,
}

impl<'a> CosineBasis<'a> {
    /// Will construct a new cosine basis with the following points as start and end points. The
    /// domain will end at the ending point, so it will not be an additional unit cell bigger.
    pub fn new(
        p0: Point3<usize>,
        p1: Point3<usize>,
        frequency: f64,
        permitivity: &'a ScalarField,
    ) -> CosineBasis<'a> {
        let normal = {
            if p0.x == p1.x {
                Direction::X
            } else if p0.y == p1.y {
                Direction::Y
            } else if p0.z == p1.z {
                Direction::Z
            } else {
                panic!("Coordinates are not on the same edge");
            }
        };
        CosineBasis {
            p0,
            p1,
            frequency,
            permitivity,
            normal,
        }
    }

    pub fn force(&self) -> Vector3<f64> {
        let (amax, bmax) = match self.normal {
            Direction::X => (self.p1.y - self.p0.y, self.p1.z - self.p0.z),
            Direction::Y => (self.p1.x - self.p0.x, self.p1.z - self.p0.z),
            Direction::Z => (self.p1.x - self.p0.x, self.p1.y - self.p0.y),
        };

        (0..3usize.min(amax))
            .into_par_iter()
            .flat_map(|na| {
                (0..3usize.min(bmax))
                    .into_par_iter()
                    .map(move |nb| {
                        let force = self.force_for_basis(na, nb);
                        println!("Force for {}, {}: {:?}", na, nb, force);
                        force
                    })
            })
            .sum()
    }

    fn get_source(&self, na: usize, nb: usize, polarization: Direction) -> VectorField {
        let mut source_field = VectorField::new(self.permitivity.size());
        match self.normal {
            Direction::X => {
                let dy = self.p1.y - self.p0.y;
                let dz = self.p1.z - self.p0.z;
                let vector = 2.0 / ((dy * dz) as f64).sqrt() * polarization.vector();
                for y in self.p0.y..self.p1.y {
                    for z in self.p0.z..self.p1.z {
                        source_field[(self.p0.x, y, z)] = vector
                            * (na as f64 * PI * y as f64 / dy as f64).cos()
                            * (nb as f64 * PI * z as f64 / dz as f64).cos();
                    }
                }
            }
            Direction::Y => {
                let dx = self.p1.x - self.p0.x;
                let dz = self.p1.z - self.p0.z;
                let vector = 2.0 / ((dx * dz) as f64).sqrt() * polarization.vector();
                for x in self.p0.x..self.p1.x {
                    for z in self.p0.z..self.p1.z {
                        source_field[(x, self.p0.y, z)] = vector
                            * (na as f64 * PI * x as f64 / dx as f64).cos()
                            * (nb as f64 * PI * z as f64 / dz as f64).cos();
                    }
                }
            }
            Direction::Z => {
                let dx = self.p1.x - self.p0.x;
                let dy = self.p1.y - self.p0.y;
                let vector = 2.0 / ((dx * dy) as f64).sqrt() * polarization.vector();
                for x in self.p0.x..self.p1.x {
                    for y in self.p0.y..self.p1.y {
                        source_field[(x, y, self.p0.z)] = vector
                            * (na as f64 * PI * x as f64 / dx as f64).cos()
                            * (nb as f64 * PI * y as f64 / dy as f64).cos();
                    }
                }
            }
        }
        source_field
    }

    pub fn force_for_basis(&self, na: usize, nb: usize) -> Vector3<f64> {
        self.stress_tensor(na, nb) * self.normal.vector()
    }

    pub fn stress_tensor(&self, na: usize, nb: usize) -> Matrix3<f64> {
        let electric_tensor = self.green_tensor(na, nb, OperatorType::Electric);
        let magnetic_tensor = self.green_tensor(na, nb, OperatorType::Magnetic);

        self.frequency * self.frequency / PI
            * ((magnetic_tensor - Matrix3::from_diagonal_element(0.5) * magnetic_tensor.trace())
                * (electric_tensor - Matrix3::from_diagonal_element(0.5) * electric_tensor.trace()))
    }

    pub fn green_tensor(&self, na: usize, nb: usize, operator_type: OperatorType) -> Matrix3<f64> {
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
    ) -> Vector3<f64> {
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

        for _ in 0..volume {
            let a_p = &a * p.clone();
            let alpha = rsold / (&p * &a_p);
            x += alpha * &p;
            r -= alpha * &a_p;
            let rsnew = &r * &r;
            if rsnew.sqrt() < 0.01 {
                break;
            }
            p = (rsnew / rsold) * p + &r;
            rsold = rsnew;
        }

        // Integrate over the surface
        let mut green = Vector3::new(0.0, 0.0, 0.0);
        match self.normal {
            Direction::X => {
                for y in self.p0.y..self.p1.y {
                    for z in self.p0.z..self.p1.z {
                        green += x[(self.p0.x, y, z)];
                    }
                }
            }
            Direction::Y => {
                for px in self.p0.x..self.p1.x {
                    for z in self.p0.z..self.p1.z {
                        green += x[(px, self.p0.y, z)];
                    }
                }
            }
            Direction::Z => {
                for px in self.p0.x..self.p1.x {
                    for y in self.p0.y..self.p1.y {
                        green += x[(px, y, self.p0.z)];
                    }
                }
            }
        }
        green
    }
}
