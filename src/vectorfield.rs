use nalgebra::*;
use std::ops::{Add, Mul, Sub};

use scalarfield::ScalarField;

#[derive(PartialEq, Clone, Debug)]
pub struct VectorField {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub vectors: DVector<Vector3<f32>>,
}

impl VectorField {
    pub fn new(nx: usize, ny: usize, nz: usize) -> VectorField {
        VectorField {
            nx,
            ny,
            nz,
            vectors: DVector::from_element(nx * ny * nz, Vector3::new(0.0, 0.0, 0.0)),
        }
    }

    /// Get the vector at location x, y and z
    pub fn at(&self, x: isize, y: isize, z: isize) -> Vector3<f32> {
        if x < 0 || x as usize >= self.nx
            || y < 0 || y as usize >= self.ny
            || z < 0 || z as usize >= self.nz
        {
            Vector3::new(0.0, 0.0, 0.0)
        } else {
            self.vectors[x as usize + self.nx * (y as usize + self.ny * (z as usize))]
        }
    }

    pub fn set(&mut self, x: isize, y: isize, z: isize, value: Vector3<f32>) {
        if x < 0 || x as usize >= self.nx
            || y < 0 || y as usize >= self.ny
            || z < 0 || z as usize >= self.nz {
            return;
        }
        self.vectors[x as usize + self.nx * (y as usize + self.ny * (z as usize))] = value;
    }

    /// Calculates the positive curl of a vector field. The new positions will be the centers of the
    /// Yee grid planes formed between four grid edges.
    pub fn curl_positive(&self) -> VectorField {
        VectorField {
            nx: self.nx,
            ny: self.ny,
            nz: self.nz,
            vectors: DVector::from_iterator(
                self.nx * self.ny * self.nz,
                (0..self.nz as isize).flat_map(|z| {
                    (0..self.ny as isize).flat_map(move |y| {
                        (0..self.nx as isize).map(move |x| {
                            Vector3::new(
                                self.at(x, y + 1, z).z - self.at(x, y, z).z - self.at(x, y, z + 1).y
                                    + self.at(x, y, z).y,
                                self.at(x, y, z + 1).x - self.at(x, y, z).x - self.at(x + 1, y, z).z
                                    + self.at(x, y, z).z,
                                self.at(x + 1, y, z).y - self.at(x, y, z).y - self.at(x, y + 1, z).x
                                    - self.at(x, y, z).x,
                            )
                        })
                    })
                }),
            ),
        }
    }

    pub fn curl_negative(&self) -> VectorField {
        VectorField {
            nx: self.nx,
            ny: self.ny,
            nz: self.nz,
            vectors: DVector::from_iterator(
                self.nx * self.ny * self.nz,
                (0..self.nz as isize).flat_map(|z| {
                    (0..self.ny as isize).flat_map(move |y| {
                        (0..self.nx as isize).map(move |x| {
                            Vector3::new(
                                self.at(x, y, z).z - self.at(x, y - 1, z).z - self.at(x, y, z).y
                                    + self.at(x, y, z - 1).y,
                                self.at(x, y, z).x - self.at(x, y, z - 1).x - self.at(x, y, z).z
                                    + self.at(x - 1, y, z).z,
                                self.at(x, y, z).y - self.at(x - 1, y, z).y - self.at(x, y, z).x
                                    - self.at(x, y - 1, z).x,
                            )
                        })
                    })
                }),
            ),
        }
    }
}

impl<'a> Add<&'a VectorField> for VectorField {
    type Output = VectorField;

    fn add(mut self, f: &'a VectorField) -> Self::Output {
        assert!(self.nx == f.nx && self.ny == f.ny && self.nz == f.nz);
        self.vectors += &f.vectors;
        self
    }
}

impl<'a> Sub<&'a VectorField> for VectorField {
    type Output = VectorField;

    fn sub(mut self, f: &'a VectorField) -> Self::Output {
        assert!(self.nx == f.nx && self.ny == f.ny && self.nz == f.nz);
        self.vectors -= &f.vectors;
        self
    }
}

impl<'a> Mul<&'a ScalarField> for VectorField {
    type Output = VectorField;

    fn mul(mut self, rhs: &'a ScalarField) -> VectorField {
        assert!(self.nx == rhs.nx && self.ny == rhs.ny && self.nz == rhs.nz);
        self.vectors = self.vectors.zip_map(&rhs.scalars, |a, b| a * b);
        self
    }
}

impl<'a, 'b> Mul<&'a VectorField> for &'b VectorField {
    type Output = f32;

    fn mul(self, rhs: &'a VectorField) -> Self::Output {
        assert!(self.nx == rhs.nx && self.ny == rhs.ny && self.nz == rhs.nz);
        self.vectors
            .iter()
            .zip(rhs.vectors.iter())
            .map(|(a, b)| a.dot(b))
            .sum()
    }
}

impl Mul<VectorField> for f32 {
    type Output = VectorField;

    fn mul(self, mut rhs: VectorField) -> Self::Output {
        rhs.vectors = rhs.vectors.map(|a| a * self);
        rhs
    }
}
