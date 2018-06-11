use nalgebra::*;
use std::ops::{Add, Mul, Sub, Index, IndexMut};

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
                                self[(x, y + 1, z)].z - self[(x, y, z)].z
                                    - self[(x, y, z + 1)].y + self[(x, y, z)].y,
                                self[(x, y, z + 1)].x - self[(x, y, z)].x
                                    - self[(x + 1, y, z)].z + self[(x, y, z)].z,
                                self[(x + 1, y, z)].y - self[(x, y, z)].y
                                    - self[(x, y + 1, z)].x + self[(x, y, z)].x,
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
                                self[(x, y, z)].z - self[(x, y - 1, z)].z
                                    - self[(x, y, z)].y + self[(x, y, z - 1)].y,
                                self[(x, y, z)].x - self[(x, y, z - 1)].x
                                    - self[(x, y, z)].z + self[(x - 1, y, z)].z,
                                self[(x, y, z)].y - self[(x - 1, y, z)].y
                                    - self[(x, y, z)].x + self[(x, y - 1, z)].x,
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
        for x in 0..self.nx {
            for y in 0..self.ny {
                for z in 0..self.nz {
                    self[(x, y, z)] *= rhs[(x, y, z)];
                }
            }
        }
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
        for x in 0..rhs.nx {
            for y in 0..rhs.ny {
                for z in 0..rhs.nz {
                    rhs[(x, y, z)] *= self;
                }
            }
        }
        rhs
    }
}

impl Index<(usize, usize, usize)> for VectorField {
    type Output = Vector3<f32>;

    fn index(&self, index: (usize, usize, usize)) -> &Vector3<f32> {
        let (x, y, z) = index;
        assert!(x < self.nx && y < self.ny && z < self.nz);
        &self.vectors[x + self.nx * (y + self.ny * z)]
    }
}

impl Index<(isize, isize, isize)> for VectorField {
    type Output = Vector3<f32>;

    fn index(&self, index: (isize, isize, isize)) -> &Vector3<f32> {
        let (x, y, z) = index;

        // Implementing periodic boundary conditions
        let nx = self.nx as isize;
        let ny = self.ny as isize;
        let nz = self.nz as isize;
        let x: usize = (((x % nx) + nx) % nx) as usize;
        let y: usize = (((y % ny) + ny) % ny) as usize;
        let z: usize = (((z % nz) + nz) % nz) as usize;
        &self[(x, y, z)]
    }
}

impl Index<Point3<usize>> for VectorField {
    type Output = Vector3<f32>;

    fn index(&self, index: Point3<usize>) -> &Vector3<f32> {
        &self[(index.x, index.y, index.z)]
    }
}

impl Index<Point3<isize>> for VectorField {
    type Output = Vector3<f32>;

    fn index(&self, index: Point3<isize>) -> &Vector3<f32> {
        &self[(index.x, index.y, index.z)]
    }
}

impl IndexMut<(usize, usize, usize)> for VectorField {
    fn index_mut<'a>(&'a mut self, index: (usize, usize, usize)) -> &'a mut Vector3<f32> {
        let (x, y, z) = index;
        assert!(x < self.nx && y < self.ny && z < self.nz);
        &mut self.vectors[x + self.nx * (y + self.ny * z)]
    }
}

impl IndexMut<Point3<usize>> for VectorField {
    fn index_mut<'a>(&'a mut self, index: Point3<usize>) -> &'a mut Vector3<f32> {
        &mut self[(index.x, index.y, index.z)]
    }
}
