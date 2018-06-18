use nalgebra::*;
use std::ops::{Add, AddAssign, Div, Index, IndexMut, Mul, Neg, Sub, SubAssign};

use scalarfield::ScalarField;

#[derive(PartialEq, Clone, Debug)]
pub struct VectorField {
    size: Vector3<usize>,
    vectors: DVector<Vector3<f32>>,
}

impl VectorField {
    pub fn new(size: Vector3<usize>) -> VectorField {
        VectorField {
            size,
            vectors: DVector::from_element(size.x * size.y * size.z, Vector3::new(0.0, 0.0, 0.0)),
        }
    }

    /// Calculates the positive curl of a vector field. The new positions will be the centers of the
    /// Yee grid planes formed between four grid edges.
    pub fn curl_positive(&self) -> VectorField {
        VectorField {
            size: self.size,
            vectors: DVector::from_iterator(
                self.vectors.len(),
                (0..self.size.z as isize).flat_map(|z| {
                    (0..self.size.y as isize).flat_map(move |y| {
                        (0..self.size.x as isize).map(move |x| {
                            Vector3::new(
                                self[(x, y + 1, z)].z - self[(x, y, z)].z - self[(x, y, z + 1)].y
                                    + self[(x, y, z)].y,
                                self[(x, y, z + 1)].x - self[(x, y, z)].x - self[(x + 1, y, z)].z
                                    + self[(x, y, z)].z,
                                self[(x + 1, y, z)].y - self[(x, y, z)].y - self[(x, y + 1, z)].x
                                    + self[(x, y, z)].x,
                            )
                        })
                    })
                }),
            ),
        }
    }

    pub fn curl_negative(&self) -> VectorField {
        VectorField {
            size: self.size,
            vectors: DVector::from_iterator(
                self.vectors.len(),
                (0..self.size.z as isize).flat_map(|z| {
                    (0..self.size.y as isize).flat_map(move |y| {
                        (0..self.size.x as isize).map(move |x| {
                            Vector3::new(
                                self[(x, y, z)].z - self[(x, y - 1, z)].z - self[(x, y, z)].y
                                    + self[(x, y, z - 1)].y,
                                self[(x, y, z)].x - self[(x, y, z - 1)].x - self[(x, y, z)].z
                                    + self[(x - 1, y, z)].z,
                                self[(x, y, z)].y - self[(x - 1, y, z)].y - self[(x, y, z)].x
                                    + self[(x, y - 1, z)].x,
                            )
                        })
                    })
                }),
            ),
        }
    }

    pub fn size(&self) -> Vector3<usize> {
        self.size
    }

    pub fn len(&self) -> usize {
        self.vectors.len()
    }
}

impl<'a> Add<&'a VectorField> for VectorField {
    type Output = VectorField;

    fn add(mut self, rhs: &'a VectorField) -> Self::Output {
        debug_assert!(self.size == rhs.size);
        self.vectors += &rhs.vectors;
        self
    }
}

impl<'a> Sub<&'a VectorField> for VectorField {
    type Output = VectorField;

    fn sub(mut self, rhs: &'a VectorField) -> Self::Output {
        debug_assert!(self.size == rhs.size);
        self.vectors -= &rhs.vectors;
        self
    }
}

impl<'a> Sub<VectorField> for &'a VectorField {
    type Output = VectorField;

    fn sub(self, mut rhs: VectorField) -> Self::Output {
        debug_assert!(self.size == rhs.size);
        rhs.vectors -= &self.vectors;
        rhs
    }
}

impl<'a> Mul<&'a ScalarField> for VectorField {
    type Output = VectorField;

    fn mul(mut self, rhs: &'a ScalarField) -> VectorField {
        debug_assert!(self.size() == rhs.size());
        for i in 0..self.vectors.len() {
            self[i] *= rhs[i];
        }
        self
    }
}

impl<'a> Div<&'a ScalarField> for VectorField {
    type Output = VectorField;

    fn div(mut self, rhs: &'a ScalarField) -> VectorField {
        debug_assert!(self.size() == rhs.size());
        for i in 0..self.vectors.len() {
            self[i] /= rhs[i];
        }
        self
    }
}

impl<'a, 'b> Mul<&'a VectorField> for &'b VectorField {
    type Output = f32;

    fn mul(self, rhs: &'a VectorField) -> Self::Output {
        let mut sum = 0.0;
        debug_assert!(self.size() == rhs.size());
        for i in 0..self.vectors.len() {
            sum += self[i].dot(&rhs[i]);
        }
        sum
    }
}

impl Mul<VectorField> for f32 {
    type Output = VectorField;

    fn mul(self, mut rhs: VectorField) -> Self::Output {
        for i in 0..rhs.vectors.len() {
            rhs[i] *= self;
        }
        rhs
    }
}

impl Index<(usize, usize, usize)> for VectorField {
    type Output = Vector3<f32>;

    fn index(&self, index: (usize, usize, usize)) -> &Vector3<f32> {
        let (x, y, z) = index;
        debug_assert!(x < self.size.x && y < self.size.y && z < self.size.z);
        &self.vectors[x + self.size.x * (y + self.size.y * z)]
    }
}

impl Index<(isize, isize, isize)> for VectorField {
    type Output = Vector3<f32>;

    fn index(&self, index: (isize, isize, isize)) -> &Vector3<f32> {
        let (x, y, z) = index;

        // Implementing periodic boundary conditions
        let nx = self.size.x as isize;
        let ny = self.size.y as isize;
        let nz = self.size.z as isize;
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
        debug_assert!(x < self.size.x && y < self.size.y && z < self.size.z);
        &mut self.vectors[x + self.size.x * (y + self.size.y * z)]
    }
}

impl IndexMut<Point3<usize>> for VectorField {
    fn index_mut<'a>(&'a mut self, index: Point3<usize>) -> &'a mut Vector3<f32> {
        &mut self[(index.x, index.y, index.z)]
    }
}

impl Index<usize> for VectorField {
    type Output = Vector3<f32>;

    fn index(&self, index: usize) -> &Vector3<f32> {
        &self.vectors[index]
    }
}

impl IndexMut<usize> for VectorField {
    fn index_mut(&mut self, index: usize) -> &mut Vector3<f32> {
        &mut self.vectors[index]
    }
}

impl Neg for VectorField {
    type Output = VectorField;
    fn neg(mut self) -> Self::Output {
        self.vectors = -self.vectors;
        self
    }
}

impl<'a> AddAssign<&'a VectorField> for VectorField {
    fn add_assign(&mut self, rhs: &VectorField) {
        debug_assert!(self.size == rhs.size);
        for i in 0..self.vectors.len() {
            self[i] += rhs[i];
        }
    }
}

impl<'a> SubAssign<&'a VectorField> for VectorField {
    fn sub_assign(&mut self, rhs: &VectorField) {
        debug_assert!(self.size == rhs.size);
        for i in 0..self.vectors.len() {
            self[i] -= rhs[i];
        }
    }
}
