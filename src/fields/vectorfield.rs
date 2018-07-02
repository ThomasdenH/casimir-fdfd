use nalgebra::*;
use std::ops::{Add, AddAssign, Div, Index, IndexMut, Mul, Neg, Sub, SubAssign};

use fields::ScalarField;

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

    pub fn curl_positive_to(&self, mut to: VectorField) -> VectorField {
        debug_assert!(self.size() == to.size());
        for x in 0..self.size.x as isize {
            for y in 0..self.size.y as isize {
                for z in 0..self.size.z as isize {
                    to[(x as usize, y as usize, z as usize)] = Vector3::new(
                        self[(x, y + 1, z)].z - self[(x, y, z)].z - self[(x, y, z + 1)].y
                            + self[(x, y, z)].y,
                        self[(x, y, z + 1)].x - self[(x, y, z)].x - self[(x + 1, y, z)].z
                            + self[(x, y, z)].z,
                        self[(x + 1, y, z)].y - self[(x, y, z)].y - self[(x, y + 1, z)].x
                            + self[(x, y, z)].x,
                    )
                }
            }
        }
        to
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

    pub fn curl_negative_to(&self, mut to: VectorField) -> VectorField {
        debug_assert!(self.size() == to.size());
        for x in 0..self.size.x as isize {
            for y in 0..self.size.y as isize {
                for z in 0..self.size.z as isize {
                    to[(x as usize, y as usize, z as usize)] = Vector3::new(
                        self[(x, y, z)].z - self[(x, y - 1, z)].z - self[(x, y, z)].y
                            + self[(x, y, z - 1)].y,
                        self[(x, y, z)].x - self[(x, y, z - 1)].x - self[(x, y, z)].z
                            + self[(x - 1, y, z)].z,
                        self[(x, y, z)].y - self[(x - 1, y, z)].y - self[(x, y, z)].x
                            + self[(x, y - 1, z)].x,
                    )
                }
            }
        }
        to
    }

    pub fn size(&self) -> Vector3<usize> {
        self.size
    }

    pub fn len(&self) -> usize {
        self.vectors.len()
    }

    pub fn clone_to(&self, mut other: VectorField) -> VectorField {
        debug_assert!(self.size() == other.size());
        for i in 0..self.vectors.len() {
            other[i] = self[i];
        }
        other
    }

    pub fn vectors(&self) -> &DVector<Vector3<f32>> {
        &self.vectors
    }

    pub fn vectors_mut(&mut self) -> &mut DVector<Vector3<f32>> {
        &mut self.vectors
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

#[test]
fn mul_vector_field_scalar_field() {
    let size = Vector3::new(2, 2, 1);
    let mut vector_field = VectorField::new(size);
    vector_field[(0, 0, 0)] = Vector3::repeat(1.0);
    vector_field[(1, 0, 0)] = Vector3::repeat(2.0);
    vector_field[(0, 1, 0)] = Vector3::repeat(3.0);
    vector_field[(1, 1, 0)] = Vector3::repeat(4.0);

    let mut scalar_field = ScalarField::ones(size);
    scalar_field[(0, 0, 0)] = 2.0;
    scalar_field[(1, 0, 0)] = 3.0;
    scalar_field[(0, 1, 0)] = 4.0;
    scalar_field[(1, 1, 0)] = 5.0;

    let result = vector_field * &scalar_field;
    assert!((result[(0usize, 0, 0)] - Vector3::repeat(2.0)).norm() < 1e-10);
    assert!((result[(1usize, 0, 0)] - Vector3::repeat(6.0)).norm() < 1e-10);
    assert!((result[(0usize, 1, 0)] - Vector3::repeat(12.0)).norm() < 1e-10);
    assert!((result[(1usize, 1, 0)] - Vector3::repeat(20.0)).norm() < 1e-10);
}

impl<'a> Mul<&'a ScalarField> for VectorField {
    type Output = VectorField;

    fn mul(mut self, rhs: &'a ScalarField) -> VectorField {
        debug_assert!(self.size() == rhs.size());
        for (element_self, element_other) in self.vectors.iter_mut().zip(rhs.scalars().iter()) {
            *element_self *= *element_other;
        }
        self
    }
}

#[test]
fn div_vector_field_scalar_field() {
    let size = Vector3::new(2, 2, 1);
    let mut vector_field = VectorField::new(size);
    vector_field[(0, 0, 0)] = Vector3::repeat(1.0);
    vector_field[(1, 0, 0)] = Vector3::repeat(2.0);
    vector_field[(0, 1, 0)] = Vector3::repeat(3.0);
    vector_field[(1, 1, 0)] = Vector3::repeat(4.0);

    let mut scalar_field = ScalarField::ones(size);
    scalar_field[(0, 0, 0)] = 2.0;
    scalar_field[(1, 0, 0)] = 3.0;
    scalar_field[(0, 1, 0)] = 4.0;
    scalar_field[(1, 1, 0)] = 5.0;

    let result = vector_field / &scalar_field;
    assert!((result[(0usize, 0, 0)] - Vector3::repeat(0.5)).norm() < 1e-10);
    assert!((result[(1usize, 0, 0)] - Vector3::repeat(2.0 / 3.0)).norm() < 1e-10);
    assert!((result[(0usize, 1, 0)] - Vector3::repeat(0.75)).norm() < 1e-10);
    assert!((result[(1usize, 1, 0)] - Vector3::repeat(0.8)).norm() < 1e-10);
}

impl<'a> Div<&'a ScalarField> for VectorField {
    type Output = VectorField;

    fn div(mut self, rhs: &'a ScalarField) -> VectorField {
        debug_assert!(self.size() == rhs.size());
        for (element_self, element_other) in self.vectors.iter_mut().zip(rhs.scalars().iter()) {
            *element_self /= *element_other;
        }
        self
    }
}

#[test]
fn mul_vector_field_vector_field() {
    let size = Vector3::new(2, 2, 1);

    let mut vector_field = VectorField::new(size);
    vector_field[(0, 0, 0)] = Vector3::repeat(1.0);
    vector_field[(1, 0, 0)] = Vector3::repeat(2.0);
    vector_field[(0, 1, 0)] = Vector3::repeat(3.0);
    vector_field[(1, 1, 0)] = Vector3::repeat(4.0);

    let mut vector_field_2 = VectorField::new(size);
    vector_field_2[(0, 0, 0)] = Vector3::repeat(4.0);
    vector_field_2[(1, 0, 0)] = Vector3::repeat(3.0);
    vector_field_2[(0, 1, 0)] = Vector3::repeat(2.0);
    vector_field_2[(1, 1, 0)] = Vector3::repeat(1.0);

    assert!(((&vector_field * &vector_field_2) - 60.0).abs() < 1e-10);
}

impl<'a, 'b> Mul<&'a VectorField> for &'b VectorField {
    type Output = f32;

    fn mul(self, rhs: &'a VectorField) -> Self::Output {
        debug_assert!(self.size() == rhs.size());
        self.vectors.iter().zip(rhs.vectors.iter())
            .map(|(element_self, element_other)| element_self.dot(element_other))
            .sum()
    }
}

#[test]
fn mul_f32_vector_field() {
    let mut field: VectorField = VectorField::new(Vector3::new(2, 2, 1));
    field[(0, 0, 0)] = Vector3::repeat(0.5);
    field[(1, 0, 0)] = Vector3::repeat(1.0);
    field[(0, 1, 0)] = Vector3::repeat(1.5);
    field[(1, 1, 0)] = Vector3::repeat(2.0);

    let result = 2.0 * field;

    assert!((result[(0usize, 0, 0)] - Vector3::repeat(1.0)).norm() < 1e-10);
    assert!((result[(1usize, 0, 0)] - Vector3::repeat(2.0)).norm() < 1e-10);
    assert!((result[(0usize, 1, 0)] - Vector3::repeat(3.0)).norm() < 1e-10);
    assert!((result[(1usize, 1, 0)] - Vector3::repeat(4.0)).norm() < 1e-10);
}

impl Mul<VectorField> for f32 {
    type Output = VectorField;

    fn mul(self, mut rhs: VectorField) -> Self::Output {
        for element in rhs.vectors.iter_mut() {
            *element *= self;
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

#[test]
fn neg_vector_field() {
    let mut field: VectorField = VectorField::new(Vector3::new(2, 2, 1));
    field[(0, 0, 0)] = Vector3::repeat(0.5);
    field[(1, 0, 0)] = Vector3::repeat(1.0);
    field[(0, 1, 0)] = Vector3::repeat(1.5);
    field[(1, 1, 0)] = Vector3::repeat(2.0);

    let result = -field;

    assert!((result[(0usize, 0, 0)] - Vector3::repeat(-0.5)).norm() < 1e-10);
    assert!((result[(1usize, 0, 0)] - Vector3::repeat(-1.0)).norm() < 1e-10);
    assert!((result[(0usize, 1, 0)] - Vector3::repeat(-1.5)).norm() < 1e-10);
    assert!((result[(1usize, 1, 0)] - Vector3::repeat(-2.0)).norm() < 1e-10);
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
