use nalgebra::*;
use std::ops::{Add, AddAssign, Div, Index, IndexMut, Mul, Neg, Sub, SubAssign};

use fields::ScalarField;

/// Represents a vector field. It has a certain size, with a vector on every coordinate on the grid.
#[derive(PartialEq, Clone, Debug)]
pub struct VectorField {
    size: Vector3<usize>,
    vectors: DVector<Vector3<f32>>,
}

impl VectorField {
    /// Create a new `VectorField` filled with (0, 0, 0) with a certain `size`.
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

    /// Same as `curl_positive`, but uses no allocations and recycles a `VectorField` instead.
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

    /// Calculates the negative curl of a vector field. The new positions will be the centers of the
    /// Yee grid planes formed between four grid edges.
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

    /// Same as `curl_negative`, but uses no allocations and recycles a `VectorField` instead.
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

    /// Returns the size in three dimensions of this field.
    pub fn size(&self) -> Vector3<usize> {
        self.size
    }

    /// Returns the number of vectors in this field.
    pub fn len(&self) -> usize {
        self.vectors.len()
    }

    /// Clone the components of this vector field to the vector field `other`.
    pub fn clone_to(&self, mut other: VectorField) -> VectorField {
        debug_assert!(self.size() == other.size());
        for i in 0..self.vectors.len() {
            other[i] = self[i];
        }
        other
    }

    /// Returns the vectors in this field in an immutable way.
    pub fn vectors(&self) -> &DVector<Vector3<f32>> {
        &self.vectors
    }

    /// Returns the vectors inside this field in a mutable way.
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

impl<'a, 'b> Mul<&'a VectorField> for &'b VectorField {
    type Output = f32;

    fn mul(self, rhs: &'a VectorField) -> Self::Output {
        debug_assert!(self.size() == rhs.size());
        self.vectors
            .iter()
            .zip(rhs.vectors.iter())
            .map(|(element_self, element_other)| element_self.dot(element_other))
            .sum()
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

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

    #[test]
    fn zero_curl_field() {
        let mut field = VectorField::new(Vector3::new(5, 5, 5));
        for vector in field.vectors_mut().iter_mut() {
            *vector = Vector3::new(3.0, 9.0, 1.0);
        }
        // Remove mutability
        let field = field;

        // Negative curl
        for vector in field.curl_negative().vectors().iter() {
            assert_approx_eq!(vector.norm(), 0.0, 1e-10);
        }

        // Positive curl
        for vector in field.curl_positive().vectors().iter() {
            assert_approx_eq!(vector.norm(), 0.0, 1e-10);
        }

        // Negative curl to
        let negative_curl = field.curl_negative_to(field.clone());
        for vector in negative_curl.vectors().iter() {
            assert_approx_eq!(vector.norm(), 0.0, 1e-10);
        }

        // Negative curl to
        let positive_curl = field.curl_positive_to(field.clone());
        for vector in positive_curl.vectors().iter() {
            assert_approx_eq!(vector.norm(), 0.0, 1e-10);
        }
    }

    #[test]
    fn constant_curl_field() {
        // Create a circling vector field
        let mut field = VectorField::new(Vector3::new(10, 10, 10));
        for x in 0..10 {
            for y in 0..10 {
                for z in 0..10 {
                    field[(x, y, z)] = Vector3::new(
                        (y as f32) - 5.0,
                        5.0 - (x as f32),
                        0.0
                    );
                }
            }
        }

        // Compute the result
        let negative_curl = field.curl_negative();
        let positive_curl = field.curl_positive();
        let negative_curl_to = field.curl_negative_to(field.clone());
        let positive_curl_to = field.curl_positive_to(field.clone());

        let expected = Vector3::new(0.0, 0.0, -2.0);
        let error = Vector3::repeat(1e-5);
        // Because of boundary conditions, avoid edges
        for x in 1usize..9 {
            for y in 1usize..9 {
                for z in 0usize..10 {
                    assert_approx_eq!(negative_curl[(x, y, z)], expected, error);
                    assert_approx_eq!(positive_curl[(x, y, z)], expected, error);
                    assert_approx_eq!(negative_curl_to[(x, y, z)], expected, error);
                    assert_approx_eq!(positive_curl_to[(x, y, z)], expected, error);
                }
            }
        }
    }

    #[test]
    fn custom_curl_field_1() {
        // Create a circling vector field
        let mut field = VectorField::new(Vector3::new(10, 10, 10));
        for x in 0..10 {
            for y in 0..10 {
                for z in 0..10 {
                    field[(x, y, z)] = Vector3::new(
                        0.0,
                        -((x as f32) - 5.0).powi(2),
                        0.0
                    );
                }
            }
        }

        // Compute the result
        let negative_curl = field.curl_negative();
        let positive_curl = field.curl_positive();
        let negative_curl_to = field.curl_negative_to(field.clone());
        let positive_curl_to = field.curl_positive_to(field.clone());

        // Error can be quite large
        let error = Vector3::repeat(1.5);
        // Because of boundary conditions, avoid edges
        for x in 1usize..9 {
            for y in 1usize..9 {
                for z in 0usize..10 {
                    let expected = Vector3::new(0.0, 0.0, -2.0 * ((x as f32) - 5.0));
                    assert_approx_eq!(negative_curl[(x, y, z)], expected, error);
                    assert_approx_eq!(positive_curl[(x, y, z)], expected, error);
                    assert_approx_eq!(negative_curl_to[(x, y, z)], expected, error);
                    assert_approx_eq!(positive_curl_to[(x, y, z)], expected, error);
                }
            }
        }
    }

    #[test]
    fn custom_curl_field_2() {
        // Create a circling vector field
        let mut field = VectorField::new(Vector3::new(10, 10, 10));
        for x in 0..10 {
            for y in 0..10 {
                for z in 0..10 {
                    field[(x, y, z)] = Vector3::new(
                        (x as f32 / 40.0 / PI).sin(),
                        (z as f32 / 40.0 / PI).cos(),
                        -(y as f32 / 40.0 / PI).cos()
                    );
                }
            }
        }

        // Compute the result
        let negative_curl = field.curl_negative();
        let positive_curl = field.curl_positive();
        let negative_curl_to = field.curl_negative_to(field.clone());
        let positive_curl_to = field.curl_positive_to(field.clone());

        // Error can be quite large
        let error = Vector3::repeat(1e-2);
        // Because of boundary conditions, avoid edges
        for x in 1usize..9 {
            for y in 1usize..9 {
                for z in 0usize..10 {
                    let expected = Vector3::new(
                        (z as f32 / 40.0 / PI).sin() * (y as f32 / 40.0 / PI).sin() / 80.0 / PI,
                        0.0,
                        0.0
                    );
                    assert_approx_eq!(negative_curl[(x, y, z)], expected, error);
                    assert_approx_eq!(positive_curl[(x, y, z)], expected, error);
                    assert_approx_eq!(negative_curl_to[(x, y, z)], expected, error);
                    assert_approx_eq!(positive_curl_to[(x, y, z)], expected, error);
                }
            }
        }
    }

    #[test]
    fn index() {
        let mut field = VectorField::new(Vector3::new(4, 3, 2));
        let error = Vector3::repeat(1e-10);

        // Check mutability and indexing in both formats
        let vec_a = Vector3::repeat(5.0);
        field[(0usize, 0, 0)] = vec_a;
        assert_approx_eq!(field[Point3::new(0usize, 0, 0)], vec_a, error);

        let vec_b = Vector3::repeat(3.0);
        field[Point3::new(1usize, 0, 0)] = vec_b;
        assert_approx_eq!(field[(1usize, 0, 0)], vec_b, error);

        // Check boundary conditions (should repeat outside domain)
        assert_approx_eq!(field[Point3::new(-3isize, 0, 0)], field[Point3::new(1isize, 0, 0)], error);
        assert_approx_eq!(field[(1000isize, 1000, 10)], field[(0usize, 1, 0)], error);

        // Check coordinate conversion
        assert_approx_eq!(field[0], field[(0usize, 0, 0)], error);
        assert_approx_eq!(field[17], field[(1usize, 1, 1)], error);
    }
}
