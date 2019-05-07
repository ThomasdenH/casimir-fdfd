use nalgebra::*;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign};

/// On a grid of a certain size, a scalar is located on every coordinate.
#[derive(PartialEq, Clone, Debug)]
pub struct ScalarField {
    size: Vector3<usize>,
    scalars: DVector<f32>,
}

impl ScalarField {
    /// Constructs a new `ScalarField` filled with the value 1.
    pub fn ones(size: Vector3<usize>) -> ScalarField {
        ScalarField {
            size,
            scalars: DVector::from_element(size.x * size.y * size.z, 1.0),
        }
    }

    /// Perform the action 1.0 / M_i on each element.
    pub fn multiplicative_invert(&mut self) {
        self.scalars = self.scalars.map(|a| 1.0 / a);
    }

    /// Get the three-dimensional size of this field.
    pub fn size(&self) -> Vector3<usize> {
        self.size
    }

    /// Get the number of scalars in this field.
    pub fn len(&self) -> usize {
        self.scalars.len()
    }

    /// Returns true if there are no scaalrs in this field.
    pub fn is_empty(&self) -> bool {
        self.scalars.is_empty()
    }

    /// Get the scalars in this field.
    pub fn scalars(&self) -> &DVector<f32> {
        &self.scalars
    }
}

impl Mul<f32> for ScalarField {
    type Output = ScalarField;

    fn mul(mut self, rhs: f32) -> Self::Output {
        self *= rhs;
        self
    }
}

impl Mul<ScalarField> for f32 {
    type Output = ScalarField;

    fn mul(self, mut rhs: ScalarField) -> Self::Output {
        rhs *= self;
        rhs
    }
}

impl<'a> Add<&'a ScalarField> for ScalarField {
    type Output = ScalarField;

    fn add(mut self, rhs: &'a ScalarField) -> ScalarField {
        debug_assert!(self.size == rhs.size);
        self += rhs;
        self
    }
}

impl<'a> Sub<&'a ScalarField> for ScalarField {
    type Output = ScalarField;

    fn sub(mut self, rhs: &'a ScalarField) -> ScalarField {
        debug_assert!(self.size == rhs.size);
        self -= rhs;
        self
    }
}

impl<'a> AddAssign<&'a ScalarField> for ScalarField {
    fn add_assign(&mut self, rhs: &'a ScalarField) {
        debug_assert!(self.size == rhs.size);
        self.scalars += &rhs.scalars;
    }
}

impl<'a> SubAssign<&'a ScalarField> for ScalarField {
    fn sub_assign(&mut self, rhs: &'a ScalarField) {
        debug_assert!(self.size == rhs.size);
        self.scalars -= &rhs.scalars;
    }
}

impl<'a> MulAssign<f32> for ScalarField {
    fn mul_assign(&mut self, rhs: f32) {
        self.scalars *= rhs;
    }
}

impl Index<(usize, usize, usize)> for ScalarField {
    type Output = f32;

    fn index(&self, index: (usize, usize, usize)) -> &f32 {
        let (x, y, z) = index;
        debug_assert!(x < self.size.x && y < self.size.y && z < self.size.z);
        &self.scalars[x as usize + self.size.x * (y as usize + self.size.y * (z as usize))]
    }
}

impl IndexMut<(usize, usize, usize)> for ScalarField {
    fn index_mut(&mut self, index: (usize, usize, usize)) -> &mut f32 {
        let (x, y, z) = index;
        debug_assert!(x < self.size.x && y < self.size.y && z < self.size.z);
        &mut self.scalars[x as usize + self.size.x * (y as usize + self.size.y * (z as usize))]
    }
}

impl IndexMut<Point3<usize>> for ScalarField {
    fn index_mut<'a>(&'a mut self, index: Point3<usize>) -> &'a mut f32 {
        &mut self[(index.x, index.y, index.z)]
    }
}

impl Index<(isize, isize, isize)> for ScalarField {
    type Output = f32;

    fn index(&self, index: (isize, isize, isize)) -> &f32 {
        let (x, y, z) = index;
        if x < 0
            || x as usize >= self.size.x
            || y < 0
            || y as usize >= self.size.y
            || z < 0
            || z as usize >= self.size.z
        {
            &1.0
        } else {
            &self.scalars[x as usize + self.size.x * (y as usize + self.size.y * (z as usize))]
        }
    }
}

impl Index<Point3<usize>> for ScalarField {
    type Output = f32;

    fn index(&self, index: Point3<usize>) -> &f32 {
        &self[(index.x, index.y, index.z)]
    }
}

impl Index<Point3<isize>> for ScalarField {
    type Output = f32;

    fn index(&self, index: Point3<isize>) -> &f32 {
        &self[(index.x, index.y, index.z)]
    }
}

impl Index<usize> for ScalarField {
    type Output = f32;

    fn index(&self, index: usize) -> &f32 {
        &self.scalars[index]
    }
}

impl IndexMut<usize> for ScalarField {
    fn index_mut(&mut self, index: usize) -> &mut f32 {
        &mut self.scalars[index]
    }
}

#[cfg(test)]
mod tests {
    use fields::ScalarField;
    use nalgebra::*;

    #[test]
    fn test_ones_and_properties() {
        let size = Vector3::new(4, 4, 4);
        let field = ScalarField::ones(size);
        assert_eq!(field.size(), size);
        assert_eq!(field.len(), 4 * 4 * 4);
        for value in field.scalars().iter() {
            assert_approx_eq!(value, 1.0, 1e-10);
        }
    }

    #[test]
    fn multiplicative_invert() {
        let mut field = ScalarField::ones(Vector3::new(2, 2, 1));
        field[(0, 0, 0)] = 2.0;
        field[(1, 0, 0)] = 10.0;
        field[(0, 1, 0)] = 0.1;
        field[(1, 1, 0)] = 0.2;
        field.multiplicative_invert();
        assert_approx_eq!(field[(0usize, 0, 0)], 0.5, 1e-10);
        assert_approx_eq!(field[(1usize, 0, 0)], 0.1, 1e-10);
        assert_approx_eq!(field[(0usize, 1, 0)], 10.0, 1e-10);
        assert_approx_eq!(field[(1usize, 1, 0)], 5.0, 1e-10);
    }

    #[test]
    fn mul_scalar_field_f32() {
        let size = Vector3::new(4, 4, 4);
        let field = ScalarField::ones(size) * 3.0;
        for value in field.scalars().iter() {
            assert_approx_eq!(value, 3.0, 1e-10);
        }
    }

    #[test]
    fn mul_assign_scalar_field_f32() {
        let size = Vector3::new(3, 2, 5);
        let mut field = ScalarField::ones(size);
        field *= 9.0;
        for value in field.scalars().iter() {
            assert_approx_eq!(value, 9.0, 1e-10);
        }
    }

    #[test]
    fn mul_f32_scalar_field() {
        let size = Vector3::new(4, 4, 4);
        let field = 4.0 * ScalarField::ones(size);
        for value in field.scalars().iter() {
            assert_approx_eq!(value, 4.0, 1e-10);
        }
    }

    #[test]
    fn add_scalar_field_scalar_field() {
        let size = Vector3::new(4, 4, 4);
        let field1 = 3.0 * ScalarField::ones(size);
        let field2 = 4.0 * ScalarField::ones(size);
        for value in (field1 + &field2).scalars().iter() {
            assert_approx_eq!(value, 7.0, 1e-10);
        }
    }

    #[test]
    #[should_panic]
    fn add_scalar_field_scalar_field_different_sizes() {
        let field1 = 3.0 * ScalarField::ones(Vector3::new(4, 3, 2));
        let field2 = 1.5 * ScalarField::ones(Vector3::new(2, 3, 4));
        let _ = field1 + &field2;
    }

    #[test]
    fn sub_scalar_field_scalar_field() {
        let size = Vector3::new(4, 3, 4);
        let field1 = 3.0 * ScalarField::ones(size);
        let field2 = 4.0 * ScalarField::ones(size);
        for value in (field1 - &field2).scalars().iter() {
            assert_approx_eq!(value, -1.0, 1e-10);
        }
    }

    #[test]
    #[should_panic]
    fn sub_scalar_field_scalar_field_different_sizes() {
        let field1 = 3.0 * ScalarField::ones(Vector3::new(4, 3, 2));
        let field2 = 1.5 * ScalarField::ones(Vector3::new(2, 3, 4));
        let _ = field1 - &field2;
    }

    #[test]
    fn add_assign_scalar_field_scalar_field() {
        let size = Vector3::new(4, 4, 4);
        let mut field1 = 3.0 * ScalarField::ones(size);
        let field2 = 4.0 * ScalarField::ones(size);
        field1 += &field2;
        for value in field1.scalars().iter() {
            assert_approx_eq!(value, 7.0, 1e-10);
        }
    }

    #[test]
    #[should_panic]
    fn add_assign_scalar_field_scalar_field_different_sizes() {
        let mut field1 = 3.0 * ScalarField::ones(Vector3::new(4, 3, 2));
        let field2 = 1.5 * ScalarField::ones(Vector3::new(2, 3, 4));
        field1 += &field2;
    }

    #[test]
    fn sub_assign_scalar_field_scalar_field() {
        let size = Vector3::new(4, 3, 4);
        let mut field1 = 3.0 * ScalarField::ones(size);
        let field2 = 4.0 * ScalarField::ones(size);
        field1 -= &field2;
        for value in field1.scalars().iter() {
            assert_approx_eq!(value, -1.0, 1e-10);
        }
    }

    #[test]
    #[should_panic]
    fn sub_assign_scalar_field_scalar_field_different_sizes() {
        let mut field1 = 3.0 * ScalarField::ones(Vector3::new(4, 3, 2));
        let field2 = 1.5 * ScalarField::ones(Vector3::new(2, 3, 4));
        field1 -= &field2;
    }

    #[test]
    fn index() {
        let mut field = 3.0 * ScalarField::ones(Vector3::new(4, 3, 2));

        // Check mutability and indexing in both formats
        field[(0usize, 0, 0)] = 5.0;
        assert_approx_eq!(field[Point3::new(0usize, 0, 0)], 5.0, 1e-10);

        field[Point3::new(1usize, 0, 0)] = 2.5;
        assert_approx_eq!(field[(1usize, 0, 0)], 2.5, 1e-10);

        // Check boundary conditions (should be 1.0 outsize domain)
        assert_approx_eq!(field[Point3::new(-3isize, 0, 0)], 1.0, 1e-10);
        assert_approx_eq!(field[(1000isize, 1000, 10)], 1.0, 1e-10);

        // Check coordinate conversion
        assert_approx_eq!(field[0], field[(0usize, 0, 0)], 1e-10);
        assert_approx_eq!(field[17], field[(1usize, 1, 1)], 1e-10);
    }
}
