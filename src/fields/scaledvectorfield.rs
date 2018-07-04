use fields::VectorField;
use nalgebra::*;
use std::ops::{AddAssign, Mul, SubAssign};

/// Represents the result of a scalar times a vector field reference. If it is added to a vector
/// field later, it requires no owned source. Internally the multiplication is postponed until the
/// next operation.
pub struct ScaledVectorField<'a> {
    field: &'a VectorField,
    scalar: f32,
}

impl<'a> ScaledVectorField<'a> {
    pub fn size(&self) -> Vector3<usize> {
        self.field.size()
    }

    pub fn len(&self) -> usize {
        self.field.len()
    }

    pub fn field(&self) -> &VectorField {
        self.field
    }
}

impl<'a> Mul<f32> for &'a VectorField {
    type Output = ScaledVectorField<'a>;

    fn mul(self, rhs: f32) -> ScaledVectorField<'a> {
        ScaledVectorField {
            field: self,
            scalar: rhs,
        }
    }
}

impl<'a> Mul<&'a VectorField> for f32 {
    type Output = ScaledVectorField<'a>;

    fn mul(self, rhs: &'a VectorField) -> ScaledVectorField<'a> {
        ScaledVectorField {
            field: rhs,
            scalar: self,
        }
    }
}

impl<'a, 'b> AddAssign<&'b ScaledVectorField<'a>> for VectorField {
    fn add_assign(&mut self, other: &'b ScaledVectorField) {
        debug_assert!(self.size() == other.size());
        for (self_element, other_element) in self.vectors_mut()
            .iter_mut()
            .zip(other.field().vectors().iter())
        {
            *self_element += other.scalar * other_element
        }
    }
}

impl<'a, 'b> SubAssign<&'b ScaledVectorField<'a>> for VectorField {
    fn sub_assign(&mut self, other: &'b ScaledVectorField) {
        debug_assert!(self.size() == other.size());
        for (self_element, other_element) in self.vectors_mut()
            .iter_mut()
            .zip(other.field().vectors().iter())
        {
            *self_element -= other.scalar * other_element
        }
    }
}

impl<'a> AddAssign<ScaledVectorField<'a>> for VectorField {
    fn add_assign(&mut self, other: ScaledVectorField) {
        *self += &other;
    }
}

impl<'a> SubAssign<ScaledVectorField<'a>> for VectorField {
    fn sub_assign(&mut self, other: ScaledVectorField) {
        *self -= &other;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mul_f32_vector_field_ref() {
        let size = Vector3::new(4, 5, 4);
        let field_a = VectorField::new(size);
        let scaled: ScaledVectorField = 3.0 * &field_a;
        assert_eq!(scaled.size(), size);
        assert_eq!(scaled.len(), field_a.len());
    }

    #[test]
    fn test_sub_assign_scaled_vector_field() {
        let mut field_a: VectorField = VectorField::new(Vector3::new(2, 2, 1));
        field_a[(0, 0, 0)] = Vector3::repeat(1.0);
        field_a[(1, 0, 0)] = Vector3::repeat(2.0);
        field_a[(0, 1, 0)] = Vector3::repeat(3.0);
        field_a[(1, 1, 0)] = Vector3::repeat(4.0);

        let cloned_a = field_a.clone();
        let field_b: ScaledVectorField = -2.0 * &cloned_a;

        field_a -= &field_b;

        assert!((field_a[(0usize, 0, 0)] - Vector3::repeat(3.0)).norm() < 1e-10);
        assert!((field_a[(1usize, 0, 0)] - Vector3::repeat(6.0)).norm() < 1e-10);
        assert!((field_a[(0usize, 1, 0)] - Vector3::repeat(9.0)).norm() < 1e-10);
        assert!((field_a[(1usize, 1, 0)] - Vector3::repeat(12.0)).norm() < 1e-10);
    }

    #[test]
    #[should_panic]
    fn test_sub_assign_scaled_vector_field_different_sizes() {
        let mut field_a: VectorField = VectorField::new(Vector3::new(2, 2, 1));
        let field_b = VectorField::new(Vector3::new(1, 2, 1));
        let scaled_field_b: ScaledVectorField = -2.0 * &field_b;
        field_a -= &scaled_field_b;
    }

    #[test]
    fn test_add_assign_scaled_vector_field() {
        let mut field_a: VectorField = VectorField::new(Vector3::new(2, 2, 1));
        field_a[(0, 0, 0)] = Vector3::repeat(1.0);
        field_a[(1, 0, 0)] = Vector3::repeat(2.0);
        field_a[(0, 1, 0)] = Vector3::repeat(3.0);
        field_a[(1, 1, 0)] = Vector3::repeat(4.0);

        let cloned_a = field_a.clone();
        let field_b: ScaledVectorField = 2.0 * &cloned_a;

        field_a += &field_b;

        assert!((field_a[(0usize, 0, 0)] - Vector3::repeat(3.0)).norm() < 1e-10);
        assert!((field_a[(1usize, 0, 0)] - Vector3::repeat(6.0)).norm() < 1e-10);
        assert!((field_a[(0usize, 1, 0)] - Vector3::repeat(9.0)).norm() < 1e-10);
        assert!((field_a[(1usize, 1, 0)] - Vector3::repeat(12.0)).norm() < 1e-10);
    }

    #[test]
    #[should_panic]
    fn test_add_assign_scaled_vector_field_different_sizes() {
        let mut field_a: VectorField = VectorField::new(Vector3::new(2, 2, 1));
        let field_b = VectorField::new(Vector3::new(1, 2, 1));
        let scaled_field_b: ScaledVectorField = -2.0 * &field_b;
        field_a += &scaled_field_b;
    }
}
