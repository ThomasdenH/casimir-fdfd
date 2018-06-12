use nalgebra::*;
use std::ops::{AddAssign, Mul, SubAssign};
use vectorfield::VectorField;

/// Represents the result of a scalar times a vector field reference. If it is added to a vector
/// field later, it requires no owned source. Internally the addition is postponed until the next
/// operation.
pub struct ScaledVectorField<'a> {
    field: &'a VectorField,
    scalar: f64,
}

impl<'a> ScaledVectorField<'a> {
    pub fn size(&self) -> Vector3<usize> {
        self.field.size()
    }
}

impl<'a> Mul<f64> for &'a VectorField {
    type Output = ScaledVectorField<'a>;

    fn mul(self, rhs: f64) -> ScaledVectorField<'a> {
        ScaledVectorField {
            field: self,
            scalar: rhs,
        }
    }
}

impl<'a> Mul<&'a VectorField> for f64 {
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
        for x in 0..self.size().x {
            for y in 0..self.size().y {
                for z in 0..self.size().z {
                    self[(x, y, z)] += other.scalar * other.field[(x, y, z)];
                }
            }
        }
    }
}

impl<'a, 'b> SubAssign<&'b ScaledVectorField<'a>> for VectorField {
    fn sub_assign(&mut self, other: &'b ScaledVectorField) {
        debug_assert!(self.size() == other.size());
        for x in 0..self.size().x {
            for y in 0..self.size().y {
                for z in 0..self.size().z {
                    self[(x, y, z)] -= other.scalar * other.field[(x, y, z)];
                }
            }
        }
    }
}

impl<'a> AddAssign<ScaledVectorField<'a>> for VectorField {
    fn add_assign(&mut self, other: ScaledVectorField) {
        debug_assert!(self.size() == other.size());
        for x in 0..self.size().x {
            for y in 0..self.size().y {
                for z in 0..self.size().z {
                    self[(x, y, z)] += other.scalar * other.field[(x, y, z)];
                }
            }
        }
    }
}

impl<'a> SubAssign<ScaledVectorField<'a>> for VectorField {
    fn sub_assign(&mut self, other: ScaledVectorField) {
        debug_assert!(self.size() == other.size());
        for x in 0..self.size().x {
            for y in 0..self.size().y {
                for z in 0..self.size().z {
                    self[(x, y, z)] -= other.scalar * other.field[(x, y, z)];
                }
            }
        }
    }
}
