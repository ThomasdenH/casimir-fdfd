use nalgebra::*;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign};

#[derive(PartialEq, Clone, Debug)]
pub struct ScalarField {
    size: Vector3<usize>,
    scalars: DVector<f32>,
}

impl ScalarField {
    pub fn ones(size: Vector3<usize>) -> ScalarField {
        ScalarField {
            size,
            scalars: DVector::from_element(size.x * size.y * size.z, 1.0),
        }
    }

    pub fn set(&mut self, x: isize, y: isize, z: isize, value: f32) {
        if x < 0
            || x as usize >= self.size.x
            || y < 0
            || y as usize >= self.size.y
            || z < 0
            || z as usize >= self.size.z
        {
            return;
        }
        self.scalars[x as usize + self.size.x * (y as usize + self.size.y * (z as usize))] = value;
    }

    /// Perform the action 1.0 / M_i on each element.
    pub fn multiplicative_invert(&mut self) {
        self.scalars = self.scalars.map(|a| 1.0 / a);
    }

    pub fn size(&self) -> Vector3<usize> {
        self.size
    }

    pub fn len(&self) -> usize {
        self.scalars.len()
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
