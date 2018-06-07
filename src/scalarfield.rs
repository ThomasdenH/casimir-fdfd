use nalgebra::*;
use std::ops::{Add, Mul};

#[derive(PartialEq, Clone, Debug)]
pub struct ScalarField {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub scalars: DVector<f32>,
}

impl ScalarField {
    pub fn ones(nx: usize, ny: usize, nz: usize) -> ScalarField {
        ScalarField {
            nx,
            ny,
            nz,
            scalars: DVector::from_element(nx * ny * nz, 1.0),
        }
    }

    pub fn at(&self, x: isize, y: isize, z: isize) -> f32 {
        if x < 0 || x as usize >= self.nx || y < 0 || y as usize >= self.ny || z < 0
            || z as usize >= self.nz
        {
            1.0
        } else {
            self.scalars[x as usize + self.nx * (y as usize + self.ny * (z as usize))]
        }
    }

    pub fn set(&mut self, x: isize, y: isize, z: isize, value: f32) {
        if x < 0 || x as usize >= self.nx || y < 0 || y as usize >= self.ny || z < 0
            || z as usize >= self.nz
        {
            return;
        }
        self.scalars[x as usize + self.nx * (y as usize + self.ny * (z as usize))] = value;
    }

    /// Perform the action 1.0 / M_i on each element.
    pub fn multiplicative_invert(&mut self) {
        self.scalars = self.scalars.map(|a| 1.0 / a);
    }
}

impl Mul<f32> for ScalarField {
    type Output = ScalarField;

    fn mul(mut self, rhs: f32) -> Self::Output {
        self.scalars *= rhs;
        self
    }
}

impl Mul<ScalarField> for f32 {
    type Output = ScalarField;

    fn mul(self, mut rhs: ScalarField) -> Self::Output {
        rhs.scalars *= self;
        rhs
    }
}

impl<'a> Add<&'a ScalarField> for ScalarField {
    type Output = ScalarField;

    fn add(mut self, rhs: &'a ScalarField) -> ScalarField {
        assert!(self.nx == rhs.nx && self.ny == rhs.ny && self.nz == rhs.nz);
        self.scalars += &rhs.scalars;
        self
    }
}