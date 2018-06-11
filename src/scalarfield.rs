use nalgebra::*;
use std::ops::{Add, Mul, Index};

#[derive(PartialEq, Clone, Debug)]
pub struct ScalarField {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub scalars: DVector<f64>,
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

    pub fn set(&mut self, x: isize, y: isize, z: isize, value: f64) {
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

impl Mul<f64> for ScalarField {
    type Output = ScalarField;

    fn mul(mut self, rhs: f64) -> Self::Output {
        self.scalars *= rhs;
        self
    }
}

impl Mul<ScalarField> for f64 {
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

impl Index<(usize, usize, usize)> for ScalarField {
    type Output = f64;

    fn index(&self, index: (usize, usize, usize)) -> &f64 {
        let (x, y, z) = index;
        assert!(x < self.nx && y < self.ny && z < self.nz);
        &self.scalars[x as usize + self.nx * (y as usize + self.ny * (z as usize))]
    }
}

impl Index<(isize, isize, isize)> for ScalarField {
    type Output = f64;

    fn index(&self, index: (isize, isize, isize)) -> &f64 {
        let (x, y, z) = index;
        if x < 0 || x as usize >= self.nx
            || y < 0 || y as usize >= self.ny
            || z < 0 || z as usize >= self.nz
            {
                &1.0
            } else {
            &self.scalars[x as usize + self.nx * (y as usize + self.ny * (z as usize))]
        }
    }
}

impl Index<Point3<usize>> for ScalarField {
    type Output = f64;

    fn index(&self, index: Point3<usize>) -> &f64 {
        &self[(index.x, index.y, index.z)]
    }
}

impl Index<Point3<isize>> for ScalarField {
    type Output = f64;

    fn index(&self, index: Point3<isize>) -> &f64 {
        &self[(index.x, index.y, index.z)]
    }
}
