use nalgebra::*;
use scalarfield::ScalarField;
use std::f64::consts::PI;
use std::ops::Mul;
use vectorfield::VectorField;

/// A representation of the green function operator
struct A<'a> {
    scalar_primary: &'a ScalarField,
    scalar_secundary: &'a ScalarField,
    frequency_2: f64,
}

impl<'b> A<'b> {
    fn new<'a>(
        frequency: f64,
        scalar_primary: &'a ScalarField,
        scalar_secundary: &'a ScalarField,
    ) -> A<'a> {
        A {
            scalar_primary,
            scalar_secundary,
            frequency_2: frequency * frequency,
        }
    }
}

impl<'a, 'b> Mul<VectorField> for &'a A<'b> {
    type Output = VectorField;

    fn mul(self, x: VectorField) -> Self::Output {
        let curl_part = (x.clone().curl_positive() * self.scalar_secundary).curl_negative();
        let scalar_part = self.frequency_2 * (x * self.scalar_primary);
        curl_part + &scalar_part
    }
}

pub fn stress_tensor(
    frequency: f64,
    point: Point3<usize>,
    electric_field: &ScalarField,
    magnetic_field: &ScalarField,
    inv_electric_field: &ScalarField,
    inv_magnetic_field: &ScalarField,
    size: Vector3<usize>,
) -> Matrix3<f64> {
    let electric_tensor = green_tensor(frequency, point, electric_field, inv_magnetic_field, size);
    let magnetic_tensor = green_tensor(frequency, point, magnetic_field, inv_electric_field, size);

    frequency * frequency / PI
        * (magnetic_field[point]
            * (magnetic_tensor - Matrix3::from_diagonal_element(0.5) * magnetic_tensor.trace())
            + electric_field[point]
                * (electric_tensor - Matrix3::from_diagonal_element(0.5) * electric_tensor.trace()))
}

pub fn green_tensor(
    frequency: f64,
    point: Point3<usize>,
    scalar_primary: &ScalarField,
    scalar_secundary: &ScalarField,
    size: Vector3<usize>,
) -> Matrix3<f64> {
    Matrix3::from_columns(&[
        green_function(
            frequency,
            point,
            Vector3::new(1.0, 0.0, 0.0),
            scalar_primary,
            scalar_secundary,
            size,
        ),
        green_function(
            frequency,
            point,
            Vector3::new(0.0, 1.0, 0.0),
            scalar_primary,
            scalar_secundary,
            size,
        ),
        green_function(
            frequency,
            point,
            Vector3::new(0.0, 0.0, 1.0),
            scalar_primary,
            scalar_secundary,
            size,
        ),
    ])
}

pub fn green_function(
    frequency: f64,
    point: Point3<usize>,
    polarization: Vector3<f64>,
    scalar_primary: &ScalarField,
    scalar_secundary: &ScalarField,
    size: Vector3<usize>,
) -> Vector3<f64> {
    // The delta function right hand side
    let mut b = VectorField::new(size.x, size.y, size.z);
    b[point] = polarization;

    // The operator
    let a = A::new(frequency, scalar_primary, scalar_secundary);

    let mut x = VectorField::new(size.x, size.y, size.z);
    let mut r = &b - &a * x.clone();
    let mut p = r.clone();
    let mut rsold = &r * &r;

    let n = size.x * size.y * size.z;

    for i in 0..n {
        let a_p = &a * p.clone();
        let alpha = rsold / (&p * &a_p);
        x += alpha * &p;
        r -= alpha * &a_p;
        let rsnew = &r * &r;
        if rsnew.sqrt() / (n as f64) < 1e-8 {
            break;
        }
        p = (rsnew / rsold) * p + &r;
        rsold = rsnew;
    }
    x[point]
}
