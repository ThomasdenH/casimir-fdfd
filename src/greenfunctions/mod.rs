use nalgebra::*;
use scalarfield::ScalarField;
use std::f64::consts::PI;
use vectorfield::VectorField;

mod operator;

use greenfunctions::operator::{OperatorType, Operator};

pub fn stress_tensor(
    frequency: f64,
    point: Point3<usize>,
    permitivity: &ScalarField,
    size: Vector3<usize>,
) -> Matrix3<f64> {
    let electric_tensor = green_tensor(frequency, point, permitivity, size, OperatorType::Electric);
    let magnetic_tensor = green_tensor(frequency, point, permitivity, size, OperatorType::Magnetic);

    frequency * frequency / PI * (
            (magnetic_tensor - Matrix3::from_diagonal_element(0.5) * magnetic_tensor.trace())
            + permitivity[point] * (electric_tensor - Matrix3::from_diagonal_element(0.5) * electric_tensor.trace())
    )
}

pub fn green_tensor(
    frequency: f64,
    point: Point3<usize>,
    permitivity: &ScalarField,
    size: Vector3<usize>,
    operator_type: OperatorType
) -> Matrix3<f64> {
    Matrix3::from_columns(&[
        green_function(
            frequency,
            point,
            Vector3::new(1.0, 0.0, 0.0),
            permitivity,
            size,
            operator_type
        ),
        green_function(
            frequency,
            point,
            Vector3::new(0.0, 1.0, 0.0),
            permitivity,
            size,
            operator_type
        ),
        green_function(
            frequency,
            point,
            Vector3::new(0.0, 0.0, 1.0),
            permitivity,
            size,
            operator_type
        ),
    ])
}

pub fn green_function(
    frequency: f64,
    point: Point3<usize>,
    polarization: Vector3<f64>,
    permitivity: &ScalarField,
    size: Vector3<usize>,
    operator_type: OperatorType
) -> Vector3<f64> {
    // The delta function right hand side
    let mut b = VectorField::new(size);
    b[point] = polarization;

    // The operator
    let a = Operator::new(frequency, permitivity, operator_type);

    let mut x = VectorField::new(size);
    let mut r = &b - &a * x.clone();
    let mut p = r.clone();
    let mut rsold = &r * &r;

    for _ in 0..100 {
        let a_p = &a * p.clone();
        let alpha = rsold / (&p * &a_p);
        x += alpha * &p;
        r -= alpha * &a_p;
        let rsnew = &r * &r;
        p = (rsnew / rsold) * p + &r;
        rsold = rsnew;
    }
    x[point]
}
