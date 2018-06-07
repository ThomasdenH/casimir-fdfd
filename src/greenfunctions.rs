use nalgebra::*;
use scalarfield::ScalarField;
use std::f32::consts::PI;
use std::ops::Mul;
use vectorfield::VectorField;

/// A representation of the green function operator
struct A<'a> {
    scalar_primary: &'a ScalarField,
    scalar_secundary: &'a ScalarField,
    frequency_2: f32,
}

impl<'b> A<'b> {
    fn new<'a>(
        frequency: f32,
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
        let scalar_part = self.frequency_2 * (x.clone() * self.scalar_primary);
        let curl_part = (x.curl_positive() * self.scalar_secundary).curl_negative();
        curl_part + &scalar_part
    }
}

pub fn stress_tensor(
    frequency: f32,
    point: Point3<usize>,
    electric_field: &ScalarField,
    magnetic_field: &ScalarField,
    inv_electric_field: &ScalarField,
    inv_magnetic_field: &ScalarField,
    size: Vector3<usize>,
) -> Matrix3<f32> {
    let electric_tensor = green_tensor(frequency, point, electric_field, inv_magnetic_field, size);
    let magnetic_tensor = green_tensor(frequency, point, magnetic_field, inv_electric_field, size);

    frequency * frequency / PI
        * (magnetic_field.at(point.x as isize, point.y as isize, point.z as isize)
            * (magnetic_tensor - Matrix3::from_diagonal_element(0.5) * magnetic_tensor.trace())
            + electric_field.at(point.x as isize, point.y as isize, point.z as isize)
                * (electric_tensor - Matrix3::from_diagonal_element(0.5) * electric_tensor.trace()))
}

pub fn green_tensor(
    frequency: f32,
    point: Point3<usize>,
    scalar_primary: &ScalarField,
    scalar_secundary: &ScalarField,
    size: Vector3<usize>,
) -> Matrix3<f32> {
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
    frequency: f32,
    point: Point3<usize>,
    polarization: Vector3<f32>,
    scalar_primary: &ScalarField,
    scalar_secundary: &ScalarField,
    size: Vector3<usize>,
) -> Vector3<f32> {
    // The delta function right hand side
    let mut b = VectorField::new(size.x, size.y, size.z);
    b.set(
        point.x as isize,
        point.y as isize,
        point.z as isize,
        polarization,
    );

    // The operator
    let a = A::new(frequency, scalar_primary, scalar_secundary);

    let mut x = VectorField::new(size.x, size.y, size.z);

    let mut r = b - &(&a * x.clone());

    for _ in 0..*size.iter().max().unwrap() {
        let a_r = &a * r.clone();
        let alpha = (&r * &r) / (&r * &a_r);
        let alpha_a_r = alpha * a_r;
        r = r - &alpha_a_r;
        x = x + &alpha_a_r;
    }

    x.at(point.x as isize, point.y as isize, point.z as isize)
}
