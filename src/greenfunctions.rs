use nalgebra::*;
use scalarfield::ScalarField;
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

pub fn green_function(
    frequency: f32,
    point: Vector3<usize>,
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

    for _ in 0..10 {
        let a_r = &a * r.clone();
        let alpha = (&r * &r) / (&r * &a_r);
        let alpha_a_r = alpha * a_r;
        r = r - &alpha_a_r;
        x = x + &alpha_a_r;
    }

    x.at(point.x as isize, point.y as isize, point.z as isize)
}
