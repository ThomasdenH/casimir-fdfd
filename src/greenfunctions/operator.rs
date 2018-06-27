use std::ops::Mul;
use fields::{ScalarField, VectorField};

/// A representation of the green function operator
#[derive(Clone, PartialEq, Debug)]
pub struct Operator<'a> {
    permitivity: &'a ScalarField,
    frequency_2: f32,
    operator_type: OperatorType,
}

/// The type of operator to use. The structure is very similar, but a scalar multiplication is
/// changed.
#[derive(Copy, Clone, Eq, PartialEq, Debug, Hash)]
pub enum OperatorType {
    Electric,
    Magnetic,
}

impl<'b> Operator<'b> {
    pub fn new(frequency: f32, permitivity: &ScalarField, operator_type: OperatorType) -> Operator {
        Operator {
            permitivity,
            frequency_2: frequency * frequency,
            operator_type,
        }
    }

    /// The operation, only without new allocations using two temporary vectors.
    pub fn mul_with_temps(
        &self,
        x: VectorField,
        temp1: VectorField,
        temp2: VectorField,
    ) -> (VectorField, VectorField, VectorField) {
        match self.operator_type {
            OperatorType::Electric => {
                let curl_first = x.curl_positive_to(temp1);
                let curl_part = curl_first.curl_negative_to(temp2);
                let scalar_part = self.frequency_2 * (x * self.permitivity);
                (curl_part + &scalar_part, scalar_part, curl_first)
            }
            OperatorType::Magnetic => {
                let curl_first = x.curl_positive_to(temp1) / self.permitivity;
                let curl_part = curl_first.curl_negative_to(temp2);
                let scalar_part = self.frequency_2 * x;
                (curl_part + &scalar_part, scalar_part, curl_first)
            }
        }
    }
}

impl<'a, 'b> Mul<VectorField> for &'a Operator<'b> {
    type Output = VectorField;

    fn mul(self, x: VectorField) -> Self::Output {
        match self.operator_type {
            OperatorType::Electric => {
                let curl_part = x.curl_positive().curl_negative();
                let scalar_part = self.frequency_2 * (x * self.permitivity);
                curl_part + &scalar_part
            }
            OperatorType::Magnetic => {
                let curl_part = (x.curl_positive() / self.permitivity).curl_negative();
                let scalar_part = self.frequency_2 * x;
                curl_part + &scalar_part
            }
        }
    }
}
