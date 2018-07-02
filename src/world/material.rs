use serde::de::{self, Deserialize, Deserializer, MapAccess, Visitor};
use std::f32::consts::PI;
use std::fmt;
use std::marker::PhantomData;
use std::str::FromStr;

/// A material is any structure that provides a permittivity.
pub trait Material {
    /// Get the permittivity for the given frequency.
    fn permittivity(&self, freq: f32) -> f32;
}

/// A `DrudeMaterial` is a material of which the permittivity is given by the Drude model.
#[derive(Debug, Deserialize, PartialEq, Copy, Clone)]
pub struct DrudeMaterial {
    /// The plasma frequency in the model.
    omega_p: f32,
    /// The relaxation frequency in the model.
    omega_tau: f32,
    /// For the frequency integration, this step size will be used.
    step: f32,
    /// The frequency integration will stop after the contribution sinks below this value.
    precision: f32,
}

/// This error indicates that a string could not be parsed to a `DrudeMaterial`.
#[derive(Debug, Fail, Clone, Eq, PartialEq, Hash)]
#[fail(display = "unknown drude material: {}", name)]
pub struct UnknownDrudeMaterialError {
    name: String,
}

impl FromStr for DrudeMaterial {
    type Err = UnknownDrudeMaterialError;

    /// Some predefined material names can be converted to a `DrudeMaterial`. If it cannot be
    /// converted, an `Err(UnknownDrudeMaterialError)` is returned.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "gold" => Ok(DrudeMaterial {
                omega_p: 7.79,
                omega_tau: 48.8,
                step: 0.1,
                precision: 0.001,
            }),
            other => Err(UnknownDrudeMaterialError {
                name: other.to_string(),
            }),
        }
    }
}

impl fmt::Display for DrudeMaterial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Drude material: ωp: {}, ωτ: {}, step size: {}, precision: {}",
            self.omega_p, self.omega_tau, self.step, self.precision
        )
    }
}

impl Material for DrudeMaterial {
    fn permittivity(&self, freq: f32) -> f32 {
        let mut total = 0.0;
        for i in 0.. {
            let omega = i as f32 * self.step;
            let added = (self.omega_p * self.omega_p * self.omega_tau)
                / (omega * omega + self.omega_tau * self.omega_tau)
                / (omega * omega + freq * freq) * self.step;
            total += added;
            if added < self.precision {
                break;
            }
        }
        1.0 + 2.0 / PI * total
    }
}

/// This function is able to deserialize a `DrudeMaterial` from both the struct as well as the
/// string representation.
pub fn string_or_struct<'de, T, D>(deserializer: D) -> Result<T, D::Error>
where
    T: Deserialize<'de> + FromStr<Err = UnknownDrudeMaterialError>,
    D: Deserializer<'de>,
{
    // This is a Visitor that forwards string types to T's `FromStr` impl and
    // forwards map types to T's `Deserialize` impl. The `PhantomData` is to
    // keep the compiler from complaining about T being an unused generic type
    // parameter. We need T in order to know the Value type for the Visitor
    // impl.
    struct StringOrStruct<T>(PhantomData<fn() -> T>);

    impl<'de, T> Visitor<'de> for StringOrStruct<T>
    where
        T: Deserialize<'de> + FromStr<Err = UnknownDrudeMaterialError>,
    {
        type Value = T;

        fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
            formatter.write_str("string or map")
        }

        fn visit_str<E>(self, value: &str) -> Result<T, E>
        where
            E: de::Error,
        {
            Ok(FromStr::from_str(value).unwrap())
        }

        fn visit_map<M>(self, visitor: M) -> Result<T, M::Error>
        where
            M: MapAccess<'de>,
        {
            // `MapAccessDeserializer` is a wrapper that turns a `MapAccess`
            // into a `Deserializer`, allowing it to be used as the input to T's
            // `Deserialize` implementation. T then deserializes itself using
            // the entries from the map visitor.
            Deserialize::deserialize(de::value::MapAccessDeserializer::new(visitor))
        }
    }

    deserializer.deserialize_any(StringOrStruct(PhantomData))
}
