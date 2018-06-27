use std::f32::consts::PI;
use std::str::FromStr;
use std::fmt;
use std::marker::PhantomData;
use serde::de::{self, Deserialize, Deserializer, Visitor, MapAccess};

pub trait Material {
    fn permitivity(&self, freq: f32) -> f32;
}

#[derive(Debug, Deserialize, PartialEq, Copy, Clone)]
pub struct DrudeMaterial {
    omega_p: f32,
    omega_tau: f32,
    step: f32,
    precision: f32
}

#[derive(Debug, Fail, Clone, Eq, PartialEq, Hash)]
#[fail(display = "unknown drude material: {}", name)]
pub struct UnknownDrudeMaterialError {
    name: String
}

impl FromStr for DrudeMaterial {
    // This implementation of `from_str` can never fail, so use the impossible
    // `Void` type as the error type.
    type Err = UnknownDrudeMaterialError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "gold" => Ok(DrudeMaterial {
                omega_p: 7.79,
                omega_tau: 48.8,
                step: 0.1,
                precision: 0.001
            }),
            other => Err(UnknownDrudeMaterialError{ name: other.to_string() })
        }
    }
}

impl fmt::Display for DrudeMaterial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Drude material: ωp: {}, ωτ: {}, step size: {}, precision: {}", self.omega_p, self.omega_tau, self.step, self.precision)
    }
}

impl Material for DrudeMaterial {
    fn permitivity(&self, freq: f32) -> f32 {
        let mut total = 0.0;
        for i in 0.. {
            let omega = i as f32 * self.step;
            let added = (self.omega_p * self.omega_p * self.omega_tau)
                / (omega * omega + self.omega_tau * self.omega_tau)
                / (omega * omega + freq * freq)
                * self.step;
            total += added;
            if added < self.precision {
                break;
            }
        }
        1.0 + 2.0 / PI * total
    }
}

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
