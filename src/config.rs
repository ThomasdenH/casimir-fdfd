use serde::Deserialize;
use std::fmt;

/// Contains the configuration related to the simulation itself. Generally, these settings determine
/// the accuracy of the result and the efficiency of the simulation.
#[derive(Copy, Clone, PartialEq, Debug, Deserialize)]
#[serde(default = "SimulationConfig::default")]
pub struct SimulationConfig {
    /// If the force from two frequencies are within this range from eachother, no subdivision will
    /// be used in the integration
    pub frequency_threshold: f32,
    /// If the norm of the vector is lower than this value, the field will be considered to have
    /// converged.
    pub fdfd_convergence: f32,
    /// After the cosine expansion adds less than this ratio of the n = (0, 0) force to the total,
    /// it will stop.
    pub cosine_cutoff: f32,
    /// The frequency range to integrate over.
    pub frequency_range: [f32; 2],
}

impl Default for SimulationConfig {
    fn default() -> SimulationConfig {
        SimulationConfig {
            frequency_threshold: 0.005,
            fdfd_convergence: 0.00001,
            cosine_cutoff: 0.005,
            frequency_range: [0.0001, 1.0],
        }
    }
}

impl fmt::Display for SimulationConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Simulation configuration:")?;
        writeln!(f, "\tFrequency threshold: {}", self.frequency_threshold)?;
        writeln!(
            f,
            "\tFrequency range: {} .. {}",
            self.frequency_range[0], self.frequency_range[1]
        )?;
        writeln!(f, "\tFDFD convergence: {}", self.fdfd_convergence)?;
        write!(f, "\tCosine cutoff: {}", self.cosine_cutoff)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;
    use serde_json::from_str;

    #[test]
    fn deserialize() {
        let config: SimulationConfig = from_str(
            r#"{
            "frequency_threshold": 0.1,
            "fdfd_convergence": 2.0,
            "cosine_cutoff": 1.0,
            "frequency_range": [0.1, 1.0]
        }"#,
        )
        .unwrap();
        assert_approx_eq!(config.frequency_threshold, 0.1, 1e-10);
        assert_approx_eq!(config.fdfd_convergence, 2.0, 1e-10);
        assert_approx_eq!(config.cosine_cutoff, 1.0, 1e-10);
        assert_approx_eq!(config.frequency_range[0], 0.1, 1e-10);
        assert_approx_eq!(config.frequency_range[1], 1.0, 1e-10);
    }
}
