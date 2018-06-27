use std::fmt;

/// Contains the configuration related to the simulation itself.
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
            frequency_threshold: 0.01,
            fdfd_convergence: 0.00001,
            cosine_cutoff: 0.01,
            frequency_range: [0.01, 1.0],
        }
    }
}

impl fmt::Display for SimulationConfig {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
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
