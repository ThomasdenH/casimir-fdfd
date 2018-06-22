use failure::Error;
use nalgebra::*;
use serde_json;
use serde_json::Value;
use std::fmt;
use std::fs::File;
use std::io::Read;
use world::World;

/// Contains the json representation of a configuration file.
#[derive(Clone, PartialEq, Debug)]
pub struct ConfigFile {
    json_value: serde_json::Value,
}

/// Contains the configuration related to the simulation itself.
#[derive(Copy, Clone, PartialEq, Debug)]
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

impl ConfigFile {
    pub fn from_file(path: &str) -> Result<ConfigFile, Error> {
        let mut file = File::open(path)?;
        let mut content = String::new();
        file.read_to_string(&mut content)?;
        let json_value = serde_json::from_str(&content)?;
        println!("Loaded file '{}'.", path);
        Ok(ConfigFile { json_value })
    }

    fn read_simulation_config(&self) -> SimulationConfig {
        let json_config = &self.json_value["simulation_config"];
        let mut config = SimulationConfig::default();
        if let Value::Object(map) = json_config {
            if let Value::Number(ref num) = map["frequency_threshold"] {
                if let Some(frequency_threshold) = num.as_f64() {
                    config.frequency_threshold = frequency_threshold as f32;
                }
            }
            if let Value::Number(ref num) = map["fdfd_convergence"] {
                if let Some(fdfd_convergence) = num.as_f64() {
                    config.fdfd_convergence = fdfd_convergence as f32;
                }
            }
            if let Value::Number(ref num) = map["cosine_cutoff"] {
                if let Some(cosine_cutoff) = num.as_f64() {
                    config.cosine_cutoff = cosine_cutoff as f32;
                }
            }
            if let Value::Array(ref range) = map["frequency_range"] {
                if let (Some(start), Some(end)) = (range[0].as_f64(), range[1].as_f64()) {
                    config.frequency_range = [start as f32, end as f32];
                }
            }
        }
        config
    }

    pub fn to_world(&self) -> Result<World, Error> {
        let size = &self.json_value["size"];
        let simulation_config = self.read_simulation_config();
        let mut world = World::new(
            Vector3::new(
                size[0].as_u64().unwrap() as usize,
                size[1].as_u64().unwrap() as usize,
                size[2].as_u64().unwrap() as usize,
            ),
            simulation_config,
        );

        for object in self.json_value["objects"].as_array().unwrap() {
            if let Value::String(ref a) = object["type"] {
                match a.as_ref() {
                    "box" => ConfigFile::read_box(&mut world, &object),
                    "sphere" => ConfigFile::read_sphere(&mut world, &object),
                    a => panic!("Unknown object type: {}", a),
                }
            }
        }
        Ok(world)
    }

    fn read_box(world: &mut World, object: &Value) {
        let p0 = &object["p0"];
        let p1 = &object["p1"];
        world.add_box(
            p0[0].as_u64().unwrap() as usize,
            p0[1].as_u64().unwrap() as usize,
            p0[2].as_u64().unwrap() as usize,
            p1[0].as_u64().unwrap() as usize,
            p1[1].as_u64().unwrap() as usize,
            p1[2].as_u64().unwrap() as usize,
        );
    }

    fn read_sphere(world: &mut World, object: &Value) {
        let p = &object["point"];
        let radius = &object["radius"];
        world.add_sphere(
            Point3::new(
                p[0].as_u64().unwrap() as usize,
                p[1].as_u64().unwrap() as usize,
                p[2].as_u64().unwrap() as usize,
            ),
            radius.as_f64().unwrap() as f32,
        );
    }
}
