extern crate failure;
extern crate nalgebra;
extern crate pbr;
extern crate rayon;
extern crate serde_json;

mod config;
mod greenfunctions;
mod scalarfield;
mod scaledvectorfield;
mod vectorfield;
mod world;

use config::ConfigFile;
use failure::Error;

fn main() -> Result<(), Error> {
    let config = ConfigFile::from_file("worlds/plates_tiny.json")?;
    let force = config.to_world()?.force_on(0);
    println!("Calculated force: ({}, {}, {})", force.x, force.y, force.z);
    Ok(())
}
