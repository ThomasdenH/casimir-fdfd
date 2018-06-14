extern crate nalgebra;
extern crate pbr;
extern crate rayon;
extern crate serde_json;
extern crate failure;

mod greenfunctions;
mod scalarfield;
mod scaledvectorfield;
mod vectorfield;
mod world;
mod config;

use config::ConfigFile;
use failure::Error;

fn main() -> Result<(), Error> {
    let config = ConfigFile::from_file("worlds/plates_medium.json")?;
    let force = config.to_world()?.force_on(0);
    println!("Calculated force: ({}, {}, {})", force.x, force.y, force.z);
    Ok(())
}
