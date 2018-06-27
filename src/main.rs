#[macro_use]
extern crate failure;
extern crate nalgebra;
extern crate pbr;
extern crate rayon;
extern crate serde_json;
#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate clap;

mod config;
mod greenfunctions;
mod fields;
mod world;

use failure::Error;
use world::World;
use std::fs::File;
use std::path::Path;
use clap::{Arg, App};

fn main() -> Result<(), Error> {
    let matches = App::new("casimir-fdfd")
        .arg(Arg::with_name("INPUT")
            .help("Sets the input file to use")
            .required(true)
            .index(1))
        .get_matches();

    let world_path = Path::new(matches.value_of("INPUT").unwrap());
    let world_file = File::open(world_path)?;
    let world: World = serde_json::from_reader(world_file)?;
    world.validate()?;

    for (i, force) in world.forces().enumerate() {
        println!("Calculated force for object {}: ({}, {}, {})", i, force.x, force.y, force.z);
    }

    Ok(())
}
