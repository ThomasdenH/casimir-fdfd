#![warn(missing_docs)]

//! `casimir-fdfd` is an implementation of a stress-tensor based FDFD method for computing Casimir
//! forces. Based on [Virtual photons in imaginary time: Computing exact Casimir forces via standard numerical-electromagnetism techniques](https://arxiv.org/abs/0705.3661).

#[macro_use]
extern crate failure;
extern crate nalgebra;
extern crate pbr;
extern crate rayon;
extern crate serde_json;
#[macro_use]
extern crate serde_derive;
extern crate clap;
extern crate serde;
#[cfg(test)]
#[macro_use]
extern crate assert_approx_eq;

/// Contains the simulation configuration options.
pub mod config;
/// Contains the fields used in the simulation.
pub mod fields;
/// Contains the functions to compute Green functions and stress tensors.
pub mod greenfunctions;
/// Geometry and materials that describe the simulated environment.
pub mod world;

use clap::{App, Arg};
use failure::Error;
use std::fs::File;
use std::path::Path;
use crate::world::World;

fn main() -> Result<(), Error> {
    let matches = App::new("casimir-fdfd")
        .arg(
            Arg::with_name("INPUT")
                .help("Sets the input file to use")
                .required(true)
                .index(1),
        )
        .arg(Arg::from_usage(
            "-p, --progressbar 'Display a progress bar while running'",
        ))
        .get_matches();

    let world_path = Path::new(matches.value_of("INPUT").unwrap());
    let world_file = File::open(world_path)?;
    let mut world: World = serde_json::from_reader(world_file)?;
    world.validate()?;

    world.set_progress_bar(matches.is_present("progressbar"));

    for (i, force) in world.forces().enumerate() {
        println!(
            "Calculated force for object {}: ({}, {}, {})",
            i, force.x, force.y, force.z
        );
    }

    Ok(())
}
