#![warn(missing_docs)]
#![forbid(unsafe_code)]

//! `casimir-fdfd` is an implementation of a stress-tensor based FDFD method for computing Casimir
//! forces. Based on [Virtual photons in imaginary time: Computing exact Casimir forces via standard numerical-electromagnetism techniques](https://arxiv.org/abs/0705.3661).

/// Contains the simulation configuration options.
pub mod config;
/// Contains the fields used in the simulation.
pub mod fields;
/// Contains the functions to compute Green functions and stress tensors.
pub mod greenfunctions;
/// Geometry and materials that describe the simulated environment.
pub mod world;

use crate::world::World;
use clap::{App, Arg};
use serde_json;
use snafu::{ResultExt, Snafu};
use std::fs::File;
use std::path::Path;

#[derive(Debug, Snafu)]
enum Error {
    #[snafu(display("could not open world file: {}", source))]
    OpenWorld {
        /// The source that caused this error.
        source: std::io::Error,
    },
    #[snafu(display("could not deserialize world: {}", source))]
    DeserializeWorld {
        /// The source that caused this error.
        source: serde_json::error::Error,
    },
    #[snafu(display("invalid world: {}", source))]
    ValidateWorld {
        /// The source that caused this error.
        source: world::WorldError,
    },
}

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
    let world_file = File::open(world_path).context(OpenWorld {})?;
    let mut world: World = serde_json::from_reader(world_file).context(DeserializeWorld {})?;
    world.validate().context(ValidateWorld {})?;

    world.set_progress_bar(matches.is_present("progressbar"));

    for (i, force) in world.forces().enumerate() {
        println!(
            "Calculated force for object {}: ({}, {}, {})",
            i, force.x, force.y, force.z
        );
    }

    Ok(())
}
