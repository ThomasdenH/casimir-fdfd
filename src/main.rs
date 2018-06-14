extern crate nalgebra;
extern crate pbr;
extern crate rayon;

mod greenfunctions;
mod scalarfield;
mod scaledvectorfield;
mod vectorfield;
mod world;

use nalgebra::*;
use world::World;

fn main() {
    let mut world = World::new(Vector3::new(60, 60, 60));
    world.add_box(1, 1, 1, 10, 10, 2);
    world.add_box(1, 1, 10, 10, 10, 11);

    let force = world.force_on(0);
    println!("Calculated force: ({}, {}, {})", force.x, force.y, force.z);
}
