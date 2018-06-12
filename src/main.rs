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
    let mut world = World::new(Vector3::new(50, 50, 50));
    world.add_box(1, 1, 1, 6, 6, 2);
    world.add_box(1, 1, 11, 11, 11, 12);

    let force = world.force_on(0);
    println!("Calculated force: ({}, {}, {})", force.x, force.y, force.z);
}
