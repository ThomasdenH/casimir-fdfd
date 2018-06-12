extern crate nalgebra;
extern crate pbr;
extern crate rayon;

mod greenfunctions;
mod scalarfield;
mod vectorfield;
mod world;
mod scaledvectorfield;

use world::World;

fn main() {
    let mut world = World::new(15, 15, 25);
    world.add_box(5, 5, 1, 10, 10, 2);
    world.add_box(5, 5, 10, 10, 10, 11);
    println!("{:?}", world.force_on(0));
}
