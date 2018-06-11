extern crate nalgebra;
extern crate pbr;
extern crate rayon;

mod greenfunctions;
mod scalarfield;
mod vectorfield;
mod world;

use world::World;

fn main() {
    let mut world = World::new(15, 15, 25);
    world.add_box(5, 5, 5, 10, 10, 10);
    world.add_box(5, 5, 15, 10, 10, 20);
    println!("{:?}", world.force_on(0));
}
