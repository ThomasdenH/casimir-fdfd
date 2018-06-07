extern crate nalgebra;
extern crate pbr;

mod greenfunctions;
mod scalarfield;
mod vectorfield;
mod world;

use world::World;

fn main() {
    let mut world = World::new(30, 30, 50);
    world.add_box(10, 10, 10, 20, 20, 20);
    world.add_box(10, 10, 30, 20, 20, 40);
    println!("{:?}", world.force_on(0));
}
