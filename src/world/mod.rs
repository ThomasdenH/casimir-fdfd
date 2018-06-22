use nalgebra::*;
use scalarfield::ScalarField;
use std::f32::consts::PI;
use std::sync::{Arc, Mutex};

mod boundingbox;

use config::SimulationConfig;
use greenfunctions::cosinebasis::{CosineBasis, Direction};
use pbr::ProgressBar;
use rayon::iter::*;
use std::fmt;
use std::io::Stdout;
use world::boundingbox::BoundingBox;

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum Shape {
    Box {
        bbox: BoundingBox,
    },
    Sphere {
        point: Point3<usize>,
        radius: f32,
        bbox: BoundingBox,
    },
}

impl Shape {
    pub fn bbox(&self) -> &BoundingBox {
        match self {
            Shape::Box { bbox } => &bbox,
            Shape::Sphere { bbox, .. } => &bbox,
        }
    }
}

impl fmt::Display for Shape {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Shape::Box { bbox } => write!(f, "{}", bbox),
            Shape::Sphere { point, radius, .. } => write!(
                f,
                "Sphere: ({}, {}, {}), with radius \
                 {}",
                point.x, point.y, point.z, radius
            ),
        }
    }
}

/// A struct representing the geometry
#[derive(PartialEq, Clone, Debug)]
pub struct World {
    size: Vector3<usize>,
    shapes: Vec<Shape>,
    simulation_config: SimulationConfig,
}

impl World {
    /// Create a new empty world.
    pub fn new(size: Vector3<usize>, simulation_config: SimulationConfig) -> World {
        World {
            size,
            shapes: Vec::new(),
            simulation_config,
        }
    }

    /// Add a box to this world.
    pub fn add_box(&mut self, x0: usize, y0: usize, z0: usize, x1: usize, y1: usize, z1: usize) {
        debug_assert!(
            x0 < x1
                && y0 < y1
                && z0 < z1
                && x1 < self.size.x
                && y1 < self.size.y
                && z1 < self.size.z
        );
        let bbox = BoundingBox::new(x0, y0, z0, x1, y1, z1);
        debug_assert!(!self.shapes.iter().any(|b| b.bbox().intersects(&bbox)));
        self.shapes.push(Shape::Box { bbox });
    }

    pub fn add_sphere(&mut self, p: Point3<usize>, radius: f32) {
        let x0 = p.x - radius.ceil() as usize - 1;
        let y0 = p.y - radius.ceil() as usize - 1;
        let z0 = p.z - radius.ceil() as usize - 1;
        let x1 = p.x + radius.ceil() as usize + 1;
        let y1 = p.y + radius.ceil() as usize + 1;
        let z1 = p.z + radius.ceil() as usize + 1;
        debug_assert!(
            x0 > 0 && y0 > 0 && z0 > 0 && x1 < self.size.x && y1 < self.size.y && z1 < self.size.z
        );
        let bbox = BoundingBox::new(x0, y0, z0, x1, y1, z1);
        debug_assert!(!self.shapes.iter().any(|b| b.bbox().intersects(&bbox)));
        self.shapes.push(Shape::Sphere {
            bbox,
            point: p,
            radius,
        })
    }

    /// Compute the force on the n'th object.
    pub fn force_on(&self, i: usize) -> Vector3<f32> {
        println!("Geometry:");
        println!(
            "\tWorld size: ({}, {}, {})",
            self.size.x, self.size.y, self.size.z
        );
        for (i, bbox) in self.shapes.iter().enumerate() {
            println!("\t{}: {}", i, bbox);
        }

        println!("{}", self.simulation_config);

        // The maximum frequency is given by (speed of light) / (grid element size)
        let start_freq = self.simulation_config.frequency_range[0];
        let end_freq = self.simulation_config.frequency_range[1];
        let start_force = self.force_on_for_freq(i, start_freq);
        let end_force = self.force_on_for_freq(i, end_freq);

        self.integrate_force_between_frequencies(
            i,
            start_freq,
            end_freq,
            start_force,
            end_force,
            start_force.norm(),
        )
    }

    pub fn integrate_force_between_frequencies(
        &self,
        i: usize,
        start_frequency: f32,
        end_frequency: f32,
        start_value: Vector3<f32>,
        end_value: Vector3<f32>,
        max: f32,
    ) -> Vector3<f32> {
        // Do a recursive integration. The function should be smooth.
        let middle_frequency = 0.5 * (start_frequency + end_frequency);
        let middle_value = self.force_on_for_freq(i, middle_frequency);

        let average = (start_value + end_value) / 2.0;

        if (average - middle_value).norm() * (end_frequency - start_frequency)
            < self.simulation_config.frequency_threshold * max
        {
            // The change in area from the middle value vs extrapolation is less than the threshold
            0.5 * (start_value + 2.0 * middle_value + end_value) * (end_frequency - start_frequency)
        } else {
            self.integrate_force_between_frequencies(
                i,
                start_frequency,
                middle_frequency,
                start_value,
                middle_value,
                max,
            )
                + self.integrate_force_between_frequencies(
                    i,
                    middle_frequency,
                    end_frequency,
                    middle_value,
                    end_value,
                    max,
                )
        }
    }

    fn force_on_for_freq(&self, i: usize, frequency: f32) -> Vector3<f32> {
        // Progress bar
        let bbox = &self.shapes[i].bbox();
        let dx = bbox.x1 - bbox.x0 + 4;
        let dy = bbox.y1 - bbox.y0 + 4;
        let dz = bbox.z1 - bbox.z0 + 4;
        let count = 2 * (dx * dy + dy * dz + dz * dx) * (1 + self.shapes.len());
        let progress_bar = Arc::new(Mutex::new(ProgressBar::new(count as u64)));
        progress_bar.lock().unwrap().format("╢▌▌░╟");
        progress_bar.lock().unwrap().tick();

        let perm_all_geom = &self.permitivity_field_all_geometry(frequency);
        let mut total_force =
            self.force_on_for_freq_and_geometry(frequency, perm_all_geom, bbox, &progress_bar);

        // Discretization gives rise to forces of an object on itself. Removing these gives more
        // accurate results.
        for other in &self.shapes {
            let perm = &self.permitivity_field(frequency, &[*other]);
            total_force -=
                self.force_on_for_freq_and_geometry(frequency, perm, bbox, &progress_bar);
        }

        progress_bar.lock().unwrap().finish_println("");
        println!(
            "Force for frequency {}: ({}, {}, {})",
            frequency, total_force.x, total_force.y, total_force.z
        );

        total_force
    }

    fn force_on_for_freq_and_geometry(
        &self,
        frequency: f32,
        perm: &ScalarField,
        bbox: &BoundingBox,
        progress_bar: &Arc<Mutex<ProgressBar<Stdout>>>,
    ) -> Vector3<f32> {
        (0..6)
            .into_par_iter()
            .map(|face| match face {
                0 => CosineBasis::new(
                    Point3::new(bbox.x0 - 2, bbox.y0 - 2, bbox.z0 - 2),
                    Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z0 - 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                    Direction::NegZ,
                ).with_progress_bar(progress_bar.clone())
                    .force(),
                1 => CosineBasis::new(
                    Point3::new(bbox.x0 - 2, bbox.y0 - 2, bbox.z1 + 2),
                    Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                    Direction::Z,
                ).with_progress_bar(progress_bar.clone())
                    .force(),
                2 => CosineBasis::new(
                    Point3::new(bbox.x0 - 2, bbox.y0 - 2, bbox.z0 - 2),
                    Point3::new(bbox.x1 + 2, bbox.y0 - 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                    Direction::NegY,
                ).with_progress_bar(progress_bar.clone())
                    .force(),
                3 => CosineBasis::new(
                    Point3::new(bbox.x0 - 2, bbox.y1 + 2, bbox.z0 - 2),
                    Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                    Direction::Y,
                ).with_progress_bar(progress_bar.clone())
                    .force(),
                4 => CosineBasis::new(
                    Point3::new(bbox.x0 - 2, bbox.y0 - 2, bbox.z0 - 2),
                    Point3::new(bbox.x0 - 2, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                    Direction::NegX,
                ).with_progress_bar(progress_bar.clone())
                    .force(),
                5 => CosineBasis::new(
                    Point3::new(bbox.x1 + 2, bbox.y0 - 2, bbox.z0 - 2),
                    Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                    Direction::X,
                ).with_progress_bar(progress_bar.clone())
                    .force(),

                i => panic!("Face index out of bounds: {}", i),
            })
            .sum()
    }

    /// Returns a scalar field representing the permitivity of a vector of bounding boxes.
    fn permitivity_field(&self, freq: f32, shapes: &[Shape]) -> ScalarField {
        let mut permitivity_field = ScalarField::ones(self.size);
        let permitivity = World::gold_permitivity(freq);
        for shape in shapes {
            match shape {
                Shape::Box { bbox } => {
                    for x in bbox.x0..bbox.x1 {
                        for y in bbox.y0..bbox.y1 {
                            for z in bbox.z0..bbox.z1 {
                                permitivity_field.set(
                                    x as isize,
                                    y as isize,
                                    z as isize,
                                    permitivity,
                                );
                            }
                        }
                    }
                }
                Shape::Sphere {
                    point,
                    radius,
                    bbox,
                } => {
                    for x in bbox.x0..bbox.x1 {
                        for y in bbox.y0..bbox.y1 {
                            for z in bbox.z0..bbox.z1 {
                                if Vector3::new(
                                    (point.x - x) as f32,
                                    (point.y - y) as f32,
                                    (point.z - z) as f32,
                                ).norm() <= *radius
                                {
                                    permitivity_field.set(
                                        x as isize,
                                        y as isize,
                                        z as isize,
                                        permitivity,
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
        permitivity_field
    }

    /// Returns a scalar field representing the permitivity of the entire geometry.
    fn permitivity_field_all_geometry(&self, freq: f32) -> ScalarField {
        self.permitivity_field(freq, &self.shapes)
    }

    /// Calculates the permitivity of gold for a particular imaginary frequency.
    fn gold_permitivity(freq: f32) -> f32 {
        let omega_p = 7.79;
        let omega_tau = 48.8;
        let mut total = 0.0;
        for i in 0.. {
            let omega = i as f32 * 0.1;
            let added = (omega_p * omega_p * omega_tau)
                / (omega * omega + omega_tau * omega_tau)
                / (omega * omega + freq * freq) * 0.1;
            total += added;
            if added < 0.001 {
                break;
            }
        }
        1.0 + 2.0 / PI * total
    }
}
