use nalgebra::*;
use scalarfield::ScalarField;
use std::f32::consts::PI;
use std::sync::{Arc, Mutex};

mod boundingbox;

use config::SimulationConfig;
use greenfunctions::cosinebasis::CosineBasis;
use pbr::ProgressBar;
use rayon::iter::*;
use std::io::Stdout;
use world::boundingbox::BoundingBox;

/// A struct representing the geometry
#[derive(PartialEq, Clone, Debug)]
pub struct World {
    size: Vector3<usize>,
    bboxes: Vec<BoundingBox>,
    simulation_config: SimulationConfig,
}

impl World {
    /// Create a new empty world.
    pub fn new(size: Vector3<usize>, simulation_config: SimulationConfig) -> World {
        World {
            size,
            bboxes: Vec::new(),
            simulation_config,
        }
    }

    /// Add a box to this world.
    pub fn add_box(&mut self, x0: usize, y0: usize, z0: usize, x1: usize, y1: usize, z1: usize) {
        debug_assert!(
            x0 < x1 && y0 < y1 && z0 < z1 && x1 < self.size.x && y1 < self.size.y
                && z1 < self.size.z
        );
        let bbox = BoundingBox::new(x0, y0, z0, x1, y1, z1);
        debug_assert!(!self.bboxes.iter().any(|b| b.intersects(&bbox)));
        self.bboxes.push(bbox);
    }

    /// Compute the force on the n'th object.
    pub fn force_on(&self, i: usize) -> Vector3<f32> {
        println!("Geometry:");
        println!(
            "\tWorld size: ({}, {}, {})",
            self.size.x, self.size.y, self.size.z
        );
        for (i, bbox) in self.bboxes.iter().enumerate() {
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
        if (start_value - end_value).norm() < self.simulation_config.frequency_threshold * max {
            // The difference is small enough to do a line approximation
            0.5 * (start_value + end_value) * (end_frequency - start_frequency)
        } else {
            let middle_frequency = 0.5 * (start_frequency + end_frequency);
            let middle_value = self.force_on_for_freq(i, middle_frequency);
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
        let bbox = &self.bboxes[i];
        let dx = (bbox.x1 - bbox.x0 + 2).min(self.simulation_config.cosine_depth);
        let dy = (bbox.y1 - bbox.y0 + 2).min(self.simulation_config.cosine_depth);
        let dz = (bbox.z1 - bbox.z0 + 2).min(self.simulation_config.cosine_depth);
        let count = 2 * (dx * dy + dy * dz + dz * dx) * (1 + self.bboxes.len());
        let progress_bar = Arc::new(Mutex::new(ProgressBar::new(count as u64)));
        progress_bar.lock().unwrap().format("╢▌▌░╟");
        progress_bar.lock().unwrap().tick();

        let perm_all_geom = &self.permitivity_field_all_geometry(frequency);
        let mut total_force = self.force_on_for_freq_and_geometry(
            frequency,
            perm_all_geom,
            &self.bboxes[i],
            &progress_bar,
        );

        // Discretization gives rise to forces of an object on itself. Removing these gives more
        // accurate results.
        for bbox in &self.bboxes {
            let perm = &self.permitivity_field(frequency, &[*bbox]);
            total_force -= self.force_on_for_freq_and_geometry(
                frequency,
                perm,
                &self.bboxes[i],
                &progress_bar,
            );
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
                0 => -CosineBasis::new(
                    Point3::new(bbox.x0 - 1, bbox.y0 - 1, bbox.z0 - 1),
                    Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z0 - 1),
                    frequency,
                    perm,
                    &self.simulation_config,
                ).with_progress_bar(progress_bar.clone())
                    .force(),
                1 => CosineBasis::new(
                    Point3::new(bbox.x0 - 1, bbox.y0 - 1, bbox.z1 + 2),
                    Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                ).with_progress_bar(progress_bar.clone())
                    .force(),
                2 => -CosineBasis::new(
                    Point3::new(bbox.x0 - 1, bbox.y0 - 1, bbox.z0 - 1),
                    Point3::new(bbox.x1 + 2, bbox.y0 - 1, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                ).with_progress_bar(progress_bar.clone())
                    .force(),
                3 => CosineBasis::new(
                    Point3::new(bbox.x0 - 1, bbox.y1 + 2, bbox.z0 - 1),
                    Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                ).with_progress_bar(progress_bar.clone())
                    .force(),
                4 => -CosineBasis::new(
                    Point3::new(bbox.x0 - 1, bbox.y0 - 1, bbox.z0 - 1),
                    Point3::new(bbox.x0 - 1, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                ).with_progress_bar(progress_bar.clone())
                    .force(),
                5 => CosineBasis::new(
                    Point3::new(bbox.x1 + 2, bbox.y0 - 1, bbox.z0 - 1),
                    Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                ).with_progress_bar(progress_bar.clone())
                    .force(),

                i => panic!("Face index out of bounds: {}", i),
            })
            .sum()
    }

    /// Returns a scalar field representing the permitivity of a vector of bounding boxes.
    fn permitivity_field(&self, freq: f32, boxes: &[BoundingBox]) -> ScalarField {
        let mut permitivity_field = ScalarField::ones(self.size);
        let permitivity = World::gold_permitivity(freq);
        for bbox in boxes {
            for x in bbox.x0..bbox.x1 {
                for y in bbox.y0..bbox.y1 {
                    for z in bbox.z0..bbox.z1 {
                        permitivity_field.set(x as isize, y as isize, z as isize, permitivity);
                    }
                }
            }
        }
        permitivity_field
    }

    /// Returns a scalar field representing the permitivity of the entire geometry.
    fn permitivity_field_all_geometry(&self, freq: f32) -> ScalarField {
        self.permitivity_field(freq, &self.bboxes)
    }

    /// Calculates the permitivity of gold for a particular imaginary frequency.
    fn gold_permitivity(freq: f32) -> f32 {
        let omega_p = 7.79;
        let omega_tau = 48.8;
        let mut total = 0.0;
        for i in 0.. {
            let omega = i as f32 * 0.1;
            let added = (omega_p * omega_p * omega_tau) / (omega * omega + omega_tau * omega_tau)
                / (omega * omega + freq * freq) * 0.1;
            total += added;
            if added < 0.001 {
                break;
            }
        }
        1.0 + 2.0 / PI * total
    }
}
