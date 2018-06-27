mod boundingbox;
mod shape;
mod material;

use nalgebra::*;
use fields::ScalarField;
use std::sync::{Arc, Mutex};
use config::SimulationConfig;
use greenfunctions::cosinebasis::{CosineBasis, Direction};
use pbr::ProgressBar;
use rayon::iter::*;
use std::io::Stdout;
use world::boundingbox::BoundingBox;
use world::shape::Shape;

/// A struct representing the geometry
#[derive(PartialEq, Clone, Debug, Deserialize)]
pub struct World {
    size: Vector3<usize>,
    objects: Vec<Shape>,
    simulation_config: SimulationConfig,
    #[serde(skip)]
    use_progress_bar: bool
}

pub struct ForceIterator<'a> {
    i: usize,
    world: &'a World
}

impl<'a> Iterator for ForceIterator<'a> {
    type Item = Vector3<f32>;

    fn next(&mut self) -> Option<Vector3<f32>> {
        if self.i < self.world.objects.len() {
            let force = self.world.force_on(self.i);
            self.i += 1;
            Some(force)
        } else {
            None
        }
    }
}

#[derive(Debug, Fail)]
pub enum WorldError {
    #[fail(display = "shape {} too close to edge", index)]
    ShapeTooCloseToEdge { index: usize },

    #[fail(display = "bounding boxes of shapes {} and {} intersect", index_1, index_2)]
    ShapesIntersect { index_1: usize, index_2: usize }
}

impl World {
    pub fn set_progress_bar(&mut self, enable: bool) {
        self.use_progress_bar = enable;
    }

    pub fn forces(&self) -> ForceIterator {
        ForceIterator {
            i: 0,
            world: &self
        }
    }

    /// Compute the force on the n'th object.
    pub fn force_on(&self, i: usize) -> Vector3<f32> {
        println!("Geometry:");
        println!(
            "\tWorld size: ({}, {}, {})",
            self.size.x, self.size.y, self.size.z
        );
        for (i, bbox) in self.objects.iter().enumerate() {
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

    /// Checks whether the bounding boxes of any of the objects intersect. If they do, the geometry
    /// might be too complicated, or the distance between objects could be too small.
    pub fn validate(&self) -> Result<(), WorldError> {
        let bbox_world = BoundingBox::new(0, 0, 0, self.size.x, self.size.y, self.size.z);

        let expanded_boxes = self.objects.iter().enumerate()
            .map(|(index, obj)| obj.bbox()
                .expanded(2)
                .map_err(|_| WorldError::ShapeTooCloseToEdge { index })
            )
            .collect::<Result<Vec<_>, _>>()?;

        for (i, bbox_1) in expanded_boxes.iter().enumerate() {
            // Check for intersection with world

            if !bbox_1.inside(&bbox_world) {
                return Err(WorldError::ShapeTooCloseToEdge {
                    index: i
                });
            }

            // Check for intersection with other objects
            for (j, bbox_2) in expanded_boxes.iter().enumerate() {
                if i < j && bbox_1.intersects(&bbox_2) {
                    return Err(WorldError::ShapesIntersect {
                        index_1: i,
                        index_2: j
                    });
                }
            }
        }
        Ok(())
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
        let bbox = &self.objects[i].bbox();
        let dx = bbox.x1 - bbox.x0 + 4;
        let dy = bbox.y1 - bbox.y0 + 4;
        let dz = bbox.z1 - bbox.z0 + 4;
        let count = 2 * (dx * dy + dy * dz + dz * dx) * (1 + self.objects.len());
        let progress_bar = if self.use_progress_bar {
            let progress_bar = Arc::new(Mutex::new(ProgressBar::new(count as u64)));
            progress_bar.lock().unwrap().format("[=>~]");
            progress_bar.lock().unwrap().tick();
            Some(progress_bar)
        } else {
            None
        };

        let perm_all_geom = &self.permitivity_field_all_geometry(frequency);
        let mut total_force =
            self.force_on_for_freq_and_geometry(frequency, perm_all_geom, bbox, &progress_bar);

        // Discretization gives rise to forces of an object on itself. Removing these gives more
        // accurate results.
        for other in &self.objects {
            let perm = &self.permitivity_field(frequency, &[*other]);
            total_force -=
                self.force_on_for_freq_and_geometry(frequency, perm, bbox, &progress_bar);
        }

        if let Some(progress_bar) = progress_bar {
            progress_bar.lock().unwrap().finish_println("");
        }
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
        progress_bar: &Option<Arc<Mutex<ProgressBar<Stdout>>>>,
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
                ).with_progress_bar(progress_bar.clone().map(|a| a.clone()))
                    .force(),
                1 => CosineBasis::new(
                    Point3::new(bbox.x0 - 2, bbox.y0 - 2, bbox.z1 + 2),
                    Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                    Direction::Z,
                ).with_progress_bar(progress_bar.clone().map(|a| a.clone()))
                    .force(),
                2 => CosineBasis::new(
                    Point3::new(bbox.x0 - 2, bbox.y0 - 2, bbox.z0 - 2),
                    Point3::new(bbox.x1 + 2, bbox.y0 - 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                    Direction::NegY,
                ).with_progress_bar(progress_bar.clone().map(|a| a.clone()))
                    .force(),
                3 => CosineBasis::new(
                    Point3::new(bbox.x0 - 2, bbox.y1 + 2, bbox.z0 - 2),
                    Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                    Direction::Y,
                ).with_progress_bar(progress_bar.clone().map(|a| a.clone()))
                    .force(),
                4 => CosineBasis::new(
                    Point3::new(bbox.x0 - 2, bbox.y0 - 2, bbox.z0 - 2),
                    Point3::new(bbox.x0 - 2, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                    Direction::NegX,
                ).with_progress_bar(progress_bar.clone().map(|a| a.clone()))
                    .force(),
                5 => CosineBasis::new(
                    Point3::new(bbox.x1 + 2, bbox.y0 - 2, bbox.z0 - 2),
                    Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                    frequency,
                    perm,
                    &self.simulation_config,
                    Direction::X,
                ).with_progress_bar(progress_bar.clone().map(|a| a.clone()))
                    .force(),

                i => panic!("Face index out of bounds: {}", i),
            })
            .sum()
    }

    /// Returns a scalar field representing the permitivity of a vector of bounding boxes.
    fn permitivity_field(&self, freq: f32, objects: &[Shape]) -> ScalarField {
        let mut permitivity_field = ScalarField::ones(self.size);
        for shape in objects {
            shape.draw_permitivity(&mut permitivity_field, freq);
        }
        permitivity_field
    }

    /// Returns a scalar field representing the permitivity of the entire geometry.
    fn permitivity_field_all_geometry(&self, freq: f32) -> ScalarField {
        self.permitivity_field(freq, &self.objects)
    }
}
