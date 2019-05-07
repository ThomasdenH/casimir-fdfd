mod boundingbox;
mod material;
mod shape;

use config::SimulationConfig;
use fields::ScalarField;
use greenfunctions::cosinebasis::{CosineBasis, Direction};
use nalgebra::*;
use pbr::ProgressBar;
use rayon::iter::*;
use std::io::Stdout;
use std::sync::{Arc, Mutex};
use world::boundingbox::BoundingBox;
use world::shape::Shape;

/// A struct representing the world.
#[derive(PartialEq, Clone, Debug, Deserialize)]
pub struct World {
    /// The grid size corresponding to this world.
    size: Vector3<usize>,
    /// The list with the different shapes.
    objects: Vec<Shape>,
    /// The simulation configuration.
    simulation_config: SimulationConfig,
    #[serde(skip)]
    use_progress_bar: bool,
}

/// A struct that iterates over the different forces of the objects.
pub struct ForceIterator<'a> {
    i: usize,
    world: &'a World,
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

/// These errors can be thrown when validating a world.
#[derive(Debug, Fail)]
pub enum WorldError {
    /// An error indicating that one of the shapes is too close too the edge, causing the bounding
    /// box to cross the boundary.
    #[fail(display = "shape {} too close to edge", index)]
    ShapeTooCloseToEdge {
        /// The index of the violating object.
        index: usize,
    },

    /// An error indicating that the bounding boxes of two objects in this world intersect and
    /// and therefore this world is invalid.
    #[fail(
        display = "bounding boxes of shapes {} and {} intersect",
        index_1, index_2
    )]
    ShapesIntersect {
        /// The index of the first object intersecting with the second.
        index_1: usize,
        /// The index of the second object intersecting with the first.
        index_2: usize,
    },
}

impl World {
    /// Enable or disable the progress bar for the simulation.
    pub fn set_progress_bar(&mut self, enable: bool) {
        self.use_progress_bar = enable;
    }

    /// Obtain a force iterator for all objects in this world.
    pub fn forces(&self) -> ForceIterator {
        ForceIterator { i: 0, world: &self }
    }

    /// Compute the force on the `i`'th object.
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

    /// This function validates the geometry of the world. The function should be called, because it
    /// guarantees overflows later on in the simulation.
    ///
    /// # Errors
    /// - If any of the shapes is too close to the world border, the simulation can't be run and a
    /// `WorldError::ShapeTooCloseToEdge` will be returned. To fix this, move the object or increase
    /// the grid size.
    ///
    /// - If any of the objects are too close too eachother, their boundingboxes might intersect and
    /// the results will be invalid. If this is the case, a `WorldError::ShapesIntersect` will be
    /// returned. The violating shape indexes will be contained within. To fix this, move one or
    /// both of the objects.
    pub fn validate(&self) -> Result<(), WorldError> {
        let bbox_world = BoundingBox::new(0, 0, 0, self.size.x, self.size.y, self.size.z);

        let expanded_boxes = self
            .objects
            .iter()
            .enumerate()
            .map(|(index, obj)| {
                obj.bbox()
                    .expanded(2)
                    .map_err(|_| WorldError::ShapeTooCloseToEdge { index })
            })
            .collect::<Result<Vec<_>, _>>()?;

        for (i, bbox_1) in expanded_boxes.iter().enumerate() {
            // Check for intersection with world

            if !bbox_1.inside(&bbox_world) {
                return Err(WorldError::ShapeTooCloseToEdge { index: i });
            }

            // Check for intersection with other objects
            for (j, bbox_2) in expanded_boxes.iter().enumerate() {
                if i < j && bbox_1.intersects(&bbox_2) {
                    return Err(WorldError::ShapesIntersect {
                        index_1: i,
                        index_2: j,
                    });
                }
            }
        }
        Ok(())
    }

    /// Performs a recursive integration between two frequencies. If the difference between the
    /// midpoint force and the linear interpolation is too large, both sides of the domain will use
    /// this function recursively to integrate the force.
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
            ) + self.integrate_force_between_frequencies(
                i,
                middle_frequency,
                end_frequency,
                middle_value,
                end_value,
                max,
            )
        }
    }

    /// Compute the force on object `i` for a certain `frequency`. This method will also subtract
    /// the error forces due to quantization by subtracting the force due to single objects.
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

        let perm_all_geom = &self.permittivity_field_all_geometry(frequency);
        let mut total_force =
            self.force_on_for_freq_and_geometry(frequency, perm_all_geom, bbox, &progress_bar);

        // Discretization gives rise to forces of an object on itself. Removing these gives more
        // accurate results.
        for other in &self.objects {
            let perm = &self.permittivity_field(frequency, &[*other]);
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

    /// Compute the force on the geometry inside `BoundingBox`, for the given permittivity field
    /// `perm` and `BoundingBox` `bbox`.
    fn force_on_for_freq_and_geometry(
        &self,
        frequency: f32,
        perm: &ScalarField,
        bbox: &BoundingBox,
        progress_bar: &Option<Arc<Mutex<ProgressBar<Stdout>>>>,
    ) -> Vector3<f32> {
        (0..6)
            .into_par_iter()
            .map(|face| {
                (match face {
                    0 => CosineBasis::new(
                        Point3::new(bbox.x0 - 2, bbox.y0 - 2, bbox.z0 - 2),
                        Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z0 - 2),
                        frequency,
                        perm,
                        &self.simulation_config,
                        Direction::NegZ,
                    ),
                    1 => CosineBasis::new(
                        Point3::new(bbox.x0 - 2, bbox.y0 - 2, bbox.z1 + 2),
                        Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                        frequency,
                        perm,
                        &self.simulation_config,
                        Direction::Z,
                    ),
                    2 => CosineBasis::new(
                        Point3::new(bbox.x0 - 2, bbox.y0 - 2, bbox.z0 - 2),
                        Point3::new(bbox.x1 + 2, bbox.y0 - 2, bbox.z1 + 2),
                        frequency,
                        perm,
                        &self.simulation_config,
                        Direction::NegY,
                    ),
                    3 => CosineBasis::new(
                        Point3::new(bbox.x0 - 2, bbox.y1 + 2, bbox.z0 - 2),
                        Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                        frequency,
                        perm,
                        &self.simulation_config,
                        Direction::Y,
                    ),
                    4 => CosineBasis::new(
                        Point3::new(bbox.x0 - 2, bbox.y0 - 2, bbox.z0 - 2),
                        Point3::new(bbox.x0 - 2, bbox.y1 + 2, bbox.z1 + 2),
                        frequency,
                        perm,
                        &self.simulation_config,
                        Direction::NegX,
                    ),
                    5 => CosineBasis::new(
                        Point3::new(bbox.x1 + 2, bbox.y0 - 2, bbox.z0 - 2),
                        Point3::new(bbox.x1 + 2, bbox.y1 + 2, bbox.z1 + 2),
                        frequency,
                        perm,
                        &self.simulation_config,
                        Direction::X,
                    ),
                    i => panic!("Face index out of bounds: {}", i),
                })
                .with_progress_bar(progress_bar.clone())
                .force()
            })
            .sum()
    }

    /// Returns a scalar field representing the permittivity of a vector of bounding boxes.
    fn permittivity_field(&self, freq: f32, objects: &[Shape]) -> ScalarField {
        let mut permittivity_field = ScalarField::ones(self.size);
        for shape in objects {
            shape.draw_permittivity(&mut permittivity_field, freq);
        }
        permittivity_field
    }

    /// Returns a scalar field representing the permittivity of the entire geometry.
    fn permittivity_field_all_geometry(&self, freq: f32) -> ScalarField {
        self.permittivity_field(freq, &self.objects)
    }
}
