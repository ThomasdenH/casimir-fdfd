use fields::ScalarField;
use nalgebra::*;
use std::fmt;
use world;
use world::boundingbox::BoundingBox;
use world::material::{DrudeMaterial, Material};

/// A `Shape` represents geometry in space. It can be drawn to a permittivity field using
#[derive(PartialEq, Copy, Clone, Debug, Deserialize)]
pub enum Shape {
    Box {
        p0: Point3<usize>,
        p1: Point3<usize>,
        #[serde(deserialize_with = "world::material::string_or_struct")]
        material: DrudeMaterial,
    },
    Sphere {
        point: Point3<usize>,
        radius: f32,
        #[serde(deserialize_with = "world::material::string_or_struct")]
        material: DrudeMaterial,
    },
}

impl fmt::Display for Shape {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Shape::Box { p0, p1, material } => write!(
                f,
                "Box: ({}, {}, {}), ({}, {}, {}) with material {}",
                p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, material
            ),
            Shape::Sphere {
                point,
                radius,
                material,
            } => write!(
                f,
                "Sphere: ({}, {}, {}), with radius \
                 {} and material {}",
                point.x, point.y, point.z, radius, material
            ),
        }
    }
}

impl Shape {
    /// Create a new box, with the lowest coordinates in `p0`, and the highest coordinates in `p1`.
    /// `material` is the material to use for this object.
    pub fn new_box(p0: Point3<usize>, p1: Point3<usize>, material: DrudeMaterial) -> Shape {
        Shape::Box { p0, p1, material }
    }

    /// Create a new sphere, centered at `point` and with a `radius`. `material` is the material to
    /// use for this object.
    pub fn new_sphere(point: Point3<usize>, radius: f32, material: DrudeMaterial) -> Shape {
        Shape::Sphere {
            point,
            radius,
            material,
        }
    }

    /// Generate the `BoundingBox` for this shape. The shape can touch the bouningbox, but it is
    /// guaranteed not to be outside it.
    pub fn bbox(&self) -> BoundingBox {
        match self {
            Shape::Box { p0, p1, .. } => BoundingBox::new(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z),
            Shape::Sphere { point, radius, .. } => BoundingBox::new(
                (point.x as f32 - radius).floor() as usize,
                (point.y as f32 - radius).floor() as usize,
                (point.z as f32 - radius).floor() as usize,
                (point.x as f32 + radius).ceil() as usize,
                (point.y as f32 + radius).ceil() as usize,
                (point.z as f32 + radius).ceil() as usize,
            ),
        }
    }

    /// Draw the permittivity corresponding to this shape on the provided `ScalarField`. Since the
    /// permittivity can be dependent on the frequency, it should be passed as well.
    pub fn draw_permittivity(&self, field: &mut ScalarField, frequency: f32) {
        match *self {
            Shape::Box { material, .. } => {
                let permittivity = material.permittivity(frequency);
                let bbox = self.bbox();
                for x in bbox.x0..bbox.x1 {
                    for y in bbox.y0..bbox.y1 {
                        for z in bbox.z0..bbox.z1 {
                            field[(x, y, z)] = permittivity;
                        }
                    }
                }
            }
            Shape::Sphere {
                point,
                radius,
                material,
            } => {
                let permittivity = material.permittivity(frequency);
                let bbox = self.bbox();
                for x in bbox.x0..bbox.x1 {
                    for y in bbox.y0..bbox.y1 {
                        for z in bbox.z0..bbox.z1 {
                            if Vector3::new(
                                (point.x - x) as f32,
                                (point.y - y) as f32,
                                (point.z - z) as f32,
                            ).norm() <= radius
                            {
                                field[(x, y, z)] = permittivity;
                            }
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Point3;
    use world::shape::Shape;

    #[test]
    fn test_sphere_bbox_integer_radius() {
        let sphere = Shape::new_sphere(Point3::new(5, 5, 5), 2.0, "gold".parse().unwrap());
        let bbox = sphere.bbox();
        assert_eq!(bbox.x0, 3);
        assert_eq!(bbox.y0, 3);
        assert_eq!(bbox.z0, 3);
        assert_eq!(bbox.x1, 7);
        assert_eq!(bbox.y1, 7);
        assert_eq!(bbox.z1, 7);
    }

    #[test]
    fn test_sphere_bbox_float_radius() {
        let sphere = Shape::new_sphere(Point3::new(5, 5, 5), 2.5, "gold".parse().unwrap());
        let bbox = sphere.bbox();
        assert_eq!(bbox.x0, 2);
        assert_eq!(bbox.y0, 2);
        assert_eq!(bbox.z0, 2);
        assert_eq!(bbox.x1, 8);
        assert_eq!(bbox.y1, 8);
        assert_eq!(bbox.z1, 8);
    }

    #[test]
    fn test_box_bbox_integer_radius() {
        let box_shape = Shape::new_box(
            Point3::new(5, 5, 5),
            Point3::new(7, 7, 7),
            "gold".parse().unwrap(),
        );
        let bbox = box_shape.bbox();
        assert_eq!(bbox.x0, 5);
        assert_eq!(bbox.y0, 5);
        assert_eq!(bbox.z0, 5);
        assert_eq!(bbox.x1, 7);
        assert_eq!(bbox.y1, 7);
        assert_eq!(bbox.z1, 7);
    }
}
