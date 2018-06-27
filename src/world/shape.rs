use std::fmt;
use fields::ScalarField;
use world::boundingbox::BoundingBox;
use nalgebra::*;

#[derive(PartialEq, Copy, Clone, Debug, Deserialize)]
pub enum Shape {
    Box {
        p0: Point3<usize>,
        p1: Point3<usize>
    },
    Sphere {
        point: Point3<usize>,
        radius: f32
    },
}

impl fmt::Display for Shape {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Shape::Box { p0, p1 } => write!(f, "Box: ({}, {}, {}), ({}, {}, {})", p0.x, p0.y, p0.z, p1.x, p1.y, p1.z),
            Shape::Sphere { point, radius, .. } => write!(
                f,
                "Sphere: ({}, {}, {}), with radius \
                 {}",
                point.x, point.y, point.z, radius
            ),
        }
    }
}

impl Shape {
    pub fn bbox(&self) -> BoundingBox {
        match self {
            Shape::Box { p0, p1 } => BoundingBox::new(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z),
            Shape::Sphere { point, radius } => BoundingBox::new(
                (point.x as f32 - radius).floor() as usize,
                (point.y as f32 - radius).floor() as usize,
                (point.z as f32 - radius).floor() as usize,
                (point.x as f32 + radius).floor() as usize,
                (point.y as f32 + radius).floor() as usize,
                (point.z as f32 + radius).floor() as usize,
            )
        }
    }

    pub fn draw_permitivity(&self, field: &mut ScalarField, permitivity: f32) {
        match self {
            &Shape::Box {..} => {
                let bbox = self.bbox();
                for x in bbox.x0..bbox.x1 {
                    for y in bbox.y0..bbox.y1 {
                        for z in bbox.z0..bbox.z1 {
                            field.set(
                                x as isize,
                                y as isize,
                                z as isize,
                                permitivity,
                            );
                        }
                    }
                }
            }
            &Shape::Sphere { point, radius } => {
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
                                    field.set(
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
}
