use nalgebra::Vector3;

mod world;
mod body;
mod constraint;

pub type Precision = f64;

pub const EPSILON: Precision = 1e-9;
pub const EPSILON_SQUARED: Precision = EPSILON * EPSILON;

pub const AXES: [Vector3<Precision>; 3] = [
    Vector3::new(1.0, 0.0, 0.0),
    Vector3::new(0.0, 1.0, 0.0),
    Vector3::new(0.0, 0.0, 1.0),
];

pub use world::*;
pub use body::*;
pub use constraint::*;