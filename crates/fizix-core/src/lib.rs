mod world;
mod body;
mod constraint;

pub type Precision = f64;

pub const EPSILON: Precision = 1e-9;
pub const EPSILON_SQUARED: Precision = EPSILON * EPSILON;

pub use world::*;
pub use body::*;
pub use constraint::*;