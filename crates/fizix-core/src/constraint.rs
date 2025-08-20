use crate::{BodySet, Precision};

pub trait Constraint {
    fn solve(&mut self, bodies: &mut BodySet, dt: Precision);
}