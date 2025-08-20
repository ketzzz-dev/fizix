use crate::{Body, Precision};

pub trait Constraint {
    fn solve(&mut self, bodies: &mut [Body], dt: Precision);
}