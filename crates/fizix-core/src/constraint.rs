use crate::BodySet;

pub trait Constraint {
    fn solve(&mut self, bodies: &mut BodySet);
}