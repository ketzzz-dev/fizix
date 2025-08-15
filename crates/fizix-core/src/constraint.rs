use crate::BodySet;

pub trait Constraint {
    fn project(&mut self, bodies: &mut BodySet);
}