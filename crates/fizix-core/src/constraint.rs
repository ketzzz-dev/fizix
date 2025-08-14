use crate::body::BodySet;

pub trait Constraint {
    fn project(&self, bodies: &mut BodySet);
}