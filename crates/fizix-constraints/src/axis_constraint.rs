use fizix_core::{BodyHandle, BodySet, Constraint, CorrectionData, Precision, EPSILON_SQUARED};
use nalgebra::Vector3;

pub struct AxisConstraint {
    body_a: BodyHandle, body_b: BodyHandle,

    local_axis_a: Vector3<Precision>,
    local_axis_b: Vector3<Precision>
}
impl AxisConstraint {
    pub fn new(
        body_a: BodyHandle, body_b: BodyHandle,

        local_axis_a: Vector3<Precision>,
        local_axis_b: Vector3<Precision>
    ) -> Self {
        Self {
            body_a, body_b,

            local_axis_a,
            local_axis_b
        }
    }
}

impl Constraint for AxisConstraint {
    fn compute_correction_data(&self, bodies: &BodySet) -> Option<CorrectionData> {
        let body_a: usize = self.body_a.0;
        let body_b: usize = self.body_b.0;

        let world_axis_a = bodies.transform[body_a].transform_vector(&self.local_axis_a);
        let world_axis_b = bodies.transform[body_b].transform_vector(&self.local_axis_b);

        let orthogonal = world_axis_a.cross(&world_axis_b);
        let magnitude_squared = orthogonal.norm_squared();

        if magnitude_squared < EPSILON_SQUARED { return None; }

        let magnitude = magnitude_squared.sqrt();
        let axis = orthogonal / magnitude;

        Some(CorrectionData::Rotational {
            error: magnitude,
            body_handles: vec![self.body_a, self.body_b],
            axes: vec![axis, -axis]
        })
    }
}