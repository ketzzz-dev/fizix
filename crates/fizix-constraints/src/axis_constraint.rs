use fizix_core::{BodyHandle, BodySet, Constraint, CorrectionData, Precision, EPSILON};
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

        let cross = world_axis_a.cross(&world_axis_b);
        
        let sin_theta = cross.norm();
        let cos_theta = world_axis_a.dot(&world_axis_b).clamp(-1.0, 1.0);

        let angle_error = sin_theta.atan2(cos_theta);

        if angle_error.abs() < EPSILON { return None; }

        let axis = cross / sin_theta;

        Some(CorrectionData::Rotational {
            error: angle_error,
            body_handles: vec![self.body_a, self.body_b],
            axes: vec![axis, -axis]
        })
    }
}