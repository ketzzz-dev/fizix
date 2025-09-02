use fizix_core::{BodyHandle, BodySet, Constraint, CorrectionData, Precision, EPSILON, EPSILON_SQUARED};
use nalgebra::{UnitVector3, Vector3};

pub struct AxisConstraint {
    pub body_a: BodyHandle,
    pub body_b: BodyHandle,

    pub local_axis_a: UnitVector3<Precision>,
    pub local_axis_b: UnitVector3<Precision>,

    pub compliance: Precision // inverse stiffness
}

impl Default for AxisConstraint {
    fn default() -> Self {
        Self {
            body_a: BodyHandle::INVALID,
            body_b: BodyHandle::INVALID,

            local_axis_a: Vector3::x_axis(),
            local_axis_b: Vector3::x_axis(),

            compliance: 0.0
        }
    }
}

impl Constraint for AxisConstraint {
    fn compute_correction(&self, bodies: &BodySet) -> Option<CorrectionData> {
        let body_a: usize = *self.body_a;
        let body_b: usize = *self.body_b;

        let orient_a = bodies.orientation[body_a];
        let orient_b = bodies.orientation[body_b];

        // world space axes
        let u_a = orient_a.transform_vector(&self.local_axis_a);
        let u_b = orient_b.transform_vector(&self.local_axis_b);

        let orthogonal = u_a.cross(&u_b);

        let sin_theta = orthogonal.norm();
        let cos_theta = u_a.dot(&u_b);
        let phi = sin_theta.atan2(cos_theta);

        if phi.abs() < EPSILON { return None; }

        let axis = UnitVector3::new_unchecked(orthogonal / sin_theta);

        Some(CorrectionData::Rotational {
            handles: vec![self.body_a, self.body_b],
            axes: vec![-axis, axis],

            error: phi,
            alpha: self.compliance
        })
    }
}