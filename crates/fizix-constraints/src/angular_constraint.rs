use fizix_core::{BodyHandle, BodySet, Constraint, CorrectionData, Precision};
use nalgebra::{UnitVector3, Vector3};

pub struct AngularConstraint {
    pub body_a: BodyHandle,
    pub body_b: BodyHandle,

    pub local_axis_a: UnitVector3<Precision>, // relative to body A
    pub local_axis_b: UnitVector3<Precision>, // relative to body B

    pub min_angle: Precision,
    pub max_angle: Precision,

    pub compliance: Precision // inverse stiffness
}

impl Default for AngularConstraint {
    fn default() -> Self {
        Self {
            body_a: BodyHandle::INVALID,
            body_b: BodyHandle::INVALID,

            local_axis_a: Vector3::x_axis(),
            local_axis_b: Vector3::x_axis(),

            min_angle: Precision::NEG_INFINITY,
            max_angle: Precision::INFINITY,

            compliance: 0.0
        }
    }
}

impl Constraint for AngularConstraint {
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

        if phi >= self.min_angle && phi <= self.max_angle { return None; }

        let error = phi.clamp(self.min_angle, self.max_angle);
        let normal = orthogonal / sin_theta;

        Some(CorrectionData::Rotational {
            handles: vec![self.body_a, self.body_b],
            axes: vec![normal, -normal],

            error,
            alpha: self.compliance
        })
    }
}