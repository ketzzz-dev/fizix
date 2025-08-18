use fizix_core::{BodyHandle, BodySet, Constraint, Precision, AXES, EPSILON};
use nalgebra::Point3;

pub struct BallJointConstraint {
    pub body_a: BodyHandle,
    pub body_b: BodyHandle,

    pub local_anchor_a: Point3<Precision>,
    pub local_anchor_b: Point3<Precision>,
}

impl BallJointConstraint {
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,

        local_anchor_a: Point3<Precision>,
        local_anchor_b: Point3<Precision>
    ) -> Self {
        Self {
            body_a,
            body_b,

            local_anchor_a,
            local_anchor_b,
        }
    }
}

impl Constraint for BallJointConstraint {
    fn solve(&mut self, bodies: &mut BodySet) {
        let body_a = *self.body_a;
        let body_b = *self.body_b;

        let world_anchor_a = bodies.position[body_a] + bodies.orientation[body_a].transform_point(&self.local_anchor_a).coords;
        let world_anchor_b = bodies.position[body_b] + bodies.orientation[body_b].transform_point(&self.local_anchor_b).coords;

        let error = world_anchor_b - world_anchor_a;

        let inverse_mass_a = bodies.inverse_mass[body_a];
        let inverse_mass_b = bodies.inverse_mass[body_b];

        let inverse_inertia_tensor_a = bodies.inverse_inertia_tensor_world[body_a];
        let inverse_inertia_tensor_b = bodies.inverse_inertia_tensor_world[body_b];

        for i in 0..3 {
            let j_xa = AXES[i];
            let j_xb = -AXES[i];

            let j_ta = -world_anchor_a.coords.cross(&j_xa);
            let j_tb = world_anchor_b.coords.cross(&j_xb);

            let w_a = inverse_mass_a * j_xa.norm_squared() + (j_ta.transpose() * inverse_inertia_tensor_a * j_ta).trace(); 
            let w_b = inverse_mass_b * j_xb.norm_squared() + (j_tb.transpose() * inverse_inertia_tensor_b * j_tb).trace();

            let denom = w_a + w_b;

            if denom < EPSILON { continue; }

            let lambda = -error[i] / denom;

            if bodies.has_finite_mass(body_a) {
                bodies.position[body_a] += bodies.inverse_mass[body_a] * j_xa * lambda;

                let d_theta = bodies.inverse_inertia_tensor_world[body_a] * j_ta * lambda;

                bodies.apply_rotation_delta(body_a, d_theta);
                bodies.update_derived_data(body_a);
            }
            if bodies.has_finite_mass(body_b) {
                bodies.position[body_b] += bodies.inverse_mass[body_b] * j_xb * lambda;

                let d_theta = bodies.inverse_inertia_tensor_world[body_b] * j_tb * lambda;

                bodies.apply_rotation_delta(body_b, d_theta);
                bodies.update_derived_data(body_b);
            }
        }
    }
}