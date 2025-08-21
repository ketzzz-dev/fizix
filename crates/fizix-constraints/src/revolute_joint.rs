use fizix_core::{BodyHandle, BodySet, Constraint, EPSILON, EPSILON_SQUARED, Precision};
use nalgebra::{Point3, UnitVector3};

pub struct RevoluteJoint {
    pub body_a: BodyHandle,
    pub body_b: BodyHandle,

    pub local_anchor_a: Point3<Precision>,
    pub local_anchor_b: Point3<Precision>,

    pub local_axis_a: UnitVector3<Precision>,
    pub local_axis_b: UnitVector3<Precision>,

    pub compliance: Precision,

    lambda_pos: Precision,
    lambda_rot: Precision
}

impl RevoluteJoint {
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,

        local_anchor_a: Point3<Precision>,
        local_anchor_b: Point3<Precision>,

        local_axis_a: UnitVector3<Precision>,
        local_axis_b: UnitVector3<Precision>,

        compliance: Precision,
    ) -> Self {
        Self {
            body_a,
            body_b,

            local_anchor_a,
            local_anchor_b,

            local_axis_a,
            local_axis_b,

            compliance,

            lambda_pos: 0.0,
            lambda_rot: 0.0
        }
    }
}

impl Constraint for RevoluteJoint {
    fn solve(&mut self, bodies: &mut BodySet, dt: Precision) {
        let body_a = *self.body_a;
        let body_b = *self.body_b;

        let pos_a = bodies.position[body_a];
        let pos_b = bodies.position[body_b];

        let orient_a = bodies.orientation[body_a];
        let orient_b = bodies.orientation[body_b];

        let inv_mass_a = bodies.inverse_mass[body_a];
        let inv_mass_b = bodies.inverse_mass[body_b];

        let inv_inertia_a = bodies.inverse_inertia_tensor_world[body_a];
        let inv_inertia_b = bodies.inverse_inertia_tensor_world[body_b];

        // relative space anchor points
        let r_a = orient_a.transform_point(&self.local_anchor_a);
        let r_b = orient_b.transform_point(&self.local_anchor_b);

        // absolute space anchor points
        let p_a = pos_a + r_a.coords;
        let p_b = pos_b + r_b.coords;

        let difference = p_b - p_a;
        let distance_sq = difference.norm_squared();

        // alpha tilde = alpha / dt^2
        let alpha_t = if self.compliance > 0.0 {
            self.compliance / (dt * dt)
        } else {
            0.0
        };

        if distance_sq > EPSILON_SQUARED {
            let distance = distance_sq.sqrt();
            let delta = difference / distance;

            let perp_a = r_a.coords.cross(&delta);
            let perp_b = r_b.coords.cross(&delta);

            // generalized inverse mass
            let w = inv_mass_a + inv_mass_b + (inv_inertia_a * perp_a).dot(&perp_a) + (inv_inertia_b * perp_b).dot(&perp_b);

            if w > EPSILON {
                // compute the impulse to apply
                let d_lambda = -(distance + alpha_t * self.lambda_pos) / (w + alpha_t);
                let d_x = delta * d_lambda;

                self.lambda_pos += d_lambda;

                if bodies.has_finite_mass(body_a) {
                    let d_theta = -r_a.coords.cross(&d_x);

                    bodies.position[body_a] -= inv_mass_a * d_x;

                    bodies.apply_rotation_delta(body_a, inv_inertia_a * d_theta);
                    bodies.update_derived_data(body_a);
                }
                if bodies.has_finite_mass(body_b) {
                    let d_theta = r_b.coords.cross(&d_x);

                    bodies.position[body_b] += inv_mass_b * d_x;

                    bodies.apply_rotation_delta(body_b, inv_inertia_b * d_theta);
                    bodies.update_derived_data(body_b);
                }
            }
        }

        // recompute orientaton and inertia tensors after position update
        let orient_a = bodies.orientation[body_a];
        let orient_b = bodies.orientation[body_b];

        let inv_inertia_a = bodies.inverse_inertia_tensor_world[body_a];
        let inv_inertia_b = bodies.inverse_inertia_tensor_world[body_b];

        // world space axis vectors
        let u_a = orient_a.transform_vector(&self.local_axis_a);
        let u_b = orient_b.transform_vector(&self.local_axis_b);

        let cross = u_a.cross(&u_b);
        let sin_sq_theta = cross.norm_squared();

        if sin_sq_theta < EPSILON_SQUARED { return; } // the axes are parallel, no rotation needed

        let sin_theta = sin_sq_theta.sqrt();
        let cos_theta = u_a.dot(&u_b);
        let angle = sin_theta.atan2(cos_theta);

        if angle.abs() > EPSILON {
            let axis = cross / sin_theta;

            // (I_a * u) . u + (I_b * u) . u = ((I_a + I_b) * u) . u
            let w = ((inv_inertia_a + inv_inertia_b) * axis).dot(&axis);

            if w > EPSILON {
                // compute the impulse to apply
                let d_lambda = -(angle + self.lambda_rot) / w;
                let d_theta = axis * d_lambda;

                self.lambda_rot += d_lambda;

                if bodies.has_finite_mass(body_a) {
                    bodies.apply_rotation_delta(body_a, inv_inertia_a * -d_theta);
                    bodies.update_derived_data(body_a);
                }
                if bodies.has_finite_mass(body_b) {
                    bodies.apply_rotation_delta(body_b, inv_inertia_b * d_theta);
                    bodies.update_derived_data(body_b);
                }
            }
        }
    }
}
