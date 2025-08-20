use fizix_core::{get_two_mut, Body, BodyHandle, Constraint, Precision, EPSILON, EPSILON_SQUARED};
use nalgebra::Point3;

pub struct SphericalJoint {
    pub body_a: BodyHandle,
    pub body_b: BodyHandle,

    pub local_anchor_a: Point3<Precision>,
    pub local_anchor_b: Point3<Precision>,

    pub compliance: Precision, // inverse stiffness
}

impl SphericalJoint {
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,

        local_anchor_a: Point3<Precision>,
        local_anchor_b: Point3<Precision>,

        compliance: Precision
    ) -> Self {
        Self {
            body_a,
            body_b,

            local_anchor_a,
            local_anchor_b,

            compliance
        }
    }
}

impl Constraint for SphericalJoint {
    fn solve(&mut self, bodies: &mut [Body], dt: Precision) {
        let (body_a, body_b) = if let Some((a, b)) = get_two_mut(bodies, *self.body_a, *self.body_b) {
            (a, b)
        } else {
            return; // Invalid handles
        };

        let pos_a = body_a.position;
        let pos_b = body_b.position;

        let orient_a = body_a.orientation;
        let orient_b = body_b.orientation;

        let inv_mass_a = body_a.inverse_mass;
        let inv_mass_b = body_b.inverse_mass;

        let inv_inertia_a = body_a.inverse_inertia_world;
        let inv_inertia_b = body_b.inverse_inertia_world;

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
            let w_a = inv_mass_a + (inv_inertia_a * perp_a).dot(&perp_a);
            let w_b = inv_mass_b + (inv_inertia_b * perp_b).dot(&perp_b);

            let denom = w_a + w_b + alpha_t;

            if denom > EPSILON {
                // compute the impulse to apply
                let d_lambda = -distance / denom;
                let d_x = delta * d_lambda;

                if body_a.has_finite_mass() {
                    let d_theta = -r_a.coords.cross(&d_x);

                    body_a.position -= inv_mass_a * d_x;

                    body_a.apply_rotation_delta(inv_inertia_a * d_theta);
                }
                if body_b.has_finite_mass() {
                    let d_theta = r_b.coords.cross(&d_x);

                    body_b.position += inv_mass_b * d_x;

                    body_b.apply_rotation_delta(inv_inertia_b * d_theta);
                }
            }
        }
    }
}