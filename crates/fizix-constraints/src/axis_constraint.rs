use fizix_core::{BodyHandle, BodySet, Constraint, Precision, EPSILON, EPSILON_SQUARED};
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
    fn project(&mut self, bodies: &mut BodySet) {
        let body_a: usize = self.body_a.0;
        let body_b: usize = self.body_b.0;

        let world_axis_a = bodies.transform[body_a].transform_vector(&self.local_axis_a);
        let world_axis_b = bodies.transform[body_b].transform_vector(&self.local_axis_b);

        let orthogonal = world_axis_a.cross(&world_axis_b);
        let magnitude_squared = orthogonal.norm_squared();

        if magnitude_squared < EPSILON_SQUARED { return; }

        let magnitude = magnitude_squared.sqrt();
        let normal = orthogonal / magnitude;

        let inverse_inertia_a = (bodies.inverse_inertia_tensor_world[body_a] * normal).dot(&normal);
        let inverse_inertia_b = (bodies.inverse_inertia_tensor_world[body_b] * normal).dot(&normal);

        let total_inverse_mass = inverse_inertia_a + inverse_inertia_b;

        if total_inverse_mass < EPSILON { return; }

        let lambda = -magnitude / total_inverse_mass;
        let rotational_correction = normal * lambda;

        if bodies.has_finite_mass(body_a) {
            bodies.apply_rotation_delta(body_a, bodies.inverse_inertia_tensor_world[body_a] * rotational_correction);
            bodies.update_derived_data(body_a);
        }
        if bodies.has_finite_mass(body_b) {
            bodies.apply_rotation_delta(body_b, bodies.inverse_inertia_tensor_world[body_b] * -rotational_correction);
            bodies.update_derived_data(body_b);
        }
    }
}