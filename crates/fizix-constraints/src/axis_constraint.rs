use fizix_core::{body::BodySet, constraint::Constraint, Precision, EPSILON};
use nalgebra::Vector3;

pub struct AxisConstraint {
    body_a: usize, body_b: usize,

    local_axis_a: Vector3<Precision>,
    local_axis_b: Vector3<Precision>
}
impl AxisConstraint {
    pub fn new(
        body_a: usize, body_b: usize,

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
    fn project(&self, bodies: &mut BodySet) {
        if !bodies.has_finite_mass(self.body_a) && !bodies.has_finite_mass(self.body_b) { return };

        let world_axis_a = bodies.transform[self.body_a].transform_vector(&self.local_axis_a);
        let world_axis_b = bodies.transform[self.body_b].transform_vector(&self.local_axis_b);

        let orthogonal = world_axis_a.cross(&world_axis_b);
        let magnitude = orthogonal.norm();

        if magnitude < EPSILON { return; }

        let normal = orthogonal / magnitude;
        let normal_transpose = normal.transpose();

        let inverse_inertia_a = (normal_transpose * bodies.inverse_inertia_tensor_world[self.body_a] * normal).x;
        let inverse_inertia_b = (normal_transpose * bodies.inverse_inertia_tensor_world[self.body_b] * normal).x;

        let total_inverse_mass = inverse_inertia_a + inverse_inertia_b;

        let lambda = -magnitude / total_inverse_mass;
        let rotational_correction = normal * lambda;

        if bodies.has_finite_mass(self.body_a) {
            bodies.apply_rotation_delta(self.body_a, bodies.inverse_inertia_tensor_world[self.body_a] * rotational_correction);
            bodies.update_derived_data(self.body_a);
        }
        if bodies.has_finite_mass(self.body_b) {
            bodies.apply_rotation_delta(self.body_b, bodies.inverse_inertia_tensor_world[self.body_b] * -rotational_correction);
            bodies.update_derived_data(self.body_b);
        }
    }
}