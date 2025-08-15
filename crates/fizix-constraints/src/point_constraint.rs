use fizix_core::{BodyHandle, BodySet, Constraint, Precision, EPSILON, EPSILON_SQUARED};
use nalgebra::Point3;

pub struct PointConstraint {
    body_a: BodyHandle, body_b: BodyHandle,

    local_point_a: Point3<Precision>,
    local_point_b: Point3<Precision>
}
impl PointConstraint {
    pub fn new(
        body_a: BodyHandle, body_b: BodyHandle,

        local_point_a: Point3<Precision>,
        local_point_b: Point3<Precision>
    ) -> Self {
        Self {
            body_a, body_b,

            local_point_a,
            local_point_b
        }
    }
}

impl Constraint for PointConstraint {
    fn project(&mut self, bodies: &mut BodySet) {
        let body_a = self.body_a.0;
        let body_b = self.body_b.0;

        let world_point_a = bodies.transform[body_a].transform_point(&self.local_point_a);
        let world_point_b = bodies.transform[body_b].transform_point(&self.local_point_b);

        let difference = world_point_a - world_point_b;
        let distance_squared = difference.norm_squared();

        if distance_squared < EPSILON_SQUARED { return; }

        let distance = distance_squared.sqrt();
        let normal = difference / distance;

        let relative_point_a = world_point_a - bodies.position[body_a];
        let relative_point_b = world_point_b - bodies.position[body_b];

        let perpendicular_a = relative_point_a.cross(&normal);
        let perpendicular_b = relative_point_b.cross(&normal);

        let inverse_inertia_a = (bodies.inverse_inertia_tensor_world[body_a] * perpendicular_a).dot(&perpendicular_a);
        let inverse_inertia_b = (bodies.inverse_inertia_tensor_world[body_b] * perpendicular_b).dot(&perpendicular_b);

        let total_inverse_mass = bodies.inverse_mass[body_a] + bodies.inverse_mass[body_b] + inverse_inertia_a + inverse_inertia_b;

        if total_inverse_mass < EPSILON { return; }

        let lambda = -distance / total_inverse_mass;
        let translational_correction = normal * lambda;

        if bodies.has_finite_mass(body_a) {
            let rotational_correction = relative_point_a.cross(&translational_correction);

            bodies.position[body_a] += bodies.inverse_mass[body_a] * translational_correction;

            bodies.apply_rotation_delta(body_a, bodies.inverse_inertia_tensor_world[body_a] * rotational_correction);
            bodies.update_derived_data(body_a);
        }
        if bodies.has_finite_mass(body_b) {
            let rotational_correction = relative_point_b.cross(&translational_correction);

            bodies.position[body_b] -= bodies.inverse_mass[body_b] * translational_correction;

            bodies.apply_rotation_delta(body_b, bodies.inverse_inertia_tensor_world[body_b] * -rotational_correction);
            bodies.update_derived_data(body_b);
        }
    }
}