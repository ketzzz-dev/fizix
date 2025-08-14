use fizix_core::{body::BodySet, constraint::Constraint, Precision, EPSILON};
use nalgebra::Point3;

pub struct PointConstraint {
    body_a: usize, body_b: usize,

    local_point_a: Point3<Precision>,
    local_point_b: Point3<Precision>
}
impl PointConstraint {
    pub fn new(
        body_a: usize, body_b: usize,

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
    fn project(&self, bodies: &mut BodySet) {
        if !bodies.has_finite_mass(self.body_a) && !bodies.has_finite_mass(self.body_b) { return };

        let world_point_a = bodies.transform[self.body_a].transform_point(&self.local_point_a);
        let world_point_b = bodies.transform[self.body_b].transform_point(&self.local_point_b);

        let difference = world_point_a - world_point_b;
        let distance = difference.norm();

        if distance < EPSILON { return; }

        let normal = difference / distance;

        let relative_point_a = world_point_a - bodies.position[self.body_a];
        let relative_point_b = world_point_b - bodies.position[self.body_b];

        let perpendicular_a = relative_point_a.cross(&normal);
        let perpendicular_b = relative_point_b.cross(&normal);

        let inverse_inertia_a = (perpendicular_a.transpose() * bodies.inverse_inertia_tensor_world[self.body_a] * perpendicular_a).x;
        let inverse_inertia_b = (perpendicular_b.transpose() * bodies.inverse_inertia_tensor_world[self.body_b] * perpendicular_b).x;

        let total_inverse_mass = bodies.inverse_mass[self.body_a] + bodies.inverse_mass[self.body_b] + inverse_inertia_a + inverse_inertia_b;

        let lambda = -distance / total_inverse_mass;
        let translational_correction = normal * lambda;

        if bodies.has_finite_mass(self.body_a) {
            let rotational_correction = relative_point_a.cross(&translational_correction);

            bodies.position[self.body_a] += bodies.inverse_mass[self.body_a] * translational_correction;

            bodies.apply_rotation_delta(self.body_a, bodies.inverse_inertia_tensor_world[self.body_a] * rotational_correction);
            bodies.update_derived_data(self.body_a);
        }
        if bodies.has_finite_mass(self.body_b) {
            let rotational_correction = relative_point_b.cross(&translational_correction);

            bodies.position[self.body_b] -= bodies.inverse_mass[self.body_b] * translational_correction;

            bodies.apply_rotation_delta(self.body_b, bodies.inverse_inertia_tensor_world[self.body_b] * -rotational_correction);
            bodies.update_derived_data(self.body_b);
        }
    }
}