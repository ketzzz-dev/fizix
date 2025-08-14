use fizix_core::{body::{BodyHandle, BodySet}, constraint::Constraint, Precision, EPSILON};
use nalgebra::{Point3, Vector3};

pub struct LineConstraint {
    body_a: BodyHandle, body_b: BodyHandle,

    local_point_a: Point3<Precision>,
    local_point_b: Point3<Precision>,

    local_axis: Vector3<Precision> // relative to B
}
impl LineConstraint {
    pub fn new(
        body_a: BodyHandle, body_b: BodyHandle,

        local_point_a: Point3<Precision>,
        local_point_b: Point3<Precision>,

        local_axis: Vector3<Precision>
    ) -> Self {
        Self {
            body_a, body_b,

            local_point_a,
            local_point_b,

            local_axis
        }
    }
}

impl Constraint for LineConstraint {
    fn project(&mut self, bodies: &mut BodySet) {
        let body_a: usize = self.body_a.0;
        let body_b = self.body_b.0;

        if !bodies.has_finite_mass(body_a) && !bodies.has_finite_mass(body_b) { return };

        let world_point_a = bodies.transform[body_a].transform_point(&self.local_point_a);
        let world_point_b = bodies.transform[body_b].transform_point(&self.local_point_b);
        let world_direction = bodies.transform[body_b].transform_vector(&self.local_axis);

        let projected_point = world_point_b + world_direction
            * world_direction.dot(&(world_point_a - world_point_b));

        let difference = world_point_a - projected_point;
        let distance = difference.norm();
        
        if distance < EPSILON { return; }

        let normal = difference / distance;

        let relative_point_a = world_point_a - bodies.position[body_a];
        let relative_point_b = projected_point - bodies.position[body_b];

        let perpendicular_a = relative_point_a.cross(&normal);
        let perpendicular_b = relative_point_b.cross(&normal);

        let inverse_inertia_a = (perpendicular_a.transpose() * bodies.inverse_inertia_tensor_world[body_a] * perpendicular_a).x;
        let inverse_inertia_b = (perpendicular_b.transpose() * bodies.inverse_inertia_tensor_world[body_b] * perpendicular_b).x;

        let total_inverse_mass = bodies.inverse_mass[body_a] + bodies.inverse_mass[body_b] + inverse_inertia_a + inverse_inertia_b;

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