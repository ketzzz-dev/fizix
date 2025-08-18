use crate::{BodyHandle, BodySet, Constraint, Precision};
use nalgebra::{Matrix3, Point3, UnitQuaternion, Vector3};

pub struct World {
    pub bodies: BodySet,
    pub constraints: Vec<Box<dyn Constraint>>,

    pub gravity: Vector3<Precision>,

    pub sub_steps: usize,
    pub constraint_iterations: usize
}

impl World {
    pub fn new(gravity: Vector3<Precision>, sub_steps: usize, constraint_iterations: usize) -> Self {
        Self {
            bodies: BodySet::new(),
            constraints: Vec::new(),

            gravity, sub_steps, constraint_iterations
        }
    }

    pub fn add_body(
        &mut self,
        position: Point3<Precision>,
        orientation: UnitQuaternion<Precision>,
        mass: Precision,
        inertia_tensor: Matrix3<Precision>
    ) -> BodyHandle {
        let rotation = orientation.to_rotation_matrix();

        let is_mass_valid = mass.is_finite() && mass > 0.0;
        let inverse_mass = if is_mass_valid { 1.0 / mass } else { 0.0 };
        let inverse_inertia_tensor = if is_mass_valid {
            inertia_tensor.try_inverse().unwrap_or(Matrix3::zeros())
        } else {
            Matrix3::zeros()
        };
        let inverse_inertia_tensor_world = rotation * inverse_inertia_tensor * rotation.transpose();

        self.bodies.position.push(position);
        self.bodies.orientation.push(orientation);

        self.bodies.last_position.push(position);
        self.bodies.last_orientation.push(orientation);

        self.bodies.linear_velocity.push(Vector3::zeros());
        self.bodies.angular_velocity.push(Vector3::zeros());

        self.bodies.force.push(Vector3::zeros());
        self.bodies.torque.push(Vector3::zeros());

        self.bodies.inverse_mass.push(inverse_mass);
        self.bodies.inverse_inertia_tensor_local.push(inverse_inertia_tensor);

        self.bodies.inverse_inertia_tensor_world.push(inverse_inertia_tensor_world);

        BodyHandle::new(self.bodies.position.len() - 1)
    }

    pub fn add_constraint<C>(&mut self, constraint: C)
    where C: Constraint + 'static {
        self.constraints.push(Box::new(constraint));
    }

    pub fn step(&mut self, delta_time: Precision) {
        let sub_delta_time = delta_time / self.sub_steps as Precision;
        let inverse_delta_time = 1.0 / sub_delta_time;

        for _ in 0..self.sub_steps {
            // pre-solve (integration)
            for i in 0..self.bodies.position.len() {
                if !self.bodies.has_finite_mass(i) { continue; }

                self.bodies.last_position[i] = self.bodies.position[i];
                self.bodies.last_orientation[i] = self.bodies.orientation[i];

                let linear_acceleration = self.gravity + self.bodies.inverse_mass[i] * self.bodies.force[i];
                let angular_acceleration = self.bodies.inverse_inertia_tensor_world[i] * self.bodies.torque[i];

                self.bodies.force[i] = Vector3::zeros();
                self.bodies.torque[i] = Vector3::zeros();
                
                self.bodies.linear_velocity[i] += linear_acceleration * sub_delta_time;
                self.bodies.angular_velocity[i] += angular_acceleration * sub_delta_time;

                self.bodies.position[i] += self.bodies.linear_velocity[i] * sub_delta_time;
                
                self.bodies.apply_rotation_delta(i, self.bodies.angular_velocity[i] * sub_delta_time);
                self.bodies.update_derived_data(i);
            }

            // solve (constraints)
            for _ in 0..self.constraint_iterations {
                for constraint in &mut self.constraints {
                    constraint.solve(&mut self.bodies);
                }
            }

            // post-solve (velocity update)
            for i in 0..self.bodies.position.len() {
                if !self.bodies.has_finite_mass(i) { continue; }

                let delta_orientation = self.bodies.orientation[i] * self.bodies.last_orientation[i].conjugate();

                self.bodies.linear_velocity[i] = (self.bodies.position[i] - self.bodies.last_position[i]) * inverse_delta_time;
                self.bodies.angular_velocity[i] = delta_orientation.scaled_axis() * inverse_delta_time;
            }
        }
    }
}