use crate::{BodyHandle, BodySet, Constraint, Precision};
use itertools::izip;
use nalgebra::{Matrix3, Point3, UnitQuaternion, Vector3};

pub struct World {
    pub bodies: BodySet,
    pub constraints: Vec<Box<dyn Constraint>>,

    gravity: Vector3<Precision>,

    sub_steps: usize,
    constraint_iterations: usize
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

    pub fn add_constraint(&mut self, constraint: impl Constraint + 'static) {
        self.constraints.push(Box::new(constraint));
    }

    pub fn step(&mut self, dt: Precision) {
        let sub_dt = dt / self.sub_steps as Precision;
        let inv_dt = 1.0 / sub_dt;

        for _ in 0..self.sub_steps {
            // integration
            for i in 0..self.bodies.position.len() {
                if !self.bodies.has_finite_mass(i) { continue; }

                self.bodies.last_position[i] = self.bodies.position[i];
                self.bodies.last_orientation[i] = self.bodies.orientation[i];

                let linear_acc = self.gravity + self.bodies.inverse_mass[i] * self.bodies.force[i];
                let angular_acc = self.bodies.inverse_inertia_tensor_world[i] * self.bodies.torque[i];

                self.bodies.force[i].fill(0.0);
                self.bodies.torque[i].fill(0.0);
                
                self.bodies.linear_velocity[i] += linear_acc * sub_dt;
                self.bodies.angular_velocity[i] += angular_acc * sub_dt;

                self.bodies.position[i] += self.bodies.linear_velocity[i] * sub_dt;
                
                self.bodies.apply_rotation_delta(i, self.bodies.angular_velocity[i] * sub_dt);
                self.bodies.update_derived_data(i);
            }

            // constraint solve
            let mut lambdas = vec![0.0; self.constraints.len()];

            for _ in 0..self.constraint_iterations {
                for (constraint, lambda) in izip!(&mut self.constraints, &mut lambdas) {
                    constraint.solve(&mut self.bodies, lambda, sub_dt);
                }
            }

            // velocity update
            for i in 0..self.bodies.position.len() {
                if !self.bodies.has_finite_mass(i) { continue; }

                let delta_q = self.bodies.orientation[i] * self.bodies.last_orientation[i].conjugate();

                self.bodies.linear_velocity[i] = (self.bodies.position[i] - self.bodies.last_position[i]) * inv_dt;
                self.bodies.angular_velocity[i] = delta_q.scaled_axis() * inv_dt;
            }
        }
    }
}