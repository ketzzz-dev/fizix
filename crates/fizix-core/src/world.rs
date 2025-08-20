use crate::{Body, BodyHandle, Constraint, Precision};
use nalgebra::{Matrix3, Point3, UnitQuaternion, Vector3};

pub struct World {
    pub bodies: Vec<Body>,
    pub constraints: Vec<Box<dyn Constraint>>,

    pub gravity: Vector3<Precision>,

    pub sub_steps: usize,
    pub constraint_iterations: usize
}

impl World {
    pub fn new(gravity: Vector3<Precision>, sub_steps: usize, constraint_iterations: usize) -> Self {
        Self {
            bodies: Vec::new(),
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
        let inverse_inertia_local = if is_mass_valid {
            inertia_tensor.try_inverse().unwrap_or(Matrix3::zeros())
        } else {
            Matrix3::zeros()
        };
        let inverse_inertia_world = rotation * inverse_inertia_local * rotation.transpose();

        self.bodies.push(Body {
            position,
            orientation,

            inverse_mass,
            inverse_inertia_local,
            inverse_inertia_world,

            ..Default::default()
        });

        BodyHandle::new(self.bodies.len() - 1)
    }

    pub fn add_constraint<C>(&mut self, constraint: C)
    where C: Constraint + 'static {
        self.constraints.push(Box::new(constraint));
    }

    pub fn step(&mut self, dt: Precision) {
        let sub_dt = dt / self.sub_steps as Precision;
        let inv_dt = 1.0 / sub_dt;

        for _ in 0..self.sub_steps {
            // integration
            for body in &mut self.bodies {
                if !body.has_finite_mass() { continue; }

                body.last_position = body.position;
                body.last_orientation = body.orientation;

                let linear_acc = self.gravity + body.inverse_mass * body.force;
                let angular_acc = body.inverse_inertia_world * body.torque;

                body.force = Vector3::zeros();
                body.torque = Vector3::zeros();

                body.linear_velocity += linear_acc * sub_dt;
                body.angular_velocity += angular_acc * sub_dt;

                body.position += body.linear_velocity * sub_dt;

                body.apply_rotation_delta(body.angular_velocity * sub_dt);
            }

            // constraint solve
            for _ in 0..self.constraint_iterations {
                for constraint in &mut self.constraints {
                    constraint.solve(&mut self.bodies, sub_dt);
                }
            }

            // velocity update
            for body in &mut self.bodies {
                if !body.has_finite_mass() { continue; }

                let delta_q = body.orientation * body.last_orientation.conjugate();

                body.linear_velocity = (body.position - body.last_position) * inv_dt;
                body.angular_velocity = delta_q.scaled_axis() * inv_dt;
            }
        }
    }
}