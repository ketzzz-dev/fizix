use nalgebra::{Point3, Vector3};

use crate::{BodyHandle, BodySet, Precision, EPSILON};

pub enum CorrectionData {
    Translational {
        handles: Vec<BodyHandle>,
        relative_points: Vec<Point3<Precision>>,
        normals: Vec<Vector3<Precision>>,

        error: Precision,
        alpha: Precision // compliance
    },
    Rotational {
        handles: Vec<BodyHandle>,
        axes: Vec<Vector3<Precision>>,

        error: Precision,
        alpha: Precision // compliance
    }
}

pub trait Constraint {
    fn compute_correction(&self, bodies: &BodySet) -> Option<CorrectionData>;

    fn solve(&mut self, bodies: &mut BodySet, dt: Precision) {
        let correction_data = match self.compute_correction(bodies) {
            Some(data) => data,
            None => return
        };

        match correction_data {
            CorrectionData::Translational { handles, relative_points, normals, error, alpha } =>
                self.apply_translational_correction(bodies, &handles, &relative_points, &normals, error, alpha, dt),
            CorrectionData::Rotational { handles, axes, error, alpha } =>
                self.apply_rotational_correction(bodies, &handles, &axes, error, alpha, dt)
        }
    }

    fn apply_translational_correction(
        &self,
        bodies: &mut BodySet,

        handles: &[BodyHandle],
        relative_points: &[Point3<Precision>],
        normals: &[Vector3<Precision>],

        error: Precision,
        alpha: Precision,

        dt: Precision
    ) {
        let mut total_inverse_mass = if alpha > 0.0 {
            1.0 / (alpha * dt * dt)
        } else {
            0.0
        };

        for (i, &handle) in handles.iter().enumerate() {
            let body = *handle;

            let perpendicular = relative_points[i].coords.cross(&normals[i]);

            total_inverse_mass += bodies.inverse_mass[body] + (bodies.inverse_inertia_tensor_world[body] * perpendicular).dot(&perpendicular);
        }

        if total_inverse_mass < EPSILON { return; }

        let lambda = -error / total_inverse_mass;

        for (i, &handle) in handles.iter().enumerate() {
            let body = *handle;

            if !bodies.has_finite_mass(body) { continue; }

            let correction_impulse = normals[i] * lambda;
            let rotational_correction = relative_points[i].coords.cross(&correction_impulse);

            bodies.position[body] += bodies.inverse_mass[body] * correction_impulse;

            bodies.apply_rotation_delta(body, bodies.inverse_inertia_tensor_world[body] * rotational_correction);
        }
    }

    fn apply_rotational_correction(
        &self,
        bodies: &mut BodySet,

        handles: &[BodyHandle],
        axes: &[Vector3<Precision>],

        error: Precision,
        alpha: Precision,

        dt: Precision
    ) {
        let mut total_inverse_mass = if alpha > 0.0 {
            1.0 / (alpha * dt * dt)
        } else {
            0.0
        };

        for (i, &handle) in handles.iter().enumerate() {
            let body = *handle;

            total_inverse_mass += (bodies.inverse_inertia_tensor_world[body] * axes[i]).dot(&axes[i]);
        }

        if total_inverse_mass < EPSILON { return; }

        let lambda = -error / total_inverse_mass;

        for (i, &handle) in handles.iter().enumerate() {
            let body = *handle;

            if !bodies.has_finite_mass(body) { continue; }

            let correction_impulse = axes[i] * lambda;

            bodies.apply_rotation_delta(body, bodies.inverse_inertia_tensor_world[body] * correction_impulse);
        }
    }
}