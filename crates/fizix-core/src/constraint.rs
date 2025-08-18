use nalgebra::{Point3, Vector3};

use crate::{BodyHandle, BodySet, Precision, EPSILON};

/// Represents the geometric and kinematic data required for constraint correction.
pub enum CorrectionData {
    Translational {
        error: Precision,

        body_handles: Vec<BodyHandle>,
        relative_points: Vec<Point3<Precision>>,
        normals: Vec<Vector3<Precision>>
    },
    Rotational {
        error: Precision,

        body_handles: Vec<BodyHandle>,
        axes: Vec<Vector3<Precision>>
    }
}

pub trait Constraint {
    fn compute_correction_data(&self, bodies: &BodySet) -> Option<CorrectionData>;

    fn apply_correction(&self, bodies: &mut BodySet) {
        let correction_data = match self.compute_correction_data(bodies) {
            Some(data) => data,
            None => return
        };

        match correction_data {
            CorrectionData::Translational { body_handles, relative_points, normals, error } => {
                self.apply_translational_correction(bodies, &body_handles, &relative_points, &normals, error);
            },
            CorrectionData::Rotational { body_handles, axes, error } => {
                self.apply_rotational_correction(bodies, &body_handles, &axes, error);
            }
        }
    }

    fn apply_translational_correction(&self, bodies: &mut BodySet, body_handles: &[BodyHandle], relative_points: &[Point3<Precision>], normals: &[Vector3<Precision>], error: Precision) {
        let mut total_inverse_mass = 0.0;

        for (i, &BodyHandle(body_index)) in body_handles.iter().enumerate() {
            let perpendicular = relative_points[i].coords.cross(&normals[i]);

            total_inverse_mass += bodies.inverse_mass[body_index] + (bodies.inverse_inertia_tensor_world[body_index] * perpendicular).dot(&perpendicular);
        }

        if total_inverse_mass < EPSILON { return; }

        let lambda = -error / total_inverse_mass;

        for (i, &BodyHandle(body_index)) in body_handles.iter().enumerate() {
            if !bodies.has_finite_mass(body_index) { continue; }

            let correction_impulse = normals[i] * lambda;
            let rotational_correction = relative_points[i].coords.cross(&correction_impulse);

            bodies.position[body_index] += bodies.inverse_mass[body_index] * correction_impulse;

            bodies.apply_rotation_delta(body_index, bodies.inverse_inertia_tensor_world[body_index] * rotational_correction);
            bodies.update_derived_data(body_index);
        }
    }

    fn apply_rotational_correction(&self, bodies: &mut BodySet, body_handles: &[BodyHandle], axes: &[Vector3<Precision>], error: Precision) {
        let mut total_inverse_inertia = 0.0;

        for (i, &BodyHandle(body_index)) in body_handles.iter().enumerate() {
            total_inverse_inertia += (bodies.inverse_inertia_tensor_world[body_index] * axes[i]).dot(&axes[i]);
        }

        if total_inverse_inertia < EPSILON { return; }

        let lambda = -error / total_inverse_inertia;

        for (i, &BodyHandle(body_index)) in body_handles.iter().enumerate() {
            if !bodies.has_finite_mass(body_index) { continue; }

            let correction_impulse = axes[i] * lambda;

            bodies.apply_rotation_delta(body_index, bodies.inverse_inertia_tensor_world[body_index] * correction_impulse);
            bodies.update_derived_data(body_index);
        }
    }
}