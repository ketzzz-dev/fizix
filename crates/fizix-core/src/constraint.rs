use itertools::izip;
use nalgebra::{Point3, UnitVector3, Vector3};
use crate::{BodyHandle, BodySet, Precision, EPSILON};

pub trait Constraint {
    fn compute_correction(&self, bodies: &BodySet) -> Option<CorrectionData>;

    fn solve(&mut self, bodies: &mut BodySet, lambda: &mut Precision, dt: Precision) {
        if let Some(correction) = self.compute_correction(bodies) {
            correction.apply_correction(bodies, lambda, dt);
        }
    }
}

pub enum CorrectionData {
    Translational {
        handles: Vec<BodyHandle>,
        relative_points: Vec<Point3<Precision>>,
        normals: Vec<UnitVector3<Precision>>,

        error: Precision,
        alpha: Precision // compliance
    },
    Rotational {
        handles: Vec<BodyHandle>,
        axes: Vec<UnitVector3<Precision>>,

        error: Precision,
        alpha: Precision // compliance
    }
}

impl CorrectionData {
    pub fn apply_correction(&self, bodies: &mut BodySet, lambda: &mut Precision, dt: Precision) {
        match self {
            CorrectionData::Translational { handles, relative_points, normals, error, alpha } =>
                Self::apply_translational_correction(bodies, lambda, handles, relative_points, normals, *error, *alpha, dt),
            CorrectionData::Rotational { handles, axes, error, alpha } =>
                Self::apply_rotational_correction(bodies, lambda, handles, axes, *error, *alpha, dt)
        }
    }

    fn apply_translational_correction(
        bodies: &mut BodySet,
        lambda: &mut Precision,

        handles: &[BodyHandle],
        relative_points: &[Point3<Precision>],
        normals: &[UnitVector3<Precision>],

        error: Precision,
        alpha: Precision,

        dt: Precision
    ) {
        let alpha_tilde = if alpha > 0.0 {
            alpha / (dt * dt)
        } else {
            0.0
        };
        let mut total_inverse_mass = alpha_tilde;

        for (handle, relative_point, normal) in izip!(handles, relative_points, normals) {
            let body = **handle;

            let perpendicular = relative_point.coords.cross(normal);

            total_inverse_mass += bodies.inverse_mass[body] + (bodies.inverse_inertia_tensor_world[body] * perpendicular).dot(&perpendicular);
        }

        if total_inverse_mass < EPSILON { return; }

        let d_lambda = -(error + *lambda * alpha_tilde) / total_inverse_mass;

        *lambda += d_lambda;

        for (handle, relative_point, normal) in izip!(handles, relative_points, normals) {
            let body = **handle;

            if !bodies.has_finite_mass(body) { continue; }

            let correction_impulse = normal.into_inner() * d_lambda;
            let rotational_correction = relative_point.coords.cross(&correction_impulse);

            bodies.position[body] += bodies.inverse_mass[body] * correction_impulse;

            bodies.apply_rotation_delta(body, bodies.inverse_inertia_tensor_world[body] * rotational_correction);
            bodies.update_derived_data(body);
        }
    }

    fn apply_rotational_correction(
        bodies: &mut BodySet,
        lambda: &mut Precision,

        handles: &[BodyHandle],
        axes: &[UnitVector3<Precision>],

        error: Precision,
        alpha: Precision,

        dt: Precision
    ) {
        let alpha_tilde = if alpha > 0.0 {
            alpha / (dt * dt)
        } else {
            0.0
        };
        let mut total_inverse_mass = alpha_tilde;

        for (handle, axis) in izip!(handles, axes) {
            let body = **handle;

            total_inverse_mass += (bodies.inverse_inertia_tensor_world[body] * axis.into_inner()).dot(&axis);
        }

        if total_inverse_mass < EPSILON { return; }

        let d_lambda = -(error + *lambda * alpha_tilde) / total_inverse_mass;

        *lambda += d_lambda;

        for (handle, axis) in izip!(handles, axes) {
            let body = **handle;

            if !bodies.has_finite_mass(body) { continue; }

            let correction_impulse = axis.into_inner() * d_lambda;

            bodies.apply_rotation_delta(body, bodies.inverse_inertia_tensor_world[body] * correction_impulse);
            bodies.update_derived_data(body);
        }
    }
}