use std::ops::Deref;

use crate::Precision;
use nalgebra::{Matrix3, Point3, UnitQuaternion, Vector3};

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct BodyHandle(usize);

impl BodyHandle {
    pub fn new(index: usize) -> Self {
        Self(index)
    }
}

impl Deref for BodyHandle {
    type Target = usize;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

pub struct Body {
    pub position: Point3<Precision>,
    pub orientation: UnitQuaternion<Precision>,

    pub last_position: Point3<Precision>,
    pub last_orientation: UnitQuaternion<Precision>,

    pub linear_velocity: Vector3<Precision>,
    pub angular_velocity: Vector3<Precision>,

    pub force: Vector3<Precision>,
    pub torque: Vector3<Precision>,
   
    pub inverse_mass: Precision,
    pub inverse_inertia_local: Matrix3<Precision>,
    pub inverse_inertia_world: Matrix3<Precision>,
}

impl Default for Body {
    fn default() -> Self {
        Self {
            position: Point3::origin(),
            orientation: UnitQuaternion::identity(),

            last_position: Point3::origin(),
            last_orientation: UnitQuaternion::identity(),

            linear_velocity: Vector3::zeros(),
            angular_velocity: Vector3::zeros(),

            force: Vector3::zeros(),
            torque: Vector3::zeros(),

            inverse_mass: 0.0,
            inverse_inertia_local: Matrix3::identity(),
            inverse_inertia_world: Matrix3::identity(),
        }
    }
}

impl Body {
    #[inline]
    pub fn has_finite_mass(&self) -> bool {
        self.inverse_mass > 0.0
    }

    pub fn apply_rotation_delta(&mut self, rotation: Vector3<Precision>) {
        let q = self.orientation;

        self.orientation = UnitQuaternion::from_scaled_axis(rotation) * q;
        self.orientation.renormalize();

        // Update the world inertia tensor after changing the orientation
        let rot = self.orientation.to_rotation_matrix();
        
        self.inverse_inertia_world = rot * self.inverse_inertia_local * rot.transpose();
    }
}