use std::ops::Deref;

use crate::Precision;
use nalgebra::{Matrix3, Point3, UnitQuaternion, Vector3};

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct BodyHandle(usize);

impl BodyHandle {
    pub const INVALID: Self = Self(usize::MAX);

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

#[derive(Default)]
pub struct BodySet {
    pub position: Vec<Point3<Precision>>,
    pub orientation: Vec<UnitQuaternion<Precision>>,

    pub last_position: Vec<Point3<Precision>>,
    pub last_orientation: Vec<UnitQuaternion<Precision>>,

    pub linear_velocity: Vec<Vector3<Precision>>,
    pub angular_velocity: Vec<Vector3<Precision>>,

    pub force: Vec<Vector3<Precision>>,
    pub torque: Vec<Vector3<Precision>>,

    pub inverse_mass: Vec<Precision>,
    pub inverse_inertia_tensor_local: Vec<Matrix3<Precision>>,

    // derived data
    pub inverse_inertia_tensor_world: Vec<Matrix3<Precision>>,
}

impl BodySet {
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn has_finite_mass(&self, i: usize) -> bool {
        self.inverse_mass[i] > 0.0
    }

    pub fn apply_rotation_delta(&mut self, i: usize, rotation: Vector3<Precision>) {
        let q = self.orientation[i];

        self.orientation[i] = UnitQuaternion::from_scaled_axis(rotation) * q;

        self.orientation[i].renormalize();
    }

    pub fn update_derived_data(&mut self, i: usize) {
        let rot= self.orientation[i].to_rotation_matrix();
        
        self.inverse_inertia_tensor_world[i] = rot * self.inverse_inertia_tensor_local[i] * rot.transpose();
    }
}