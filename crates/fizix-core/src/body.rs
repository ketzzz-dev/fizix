use crate::Precision;
use nalgebra::{Isometry3, Matrix3, Point3, UnitQuaternion, Vector3};

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct BodyHandle(pub usize);


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
    pub transform: Vec<Isometry3<Precision>>,
}

impl BodySet {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn has_finite_mass(&self, i: usize) -> bool {
        self.inverse_mass[i] > 0.0
    }

    pub fn apply_rotation_delta(&mut self, i: usize, rotation: Vector3<Precision>) {
        let current_orientation = self.orientation[i];

        self.orientation[i] = UnitQuaternion::from_scaled_axis(rotation) * current_orientation;

        self.orientation[i].renormalize();
    }

    pub fn update_derived_data(&mut self, i: usize) {
        self.transform[i] = Isometry3::from_parts(self.position[i].into(), self.orientation[i]);
        
        let rotation_matrix = self.transform[i].rotation.to_rotation_matrix();

        self.inverse_inertia_tensor_world[i] = rotation_matrix * self.inverse_inertia_tensor_local[i] * rotation_matrix.transpose();
    }
}