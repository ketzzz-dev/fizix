use fizix_core::{BodyHandle, BodySet, Constraint, CorrectionData, Precision, EPSILON, EPSILON_SQUARED};
use nalgebra::Point3;

pub struct DistanceConstraint {
    pub body_a: BodyHandle,
    pub body_b: BodyHandle,

    pub local_point_a: Point3<Precision>,
    pub local_point_b: Point3<Precision>,

    pub rest_length: Precision,

    pub compliance: Precision, // inverse stiffness
}

impl Default for DistanceConstraint {
    fn default() -> Self {
        Self {
            body_a: BodyHandle::INVALID,
            body_b: BodyHandle::INVALID,

            local_point_a: Point3::origin(),
            local_point_b: Point3::origin(),

            rest_length: 0.0,

            compliance: 0.0
        }
    }
}

impl Constraint for DistanceConstraint {
    fn compute_correction(&self, bodies: &BodySet) -> Option<CorrectionData> {
        let body_a = *self.body_a;
        let body_b = *self.body_b;

        let pos_a = bodies.position[body_a];
        let pos_b = bodies.position[body_b];

        let orient_a = bodies.orientation[body_a];
        let orient_b = bodies.orientation[body_b];

        // relative space points
        let r_a = orient_a.transform_point(&self.local_point_a);
        let r_b = orient_b.transform_point(&self.local_point_b);

        // absolute space points
        let p_a = pos_a + r_a.coords;
        let p_b = pos_b + r_b.coords;

        let difference = p_b - p_a;
        let distance_sq = difference.norm_squared();

        if distance_sq < EPSILON_SQUARED { return None; }

        let distance = distance_sq.sqrt();
        let error = distance - self.rest_length;

        if error.abs() < EPSILON { return None; }

        let normal = difference / distance;

        Some(CorrectionData::Translational {
            handles: vec![self.body_a, self.body_b],
            relative_points: vec![r_a, r_b],
            normals: vec![normal, -normal],
            
            error, alpha: self.compliance
        })
    }
}