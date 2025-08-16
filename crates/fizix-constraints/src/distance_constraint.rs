use fizix_core::{BodyHandle, BodySet, Constraint, CorrectionData, Precision, EPSILON_SQUARED};
use nalgebra::Point3;

pub struct DistanceConstraint {
    body_a: BodyHandle, body_b: BodyHandle,

    local_point_a: Point3<Precision>,
    local_point_b: Point3<Precision>,

    distance: Precision
}

impl Constraint for DistanceConstraint {
    fn compute_correction_data(&self, bodies: &BodySet) -> Option<CorrectionData> {
        let body_a = self.body_a.0;
        let body_b = self.body_b.0;

        let world_point_a = bodies.transform[body_a].transform_point(&self.local_point_a);
        let world_point_b = bodies.transform[body_b].transform_point(&self.local_point_b);

        let difference = world_point_a - world_point_b;
        let distance_squared = difference.norm_squared();

        if distance_squared < EPSILON_SQUARED { return None; }

        let distance = distance_squared.sqrt();
        let normal = difference / distance;

        let relative_point_a = world_point_a - bodies.position[body_a];
        let relative_point_b = world_point_b - bodies.position[body_b];

        Some(CorrectionData::Translational {
            error: distance - self.distance,
            body_handles: vec![self.body_a, self.body_b],
            relative_points: vec![relative_point_a.into(), relative_point_b.into()], 
            normals: vec![normal, -normal]
        })
    }
}