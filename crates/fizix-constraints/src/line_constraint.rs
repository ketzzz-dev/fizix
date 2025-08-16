use fizix_core::{BodyHandle, BodySet, Constraint, CorrectionData, Precision, EPSILON_SQUARED};
use nalgebra::{Point3, UnitVector3};

pub struct LineConstraint {
    body_a: BodyHandle, body_b: BodyHandle,

    local_point_a: Point3<Precision>,
    local_point_b: Point3<Precision>,

    local_axis: UnitVector3<Precision> // relative to B
}
impl LineConstraint {
    pub fn new(
        body_a: BodyHandle, body_b: BodyHandle,

        local_point_a: Point3<Precision>,
        local_point_b: Point3<Precision>,

        local_axis: UnitVector3<Precision>
    ) -> Self {
        Self {
            body_a, body_b,

            local_point_a,
            local_point_b,

            local_axis
        }
    }
}

impl Constraint for LineConstraint {
    fn compute_correction_data(&self, bodies: &BodySet) -> Option<CorrectionData> {
        let body_a: usize = self.body_a.0;
        let body_b = self.body_b.0;

        let world_point_a = bodies.transform[body_a].transform_point(&self.local_point_a);
        let world_point_b = bodies.transform[body_b].transform_point(&self.local_point_b);
        let world_direction = bodies.transform[body_b].transform_vector(&self.local_axis);

        let projected_point = world_point_b + world_direction
            * world_direction.dot(&(world_point_a - world_point_b));

        let difference = world_point_a - projected_point;
        let distance_squared = difference.norm_squared();
        
        if distance_squared < EPSILON_SQUARED { return None; }

        let distance = distance_squared.sqrt();
        let normal = difference / distance;

        let relative_point_a = world_point_a - bodies.position[body_a];
        let relative_point_b = projected_point - bodies.position[body_b];

        Some(CorrectionData::Translational {
            error: distance,
            body_handles: vec![self.body_a, self.body_b],
            relative_points: vec![relative_point_a.into(), relative_point_b.into()], 
            normals: vec![normal, -normal]
        })
    }
}