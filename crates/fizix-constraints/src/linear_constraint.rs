use fizix_core::{BodyHandle, BodySet, Constraint, CorrectionData, Precision, EPSILON_SQUARED};
use nalgebra::{Point3, UnitVector3, Vector3};

pub struct LinearConstraint {
    body_a: BodyHandle,
    body_b: BodyHandle,

    local_point_a: Point3<Precision>,
    local_point_b: Point3<Precision>,

    local_axis: UnitVector3<Precision>, // relative to body A

    // distance along the axis should be clamped to
    min_distance: Precision,
    max_distance: Precision,

    compliance: Precision // inverse stiffness
}

impl Default for LinearConstraint {
    fn default() -> Self {
        Self {
            body_a: BodyHandle::INVALID,
            body_b: BodyHandle::INVALID,

            local_point_a: Point3::origin(),
            local_point_b: Point3::origin(),

            local_axis: Vector3::x_axis(),

            min_distance: Precision::NEG_INFINITY,
            max_distance: Precision::INFINITY,

            compliance: 0.0
        }
    }
}

impl Constraint for LinearConstraint {
    fn compute_correction(&self, bodies: &BodySet) -> Option<CorrectionData> {
        let body_a: usize = *self.body_a;
        let body_b: usize = *self.body_b;

        let pos_a = bodies.position[body_a];
        let pos_b = bodies.position[body_b];

        let orient_a = bodies.orientation[body_a];
        let orient_b = bodies.orientation[body_b];

        // relative space points
        let r_a = orient_a.transform_point(&self.local_point_a);
        let r_b = orient_b.transform_point(&self.local_point_b);

        let axis = orient_a.transform_vector(&self.local_axis);

        // absolute space points
        let p_a = pos_a + r_a.coords;
        let p_b = pos_b + r_b.coords;

        let difference = p_b - p_a;
        let d_dot = difference.dot(&axis).clamp(self.min_distance, self.max_distance);

        // project difference onto axis and subtract to get shortest distance
        let p_prime = difference - axis * d_dot;
        let distance_sq = p_prime.norm_squared();

        if distance_sq < EPSILON_SQUARED { return None; }   

        let distance = distance_sq.sqrt();
        let normal = p_prime / distance;

        Some(CorrectionData::Translational {
            handles: vec![self.body_a, self.body_b],
            relative_points: vec![r_a, r_b],
            normals: vec![normal, -normal],
            
            error: distance,
            alpha: self.compliance
        })
    }
}