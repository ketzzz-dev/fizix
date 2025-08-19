use std::f64::consts::{FRAC_PI_2};
use std::time::Instant;
use fizix_constraints::{RevoluteJoint, SphericalJoint};
use fizix_core::{Precision, World};
use kiss3d::light::Light;
use kiss3d::text::Font;
use kiss3d::window::{Window};
use nalgebra::{Matrix3, Point3, UnitQuaternion, UnitVector3, Vector3};

const ORANGE: (f32, f32, f32) = (244.0 / 255.0, 115.9 / 255.0, 51.0 / 255.0); // primary color
const LIGHT_GRAY: (f32, f32, f32) = (108.0 / 255.0, 112.0 / 255.0, 134.0 / 255.0); // secondary color
const DARK_GRAY: (f32, f32, f32) = (49.0 / 255.0, 50.0 / 255.0, 68.0 / 255.0); // static color
const WHITE: (f32, f32, f32) = (204.0 / 255.0, 214.0 / 255.0, 244.0 / 255.0); // text color

fn main() {
    let mut window = Window::new_with_size("Fizix", 1280, 720);
    let mut world = World::new(Vector3::new(0.0, -9.81 * 2.0, 0.0), 8, 2);

    window.set_light(Light::StickToCamera);
    window.set_background_color(24.0  / 255.0, 24.0 / 255.0, 37.0 / 255.0);

    let sphere = world.add_body(
        Point3::new(0.0, 2.5, 20.0), 
        UnitQuaternion::identity(),
        0.0,
        Matrix3::zeros()
    );
    let arm1 = world.add_body(
        Point3::new(0.0, 2.5, 22.5),
        UnitQuaternion::from_axis_angle(&Vector3::x_axis(), -FRAC_PI_2),
        1.0,
        cylinder_inertia_tensor(0.25, 5.0, 1.0)
    );
    let arm2 = world.add_body(
        Point3::new(0.0, 2.5, 25.0),
        UnitQuaternion::from_axis_angle(&Vector3::x_axis(), -FRAC_PI_2),
        1.0,
        cylinder_inertia_tensor(0.25, 5.0, 1.0)
    );

    let joint1 = RevoluteJoint::new(
        sphere, 
        arm1, 
        Point3::origin(),
        Point3::new(0.0, 2.5, 0.0),
        UnitVector3::new_normalize(Vector3::new(1.0, 1.0, 0.0)),
        Vector3::x_axis(),
        0.0
    );
    let joint2 = SphericalJoint::new(
        arm1, 
        arm2, 
        Point3::new(0.0, -2.5, 0.0),
        Point3::new(0.0, 2.5, 0.0),
        0.0
    );

    world.add_constraint(joint1);
    world.add_constraint(joint2);

    let mut sphere_node = window.add_sphere(0.5);
    let mut arm1_node = window.add_capsule(0.25, 5.0);
    let mut arm2_node = window.add_capsule(0.25, 5.0);

    sphere_node.set_color(DARK_GRAY.0, DARK_GRAY.1, DARK_GRAY.2);
    arm1_node.set_color(ORANGE.0, ORANGE.1, ORANGE.2);
    arm2_node.set_color(ORANGE.0, ORANGE.1, ORANGE.2);

    let mut nodes = vec![sphere_node, arm1_node, arm2_node];

    let mut last_time = Instant::now();

    let mut fps_samples = [1.0 / 75.0; 75];
    let mut sample_index = 0;

    while window.render() {
        let elapsed = last_time.elapsed().as_secs_f64();

        last_time = Instant::now();
        fps_samples[sample_index] = elapsed;
        sample_index = (sample_index + 1) % fps_samples.len();

        world.step(elapsed);

        let fps = fps_samples.len() as f64 / fps_samples.iter().sum::<f64>();

        window.draw_text(
            &format!("Body Count: {}\nConstraint Count: {}\nFPS: {:.0}", world.bodies.position.len(), world.constraints.len(), fps),
            &kiss3d::nalgebra::Point2::new(10.0, 10.0),
            42.0,
            &Font::default(),
            &kiss3d::nalgebra::Point3::new(WHITE.0, WHITE.1, WHITE.2),
        );
        
        for (i, node) in nodes.iter_mut().enumerate() {
            let position = world.bodies.position[i];
            let orientation = world.bodies.orientation[i];
            let orientation = kiss3d::nalgebra::Quaternion::new(
                orientation.w as f32,
                orientation.i as f32,
                orientation.j as f32,
                orientation.k as f32
            );

            node.set_local_translation(kiss3d::nalgebra::Translation3::new(position.x as f32, position.y as f32, position.z as f32));
            node.set_local_rotation(kiss3d::nalgebra::UnitQuaternion::from_quaternion(orientation));
        }
    }
}

pub fn cuboid_inertia_tensor(width: Precision, height: Precision, length: Precision, mass: Precision) -> Matrix3<Precision> {
    let x2 = width * width;
    let y2 = height * height;
    let z2 = length * length;

    let i = mass / 12.0;

    Matrix3::new(
        i * (y2 + z2), 0.0, 0.0,
        0.0, i * (x2 + z2), 0.0,
        0.0, 0.0, i * (x2 + y2)
    )
}

pub fn sphere_inertia_tensor(radius: Precision, mass: Precision) -> Matrix3<Precision> {
    let i = 2.0 * mass * radius * radius / 5.0;

    Matrix3::new(
        i, 0.0, 0.0,
        0.0, i, 0.0,
        0.0, 0.0, i
    )
}

pub fn cylinder_inertia_tensor(radius: Precision, height: Precision, mass: Precision) -> Matrix3<Precision> {
    let i_xz = mass * (3.0 * radius * radius + height * height) / 12.0;

    Matrix3::new(
        i_xz, 0.0, 0.0,
        0.0, 0.5 * mass * radius * radius, 0.0,
        0.0, 0.0, i_xz
    )
}