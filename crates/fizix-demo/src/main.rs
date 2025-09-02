use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};
use std::time::Instant;
use fizix_constraints::{AngularConstraint, AxisConstraint, DistanceConstraint};
use fizix_core::{Precision, World};
use kiss3d::light::Light;
use kiss3d::text::Font;
use kiss3d::window::{Window};
use nalgebra::{Matrix3, Point3, UnitQuaternion, Vector3};

const ORANGE: (f32, f32, f32) = (244.0 / 255.0, 115.9 / 255.0, 51.0 / 255.0); // primary color
const LIGHT_GRAY: (f32, f32, f32) = (108.0 / 255.0, 112.0 / 255.0, 134.0 / 255.0); // secondary color
const DARK_GRAY: (f32, f32, f32) = (49.0 / 255.0, 50.0 / 255.0, 68.0 / 255.0); // static color
const WHITE: (f32, f32, f32) = (204.0 / 255.0, 214.0 / 255.0, 244.0 / 255.0); // text color

fn main() {
    let mut window = Window::new_with_size("Fizix", 1280, 720);
    let mut world = World::new(Vector3::new(0.0, -9.81 * 2.0, 0.0), 8, 2);

    window.set_light(Light::StickToCamera);
    window.set_background_color(24.0  / 255.0, 24.0 / 255.0, 37.0 / 255.0);

    let mut anchor_node = window.add_cylinder(0.25, 0.25);
    let mut arm_1_node = window.add_cube(0.5, 0.25, 2.5);
    let mut arm_2_node = window.add_cube(0.5, 0.25, 2.5);

    anchor_node.set_color(DARK_GRAY.0, DARK_GRAY.1, DARK_GRAY.2);
    arm_1_node.set_color(LIGHT_GRAY.0, LIGHT_GRAY.1, LIGHT_GRAY.2);
    arm_2_node.set_color(LIGHT_GRAY.0, LIGHT_GRAY.1, LIGHT_GRAY.2);

    let mut nodes = vec![anchor_node, arm_1_node, arm_2_node] as Vec<kiss3d::scene::SceneNode>;

    let anchor = world.add_body(
        Point3::new(0.0, 2.5, 10.0),
        UnitQuaternion::from_euler_angles(-FRAC_PI_2, 0.0, 0.0),
        0.0,
        cylinder_inertia_tensor(0.25, 0.25, 0.0)
    );
    let arm_1 = world.add_body(
        Point3::new(1.25, 2.5, 9.75),
        UnitQuaternion::from_euler_angles(-FRAC_PI_2, 0.0, FRAC_PI_2),
        1.0,
        cuboid_inertia_tensor(1.0, 0.125, 2.5, 1.0)
    );
    let arm_2 = world.add_body(
        Point3::new(3.75, 2.5, 10.0),
        UnitQuaternion::from_euler_angles(-FRAC_PI_2, 0.0, FRAC_PI_2),
        1.0,
        cuboid_inertia_tensor(1.0, 0.125, 2.5, 1.0)
    );

    world.add_constraint(DistanceConstraint {
        body_a: anchor,
        body_b: arm_1,

        local_point_a: Point3::new(0.0, 0.125, 0.0),
        local_point_b: Point3::new(0.0, -0.125, 1.25),

        ..Default::default()
    });
    world.add_constraint(AxisConstraint {
        body_a: anchor,
        body_b: arm_1,

        local_axis_a: Vector3::y_axis(),
        local_axis_b: Vector3::y_axis(),

        ..Default::default()
    });

    world.add_constraint(DistanceConstraint {
        body_a: arm_1,
        body_b: arm_2,

        local_point_a: Point3::new(0.0, -0.125, -1.25),
        local_point_b: Point3::new(0.0, 0.125, 1.25),

        ..Default::default()
    });
    world.add_constraint(AxisConstraint {
        body_a: arm_1,
        body_b: arm_2,

        local_axis_a: Vector3::y_axis(),
        local_axis_b: Vector3::y_axis(),

        ..Default::default()
    });

    let mut last_time = Instant::now();

    let mut fps_samples = [1.0 / 75.0; 75];
    let mut sample_index = 0;

    while window.render() {
        let elapsed = last_time.elapsed().as_secs_f64();

        last_time = Instant::now();
        fps_samples[sample_index] = elapsed;
        sample_index = (sample_index + 1) % fps_samples.len();

        world.step(elapsed);

        let fps = fps_samples.len() as Precision / fps_samples.iter().sum::<Precision>();

        window.draw_text(
            &format!("Body Count: {}\nConstraint Count: {}\nFPS: {:.0}", world.bodies.position.len(), world.constraints.len(), fps),
            &kiss3d::nalgebra::Point2::new(10.0, 10.0),
            42.0,
            &Font::default(),
            &kiss3d::nalgebra::Point3::new(WHITE.0, WHITE.1, WHITE.2),
        );
        
        for (i, node) in nodes.iter_mut().enumerate() {
            if i >= world.bodies.position.len() { break; }

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