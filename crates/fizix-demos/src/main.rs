use std::time::Instant;
use fizix_core::{Precision, World};
use kiss3d::light::Light;
use kiss3d::text::Font;
use kiss3d::window::{Window};
use nalgebra::{Matrix3, Vector3};

const ORANGE: (f32, f32, f32) = (244.0 / 255.0, 115.9 / 255.0, 51.0 / 255.0); // primary color
const LIGHT_GRAY: (f32, f32, f32) = (108.0 / 255.0, 112.0 / 255.0, 134.0 / 255.0); // secondary color
const DARK_GRAY: (f32, f32, f32) = (49.0 / 255.0, 50.0 / 255.0, 68.0 / 255.0); // static color
const WHITE: (f32, f32, f32) = (204.0 / 255.0, 214.0 / 255.0, 244.0 / 255.0); // text color

fn main() {
    let mut window = Window::new_with_size("Fizix", 1280, 720);
    let mut world = World::new(Vector3::new(0.0, -9.81 * 2.0, 0.0), 8, 2);

    window.set_light(Light::StickToCamera);
    window.set_background_color(24.0  / 255.0, 24.0 / 255.0, 37.0 / 255.0);

    let mut nodes = vec![] as Vec<kiss3d::scene::SceneNode>;

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