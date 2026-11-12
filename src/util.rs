use std::f32::consts::PI;

use glam::{Vec3, vec3};
use rand::{RngExt, rngs::ThreadRng};
use rand_distr::{Distribution, Normal};

pub fn rand_from_unit_cube(rng: &mut ThreadRng) -> Vec3 {
    let x = rng.random();
    let y = rng.random();
    let z = rng.random();

    vec3(x, y, z)
}

fn plane_sdf(o: Vec3, n: Vec3, p: Vec3) -> f32 {
    (p - o).dot(n)
}

pub fn under_plane(o: Vec3, n: Vec3, p: Vec3) -> bool {
    plane_sdf(o, n, p) < 0.0
}

pub fn rand_normalized(rng: &mut ThreadRng) -> Vec3 {
    let x = rng.random::<f32>() * 2.0 - 1.0;
    let y = rng.random::<f32>() * 2.0 - 1.0;
    let z = rng.random::<f32>() * 2.0 - 1.0;

    vec3(x, y, z).normalize()
}

pub fn rand_normalized_gaussian(rng: &mut ThreadRng) -> Vec3 {
    let normal = Normal::new(0.0, 1.0).unwrap();

    let x = normal.sample(rng);
    let y = normal.sample(rng);
    let z = normal.sample(rng);

    vec3(x, y, z).normalize()
}

pub fn rand_normalized_sphere(rng: &mut ThreadRng) -> Vec3 {
    let z = rng.random::<f32>() * 2.0 - 1.0;
    let phi = rng.random::<f32>() * 2.0 * PI;

    let r = (1.0 - z * z).sqrt();

    let x = r * phi.cos();
    let y = r * phi.sin();

    vec3(x, y, z)
}

// Point on the plane P at height h above the lowest vertex of the cube (along g).
pub fn get_plane_origin(g: Vec3, h: f32) -> Vec3 {
    let v = 0.5 - g.signum() * 0.5;
    v + g * h
}

pub fn max_height(g: Vec3) -> f32 {
    let v = g.signum() * 0.5 + 0.5;

    let diag = v - vec3(1.0 - v.x, 1.0 - v.y, 1.0 - v.z);
    diag.dot(g)
}
