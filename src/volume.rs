use glam::{Vec2, Vec3, vec2};

use crate::util::{get_plane_origin, rand_from_unit_cube, under_plane};

const MONTE_CARLO_ITERATIONS: i32 = 10_000;
const ZERO_ONE: [i32; 2] = [0, 1];
const EPS: f32 = 1e-6;

pub fn monte_carlo(g: Vec3, h: f32) -> f32 {
    let mut rng = rand::rng();
    let mut under: i32 = 0;
    let o = get_plane_origin(g, h);

    for _ in 0..MONTE_CARLO_ITERATIONS {
        let p = rand_from_unit_cube(&mut rng);

        if under_plane(o, g, p) {
            under += 1;
        }
    }

    under as f32 / MONTE_CARLO_ITERATIONS as f32
}

pub fn monte_carlo_abs(g: Vec3, h: f32) -> f32 {
    let g = g.abs();
    let mut rng = rand::rng();
    let mut under: i32 = 0;
    let o = get_plane_origin(g, h);

    for _ in 0..MONTE_CARLO_ITERATIONS {
        let p = rand_from_unit_cube(&mut rng);

        if under_plane(o, g, p) {
            under += 1;
        }
    }

    under as f32 / MONTE_CARLO_ITERATIONS as f32
}

pub fn formula(g: Vec3, h: f32) -> f32 {
    let g = g.abs();

    let nx = (g.x > EPS) as i32;
    let ny = (g.y > EPS) as i32;
    let nz = (g.z > EPS) as i32;

    match nx + ny + nz {
        3 => formula3(g, h),
        2 => formula2(g, h),
        1 => formula1(h),
        _ => 0.0,
    }
}

pub fn formula3(g: Vec3, h: f32) -> f32 {
    let v = g.signum() * 0.5 + 0.5;
    let o = v - g * h;
    let d = g.dot(o);

    let mut sum = 0.0;
    for i in ZERO_ONE {
        for j in ZERO_ONE {
            for k in ZERO_ONE {
                let f = d - i as f32 * g.x - j as f32 * g.y - k as f32 * g.z;
                sum += (1 - (i + j + k) % 2 * 2) as f32 * (f * f * f).max(0.0);
            }
        }
    }

    sum / (6.0 * g.x * g.y * g.z)
}

fn formula2(g: Vec3, h: f32) -> f32 {
    let r = if g.x < EPS {
        vec2(g.y, g.z)
    }
    else if g.y < EPS {
        vec2(g.x, g.z)
    }
    else {
        vec2(g.x, g.y)
    };

    formula_square(r, h)
}

pub fn formula_square(g: Vec2, h: f32) -> f32 {
    let v = g.signum() * 0.5 + 0.5;
    let o = v - g * h;
    let d = g.dot(o);

    let mut sum = 0.0;
    for i in ZERO_ONE {
        for j in ZERO_ONE {
            let f = d - i as f32 * g.x - j as f32 * g.y;
            sum += (1 - (i + j) % 2 * 2) as f32 * (f.max(0.0) * f.max(0.0));
        }
    }

    sum / (2.0 * g.x * g.y)
}

fn formula1(h: f32) -> f32 {
    1.0 - h
}
