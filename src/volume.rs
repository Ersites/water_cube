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

pub fn formula3ex(g: Vec3, h: f32) -> f32 {
    let v = g.signum() * 0.5 + 0.5;
    let o = v - g * h;
    let d = g.dot(o);

    let mut sum = 0.0;

    sum += (d).powi(3).max(0.0);

    sum -= (d - g.x).powi(3).max(0.0);
    sum -= (d - g.y).powi(3).max(0.0);
    sum -= (d - g.z).powi(3).max(0.0);

    sum += (d - g.x - g.y).powi(3).max(0.0);
    sum += (d - g.x - g.z).powi(3).max(0.0);
    sum += (d - g.y - g.z).powi(3).max(0.0);

    sum -= (d - g.x - g.y - g.z).powi(3).max(0.0);

    sum / (6.0 * g.x * g.y * g.z)
}

pub fn formula3an(g: Vec3, h: f32) -> f32 {
    let g = g.abs();

    let (mut x, mut y, mut z) = (g.x, g.y, g.z);

    if x > y {
        std::mem::swap(&mut x, &mut y);
    }
    if y > z {
        std::mem::swap(&mut y, &mut z);
    }
    if x > y {
        std::mem::swap(&mut x, &mut y);
    }

    let g = Vec3::new(x, y, z);

    let o = 1.0 - g * h;
    let d = g.dot(o);

    let (intervals, natural) = if x + y > z {
        ([x, y, z, x + y, x + z, y + z, x + y + z], true)
    } else {
        ([x, y, x + y, z, x + z, y + z, x + y + z], false)
    };

    let i = find_interval(d, &intervals, natural);
    // dbg!(&i);

    let sum = match i {
        Interval::None => d * d * d,
        Interval::X => x * ((x - 3.0 * d) * x + 3.0 * d * d),
        Interval::Y => {
            let s = x + y;
            (s - d).powi(3) + 3.0 * x * y * (2.0 * d - s)
        }
        Interval::Z => {
            let s = x + y + z;
            let q = x * x + y * y + z * z;
            let r = x * x * x + y * y * y + z * z * z;

            ((-2.0 * d + 3.0 * s) * d - 3.0 * q) * d + r
        }
        Interval::Xy => 3.0 * x * y * (2.0 * d - x - y),
        Interval::Zxy => {
            let t = z - d;
            t * t * t + 3.0 * x * y * (2.0 * d - x - y)
        }
        Interval::Xz => {
            let s = g.y * g.z;
            let q = d - g.y - g.z;

            x * ((-x + 3.0 * q) * x + 3.0 * (2.0 * s - q * q))
        }
        Interval::Yz => {
            let s = g.x + g.y + g.z;
            let q = s * s;

            ((d - 3.0 * s) * d + 3.0 * q) * d - q * s + 6.0 * g.x * g.y * g.z
        }
        Interval::Xyz => 6.0 * x * y * z,
    };

    sum / (6.0 * x * y * z)
}

#[derive(Debug)]
enum Interval {
    None,
    X,
    Y,
    Z,
    Xy,
    Zxy,
    Xz,
    Yz,
    Xyz,
}

fn find_interval(s: f32, p: &[f32; 7], natural: bool) -> Interval {
    if s < p[0] {
        return Interval::None;
    }
    if s < p[1] {
        return Interval::X;
    }
    if s < p[2] {
        return Interval::Y;
    }
    if s < p[3] {
        if natural {
            return Interval::Z;
        } else {
            return Interval::Xy;
        }
    }
    if s < p[4] {
        return Interval::Zxy;
    }
    if s < p[5] {
        return Interval::Xz;
    }
    if s < p[6] {
        return Interval::Yz;
    }
    Interval::Xyz
}

fn formula2(g: Vec3, h: f32) -> f32 {
    let r = if g.x < EPS {
        vec2(g.y, g.z)
    } else if g.y < EPS {
        vec2(g.x, g.z)
    } else {
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
