//! Functions used to calculate the volume of the unit cube under the plane P
//! with normal `g` at height `h` above the lowest vertex (see README).
//!
//! All functions return a volume in [0, 1], increasing with `h` from 0 at
//! `h = 0` to 1 at `h = max_height(g)`.

use glam::{Vec2, Vec3, vec2};

use crate::util::{get_plane_origin, rand_from_unit_cube, under_plane};

const MONTE_CARLO_ITERATIONS: i32 = 10_000;
const ZERO_ONE: [i32; 2] = [0, 1];
pub(crate) const EPS: f32 = 1e-6;

/// Brute-force reference: estimates the volume by sampling random points in the
/// cube. Standard error is ~0.005 with the current iteration count.
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

/// Same as [`monte_carlo`] but with `g` replaced by `|g|`. Used to check that
/// mirroring `g` into the positive octant doesn't change the volume.
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

/// Exact volume for any `g`. Mirrors `g` into the positive octant, then picks
/// the 3D, 2D or 1D formula depending on how many components are non-zero
/// (a zero component means the problem is constant along that axis).
pub fn exact(g: Vec3, h: f32) -> f32 {
    let g = g.abs();

    let nx = (g.x > EPS) as i32;
    let ny = (g.y > EPS) as i32;
    let nz = (g.z > EPS) as i32;

    match nx + ny + nz {
        3 => inclusion_exclusion_3d(g, h),
        2 => exact_2d(g, h),
        1 => exact_1d(h),
        _ => 0.0,
    }
}

/// Exact volume via inclusion-exclusion over the 8 cube vertices:
/// sum of (-1)^(i+j+k) * max(d - g·(i,j,k), 0)^3 / (6 gx gy gz).
/// Works for signed `g`, but all components must be non-zero. Loses precision
/// in f32 when a component is small.
pub fn inclusion_exclusion_3d(g: Vec3, h: f32) -> f32 {
    let o = get_plane_origin(g, h);
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

/// Same as [`inclusion_exclusion_3d`] with the vertex loop written out term by
/// term.
pub fn inclusion_exclusion_3d_unrolled(g: Vec3, h: f32) -> f32 {
    let o = get_plane_origin(g, h);
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

/// Exact volume as a piecewise polynomial in `h`. Sorts |g| so x <= y <= z,
/// finds which vertices of the cube lie under the plane, and evaluates the
/// expanded polynomial for that interval. Avoids some of the cancellation of
/// [`inclusion_exclusion_3d`]. All components must be non-zero.
pub fn piecewise_3d(g: Vec3, h: f32) -> f32 {
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

    // g is non-negative, so the lowest vertex is the origin and d = h
    let d = h;

    // Vertex heights along g, sorted. With x <= y <= z the only pair whose order
    // isn't fixed is z vs x + y: the plane reaches vertex z before x + y only
    // when x + y > z (g is closer to the diagonal than to the z axis).
    let (intervals, z_before_xy) = if x + y > z {
        ([x, y, z, x + y, x + z, y + z, x + y + z], true)
    } else {
        ([x, y, x + y, z, x + z, y + z, x + y + z], false)
    };

    let i = find_interval(d, &intervals, z_before_xy);
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

/// Exact volume as a piecewise polynomial in `h`, like [`piecewise_3d`], but
/// using the symmetry V(h) = 1 - V(s - h) with s = x + y + z = `max_height(g)`.
///
/// The vertex heights are symmetric around s / 2 (x pairs with y + z, y with
/// x + z, z with x + y), so after reflecting `h` into the lower half only the
/// first 4 intervals are possible. This avoids the upper polynomials, which
/// lose precision for small components, and the `Z` polynomial is expanded
/// around s / 2 for the same reason: the error stays below ~3e-7 for any `g`
/// with non-zero components. All components must be non-zero.
pub fn piecewise_3d_symmetric(g: Vec3, h: f32) -> f32 {
    let g = g.abs();

    // sort so that x <= y <= z
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

    let s = x + y + z;
    let flip = h > 0.5 * s;
    let d = if flip { s - h } else { h };

    // 6xyz times the volume under the plane at height d <= s / 2
    let sum = if d < x {
        // only the lowest vertex is under the plane
        d * d * d
    } else if d < y {
        // vertex x is under the plane
        x * ((x - 3.0 * d) * x + 3.0 * d * d)
    } else if d < z && d < x + y {
        // vertices x and y are under the plane
        let t = x + y;
        (t - d).powi(3) + 3.0 * x * y * (2.0 * d - t)
    } else if x + y > z {
        // vertex z is under the plane (it comes before x + y). Expanded around
        // s / 2, where the volume is 1/2, to avoid cancellation: with
        // u = d - s / 2 this is 3xyz + u(c - 2u^2), where
        // c = 3/2 (2(xy + yz + zx) - (x^2 + y^2 + z^2)), written so that it
        // has no cancellation either.
        let u = d - 0.5 * s;
        let c = 1.5 * (x * (2.0 * (y + z) - x) - (z - y) * (z - y));
        3.0 * x * y * z + u * (c - 2.0 * u * u)
    } else {
        // vertex x + y is under the plane (it comes before z)
        3.0 * x * y * (2.0 * d - x - y)
    };

    let v = sum / (6.0 * x * y * z);
    if flip { 1.0 - v } else { v }
}

/// Range of `h` between two consecutive vertex heights, named after the highest
/// vertex under the plane (`X` = vertex at height x, `Xy` = at height x + y, ...).
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

/// Finds the interval of the sorted vertex heights `p` that contains `s`.
fn find_interval(s: f32, p: &[f32; 7], z_before_xy: bool) -> Interval {
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
        if z_before_xy {
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

/// Exact volume when exactly one component of `g` is zero: the cube is a square
/// extruded along that axis, so the volume equals the area under the line.
fn exact_2d(g: Vec3, h: f32) -> f32 {
    let r = if g.x < EPS {
        vec2(g.y, g.z)
    } else if g.y < EPS {
        vec2(g.x, g.z)
    } else {
        vec2(g.x, g.y)
    };

    inclusion_exclusion_2d(r, h)
}

/// Area of the unit square under the line with normal `g` at height `h`, via
/// inclusion-exclusion over the 4 square vertices. Both components must be
/// non-zero.
pub fn inclusion_exclusion_2d(g: Vec2, h: f32) -> f32 {
    let v = 0.5 - g.signum() * 0.5;
    let o = v + g * h;
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

/// Exact volume when `g` is aligned with an axis: the water fills a slab of
/// height `h`.
fn exact_1d(h: f32) -> f32 {
    h
}

/// Derivative of [`exact`] with respect to `h`. Because `g` is normalized, this
/// is the area of the water surface (the plane P clipped to the cube). Uses the
/// same dispatch as [`exact`], with the inclusion-exclusion sums differentiated
/// term by term, so it has the same f32 precision limits.
pub fn surface_area(g: Vec3, h: f32) -> f32 {
    let g = g.abs();

    let nx = (g.x > EPS) as i32;
    let ny = (g.y > EPS) as i32;
    let nz = (g.z > EPS) as i32;

    match nx + ny + nz {
        3 => {
            let mut sum = 0.0;
            for i in ZERO_ONE {
                for j in ZERO_ONE {
                    for k in ZERO_ONE {
                        let f = h - i as f32 * g.x - j as f32 * g.y - k as f32 * g.z;
                        let f = f.max(0.0);
                        sum += (1 - (i + j + k) % 2 * 2) as f32 * f * f;
                    }
                }
            }
            sum / (2.0 * g.x * g.y * g.z)
        }
        2 => {
            // length of the water line in the square, extruded by 1
            let r = if g.x < EPS {
                vec2(g.y, g.z)
            } else if g.y < EPS {
                vec2(g.x, g.z)
            } else {
                vec2(g.x, g.y)
            };

            let mut sum = 0.0;
            for i in ZERO_ONE {
                for j in ZERO_ONE {
                    let f = h - i as f32 * r.x - j as f32 * r.y;
                    sum += (1 - (i + j) % 2 * 2) as f32 * f.max(0.0);
                }
            }
            sum / (r.x * r.y)
        }
        1 => 1.0,
        _ => 0.0,
    }
}

/// Inclusion-exclusion in f64 on the same f32 inputs, accurate enough to serve
/// as a reference in tests even for small components. All components of `g`
/// must be non-zero.
#[cfg(test)]
pub(crate) fn reference(g: Vec3, h: f32) -> f64 {
    let g = glam::DVec3::new(g.x as f64, g.y as f64, g.z as f64).abs();
    let h = h as f64;
    let mut sum = 0.0;
    for i in ZERO_ONE {
        for j in ZERO_ONE {
            for k in ZERO_ONE {
                let f = h - i as f64 * g.x - j as f64 * g.y - k as f64 * g.z;
                let sign = if (i + j + k) % 2 == 0 { 1.0 } else { -1.0 };
                sum += sign * f.max(0.0).powi(3);
            }
        }
    }
    sum / (6.0 * g.x * g.y * g.z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::{max_height, rand_normalized_gaussian};
    use glam::vec3;
    use rand::RngExt;

    #[test]
    fn piecewise_3d_symmetric_matches_exact() {
        let mut rng = rand::rng();
        for _ in 0..10_000 {
            let g = rand_normalized_gaussian(&mut rng);
            // `exact` is only accurate when no component is small
            if g.abs().min_element() < 0.1 {
                continue;
            }
            let h = rng.random::<f32>() * max_height(g);
            let (a, b) = (piecewise_3d_symmetric(g, h), exact(g, h));
            assert!((a - b).abs() < 1e-4, "g = {g}, h = {h}: {a} vs {b}");
        }
    }

    #[test]
    fn piecewise_3d_symmetric_is_accurate_for_small_components() {
        let mut rng = rand::rng();
        for _ in 0..10_000 {
            let mut g = rand_normalized_gaussian(&mut rng);
            // smallest component log-uniform in [1e-6, 1e-2)
            let i = g.abs().min_position();
            g[i] = 10f32.powf(-6.0 + 4.0 * rng.random::<f32>());
            let g = g.normalize();
            let h = rng.random::<f32>() * max_height(g);

            let v = piecewise_3d_symmetric(g, h);
            let r = reference(g, h);
            assert!((v as f64 - r).abs() < 1e-6, "g = {g}, h = {h}: {v} vs {r}");
        }
    }

    #[test]
    fn piecewise_3d_symmetric_endpoints() {
        let g = vec3(0.2, -0.5, 0.9).normalize();
        assert_eq!(piecewise_3d_symmetric(g, 0.0), 0.0);
        assert!((piecewise_3d_symmetric(g, max_height(g)) - 1.0).abs() < 1e-6);
        // the plane through the cube's center halves it
        assert!((piecewise_3d_symmetric(g, 0.5 * max_height(g)) - 0.5).abs() < 1e-6);
    }
}
