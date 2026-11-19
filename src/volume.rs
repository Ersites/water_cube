//! Functions used to calculate the volume of the unit cube under the plane P
//! with normal `g` at height `h` above the lowest vertex (see README).
//!
//! All functions return a volume in [0, 1], increasing with `h` from 0 at
//! `h = 0` to 1 at `h = max_height(g)`.

use glam::Vec3;

use crate::util::{get_plane_origin, rand_from_unit_cube, under_plane};

const MONTE_CARLO_ITERATIONS: i32 = 10_000;
const ZERO_ONE: [i32; 2] = [0, 1];

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

/// Exact volume for any non-zero `g`, including components that are zero or
/// close to zero, with no threshold.
///
/// Uses the lower-half polynomials of [`piecewise_3d_symmetric`] (reflecting
/// `h` with V(h) = 1 - V(s - h)), but instead of dividing the result by 6xyz,
/// each interval divides only by components that its own range guarantees to
/// be non-zero, with the other factors cancelled by hand. So the 3D formulas
/// turn into the 2D and 1D ones as components go to zero, and the error stays
/// at f32 precision (~3e-7) for every `g`.
pub fn exact(g: Vec3, h: f32) -> f32 {
    let (x, y, z) = sorted_components(g);
    let s = x + y + z;

    // V(h) = 1 - V(s - h): evaluate the lower half only
    let flip = h > 0.5 * s;
    let d = if flip { s - h } else { h };

    let v = lower_volume(x, y, z, d);
    if flip { 1.0 - v } else { v }
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

/// Volume at height `d <= s / 2`, for sorted x <= y <= z, used by [`exact`].
pub(crate) fn lower_volume(x: f32, y: f32, z: f32, d: f32) -> f32 {
    let t = x + y;
    match LowerInterval::of_height(x, y, z, d) {
        // d < x, so x > 0, and d^3 / x stays small
        LowerInterval::None => d * d * d / (6.0 * x * y * z),
        // d < y, so y > 0; x(x^2 - 3dx + 3d^2) / (6xyz) with x cancelled
        LowerInterval::X => (x * x - 3.0 * d * x + 3.0 * d * d) / (6.0 * y * z),
        // y <= d < x + y, so x > 0, and u = x + y - d <= x keeps u^3 / x small;
        // (u^3 + 3xy(2d - t)) / (6xyz) with xy cancelled in the second term
        LowerInterval::Y => {
            let u = t - d;
            u * u * u / (6.0 * x * y * z) + (2.0 * d - t) / (2.0 * z)
        }
        // x + y > z >= y, so x > 0; expanded around s / 2 as in
        // `piecewise_3d_symmetric`
        LowerInterval::Z => {
            let u = d - 0.5 * (t + z);
            let c = 1.5 * (x * (2.0 * (y + z) - x) - (z - y) * (z - y));
            0.5 + u * (c - 2.0 * u * u) / (6.0 * x * y * z)
        }
        // 3xy(2d - t) / (6xyz) with xy cancelled; z is the largest component
        LowerInterval::Xy => (2.0 * d - t) / (2.0 * z),
    }
}

/// The absolute components of `g`, sorted so that x <= y <= z.
pub(crate) fn sorted_components(g: Vec3) -> (f32, f32, f32) {
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
    (x, y, z)
}

/// Interval of the lower half (d <= s / 2) of the piecewise volume, named
/// after the highest vertex under the plane, as in [`piecewise_3d`].
/// Above s / 2 the symmetry V(d) = 1 - V(s - d) maps back into these.
/// Intervals whose range is empty (for example `X` when x = y) are skipped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LowerInterval {
    /// d in [0, x): only the lowest vertex is under the plane
    None,
    /// d in [x, y)
    X,
    /// d in [y, min(z, x + y))
    Y,
    /// d in [z, s / 2], only when x + y > z (vertex z comes before x + y)
    Z,
    /// d in [x + y, s / 2], only when x + y <= z (vertex x + y comes before z)
    Xy,
}

impl LowerInterval {
    /// The interval containing height `d <= s / 2`, for sorted x <= y <= z.
    pub(crate) fn of_height(x: f32, y: f32, z: f32, d: f32) -> Self {
        if d < x {
            Self::None
        } else if d < y {
            Self::X
        } else if d < z && d < x + y {
            Self::Y
        } else if x + y > z {
            Self::Z
        } else {
            Self::Xy
        }
    }
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

/// Derivative of [`exact`] with respect to `h`. Because `g` is normalized, this
/// is the area of the water surface (the plane P clipped to the cube). Works
/// for any non-zero `g` with no threshold: each interval is the derivative of
/// the corresponding piece of [`exact`] and divides only by components its
/// range guarantees to be non-zero.
pub fn surface_area(g: Vec3, h: f32) -> f32 {
    let (x, y, z) = sorted_components(g);
    let t = x + y;
    let s = t + z;

    // V(h) = 1 - V(s - h), so dV/dh at h equals dV/dd of the lower half at
    // d = s - h: the area is symmetric around s / 2
    let d = if h > 0.5 * s { s - h } else { h };

    match LowerInterval::of_height(x, y, z, d) {
        LowerInterval::None => d * d / (2.0 * x * y * z),
        LowerInterval::X => (2.0 * d - x) / (2.0 * y * z),
        LowerInterval::Y => {
            let u = t - d;
            1.0 / z - u * u / (2.0 * x * y * z)
        }
        LowerInterval::Z => {
            let u = d - 0.5 * s;
            let c = 1.5 * (x * (2.0 * (y + z) - x) - (z - y) * (z - y));
            (c - 6.0 * u * u) / (6.0 * x * y * z)
        }
        LowerInterval::Xy => 1.0 / z,
    }
}

/// Inclusion-exclusion in f64 on the same f32 inputs, accurate enough to serve
/// as a reference in tests even for small components. Zero components are
/// dropped, giving the 2D or 1D formula.
#[cfg(test)]
pub(crate) fn reference(g: Vec3, h: f32) -> f64 {
    reference_inclusion_exclusion(g, h, 0)
}

/// The derivative of [`reference`] with respect to `h`, i.e. the reference
/// for the surface area.
#[cfg(test)]
pub(crate) fn reference_surface_area(g: Vec3, h: f32) -> f64 {
    reference_inclusion_exclusion(g, h, 1)
}

/// `derivative`-th derivative of the volume in f64, for the n non-zero
/// components c of `g`: the sum over subsets S of the components of
/// (-1)^|S| max(h - sum(S), 0)^(n - derivative) / ((n - derivative)! prod(c)).
#[cfg(test)]
fn reference_inclusion_exclusion(g: Vec3, h: f32, derivative: i32) -> f64 {
    let c: Vec<f64> = g
        .abs()
        .to_array()
        .iter()
        .filter(|&&c| c > 0.0)
        .map(|&c| c as f64)
        .collect();
    let n = c.len() as i32;
    let power = n - derivative;
    let factorial: f64 = (1..=power).map(|k| k as f64).product();
    let mut sum = 0.0;
    for subset in 0..(1 << n) {
        let f = h as f64
            - (0..n as usize)
                .filter(|i| subset >> i & 1 == 1)
                .map(|i| c[i])
                .sum::<f64>();
        let sign = if (subset as u32).count_ones().is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        // max(f, 0)^0 is the step function
        sum += sign * if f > 0.0 { f.powi(power) } else { 0.0 };
    }
    sum / (factorial * c.iter().product::<f64>())
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
            let h = rng.random::<f32>() * max_height(g);
            let (a, b) = (piecewise_3d_symmetric(g, h), exact(g, h));
            // both are within ~3e-7 of the f64 reference
            assert!((a - b).abs() < 1e-6, "g = {g}, h = {h}: {a} vs {b}");
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

    /// A random normal whose smallest component is set to `smallest` (with a
    /// random sign) before normalizing; `None` leaves it unchanged.
    fn normal_with_smallest(rng: &mut rand::rngs::ThreadRng, smallest: Option<f32>) -> Vec3 {
        let mut g = rand_normalized_gaussian(rng);
        if let Some(c) = smallest {
            let i = g.abs().min_position();
            g[i] = if rng.random::<bool>() { c } else { -c };
        }
        g.normalize()
    }

    /// Smallest components to test: none, log-uniform in [1e-8, 1e-2), and 0.
    fn smallest_components(rng: &mut rand::rngs::ThreadRng) -> [Option<f32>; 3] {
        [
            None,
            Some(10f32.powf(-8.0 + 6.0 * rng.random::<f32>())),
            Some(0.0),
        ]
    }

    #[test]
    fn exact_is_accurate_for_any_g() {
        let mut rng = rand::rng();
        for _ in 0..5_000 {
            for smallest in smallest_components(&mut rng) {
                let g = normal_with_smallest(&mut rng, smallest);
                let h = rng.random::<f32>() * max_height(g);
                let (v, r) = (exact(g, h), reference(g, h));
                assert!((v as f64 - r).abs() < 1e-6, "g = {g}, h = {h}: {v} vs {r}");
            }
        }
    }

    #[test]
    fn exact_axis_aligned() {
        for g in [
            vec3(0.0, 0.0, 1.0),
            vec3(-1.0, 0.0, 0.0),
            vec3(0.0, 1.0, 0.0),
        ] {
            for i in 0..=100 {
                let h = i as f32 / 100.0;
                assert!((exact(g, h) - h).abs() < 1e-6, "g = {g}, h = {h}");
            }
        }
    }

    #[test]
    fn exact_is_continuous_as_a_component_goes_to_zero() {
        let mut rng = rand::rng();
        for _ in 0..10_000 {
            let g = rand_normalized_gaussian(&mut rng);
            let i = g.abs().min_position();
            let (mut zero, mut tiny) = (g, g);
            zero[i] = 0.0;
            tiny[i] = 1e-7;
            let (zero, tiny) = (zero.normalize(), tiny.normalize());
            let h = rng.random::<f32>() * max_height(zero);
            let jump = (exact(zero, h) - exact(tiny, h)).abs();
            assert!(jump < 1e-6, "g = {zero}, h = {h}: jump {jump}");
        }
    }

    #[test]
    fn surface_area_is_accurate_for_any_g() {
        let mut rng = rand::rng();
        for _ in 0..5_000 {
            for smallest in smallest_components(&mut rng) {
                let g = normal_with_smallest(&mut rng, smallest);
                let h = rng.random::<f32>() * max_height(g);
                let (a, r) = (surface_area(g, h), reference_surface_area(g, h));
                // relative: the area is up to ~1.7, and grows like 1/z
                assert!(
                    (a as f64 - r).abs() < 1e-5 * r.max(1.0),
                    "g = {g}, h = {h}: {a} vs {r}"
                );
            }
        }
    }
}
