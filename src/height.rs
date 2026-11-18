//! Functions used to calculate the height of the plane under which there is a
//! given volume, i.e. the inverse of the functions in [`crate::volume`].
//!
//! All functions take a normal `g` and a volume `v` in [0, 1] and return the
//! height `h` in [0, `max_height(g)`] measured from the lowest vertex (see
//! README).

use glam::Vec3;

use crate::util::max_height;
use crate::volume::{EPS, exact, surface_area};

const BISECTION_MAX_ITERATIONS: i32 = 64;
const BISECTION_TOLERANCE: f32 = 1e-6;
const NEWTON_MAX_ITERATIONS: i32 = 32;
const NEWTON_TOLERANCE: f32 = 1e-6;

/// Finds `h` such that `exact(g, h) == v` by bisection on [0, `max_height(g)`].
/// Always converges because the volume is continuous and increasing in `h`,
/// but needs ~20 evaluations of [`exact`] to get within the tolerance.
pub fn bisection(g: Vec3, v: f32) -> f32 {
    let max_h = max_height(g);

    // The volume is flat near both ends (cubic in the distance to the
    // vertex), so tiny f32 errors in `exact` would move the solution a lot.
    // Return the endpoints exactly instead.
    if v <= 0.0 {
        return 0.0;
    }
    if v >= 1.0 {
        return max_h;
    }

    let mut lo = 0.0;
    let mut hi = max_h;

    for _ in 0..BISECTION_MAX_ITERATIONS {
        if hi - lo < BISECTION_TOLERANCE {
            break;
        }

        let mid = 0.5 * (lo + hi);
        if exact(g, mid) < v {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    0.5 * (lo + hi)
}

/// Finds `h` such that `exact(g, h) == v` by Newton's method, using
/// [`surface_area`] as dV/dh. The volume is flat at both ends (convex then
/// concave), so plain Newton can overshoot; a bisection bracket is kept and a
/// bisection step is taken whenever the Newton step would leave it.
pub fn newton(g: Vec3, v: f32) -> f32 {
    let max_h = max_height(g);

    // same reasoning as in `bisection`
    if v <= 0.0 {
        return 0.0;
    }
    if v >= 1.0 {
        return max_h;
    }

    let mut lo = 0.0;
    let mut hi = max_h;
    let mut h = v * max_h;

    for _ in 0..NEWTON_MAX_ITERATIONS {
        let f = exact(g, h) - v;

        if f < 0.0 {
            lo = h;
        } else {
            hi = h;
        }

        let step = f / surface_area(g, h);

        // Check convergence before the bracket test: once converged, h is
        // one of the bracket ends and a tiny (or zero) step lands on it,
        // which would otherwise be rejected in favour of a bisection step.
        // A NaN step (zero derivative) fails this test and gets bisected.
        if step.abs() < NEWTON_TOLERANCE {
            return h - step;
        }

        let mut next = h - step;
        if !(next > lo && next < hi) {
            next = 0.5 * (lo + hi);
        }
        h = next;

        if hi - lo < NEWTON_TOLERANCE {
            break;
        }
    }

    h
}

/// Interval of the lower half (d <= s / 2) of the piecewise volume, named
/// after the highest vertex under the plane, as in `volume::piecewise_3d`.
/// Above s / 2 the symmetry V(d) = 1 - V(s - d) maps back into these.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LowerInterval {
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

/// The absolute components of `g`, sorted so that x <= y <= z.
fn sorted_components(g: Vec3) -> (f32, f32, f32) {
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

/// Finds the interval containing the height whose scaled volume is `w`, where
/// `w = 6xyz * v` and `v <= 0.5`. The volume increases with the height, so the
/// interval is found by comparing `w` with the scaled volume at each
/// breakpoint, using the same polynomials as `volume::piecewise_3d_symmetric`.
fn lower_interval(x: f32, y: f32, z: f32, w: f32) -> LowerInterval {
    // at d = x (end of `None`): d^3
    let w_x = x * x * x;
    if w < w_x {
        return LowerInterval::None;
    }

    // at d = y (end of `X`): x((x - 3d)x + 3d^2)
    let w_y = x * ((x - 3.0 * y) * x + 3.0 * y * y);
    if w < w_y {
        return LowerInterval::X;
    }

    // `Y` ends at whichever of the vertices z and x + y comes first
    let t = x + y;
    let z_before_xy = t > z;
    let w_end_of_y = if z_before_xy {
        // at d = z: (t - d)^3 + 3xy(2d - t)
        (t - z).powi(3) + 3.0 * x * y * (2.0 * z - t)
    } else {
        // at d = x + y the cubic term vanishes: 3xy(2t - t)
        3.0 * x * y * t
    };
    if w < w_end_of_y {
        return LowerInterval::Y;
    }

    if z_before_xy {
        LowerInterval::Z
    } else {
        LowerInterval::Xy
    }
}

/// Finds `h` such that `exact(g, h) == v` in closed form: finds the interval
/// of the piecewise volume containing the answer, then solves that piece's
/// polynomial directly (cube root, quadratic, linear or cubic formula). Takes
/// the same steps for every input, with no iteration.
pub fn analytic(g: Vec3, v: f32) -> f32 {
    // same reasoning as in `bisection`
    if v <= 0.0 {
        return 0.0;
    }
    if v >= 1.0 {
        return max_height(g);
    }

    let g = g.abs();
    let nx = (g.x > EPS) as i32;
    let ny = (g.y > EPS) as i32;
    let nz = (g.z > EPS) as i32;

    match nx + ny + nz {
        3 => analytic_3d(g, v),
        2 => analytic_2d(g, v),
        // the water fills a slab of height h
        1 => v,
        _ => 0.0,
    }
}

/// [`analytic`] for `g` with all components non-zero.
fn analytic_3d(g: Vec3, v: f32) -> f32 {
    let (x, y, z) = sorted_components(g);

    // V(h) = 1 - V(s - h): solve in the lower half only
    if v > 0.5 {
        x + y + z - lower_height_3d(x, y, z, 1.0 - v)
    } else {
        lower_height_3d(x, y, z, v)
    }
}

/// Height `d <= s / 2` with volume `v <= 0.5`, for sorted x <= y <= z.
fn lower_height_3d(x: f32, y: f32, z: f32, v: f32) -> f32 {
    let t = x + y;
    let s = t + z;
    let w = 6.0 * x * y * z * v;

    match lower_interval(x, y, z, w) {
        // d^3 = w
        LowerInterval::None => w.cbrt(),
        // 3x d^2 - 3x^2 d + x^3 = w, the larger root of the quadratic
        LowerInterval::X => 0.5 * x + (2.0 * y * z * v - x * x / 12.0).max(0.0).sqrt(),
        // (t - d)^3 + 3xy(2d - t) = w. With u = t - d this is
        // u^3 - 6xy u + 3xy(t - 2zv) = 0, and u is in [0, x].
        LowerInterval::Y => {
            let a = (2.0 * x * y).sqrt();
            t - middle_cubic_root(a, 3.0 * x * y * (t - 2.0 * z * v))
        }
        // -2d^3 + 3s d^2 - 3(x^2 + y^2 + z^2) d + x^3 + y^3 + z^3 = w. The
        // volume is symmetric around s / 2, so with d = s / 2 + u the square
        // term vanishes: u^3 - 3a^2 u + 3xyz(v - 1/2) = 0, with
        // 4a^2 = 2(xy + yz + zx) - (x^2 + y^2 + z^2), and u is in [z - s/2, 0].
        LowerInterval::Z => {
            // written to avoid cancellation; positive because x + y > z
            let a = (0.25 * (x * (2.0 * (y + z) - x) - (z - y) * (z - y))).sqrt();
            0.5 * s + middle_cubic_root(a, 3.0 * x * y * z * (v - 0.5))
        }
        // 3xy(2d - t) = w
        LowerInterval::Xy => 0.5 * t + z * v,
    }
}

/// Root of `u^3 - 3a^2 u + q = 0` with `|u| <= a`, for `a > 0`.
///
/// On [-a, a] the cubic is decreasing, so this root exists and is unique
/// whenever `|q| <= 2a^3`. The cubic usually has three real roots there, where
/// Cardano's formula needs complex numbers, so this uses the trigonometric
/// form: u_k = 2a cos(acos(-q / (2a^3)) / 3 - 2πk/3). The root in [-a, a] is
/// always the middle one, k = 1.
fn middle_cubic_root(a: f32, q: f32) -> f32 {
    let c = (-q / (2.0 * a * a * a)).clamp(-1.0, 1.0);
    2.0 * a * (c.acos() / 3.0 - 2.0 * std::f32::consts::FRAC_PI_3).cos()
}

/// [`analytic`] for `g` with exactly one zero component: the inverse of the
/// area under a line in the unit square, which is piecewise quadratic.
fn analytic_2d(g: Vec3, v: f32) -> f32 {
    let (a, b) = if g.x < EPS {
        (g.y, g.z)
    } else if g.y < EPS {
        (g.x, g.z)
    } else {
        (g.x, g.y)
    };
    let (a, b) = (a.min(b), a.max(b));

    // A(h) = 1 - A(a + b - h): solve in the lower half only
    let lower = |v: f32| {
        if v < a / (2.0 * b) {
            // triangle: d^2 / (2ab) = v
            (2.0 * a * b * v).sqrt()
        } else {
            // trapezoid: (2d - a) / (2b) = v
            b * v + 0.5 * a
        }
    };

    if v > 0.5 {
        a + b - lower(1.0 - v)
    } else {
        lower(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::rand_normalized_gaussian;
    use crate::volume::reference;
    use glam::vec3;
    use rand::RngExt;

    /// The interval containing height `d <= s / 2`, found directly from the
    /// vertex heights. This is the reference for `lower_interval`.
    fn lower_interval_of_height(x: f32, y: f32, z: f32, d: f32) -> LowerInterval {
        if d < x {
            LowerInterval::None
        } else if d < y {
            LowerInterval::X
        } else if d < z && d < x + y {
            LowerInterval::Y
        } else if x + y > z {
            LowerInterval::Z
        } else {
            LowerInterval::Xy
        }
    }

    /// Scaled volume `6xyz * V(d)` for `d <= s / 2`, written out per interval.
    fn scaled_volume(x: f32, y: f32, z: f32, d: f32) -> f32 {
        let (t, s) = (x + y, x + y + z);
        match lower_interval_of_height(x, y, z, d) {
            LowerInterval::None => d * d * d,
            LowerInterval::X => x * ((x - 3.0 * d) * x + 3.0 * d * d),
            LowerInterval::Y => (t - d).powi(3) + 3.0 * x * y * (2.0 * d - t),
            LowerInterval::Z => {
                let q = x * x + y * y + z * z;
                let r = x * x * x + y * y * y + z * z * z;
                ((-2.0 * d + 3.0 * s) * d - 3.0 * q) * d + r
            }
            LowerInterval::Xy => 3.0 * x * y * (2.0 * d - t),
        }
    }

    /// Distance from `d` to the nearest breakpoint of the lower half.
    fn distance_to_breakpoint(x: f32, y: f32, z: f32, d: f32) -> f32 {
        [x, y, z, x + y]
            .iter()
            .map(|b| (d - b).abs())
            .fold(f32::MAX, f32::min)
    }

    /// Checks `lower_interval` against the height-based reference for many
    /// heights in [0, s / 2] and returns how often each interval was hit.
    fn check_lower_interval(x: f32, y: f32, z: f32, steps: usize) -> [usize; 5] {
        let half = 0.5 * (x + y + z);
        let mut hits = [0; 5];
        for i in 0..=steps {
            let d = half * i as f32 / steps as f32;
            let expected = lower_interval_of_height(x, y, z, d);
            let found = lower_interval(x, y, z, scaled_volume(x, y, z, d));
            // right at a breakpoint, rounding may pick the neighbouring
            // interval, where both polynomials agree anyway
            if found != expected {
                assert!(
                    distance_to_breakpoint(x, y, z, d) < 1e-5,
                    "x = {x}, y = {y}, z = {z}, d = {d}: found {found:?}, expected {expected:?}"
                );
            }
            hits[expected as usize] += 1;
        }
        hits
    }

    #[test]
    fn lower_interval_z_before_xy() {
        // x + y > z: Y ends at z, then Z up to s / 2
        let (x, y, z) = sorted_components(vec3(0.5, 0.6, 0.7).normalize());
        assert!(x + y > z);
        let hits = check_lower_interval(x, y, z, 10_000);
        assert_eq!(hits[LowerInterval::Xy as usize], 0);
        assert!(hits[LowerInterval::Z as usize] > 0);
        assert!(hits[LowerInterval::None as usize] > 0);
        assert!(hits[LowerInterval::X as usize] > 0);
        assert!(hits[LowerInterval::Y as usize] > 0);
    }

    #[test]
    fn lower_interval_xy_before_z() {
        // x + y <= z: Y ends at x + y, then Xy up to s / 2
        let (x, y, z) = sorted_components(vec3(0.2, 0.3, 0.9).normalize());
        assert!(x + y <= z);
        let hits = check_lower_interval(x, y, z, 10_000);
        assert_eq!(hits[LowerInterval::Z as usize], 0);
        assert!(hits[LowerInterval::Xy as usize] > 0);
        assert!(hits[LowerInterval::None as usize] > 0);
        assert!(hits[LowerInterval::X as usize] > 0);
        assert!(hits[LowerInterval::Y as usize] > 0);
    }

    #[test]
    fn lower_interval_degenerate_cases() {
        // x + y == z exactly: the Z and Xy intervals have zero width
        check_lower_interval(0.25, 0.5, 0.75, 10_000);
        // x == y: the X interval has zero width
        check_lower_interval(0.4, 0.4, 0.7, 10_000);
        // x == y == z
        check_lower_interval(0.5, 0.5, 0.5, 10_000);
        // a tiny component
        check_lower_interval(1e-4, 0.6, 0.8, 10_000);
    }

    #[test]
    fn lower_interval_random() {
        let mut rng = rand::rng();
        let mut both = [0; 2];
        for _ in 0..1000 {
            let (x, y, z) = sorted_components(rand_normalized_gaussian(&mut rng));
            if x < 1e-3 {
                continue;
            }
            check_lower_interval(x, y, z, 200);
            both[(x + y > z) as usize] += 1;
        }
        // both orderings of z and x + y were exercised
        assert!(both[0] > 0 && both[1] > 0, "{both:?}");
    }

    #[test]
    fn bisection_inverts_volume() {
        let mut rng = rand::rng();
        for _ in 0..1000 {
            let g = rand_normalized_gaussian(&mut rng);
            // `exact` loses precision in f32 when a component of g is small
            // (residuals up to ~1e-2 for components below 1e-3), which would
            // make this test measure `exact` rather than the bisection.
            if g.abs().min_element() < 0.1 {
                continue;
            }
            let v = rng.random::<f32>();
            let h = bisection(g, v);

            assert!((0.0..=max_height(g)).contains(&h));
            // compare volumes, not heights: near h = 0 and h = max_height the
            // volume is flat, so a tiny volume error is a large height error
            assert!((exact(g, h) - v).abs() < 1e-5, "g = {g}, v = {v}, h = {h}");
        }
    }

    #[test]
    fn bisection_edge_cases() {
        let g = vec3(1.0, 2.0, 3.0).normalize();
        assert!(bisection(g, 0.0) < 1e-5);
        assert!((bisection(g, 1.0) - max_height(g)).abs() < 1e-5);

        // axis-aligned g: h == v
        let g = vec3(0.0, 0.0, -1.0);
        assert!((bisection(g, 0.3) - 0.3).abs() < 1e-5);
    }

    #[test]
    fn newton_inverts_volume() {
        let mut rng = rand::rng();
        for _ in 0..1000 {
            let g = rand_normalized_gaussian(&mut rng);
            // see `bisection_inverts_volume`
            if g.abs().min_element() < 0.1 {
                continue;
            }
            let v = rng.random::<f32>();
            let h = newton(g, v);

            assert!((0.0..=max_height(g)).contains(&h));
            assert!((exact(g, h) - v).abs() < 1e-5, "g = {g}, v = {v}, h = {h}");
        }
    }

    #[test]
    fn newton_edge_cases() {
        let g = vec3(1.0, 2.0, 3.0).normalize();
        assert_eq!(newton(g, 0.0), 0.0);
        assert_eq!(newton(g, 1.0), max_height(g));

        // axis-aligned and 2D g
        let g = vec3(0.0, 0.0, -1.0);
        assert!((newton(g, 0.3) - 0.3).abs() < 1e-5);
        let g = vec3(1.0, 0.0, 1.0).normalize();
        assert!((exact(g, newton(g, 0.3)) - 0.3).abs() < 1e-5);
    }

    #[test]
    fn surface_area_is_volume_derivative() {
        let mut rng = rand::rng();
        let dh = 1e-2;
        for _ in 0..1000 {
            let g = rand_normalized_gaussian(&mut rng);
            if g.abs().min_element() < 0.1 {
                continue;
            }
            let h = dh + rng.random::<f32>() * (max_height(g) - 2.0 * dh);
            let fd = (exact(g, h + dh) - exact(g, h - dh)) / (2.0 * dh);
            let a = surface_area(g, h);
            assert!((a - fd).abs() < 1e-2, "g = {g}, h = {h}: {a} vs {fd}");
        }
    }

    /// Checks that `analytic(g, v)` has volume `v`, using the f64 reference
    /// volume for 3D `g` and `exact` otherwise. The f64 reference keeps the
    /// check independent of the f32 rounding of the other volume functions.
    fn assert_analytic_inverts(g: Vec3, v: f32, tolerance: f32) {
        let h = analytic(g, v);
        assert!(
            (0.0..=max_height(g)).contains(&h),
            "g = {g}, v = {v}, h = {h}"
        );
        let volume = if g.abs().min_element() > EPS {
            reference(g, h)
        } else {
            exact(g, h) as f64
        };
        assert!(
            (volume - v as f64).abs() < tolerance as f64,
            "g = {g}, v = {v}, h = {h}: volume {volume}"
        );
    }

    #[test]
    fn analytic_inverts_volume() {
        let mut rng = rand::rng();
        for _ in 0..10_000 {
            let g = rand_normalized_gaussian(&mut rng);
            assert_analytic_inverts(g, rng.random::<f32>(), 1e-6);
        }
    }

    #[test]
    fn analytic_inverts_volume_for_small_components() {
        let mut rng = rand::rng();
        for _ in 0..10_000 {
            let mut g = rand_normalized_gaussian(&mut rng);
            // smallest component log-uniform in [1e-6, 1e-2)
            let i = g.abs().min_position();
            g[i] = 10f32.powf(-6.0 + 4.0 * rng.random::<f32>());
            assert_analytic_inverts(g.normalize(), rng.random::<f32>(), 1e-6);
        }
    }

    #[test]
    fn analytic_both_orderings_of_z_and_xy() {
        for g in [
            vec3(0.5, 0.6, 0.7),   // x + y > z
            vec3(0.2, 0.3, 0.9),   // x + y < z
            vec3(0.25, 0.5, 0.75), // x + y == z
            vec3(0.4, 0.4, 0.7),   // x == y
            vec3(1.0, 1.0, 1.0),   // diagonal
        ] {
            let g = g.normalize();
            for i in 1..10_000 {
                assert_analytic_inverts(g, i as f32 / 10_000.0, 1e-6);
            }
        }
    }

    #[test]
    fn analytic_matches_bisection() {
        let mut rng = rand::rng();
        for _ in 0..1000 {
            let g = rand_normalized_gaussian(&mut rng);
            // see `bisection_inverts_volume`
            if g.abs().min_element() < 0.1 {
                continue;
            }
            // away from the flat ends, where heights are well defined
            let v = 0.05 + 0.9 * rng.random::<f32>();
            let (a, b) = (analytic(g, v), bisection(g, v));
            assert!((a - b).abs() < 1e-4, "g = {g}, v = {v}: {a} vs {b}");
        }
    }

    #[test]
    fn analytic_2d_and_1d() {
        let mut rng = rand::rng();
        for _ in 0..1000 {
            let a = rng.random::<f32>() * std::f32::consts::TAU;
            let g = vec3(a.cos(), 0.0, a.sin());
            if g.abs().min_element() < EPS || g.x.abs().min(g.z.abs()) < 0.1 {
                continue;
            }
            assert_analytic_inverts(g, rng.random::<f32>(), 1e-5);
        }
        let g = vec3(0.0, -1.0, 0.0);
        assert_eq!(analytic(g, 0.3), 0.3);
    }

    #[test]
    fn analytic_edge_cases() {
        let g = vec3(1.0, 2.0, 3.0).normalize();
        assert_eq!(analytic(g, 0.0), 0.0);
        assert_eq!(analytic(g, 1.0), max_height(g));
        assert!((analytic(g, 0.5) - 0.5 * max_height(g)).abs() < 1e-6);
    }
}
