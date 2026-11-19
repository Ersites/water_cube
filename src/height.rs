//! Functions used to calculate the height of the plane under which there is a
//! given volume, i.e. the inverse of the functions in [`crate::volume`].
//!
//! All functions take a normal `g` and a volume `v` in [0, 1] and return the
//! height `h` in [0, `max_height(g)`] measured from the lowest vertex (see
//! README).

use glam::Vec3;

use crate::util::max_height;
use crate::volume::{LowerInterval, exact, sorted_components, surface_area};

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

/// Finds `h` such that `exact(g, h) == v` in closed form, for any non-zero
/// `g`: finds the interval of the piecewise volume containing the answer, then
/// solves that piece's polynomial directly (cube root, quadratic, linear or
/// cubic formula). Takes the same steps for every input, with no iteration.
///
/// Like [`exact`], it never divides by a component that can be zero, so it
/// needs no threshold for small components: with x = 0 the `X` and `Xy`
/// solutions become the 2D ones, and with x = y = 0 only `Xy` remains, giving
/// h = v.
pub fn analytic(g: Vec3, v: f32) -> f32 {
    let (x, y, z) = sorted_components(g);
    let s = x + y + z;

    // same reasoning as in `bisection`
    if v <= 0.0 {
        return 0.0;
    }
    if v >= 1.0 {
        // s equals `max_height(g)` up to rounding. Calling `max_height(g)`
        // here instead made every call ~14 ns slower, even though this branch
        // is rarely taken.
        return s;
    }

    // V(h) = 1 - V(s - h): solve in the lower half only
    if v > 0.5 {
        s - lower_height(x, y, z, 1.0 - v)
    } else {
        lower_height(x, y, z, v)
    }
}

/// Finds the interval containing the height with volume `v <= 0.5`, for
/// sorted x <= y <= z. The volume increases with the height, so this compares
/// `v` with the volume at the end of each interval, using the polynomials of
/// [`exact`] written so that no component that can be zero is divided by. The
/// end of an empty interval (x = 0 or y = 0) gets volume 0, so it is skipped.
fn lower_interval(x: f32, y: f32, z: f32, v: f32) -> LowerInterval {
    // end of `None` at d = x: x^3 / (6xyz) with x cancelled
    let v_x = if x > 0.0 { x * x / (6.0 * y * z) } else { 0.0 };
    if v < v_x {
        return LowerInterval::None;
    }

    // end of `X` at d = y: (x^2 - 3xy + 3y^2) / (6yz)
    let v_y = if y > 0.0 {
        (x * x - 3.0 * x * y + 3.0 * y * y) / (6.0 * y * z)
    } else {
        0.0
    };
    if v < v_y {
        return LowerInterval::X;
    }

    // `Y` ends at whichever of the vertices z and x + y comes first
    let t = x + y;
    let z_before_xy = t > z;
    let v_end_of_y = if z_before_xy {
        // at d = z; x + y > z >= y means x > 0 here
        let u = t - z;
        u * u * u / (6.0 * x * y * z) + (2.0 * z - t) / (2.0 * z)
    } else {
        // at d = x + y the cubic term vanishes
        t / (2.0 * z)
    };
    if v < v_end_of_y {
        return LowerInterval::Y;
    }

    if z_before_xy {
        LowerInterval::Z
    } else {
        LowerInterval::Xy
    }
}

/// Height `d <= s / 2` with volume `v <= 0.5`, for sorted x <= y <= z. Each
/// interval only divides by components that its range guarantees are
/// non-zero.
fn lower_height(x: f32, y: f32, z: f32, v: f32) -> f32 {
    let t = x + y;

    match lower_interval(x, y, z, v) {
        // d^3 / (6xyz) = v; only reached when x > 0
        LowerInterval::None => (6.0 * x * y * z * v).cbrt(),
        // (x^2 - 3dx + 3d^2) / (6yz) = v, the larger root of the quadratic
        LowerInterval::X => 0.5 * x + (2.0 * y * z * v - x * x / 12.0).max(0.0).sqrt(),
        // With u = t - d: u^3 - 6xy u + 3xy(t - 2zv) = 0, and u is in [0, x].
        // Only reached when x > 0.
        LowerInterval::Y => {
            let a = (2.0 * x * y).sqrt();
            t - middle_cubic_root(a, 3.0 * x * y * (t - 2.0 * z * v))
        }
        // The volume is symmetric around s / 2, so with d = s / 2 + u the
        // square term vanishes: u^3 - 3a^2 u + 3xyz(v - 1/2) = 0, with
        // 4a^2 = 2(xy + yz + zx) - (x^2 + y^2 + z^2), and u is in [z - s/2, 0].
        // Only reached when x + y > z, so x > 0.
        LowerInterval::Z => {
            // written to avoid cancellation; positive because x + y > z
            let a = (0.25 * (x * (2.0 * (y + z) - x) - (z - y) * (z - y))).sqrt();
            0.5 * (t + z) + middle_cubic_root(a, 3.0 * x * y * z * (v - 0.5))
        }
        // (2d - t) / (2z) = v
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::rand_normalized_gaussian;
    use crate::volume::{lower_volume, reference};
    use glam::vec3;
    use rand::RngExt;

    /// A random normal with its smallest component set to `smallest`, if
    /// given, before normalizing.
    fn normal_with_smallest(rng: &mut rand::rngs::ThreadRng, smallest: Option<f32>) -> Vec3 {
        let mut g = rand_normalized_gaussian(rng);
        if let Some(c) = smallest {
            let i = g.abs().min_position();
            g[i] = c;
        }
        g.normalize()
    }

    /// Smallest components to test: unchanged, log-uniform in [1e-8, 1e-2),
    /// and exactly 0.
    fn smallest_components(rng: &mut rand::rngs::ThreadRng) -> [Option<f32>; 3] {
        [
            None,
            Some(10f32.powf(-8.0 + 6.0 * rng.random::<f32>())),
            Some(0.0),
        ]
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
            let expected = LowerInterval::of_height(x, y, z, d);
            let found = lower_interval(x, y, z, lower_volume(x, y, z, d));
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
        check_lower_interval(1e-7, 0.6, 0.8, 10_000);
        // one zero component: only X and Xy are left
        let hits = check_lower_interval(0.0, 0.6, 0.8, 10_000);
        assert_eq!(
            hits[LowerInterval::X as usize] + hits[LowerInterval::Xy as usize],
            10_001
        );
        // two zero components: only Xy is left
        let hits = check_lower_interval(0.0, 0.0, 1.0, 10_000);
        assert_eq!(hits[LowerInterval::Xy as usize], 10_001);
    }

    #[test]
    fn lower_interval_random() {
        let mut rng = rand::rng();
        let mut both = [0; 2];
        for _ in 0..1000 {
            let (x, y, z) = sorted_components(rand_normalized_gaussian(&mut rng));
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
            for smallest in smallest_components(&mut rng) {
                let g = normal_with_smallest(&mut rng, smallest);
                let v = rng.random::<f32>();
                let h = bisection(g, v);

                assert!((0.0..=max_height(g)).contains(&h));
                // compare volumes, not heights: near h = 0 and h = max_height
                // the volume is flat, so a tiny volume error is a large height
                // error
                assert!((exact(g, h) - v).abs() < 1e-5, "g = {g}, v = {v}, h = {h}");
            }
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
            for smallest in smallest_components(&mut rng) {
                let g = normal_with_smallest(&mut rng, smallest);
                let v = rng.random::<f32>();
                let h = newton(g, v);

                assert!((0.0..=max_height(g)).contains(&h));
                assert!((exact(g, h) - v).abs() < 1e-5, "g = {g}, v = {v}, h = {h}");
            }
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
        let dh = 1e-3;
        for _ in 0..1000 {
            for smallest in smallest_components(&mut rng) {
                let g = normal_with_smallest(&mut rng, smallest);
                // The area has kinks with slopes up to ~1/(yz), with y the
                // second smallest component (with a zero component the volume
                // is only once differentiable). A central difference across a
                // kink is off by ~dh / (4yz), so this needs y not too small.
                // `surface_area` itself is checked for every g against the
                // exact f64 derivative in
                // `volume::tests::surface_area_is_accurate_for_any_g`.
                if sorted_components(g).1 < 0.2 {
                    continue;
                }
                let h = dh + rng.random::<f32>() * (max_height(g) - 2.0 * dh);
                let fd = (exact(g, h + dh) - exact(g, h - dh)) / (2.0 * dh);
                let a = surface_area(g, h);
                assert!((a - fd).abs() < 1e-2, "g = {g}, h = {h}: {a} vs {fd}");
            }
        }
    }

    /// Checks that `analytic(g, v)` has volume `v` against the f64 reference,
    /// which handles zero components and is independent of the f32 rounding
    /// of the volume functions.
    fn assert_analytic_inverts(g: Vec3, v: f32) {
        let h = analytic(g, v);
        assert!(
            (0.0..=max_height(g) + 1e-6).contains(&h),
            "g = {g}, v = {v}, h = {h}"
        );
        let volume = reference(g, h);
        assert!(
            (volume - v as f64).abs() < 1e-6,
            "g = {g}, v = {v}, h = {h}: volume {volume}"
        );
    }

    #[test]
    fn analytic_inverts_volume_for_any_g() {
        let mut rng = rand::rng();
        for _ in 0..5_000 {
            for smallest in smallest_components(&mut rng) {
                let g = normal_with_smallest(&mut rng, smallest);
                assert_analytic_inverts(g, rng.random::<f32>());
            }
        }
    }

    #[test]
    fn analytic_special_cases() {
        for g in [
            vec3(0.5, 0.6, 0.7),   // x + y > z
            vec3(0.2, 0.3, 0.9),   // x + y < z
            vec3(0.25, 0.5, 0.75), // x + y == z
            vec3(0.4, 0.4, 0.7),   // x == y
            vec3(1.0, 1.0, 1.0),   // diagonal
            vec3(0.0, 0.6, 0.8),   // one zero component
            vec3(0.0, 1.0, 1.0),   // one zero, the others equal
        ] {
            let g = g.normalize();
            for i in 1..10_000 {
                assert_analytic_inverts(g, i as f32 / 10_000.0);
            }
        }
        // two zero components: h == v
        let g = vec3(0.0, -1.0, 0.0);
        for i in 0..=100 {
            let v = i as f32 / 100.0;
            assert!((analytic(g, v) - v).abs() < 1e-6);
        }
    }

    #[test]
    fn analytic_is_continuous_as_a_component_goes_to_zero() {
        let mut rng = rand::rng();
        for _ in 0..10_000 {
            let g = rand_normalized_gaussian(&mut rng);
            let i = g.abs().min_position();
            let (mut zero, mut tiny) = (g, g);
            zero[i] = 0.0;
            tiny[i] = 1e-7;
            let v = rng.random::<f32>();
            let jump = (analytic(zero.normalize(), v) - analytic(tiny.normalize(), v)).abs();
            assert!(jump < 1e-6, "g = {zero}, v = {v}: jump {jump}");
        }
    }

    #[test]
    fn analytic_matches_bisection() {
        let mut rng = rand::rng();
        for _ in 0..1000 {
            for smallest in smallest_components(&mut rng) {
                let g = normal_with_smallest(&mut rng, smallest);
                // away from the flat ends, where heights are well defined
                let v = 0.05 + 0.9 * rng.random::<f32>();
                let (a, b) = (analytic(g, v), bisection(g, v));
                assert!((a - b).abs() < 1e-4, "g = {g}, v = {v}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn analytic_edge_cases() {
        let g = vec3(1.0, 2.0, 3.0).normalize();
        assert_eq!(analytic(g, 0.0), 0.0);
        assert!((analytic(g, 1.0) - max_height(g)).abs() < 1e-6);
        assert!((analytic(g, 0.5) - 0.5 * max_height(g)).abs() < 1e-6);
    }
}
