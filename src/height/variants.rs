//! Iterative solvers for the height, kept for comparison and testing. Both
//! call [`exact`] repeatedly, so they are slower than [`super::analytic`],
//! and only as accurate as [`exact`].

use glam::Vec3;

use crate::util::max_height;
use crate::volume::{exact, surface_area};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::volume::{normal_with_smallest, smallest_components};
    use glam::vec3;
    use rand::RngExt;

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
}
