//! Functions used to calculate the volume of the unit cube under the plane P
//! with normal `g` at height `h` above the lowest vertex (see README).
//!
//! Use [`exact`] for the volume and [`surface_area`] for its derivative. The
//! [`variants`] module keeps the alternative implementations for comparison,
//! testing and benchmarks.
//!
//! # Inputs
//!
//! - `g` must be non-zero and should be normalized. Its signs and the order
//!   of its components don't matter. The volume only depends on `g` and `h`
//!   together, `exact(k * g, k * h) == exact(g, h)`, so with an unnormalized
//!   `g` the height is scaled by `|g|`, and [`surface_area`] is no longer
//!   dV/dh.
//! - `h` should be in [0, `max_height(g)`]. Outside that range the cube is
//!   empty (volume 0) or full (volume 1) and the surface area is 0. A NaN
//!   height gives NaN.
//!
//! The volume increases with `h`, from 0 at `h = 0` to 1 at
//! `h = max_height(g)`.

pub mod variants;

use glam::Vec3;

/// Exact volume for any non-zero `g`, including components that are zero or
/// close to zero, with no threshold.
///
/// Uses the lower-half polynomials of [`variants::piecewise_3d_symmetric`] (reflecting
/// `h` with V(h) = 1 - V(s - h)), but instead of dividing the result by 6xyz,
/// each interval divides only by components that its own range guarantees to
/// be non-zero, with the other factors cancelled by hand. So the 3D formulas
/// turn into the 2D and 1D ones as components go to zero, and the error stays
/// at f32 precision (~3e-7) for every `g`.
pub fn exact(g: Vec3, h: f32) -> f32 {
    let (x, y, z) = sorted_components(g);
    let s = x + y + z;

    // Heights outside [0, s] (for example from rounding) would reflect to a
    // negative d, which the polynomials don't handle: with a zero component
    // they divide by it. The cube is empty below 0 and full above s.
    // Written so that a NaN height still gives NaN.
    if h <= 0.0 {
        return 0.0;
    }
    if h >= s {
        return 1.0;
    }

    // V(h) = 1 - V(s - h): evaluate the lower half only
    let flip = h > 0.5 * s;
    let d = if flip { s - h } else { h };

    let v = lower_volume(x, y, z, d);
    if flip { 1.0 - v } else { v }
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

    // A NaN height would land in `Xy`, whose area doesn't depend on the
    // height, so pass it through explicitly.
    if h.is_nan() {
        return h;
    }
    // The plane doesn't cut the cube outside [0, s]; see `exact`.
    if h < 0.0 || h > s {
        return 0.0;
    }

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
        // `variants::piecewise_3d_symmetric`
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
/// after the highest vertex under the plane, as in [`variants::piecewise_3d`].
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

/// A random normal whose smallest component is set to `smallest` (with a
/// random sign) before normalizing; `None` leaves it unchanged.
#[cfg(test)]
pub(crate) fn normal_with_smallest(rng: &mut rand::rngs::ThreadRng, smallest: Option<f32>) -> Vec3 {
    let mut g = crate::util::rand_normalized_gaussian(rng);
    if let Some(c) = smallest {
        let i = g.abs().min_position();
        g[i] = if rand::RngExt::random::<bool>(rng) {
            c
        } else {
            -c
        };
    }
    g.normalize()
}

/// Smallest components to test: none, log-uniform in [1e-8, 1e-2), and 0.
#[cfg(test)]
pub(crate) fn smallest_components(rng: &mut rand::rngs::ThreadRng) -> [Option<f32>; 3] {
    [
        None,
        Some(10f32.powf(-8.0 + 6.0 * rand::RngExt::random::<f32>(rng))),
        Some(0.0),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::{max_height, rand_normalized_gaussian};
    use glam::vec3;
    use rand::RngExt;

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

    #[test]
    fn exact_and_surface_area_outside_the_cube() {
        for g in [
            vec3(0.3, 0.5, 0.8).normalize(), // 3D
            vec3(0.0, -0.6, 0.8),            // one zero component
            vec3(0.0, 0.0, 1.0),             // two zero components
        ] {
            let s = max_height(g);
            for h in [-1e-3, -1e-7, -f32::INFINITY] {
                assert_eq!(exact(g, h), 0.0, "g = {g}, h = {h}");
                assert_eq!(surface_area(g, h), 0.0, "g = {g}, h = {h}");
            }
            for h in [s + 1e-6, s + 1e-3, f32::INFINITY] {
                assert_eq!(exact(g, h), 1.0, "g = {g}, h = {h}");
                assert_eq!(surface_area(g, h), 0.0, "g = {g}, h = {h}");
            }
            assert!(exact(g, f32::NAN).is_nan());
            assert!(surface_area(g, f32::NAN).is_nan());
        }
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
}
