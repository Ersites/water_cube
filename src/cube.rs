//! [`WaterCube`]: volume, surface area and height for one fixed normal `g`,
//! with everything that depends only on `g` computed once.

use glam::Vec3;

use crate::volume::{LowerInterval, sorted_components};

/// Volume, surface area and height computations for a fixed normal `g`.
///
/// Gives the same results as [`crate::volume::exact`],
/// [`crate::volume::surface_area`] and [`crate::height::analytic`] (up to a
/// rounding step or two), but sorts the components and computes every
/// reciprocal, coefficient and breakpoint once in [`WaterCube::new`] /
/// [`WaterCube::set_g`]. Use it when many volumes or heights are needed for
/// the same `g`.
///
/// `g` must be non-zero and should be normalized; it is not normalized here.
/// See [`crate::volume`] for the conventions.
#[derive(Debug, Clone, Copy)]
pub struct WaterCube {
    g: Vec3,

    // |g| sorted so that x <= y <= z, with t = x + y and s = x + y + z, the
    // height of the cube along g
    x: f32,
    y: f32,
    z: f32,
    t: f32,
    s: f32,
    /// Whether vertex z is reached before vertex x + y (x + y > z), which
    /// decides whether the lower half ends in the `Z` or the `Xy` interval.
    z_before_xy: bool,

    // Factors of the volume and surface area polynomials. A reciprocal of a
    // component that is zero is infinite, but only intervals whose range
    // guarantees the component is non-zero use it.
    inv_6xyz: f32,
    inv_2xyz: f32,
    inv_6yz: f32,
    inv_2yz: f32,
    inv_2z: f32,
    inv_z: f32,
    /// Coefficient of the `Z` piece expanded around s / 2:
    /// 3/2 (2(xy + yz + zx) - (x^2 + y^2 + z^2)), written without cancellation.
    c: f32,

    // Volume at the end of the `None`, `X` and `Y` intervals, used by
    // `height` to find the interval; 0 for an empty interval.
    v_x: f32,
    v_y: f32,
    v_end_of_y: f32,

    // Constants of the per-interval solutions in `height`
    six_xyz: f32,
    two_yz: f32,
    x2_over_12: f32,
    three_xy: f32,
    two_z: f32,
    three_xyz: f32,
    /// The `Y` cubic u^3 - 3a^2 u + q = 0 has a^2 = 2xy; stores 2a and
    /// 1 / (2a^3).
    two_a_y: f32,
    inv_2a3_y: f32,
    /// The `Z` cubic has 4a^2 = 2(xy + yz + zx) - (x^2 + y^2 + z^2); stores 2a
    /// and 1 / (2a^3). Only meaningful when x + y > z.
    two_a_z: f32,
    inv_2a3_z: f32,
}

impl WaterCube {
    /// Precomputes everything that depends only on `g`.
    pub fn new(g: Vec3) -> Self {
        let (x, y, z) = sorted_components(g);
        let t = x + y;
        let s = t + z;
        let z_before_xy = t > z;

        let inv_6xyz = 1.0 / (6.0 * x * y * z);
        let inv_6yz = 1.0 / (6.0 * y * z);
        let inv_2z = 1.0 / (2.0 * z);
        let c = 1.5 * (x * (2.0 * (y + z) - x) - (z - y) * (z - y));

        // breakpoint volumes, written as in `volume::exact` so that no
        // component that can be zero is divided by
        let v_x = if x > 0.0 { x * x * inv_6yz } else { 0.0 };
        let v_y = if y > 0.0 {
            (x * x - 3.0 * x * y + 3.0 * y * y) * inv_6yz
        } else {
            0.0
        };
        let v_end_of_y = if z_before_xy {
            // at d = z; x + y > z >= y means x > 0 here
            let u = t - z;
            u * u * u * inv_6xyz + (2.0 * z - t) * inv_2z
        } else {
            // at d = x + y
            t * inv_2z
        };

        let a_y = (2.0 * x * y).sqrt();
        let a_z = (c / 6.0).sqrt();

        Self {
            g,
            x,
            y,
            z,
            t,
            s,
            z_before_xy,
            inv_6xyz,
            inv_2xyz: 3.0 * inv_6xyz,
            inv_6yz,
            inv_2yz: 3.0 * inv_6yz,
            inv_2z,
            inv_z: 2.0 * inv_2z,
            c,
            v_x,
            v_y,
            v_end_of_y,
            six_xyz: 6.0 * x * y * z,
            two_yz: 2.0 * y * z,
            x2_over_12: x * x / 12.0,
            three_xy: 3.0 * x * y,
            two_z: 2.0 * z,
            three_xyz: 3.0 * x * y * z,
            two_a_y: 2.0 * a_y,
            inv_2a3_y: 1.0 / (2.0 * a_y * a_y * a_y),
            two_a_z: 2.0 * a_z,
            inv_2a3_z: 1.0 / (2.0 * a_z * a_z * a_z),
        }
    }

    /// Replaces `g`, recomputing everything that depends on it.
    pub fn set_g(&mut self, g: Vec3) {
        *self = Self::new(g);
    }

    /// The normal `g`.
    pub fn g(&self) -> Vec3 {
        self.g
    }

    /// Height of the cube along `g`: the height at which it is full.
    pub fn max_height(&self) -> f32 {
        self.s
    }

    /// Volume under the plane at height `h`, as [`crate::volume::exact`].
    ///
    /// Has no data-dependent branches: it evaluates all five lower-half pieces
    /// and selects one, so random heights cause no branch mispredictions
    /// (about 1.8x fewer cycles per call than choosing the piece with a
    /// branch). Pieces of empty intervals may be infinite or NaN when a
    /// component is zero, but they are never selected.
    pub fn volume(&self, h: f32) -> f32 {
        // empty below 0, full above s; a NaN height gives NaN. These are
        // predictable for valid heights, so they stay branches.
        if h <= 0.0 {
            return 0.0;
        }
        if h >= self.s {
            return 1.0;
        }

        // V(h) = 1 - V(s - h): evaluate the lower half only
        let flip = h > 0.5 * self.s;
        let d = if flip { self.s - h } else { h };
        let (x, y, z, t) = (self.x, self.y, self.z, self.t);

        // the pieces of `volume::exact`, with the stored reciprocals
        let none = d * d * d * self.inv_6xyz;
        let in_x = (x * x - 3.0 * d * x + 3.0 * d * d) * self.inv_6yz;
        let u = t - d;
        let in_y = u * u * u * self.inv_6xyz + (2.0 * d - t) * self.inv_2z;
        let w = d - 0.5 * self.s;
        let in_z = 0.5 + w * (self.c - 2.0 * w * w) * self.inv_6xyz;
        let in_xy = (2.0 * d - t) * self.inv_2z;

        // same order as `LowerInterval::of_height`, applied from the last
        // interval back to the first (`z_before_xy` is fixed for a given g)
        let mut v = if self.z_before_xy { in_z } else { in_xy };
        // `&` rather than `&&`: both are branch-free, `&` measured slightly
        // faster
        v = if (d < z) & (d < t) { in_y } else { v };
        v = if d < y { in_x } else { v };
        v = if d < x { none } else { v };

        // `if flip { 1 - v } else { v }` without a comparison: every form of
        // that select tried (an `if`, `std::hint::select_unpredictable`,
        // f + (1 - 2f)v) compiled to a branch, which mispredicts half the time
        // for random heights. Here the choice comes from the sign bit of
        // s/2 - h instead: sign = +1 gives 0 + v, sign = -1 gives 1 - v, both
        // exactly.
        let sign = 1.0f32.copysign(0.5 * self.s - h);
        (0.5 - 0.5 * sign) + sign * v
    }

    /// Area of the water surface at height `h`, i.e. dV/dh, as
    /// [`crate::volume::surface_area`].
    ///
    /// Branch-free in the same way as [`WaterCube::volume`]: evaluates all
    /// five lower-half pieces and selects one.
    pub fn surface_area(&self, h: f32) -> f32 {
        // a NaN height would select `Xy`, whose area doesn't depend on it
        if h.is_nan() {
            return h;
        }
        // the plane doesn't cut the cube outside [0, s]
        if h < 0.0 || h > self.s {
            return 0.0;
        }

        // the area is symmetric around s / 2, so unlike the volume it needs
        // no flip afterwards
        let d = if h > 0.5 * self.s { self.s - h } else { h };
        let (x, y, z, t) = (self.x, self.y, self.z, self.t);

        // the derivatives of the pieces of `volume`
        let none = d * d * self.inv_2xyz;
        let in_x = (2.0 * d - x) * self.inv_2yz;
        let u = t - d;
        let in_y = self.inv_z - u * u * self.inv_2xyz;
        let w = d - 0.5 * self.s;
        let in_z = (self.c - 6.0 * w * w) * self.inv_6xyz;
        let in_xy = self.inv_z;

        // The intervals start at 0 <= x <= y <= min(z, x + y), so each piece
        // can overwrite the previous one once d reaches its start; this
        // selects the same piece as `LowerInterval::of_height`. (Selecting
        // from the last interval back, as in `volume`, compiled the final
        // select to a branch here: 0.3 misses per call.)
        let rest = if self.z_before_xy { in_z } else { in_xy };
        let mut a = none;
        a = if d >= x { in_x } else { a };
        a = if d >= y { in_y } else { a };
        // `|` rather than `||`, as with `&` in `volume`
        a = if (d >= z) | (d >= t) { rest } else { a };
        // Assign and return `a` rather than returning the `if` directly: the
        // direct form compiles to a branch (0.19 misses per call, ~15% slower)
        a
    }

    /// Height at which the volume under the plane is `v`, as
    /// [`crate::height::analytic`].
    pub fn height(&self, v: f32) -> f32 {
        // a NaN volume gives NaN
        if v <= 0.0 {
            return 0.0;
        }
        if v >= 1.0 {
            return self.s;
        }

        // V(h) = 1 - V(s - h): solve in the lower half only
        if v > 0.5 {
            self.s - self.lower_height(1.0 - v)
        } else {
            self.lower_height(v)
        }
    }

    /// Interval containing the height with volume `v <= 0.5`, from the stored
    /// breakpoint volumes.
    fn interval_of_volume(&self, v: f32) -> LowerInterval {
        if v < self.v_x {
            LowerInterval::None
        } else if v < self.v_y {
            LowerInterval::X
        } else if v < self.v_end_of_y {
            LowerInterval::Y
        } else if self.z_before_xy {
            LowerInterval::Z
        } else {
            LowerInterval::Xy
        }
    }

    /// Height `d <= s / 2` with volume `v <= 0.5`. Same solutions as
    /// `height::analytic`, with the constants precomputed.
    fn lower_height(&self, v: f32) -> f32 {
        let (x, z, t) = (self.x, self.z, self.t);
        match self.interval_of_volume(v) {
            // d^3 / (6xyz) = v
            LowerInterval::None => (self.six_xyz * v).cbrt(),
            // the larger root of the quadratic
            LowerInterval::X => 0.5 * x + (self.two_yz * v - self.x2_over_12).max(0.0).sqrt(),
            // u = t - d solves u^3 - 6xy u + 3xy(t - 2zv) = 0
            LowerInterval::Y => {
                let q = self.three_xy * (t - self.two_z * v);
                t - middle_cubic_root(self.two_a_y, self.inv_2a3_y, q)
            }
            // u = d - s / 2 solves u^3 - 3a^2 u + 3xyz(v - 1/2) = 0
            LowerInterval::Z => {
                let q = self.three_xyz * (v - 0.5);
                0.5 * self.s + middle_cubic_root(self.two_a_z, self.inv_2a3_z, q)
            }
            // (2d - t) / (2z) = v
            LowerInterval::Xy => 0.5 * t + z * v,
        }
    }
}

/// Root of `u^3 - 3a^2 u + q = 0` with `|u| <= a`, given `2a` and
/// `1 / (2a^3)`. Same as the helper in `height`: the trigonometric formula
/// for the middle of the three real roots.
fn middle_cubic_root(two_a: f32, inv_2a3: f32, q: f32) -> f32 {
    let c = (-q * inv_2a3).clamp(-1.0, 1.0);
    two_a * (c.acos() / 3.0 - 2.0 * std::f32::consts::FRAC_PI_3).cos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::height::analytic;
    use crate::util::max_height;
    use crate::volume::{
        exact, normal_with_smallest, reference, smallest_components, surface_area,
    };
    use glam::vec3;
    use rand::RngExt;

    #[test]
    fn matches_the_free_functions() {
        let mut rng = rand::rng();
        for _ in 0..2_000 {
            for smallest in smallest_components(&mut rng) {
                let g = normal_with_smallest(&mut rng, smallest);
                let cube = WaterCube::new(g);
                assert!((cube.max_height() - max_height(g)).abs() < 1e-6);

                for _ in 0..10 {
                    let h = rng.random::<f32>() * max_height(g);
                    let (a, b) = (cube.volume(h), exact(g, h));
                    assert!((a - b).abs() < 1e-6, "volume, g = {g}, h = {h}: {a} vs {b}");
                    let (a, b) = (cube.surface_area(h), surface_area(g, h));
                    assert!(
                        (a - b).abs() < 1e-5 * b.max(1.0),
                        "surface area, g = {g}, h = {h}: {a} vs {b}"
                    );

                    // away from the flat ends, where heights are well defined
                    let v = 0.01 + 0.98 * rng.random::<f32>();
                    let (a, b) = (cube.height(v), analytic(g, v));
                    assert!((a - b).abs() < 1e-5, "height, g = {g}, v = {v}: {a} vs {b}");
                }
            }
        }
    }

    #[test]
    fn is_accurate_for_any_g() {
        let mut rng = rand::rng();
        for _ in 0..2_000 {
            for smallest in smallest_components(&mut rng) {
                let g = normal_with_smallest(&mut rng, smallest);
                let cube = WaterCube::new(g);
                for _ in 0..10 {
                    let h = rng.random::<f32>() * cube.max_height();
                    let r = reference(g, h);
                    assert!(
                        (cube.volume(h) as f64 - r).abs() < 1e-6,
                        "volume, g = {g}, h = {h}"
                    );

                    let v = rng.random::<f32>();
                    let volume = reference(g, cube.height(v));
                    assert!(
                        (volume - v as f64).abs() < 1e-6,
                        "height, g = {g}, v = {v}: volume {volume}"
                    );
                }
            }
        }
    }

    #[test]
    fn both_orderings_and_zero_components() {
        for g in [
            vec3(0.5, 0.6, 0.7),   // x + y > z
            vec3(0.2, 0.3, 0.9),   // x + y < z
            vec3(0.25, 0.5, 0.75), // x + y == z
            vec3(0.4, 0.4, 0.7),   // x == y
            vec3(1.0, 1.0, 1.0),   // diagonal
            vec3(0.0, 0.6, 0.8),   // one zero component
            vec3(0.0, 1.0, 1.0),   // one zero, the others equal
            vec3(0.0, -1.0, 0.0),  // two zero components
        ] {
            let cube = WaterCube::new(g.normalize());
            for i in 1..10_000 {
                let v = i as f32 / 10_000.0;
                let volume = reference(cube.g(), cube.height(v));
                assert!((volume - v as f64).abs() < 1e-6, "g = {g}, v = {v}");
            }
        }
    }

    #[test]
    fn outside_the_cube_and_nan() {
        for g in [
            vec3(0.3, 0.5, 0.8).normalize(),
            vec3(0.0, -0.6, 0.8),
            vec3(0.0, 0.0, 1.0),
        ] {
            let cube = WaterCube::new(g);
            let s = cube.max_height();
            for h in [-1e-3, -1e-7, -f32::INFINITY] {
                assert_eq!(cube.volume(h), 0.0);
                assert_eq!(cube.surface_area(h), 0.0);
            }
            for h in [s + 1e-6, s + 1e-3, f32::INFINITY] {
                assert_eq!(cube.volume(h), 1.0);
                assert_eq!(cube.surface_area(h), 0.0);
            }
            assert_eq!(cube.height(-0.1), 0.0);
            assert_eq!(cube.height(1.1), s);
            assert!(cube.volume(f32::NAN).is_nan());
            assert!(cube.surface_area(f32::NAN).is_nan());
            assert!(cube.height(f32::NAN).is_nan());
        }
    }

    #[test]
    fn set_g_replaces_everything() {
        let (a, b) = (
            vec3(0.2, 0.3, 0.9).normalize(),
            vec3(0.6, 0.6, 0.5).normalize(),
        );
        let mut cube = WaterCube::new(a);
        cube.set_g(b);
        assert_eq!(cube.g(), b);
        let fresh = WaterCube::new(b);
        for i in 0..=100 {
            let h = cube.max_height() * i as f32 / 100.0;
            assert_eq!(cube.volume(h), fresh.volume(h));
            let v = i as f32 / 100.0;
            assert_eq!(cube.height(v), fresh.height(v));
        }
    }
}
