//! The intended use: many cells holding different amounts of water under the
//! same gravity, with gravity changing now and then. For each direction, one
//! `WaterCube` is set up and the water height of every cell is computed.
//!
//! Run with `cargo run --release --example many_cells`.

use std::time::Instant;

use glam::Vec3;
use rand::{RngExt, SeedableRng, rngs::StdRng};
use water_cube::WaterCube;
use water_cube::height::analytic;

const CELLS: usize = 1_000_000;

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    let fill: Vec<f32> = (0..CELLS).map(|_| rng.random::<f32>()).collect();
    let mut heights = vec![0.0; CELLS];

    // gravity tilting step by step from straight down towards a diagonal
    let directions: Vec<Vec3> = (0..=4)
        .map(|i| {
            let t = i as f32 / 4.0;
            Vec3::new(0.6 * t, 0.3 * t, 1.0).normalize()
        })
        .collect();

    let mut cube = WaterCube::new(directions[0]);
    for &g in &directions {
        // recompute the precomputed values for the new direction
        cube.set_g(g);

        let start = Instant::now();
        for (h, &v) in heights.iter_mut().zip(&fill) {
            *h = cube.height(v);
        }
        let with_cube = start.elapsed();

        // the same with the free function, for comparison
        let start = Instant::now();
        let mut check = 0.0;
        for &v in &fill {
            check += analytic(g, v);
        }
        let with_free = start.elapsed();
        std::hint::black_box(check);

        let mean = heights.iter().sum::<f32>() / CELLS as f32;
        println!(
            "g = {g:.3}: mean height {mean:.4} of {:.4}, {:.1} ns per cell with WaterCube, {:.1} ns with height::analytic",
            cube.max_height(),
            with_cube.as_nanos() as f64 / CELLS as f64,
            with_free.as_nanos() as f64 / CELLS as f64,
        );
    }
}
