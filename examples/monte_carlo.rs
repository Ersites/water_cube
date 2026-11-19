//! Checks `volume::exact` against a brute-force Monte Carlo estimate, which
//! samples 10 000 random points in the cube. The two should agree to within
//! the Monte Carlo noise of about 0.005.
//!
//! Run with `cargo run --release --example monte_carlo`.

use rand::RngExt;
use water_cube::util::{max_height, rand_normalized_gaussian};
use water_cube::volume::exact;
use water_cube::volume::variants::monte_carlo;

fn main() {
    let mut rng = rand::rng();
    println!(
        "{:>28} {:>7} {:>8} {:>11} {:>10}",
        "up direction g", "height", "exact", "monte carlo", "difference"
    );
    for _ in 0..10 {
        let g = rand_normalized_gaussian(&mut rng);
        let h = rng.random::<f32>() * max_height(g);
        let (e, m) = (exact(g, h), monte_carlo(g, h));
        println!(
            "{:>28} {h:>7.3} {e:>8.4} {m:>11.4} {:>+10.4}",
            format!("{g:.3}"),
            e - m
        );
    }
}
