//! The basic API: volume, water height and surface area for one direction of
//! gravity.
//!
//! Run with `cargo run --example basic`.

use glam::Vec3;
use water_cube::WaterCube;

fn main() {
    // "up": the direction opposite to gravity, normalized
    let g = Vec3::new(0.3, 0.5, 0.8).normalize();

    // precompute everything that depends on g once
    let cube = WaterCube::new(g);
    println!("up direction g = {g}");
    println!("the cube is full at height {:.4}", cube.max_height());

    let h = 0.4;
    println!("volume below height {h}: {:.4}", cube.volume(h));

    let v = 0.25;
    let h = cube.height(v);
    println!(
        "with volume {v}, the water surface is at height {h:.4} and has area {:.4}",
        cube.surface_area(h)
    );

    // the free functions give the same results for one-off calls
    println!(
        "free functions: volume {:.4}, height {:.4}",
        water_cube::volume::exact(g, 0.4),
        water_cube::height::analytic(g, v)
    );
}
