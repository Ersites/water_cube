//! Prints how the volume and the water surface area change with the height
//! of the water, for one direction of gravity.
//!
//! Run with `cargo run --example fill_curve`, or pass the up direction as
//! three numbers (it is normalized here): `cargo run --example fill_curve -- 0 1 1`.

use glam::Vec3;
use water_cube::WaterCube;

fn main() {
    let args: Vec<f32> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("expected three numbers, e.g. 0 1 1"))
        .collect();
    let g = match args[..] {
        [x, y, z] => Vec3::new(x, y, z),
        [] => Vec3::new(0.3, 0.5, 0.8),
        _ => panic!("expected three numbers, e.g. 0 1 1"),
    };
    assert!(g.length() > 0.0, "the up direction must be non-zero");
    let g = g.normalize();

    let cube = WaterCube::new(g);
    let steps = 20;
    let width = 40;

    println!(
        "up direction g = {g}, full at height {:.4}\n",
        cube.max_height()
    );
    println!("{:>7} {:>7} {:>7}  volume", "height", "volume", "area");
    for i in 0..=steps {
        let h = cube.max_height() * i as f32 / steps as f32;
        let v = cube.volume(h);
        let bar = "#".repeat((v * width as f32).round() as usize);
        println!("{h:>7.3} {v:>7.3} {:>7.3}  {bar}", cube.surface_area(h));
    }

    println!("\nand back: the height for a given volume");
    for v in [0.1, 0.25, 0.5, 0.75, 0.9] {
        let h = cube.height(v);
        println!(
            "  volume {v:>4} -> height {h:.4} -> volume {:.4}",
            cube.volume(h)
        );
    }
}
