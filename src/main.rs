use rand::RngExt;
use water_cube::util::{max_height, rand_normalized_gaussian};
use water_cube::volume::{formula, monte_carlo};

fn main() {
    let mut rng = rand::rng();
    let mut max_diff: f32 = 0.0;
    for _ in 0..100 {
        let g = rand_normalized_gaussian(&mut rng);
        let max_h = max_height(g);
        let h = rng.random::<f32>() * max_h;
        let v = monte_carlo(g, h);
        let vf = formula(g, h);
        max_diff = max_diff.max((v-vf).abs())
    }
    println!("{}", max_diff)
}
