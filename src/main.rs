use rand::RngExt;
use water_cube::util::{max_height, rand_normalized_gaussian};
use water_cube::volume::{exact, monte_carlo};

fn main() {
    let mut rng = rand::rng();
    for _ in 0..10 {
        let g = rand_normalized_gaussian(&mut rng);
        let max_h = max_height(g);
        let h = rng.random::<f32>() * max_h;
        // println!("({}, {}, {}) at {h}", g.x, g.y, g.z);
        println!("{}", exact(g, h) - monte_carlo(g, h));
    }
}
