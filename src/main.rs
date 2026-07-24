use rand::RngExt;
use water_cube::util::{max_height, rand_normalized_gaussian};
use water_cube::volume::{formula, formula3, formula3an, formula3ex};

fn main() {
    let mut rng = rand::rng();
    for i in 0..100 {
        let g = rand_normalized_gaussian(&mut rng);
        let max_h = max_height(g);
        let h = rng.random::<f32>() * max_h;
        // println!("({}, {}, {}) at {h}", g.x, g.y, g.z);
        println!("{}", formula(g, h) - formula3an(g, h));
    }
}
