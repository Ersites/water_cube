//! Benchmarks of the height-from-volume functions in `water_cube::height`.
//!
//! Each benchmark runs a function over a fixed, seeded batch of inputs and
//! reports throughput in solves per second.

use std::hint::black_box;
use std::time::Duration;

use criterion::measurement::WallTime;
use criterion::{
    BenchmarkGroup, BenchmarkId, Criterion, SamplingMode, Throughput, criterion_group, criterion_main,
};
use glam::{Vec3, vec3};
use rand::{RngExt, SeedableRng, rngs::StdRng};
use rand_distr::{Distribution, Normal};
use water_cube::height::{bisection, newton};
use water_cube::util::max_height;
use water_cube::volume::exact;

/// Number of distinct inputs per batch. Must be large enough that the branch
/// predictor can't memorize the sequence of inputs: with ~4k repeated inputs
/// branchy functions like `piecewise_3d` measured ~3x too fast. 2^18 inputs
/// (4 MB) still fit in L3, so memory bandwidth doesn't affect the results.
const BATCH: usize = 1 << 18;
const SAMPLES: usize = 20;
const SEED: u64 = 0x5eed;

type HeightFn = fn(Vec3, f32) -> f32;
const SOLVERS: [(&str, HeightFn); 2] = [("bisection", bisection), ("newton", newton)];

fn random_g(rng: &mut StdRng) -> Vec3 {
    let normal = Normal::new(0.0, 1.0).unwrap();
    vec3(normal.sample(rng), normal.sample(rng), normal.sample(rng)).normalize()
}

/// Random g whose smallest absolute component lies in [lo, hi).
fn random_g_with_min_component(rng: &mut StdRng, lo: f32, hi: f32) -> Vec3 {
    loop {
        let g = random_g(rng);
        let m = g.abs().min_element();
        if m >= lo && m < hi {
            return g;
        }
    }
}

/// One batch takes milliseconds, so use a fixed number of iterations per
/// sample instead of criterion's default linearly increasing ones.
fn configure(group: &mut BenchmarkGroup<'_, WallTime>) {
    group.sampling_mode(SamplingMode::Flat);
    group.sample_size(SAMPLES);
    group.measurement_time(Duration::from_secs(5));
}

fn run_batch(f: HeightFn, inputs: &[(Vec3, f32)]) -> f32 {
    let mut sum = 0.0;
    for &(g, v) in inputs {
        sum += f(black_box(g), black_box(v));
    }
    sum
}

/// Uniformly random g and v: the expected cost in a simulation.
fn random_inputs(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(SEED);
    let inputs: Vec<(Vec3, f32)> = (0..BATCH)
        .map(|_| (random_g(&mut rng), rng.random::<f32>()))
        .collect();

    let mut group = c.benchmark_group("random");
    group.throughput(Throughput::Elements(BATCH as u64));
    configure(&mut group);
    for (name, f) in SOLVERS {
        group.bench_function(name, |b| b.iter(|| run_batch(f, &inputs)));
    }

    // cost of a single volume evaluation, for scale: each solver calls
    // `exact` repeatedly
    let volume_inputs: Vec<(Vec3, f32)> = inputs
        .iter()
        .map(|&(g, v)| (g, v * max_height(g)))
        .collect();
    group.bench_function("exact (one volume evaluation)", |b| {
        b.iter(|| run_batch(exact, &volume_inputs))
    });
    group.finish();
}

/// Fixed v, random g: shows how the cost depends on the fill level.
fn by_volume(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(SEED);
    let gs: Vec<Vec3> = (0..BATCH).map(|_| random_g(&mut rng)).collect();

    let mut group = c.benchmark_group("by_volume");
    group.throughput(Throughput::Elements(BATCH as u64));
    configure(&mut group);
    for v in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99] {
        let inputs: Vec<(Vec3, f32)> = gs.iter().map(|&g| (g, v)).collect();
        for (name, f) in SOLVERS {
            group.bench_with_input(BenchmarkId::new(name, v), &inputs, |b, inputs| {
                b.iter(|| run_batch(f, inputs))
            });
        }
    }
    group.finish();
}

/// Random v, g restricted by its smallest component: small components make
/// `exact` noisy in f32, which can slow down convergence.
fn by_min_component(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(SEED);

    let mut group = c.benchmark_group("by_min_component");
    group.throughput(Throughput::Elements(BATCH as u64));
    configure(&mut group);
    for (label, lo, hi) in [
        ("1e-4..1e-3", 1e-4, 1e-3),
        ("1e-3..1e-2", 1e-3, 1e-2),
        ("1e-2..1e-1", 1e-2, 1e-1),
        ("1e-1..1", 1e-1, 1.0),
    ] {
        let inputs: Vec<(Vec3, f32)> = (0..BATCH)
            .map(|_| {
                let g = random_g_with_min_component(&mut rng, lo, hi);
                (g, rng.random::<f32>())
            })
            .collect();
        for (name, f) in SOLVERS {
            group.bench_with_input(BenchmarkId::new(name, label), &inputs, |b, inputs| {
                b.iter(|| run_batch(f, inputs))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, random_inputs, by_volume, by_min_component);
criterion_main!(benches);
