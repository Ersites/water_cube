//! Benchmarks of the volume-from-height functions in `water_cube::volume`.
//!
//! Each benchmark runs a function over a fixed, seeded batch of inputs and
//! reports throughput in evaluations per second.

use std::hint::black_box;
use std::time::Duration;

use criterion::measurement::WallTime;
use criterion::{
    BenchmarkGroup, BenchmarkId, Criterion, SamplingMode, Throughput, criterion_group,
    criterion_main,
};
use glam::{Vec3, vec3};
use rand::{RngExt, SeedableRng, rngs::StdRng};
use rand_distr::{Distribution, Normal};
use water_cube::util::max_height;
use water_cube::volume::{
    exact, inclusion_exclusion_3d, inclusion_exclusion_3d_unrolled, monte_carlo, piecewise_3d,
    piecewise_3d_symmetric, surface_area,
};

/// Number of distinct inputs per batch. Must be large enough that the branch
/// predictor can't memorize the sequence of inputs: with ~4k repeated inputs
/// branchy functions like `piecewise_3d` measured ~3x too fast. 2^18 inputs
/// (4 MB) still fit in L3, so memory bandwidth doesn't affect the results.
const BATCH: usize = 1 << 18;
const SAMPLES: usize = 30;
const SEED: u64 = 0x5eed;

type VolumeFn = fn(Vec3, f32) -> f32;
const EXACT_3D: [(&str, VolumeFn); 5] = [
    ("exact", exact),
    ("inclusion_exclusion_3d", inclusion_exclusion_3d),
    (
        "inclusion_exclusion_3d_unrolled",
        inclusion_exclusion_3d_unrolled,
    ),
    ("piecewise_3d", piecewise_3d),
    ("piecewise_3d_symmetric", piecewise_3d_symmetric),
];

fn random_g(rng: &mut StdRng) -> Vec3 {
    let normal = Normal::new(0.0, 1.0).unwrap();
    vec3(normal.sample(rng), normal.sample(rng), normal.sample(rng)).normalize()
}

/// One batch takes milliseconds, so use a fixed number of iterations per
/// sample instead of criterion's default linearly increasing ones.
fn configure(group: &mut BenchmarkGroup<'_, WallTime>) {
    group.sampling_mode(SamplingMode::Flat);
    group.sample_size(SAMPLES);
    group.measurement_time(Duration::from_secs(5));
}

fn run_batch(f: VolumeFn, inputs: &[(Vec3, f32)]) -> f32 {
    let mut sum = 0.0;
    for &(g, h) in inputs {
        sum += f(black_box(g), black_box(h));
    }
    sum
}

/// Uniformly random g and h in [0, max_height(g)].
fn random_inputs(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(SEED);
    let inputs: Vec<(Vec3, f32)> = (0..BATCH)
        .map(|_| {
            let g = random_g(&mut rng);
            (g, rng.random::<f32>() * max_height(g))
        })
        .collect();

    let mut group = c.benchmark_group("random");
    group.throughput(Throughput::Elements(BATCH as u64));
    configure(&mut group);
    for (name, f) in EXACT_3D {
        group.bench_function(name, |b| b.iter(|| run_batch(f, &inputs)));
    }
    group.bench_function("surface_area", |b| {
        b.iter(|| run_batch(surface_area, &inputs))
    });
    group.finish();
}

/// Fixed fraction of max_height, random g. `piecewise_3d` takes a different
/// branch per interval; with a fixed fraction the branches are more
/// predictable than with random h.
fn by_height(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(SEED);
    let gs: Vec<Vec3> = (0..BATCH).map(|_| random_g(&mut rng)).collect();

    let mut group = c.benchmark_group("by_height");
    group.throughput(Throughput::Elements(BATCH as u64));
    configure(&mut group);
    for t in [0.05, 0.25, 0.5, 0.75, 0.95] {
        let inputs: Vec<(Vec3, f32)> = gs.iter().map(|&g| (g, t * max_height(g))).collect();
        for (name, f) in EXACT_3D {
            group.bench_with_input(BenchmarkId::new(name, t), &inputs, |b, inputs| {
                b.iter(|| run_batch(f, inputs))
            });
        }
    }
    group.finish();
}

/// `exact` for g with 3, 2 and 1 non-zero components (3D, 2D and 1D paths).
fn by_dimension(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(SEED);

    let mut group = c.benchmark_group("by_dimension");
    group.throughput(Throughput::Elements(BATCH as u64));
    configure(&mut group);
    for (label, mask) in [
        ("3d", vec3(1.0, 1.0, 1.0)),
        ("2d", vec3(1.0, 0.0, 1.0)),
        ("1d", vec3(0.0, 0.0, 1.0)),
    ] {
        let inputs: Vec<(Vec3, f32)> = (0..BATCH)
            .map(|_| {
                let g = loop {
                    let g = random_g(&mut rng) * mask;
                    if g.length() > 1e-3 {
                        break g.normalize();
                    }
                };
                (g, rng.random::<f32>() * max_height(g))
            })
            .collect();
        group.bench_with_input(BenchmarkId::new("exact", label), &inputs, |b, inputs| {
            b.iter(|| run_batch(exact, inputs))
        });
    }
    group.finish();
}

/// The Monte Carlo reference is ~1000x slower, so it gets a small batch and
/// fewer samples.
fn monte_carlo_reference(c: &mut Criterion) {
    const MC_BATCH: usize = 16;
    let mut rng = StdRng::seed_from_u64(SEED);
    let inputs: Vec<(Vec3, f32)> = (0..MC_BATCH)
        .map(|_| {
            let g = random_g(&mut rng);
            (g, rng.random::<f32>() * max_height(g))
        })
        .collect();

    let mut group = c.benchmark_group("random");
    group.throughput(Throughput::Elements(MC_BATCH as u64));
    group.sample_size(20);
    group.measurement_time(Duration::from_secs(5));
    group.bench_function("monte_carlo", |b| {
        b.iter(|| run_batch(monte_carlo, &inputs))
    });
    group.finish();
}

criterion_group!(
    benches,
    random_inputs,
    by_height,
    by_dimension,
    monte_carlo_reference
);
criterion_main!(benches);
