# Water cube

A small Rust library for voxel water simulations. Given the direction of
gravity, it computes how much of a unit cube lies below a water surface at a
given height, and the inverse: where the water surface settles for a given
amount of water.

Both directions are computed in closed form, accurate to f32 precision for
every direction, including directions parallel to the cube's faces or edges.
A precomputed `WaterCube` makes repeated calls for the same direction fast.

## The problem

Take a unit cube between (0, 0, 0) and (1, 1, 1), filled with a fraction V of
water (or any fluid). The water settles with a flat surface, the plane P,
perpendicular to gravity.

- **G** is the normal of P: a normalized vector pointing opposite to gravity
  ("up"). Points "below the plane" are on the side of P opposite to G.
- **H** is the height of the water: the distance from the cube's lowest
  vertex (along G) to the plane P. Equivalently, project the cube and P onto
  G and measure from the lowest projected vertex to P.
- **H_max** is the height of the cube along G, the height at which the cube
  is full. It is the sum of the absolute components of G, between 1 (G along
  an axis) and √3 (G along a diagonal).

For a fixed G this gives two functions:

- **V(H)**, the volume below height H, increasing from 0 at H = 0 to 1 at
  H = H_max;
- **H(V)**, its inverse: the height of the water surface when the cube holds
  volume V.

The library also provides **dV/dH**, which equals the area of the water
surface inside the cube.

## Usage

```rust
use glam::Vec3;
use water_cube::WaterCube;

// "up" direction: opposite to gravity, normalized
let g = Vec3::new(0.3, 0.5, 0.8).normalize();

// precompute once per direction, then call as often as needed
let cube = WaterCube::new(g);

let v = cube.volume(0.4); // fraction of the cube below height 0.4
let h = cube.height(0.25); // water level when the cube is a quarter full
let a = cube.surface_area(h); // area of the water surface at that level
let h_max = cube.max_height(); // height at which the cube is full

// one-off calls without precomputation
let v2 = water_cube::volume::exact(g, 0.4);
let h2 = water_cube::height::analytic(g, 0.25);
```

| Use | For |
|---|---|
| `WaterCube` (`new`, `set_g`, `volume`, `height`, `surface_area`) | many calls with the same G |
| `volume::exact(g, h)`, `volume::surface_area(g, h)` | single volume or surface-area calls |
| `height::analytic(g, v)` | single height calls |
| `WaterCube::surface_polygon`, `WaterCube::water_polygons` | the water surface, or the whole body of water, as polygons for rendering |
| `TriangleMesh` | collecting those polygons from many cells into one triangle mesh |
| `volume::variants`, `height::variants` | alternative implementations, kept for comparison, tests and benchmarks |

### Meshes for rendering

`surface_polygon(h)` gives the water surface in a cell: the plane clipped to
the cube, a convex polygon with 3 to 6 vertices. `water_polygons(h)` gives the
whole body of water: the surface plus the wetted parts of the cube's faces,
up to 7 polygons that form a closed mesh. Coordinates are local to the cell,
in [0, 1]³. `TriangleMesh` collects the polygons of many cells, each moved to
its place in the grid, into positions, normals and triangle indices, the
layout renderers such as Bevy or rerun expect:

```rust
use water_cube::TriangleMesh;

let mut mesh = TriangleMesh::new();
for (cell_position, fill) in cells {
    let h = cube.height(fill);
    mesh.add_polygon(&cube.surface_polygon(h), cell_position); // surface only
    // or: mesh.add_water(&cube.water_polygons(h), cell_position);
}
```

### Examples

```sh
cargo run --example basic                    # the API above
cargo run --example fill_curve               # volume and surface area against height
cargo run --example fill_curve -- 0 1 1      # ... for a given up direction
cargo run --release --example many_cells     # a million cells under changing gravity
cargo run --release --example monte_carlo    # exact volumes against a Monte Carlo estimate
cargo run --example mesh_obj                 # water meshes for a grid of cells, as OBJ files
```

### Inputs

- G must be non-zero and should be normalized; it is not normalized for you.
  The volume only depends on G and H together, `exact(k·g, k·h) == exact(g, h)`,
  so with an unnormalized G the height is scaled by |G|.
- Heights outside [0, H_max] give a volume of 0 or 1 and a surface area of 0.
  Volumes outside [0, 1] give a height of 0 or H_max. NaN gives NaN.

## Accuracy and speed

The error is at most about 2·10⁻⁷ for the volume and 6·10⁻⁷ for the height,
measured against an f64 reference for directions whose smallest component
ranges from 0.58 down to exactly zero.

Time per call with random inputs, on an Intel Core i7-11390H:

| | Free function | `WaterCube` (one G, many inputs) |
|---|---|---|
| volume | 14 ns | 5 ns |
| surface area | 11 ns | 4 ns |
| height | 22 ns | 14 ns |

## How it works

- **The volume is a piecewise cubic polynomial in H.** It changes formula each
  time the plane passes a vertex of the cube, so there are up to 7 pieces
  between 0 and H_max.
- **Symmetry halves the work.** V(H) = 1 − V(H_max − H), so only heights up
  to H_max / 2 are evaluated, where at most 4 pieces can occur.
- **No special case for small components.** The usual inclusion–exclusion
  formula divides by the product of G's components and loses all precision as
  one of them approaches zero. Here each piece only divides by components its
  range guarantees to be non-zero, with the other factors cancelled by hand.
  So the 3D formulas turn into the 2D and 1D ones as components go to zero,
  without any threshold.
- **The inverse is solved directly.** After finding which piece contains the
  answer, each piece is solved with a cube root, the quadratic formula, a
  linear equation, or, for the two cubic pieces, the trigonometric form of the
  cubic formula. There is no iteration.
- **`WaterCube`** computes everything that depends only on G once, and
  evaluates the volume and surface area without data-dependent branches, so
  random inputs don't cause branch mispredictions.

## Development

```sh
cargo test                    # unit tests, including accuracy against an f64 reference
cargo bench --bench volume    # criterion benchmarks (several minutes each)
cargo bench --bench height
```

The code is in `src/`:

- `volume.rs`: `exact` and `surface_area`; `volume/variants.rs`: the
  alternatives (Monte Carlo, inclusion–exclusion, earlier piecewise versions)
- `height.rs`: `analytic`; `height/variants.rs`: bisection and Newton's method
- `cube.rs`: `WaterCube`
- `mesh.rs`: water surface and body polygons, `TriangleMesh`
- `util.rs`: helpers such as `max_height` and random directions for tests

## License

This project is licensed under the MIT License ([LICENSE](LICENSE)).
