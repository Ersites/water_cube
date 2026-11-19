//! Builds the water meshes for a small grid of cells and writes them as OBJ
//! files, which any 3D viewer (Blender, MeshLab, ...) can open:
//!
//! - `water_surface.obj`: only the water surface in each cell
//! - `water_body.obj`: the closed body of water in each cell (surface, sides
//!   and bottom)
//!
//! Run with `cargo run --example mesh_obj`. The files are written to the
//! current directory.

use std::fmt::Write as _;
use std::fs;

use glam::Vec3;
use water_cube::{TriangleMesh, WaterCube};

fn main() -> std::io::Result<()> {
    // gravity tilted away from -z, so the water surfaces are slanted
    let up = Vec3::new(0.35, 0.2, 1.0).normalize();
    let cube = WaterCube::new(up);

    let mut surface = TriangleMesh::new();
    let mut body = TriangleMesh::new();

    // an 8 x 8 grid of cells, filled from empty in one corner to full in the
    // opposite one, with a gap of 0.2 between cells so they stay distinct
    let n = 8;
    for i in 0..n {
        for j in 0..n {
            let fill = (i + j) as f32 / (2 * (n - 1)) as f32;
            let h = cube.height(fill);
            let offset = Vec3::new(i as f32 * 1.2, j as f32 * 1.2, 0.0);
            surface.add_polygon(&cube.surface_polygon(h), offset);
            body.add_water(&cube.water_polygons(h), offset);
        }
    }

    write_obj("water_surface.obj", &surface)?;
    write_obj("water_body.obj", &body)?;
    println!(
        "wrote water_surface.obj ({} triangles) and water_body.obj ({} triangles)",
        surface.triangle_count(),
        body.triangle_count()
    );
    Ok(())
}

/// Writes a mesh in the Wavefront OBJ format, with per-vertex normals.
fn write_obj(path: &str, mesh: &TriangleMesh) -> std::io::Result<()> {
    let mut obj = String::new();
    for p in &mesh.positions {
        writeln!(obj, "v {} {} {}", p.x, p.y, p.z).unwrap();
    }
    for n in &mesh.normals {
        writeln!(obj, "vn {} {} {}", n.x, n.y, n.z).unwrap();
    }
    // OBJ indices start at 1; positions and normals share the same index
    for t in mesh.indices.chunks(3) {
        let [a, b, c] = [t[0] + 1, t[1] + 1, t[2] + 1];
        writeln!(obj, "f {a}//{a} {b}//{b} {c}//{c}").unwrap();
    }
    fs::write(path, obj)
}
