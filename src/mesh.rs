//! Geometry of the water in a cube, for rendering: the water surface polygon,
//! the faces of the whole body of water, and a simple triangle mesh to
//! collect them for many cells.
//!
//! All coordinates are in the cube's local space, [0, 1]^3. Add the cell's
//! position in the grid to place them (see [`TriangleMesh::add_polygon`]).

use glam::Vec3;

use crate::WaterCube;

/// A planar convex polygon with up to 6 vertices, in counter-clockwise order
/// when seen from the side its `normal` points to.
///
/// The water surface has at most 6 vertices (the plane crosses at most 6
/// edges of the cube) and a cube face clipped by the plane at most 5.
#[derive(Debug, Clone, Copy)]
pub struct Polygon {
    vertices: [Vec3; 6],
    len: usize,
    normal: Vec3,
}

impl Polygon {
    fn empty(normal: Vec3) -> Self {
        Self {
            vertices: [Vec3::ZERO; 6],
            len: 0,
            normal,
        }
    }

    /// The vertices, in counter-clockwise order around [`Polygon::normal`].
    /// Empty when there is no polygon (for example no water surface).
    pub fn vertices(&self) -> &[Vec3] {
        &self.vertices[..self.len]
    }

    /// The unit normal: `g` for the water surface, the outward axis
    /// direction for a cube face.
    pub fn normal(&self) -> Vec3 {
        self.normal
    }

    /// Whether the polygon has no area (fewer than 3 vertices).
    pub fn is_empty(&self) -> bool {
        self.len < 3
    }

    /// Area of the polygon.
    pub fn area(&self) -> f32 {
        // a fan of triangles from the first vertex; measuring from a vertex
        // rather than the origin keeps tiny polygons from cancelling to 0
        let v = self.vertices();
        let mut twice_area = Vec3::ZERO;
        for i in 1..v.len().saturating_sub(1) {
            twice_area += (v[i] - v[0]).cross(v[i + 1] - v[0]);
        }
        0.5 * twice_area.dot(self.normal)
    }

    /// Adds a vertex, skipping it if it repeats the previous one (which
    /// happens when the plane passes exactly through a vertex of the cube).
    fn push(&mut self, p: Vec3) {
        if self.len > 0 && self.vertices[self.len - 1].distance_squared(p) < 1e-12 {
            return;
        }
        // a convex polygon cut from a cube never exceeds 6 vertices
        debug_assert!(self.len < 6, "polygon has more than 6 vertices");
        if self.len < 6 {
            self.vertices[self.len] = p;
            self.len += 1;
        }
    }

    /// Drops a last vertex that repeats the first one, and empties polygons
    /// with no area.
    fn finish(mut self) -> Self {
        if self.len > 1 && self.vertices[self.len - 1].distance_squared(self.vertices[0]) < 1e-12 {
            self.len -= 1;
        }
        if self.len < 3 {
            self.len = 0;
        }
        self
    }
}

/// The faces of the body of water in a cube: the water surface and the parts
/// of the cube's faces below it, up to 7 polygons with outward normals.
/// Together they form a closed mesh.
#[derive(Debug, Clone, Copy)]
pub struct WaterPolygons {
    faces: [Polygon; 7],
    len: usize,
}

impl WaterPolygons {
    /// The non-empty faces: the water surface first (if there is one), then
    /// the wetted parts of the cube faces.
    pub fn faces(&self) -> &[Polygon] {
        &self.faces[..self.len]
    }
}

/// The 8 corners of the unit cube; bit i of the index is coordinate i.
const CORNERS: [Vec3; 8] = [
    Vec3::new(0.0, 0.0, 0.0),
    Vec3::new(1.0, 0.0, 0.0),
    Vec3::new(0.0, 1.0, 0.0),
    Vec3::new(1.0, 1.0, 0.0),
    Vec3::new(0.0, 0.0, 1.0),
    Vec3::new(1.0, 0.0, 1.0),
    Vec3::new(0.0, 1.0, 1.0),
    Vec3::new(1.0, 1.0, 1.0),
];

/// The 12 edges of the cube, as pairs of corner indices differing in one bit.
const EDGES: [(usize, usize); 12] = [
    (0, 1),
    (2, 3),
    (4, 5),
    (6, 7),
    (0, 2),
    (1, 3),
    (4, 6),
    (5, 7),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7),
];

/// The 6 faces of the cube as corner indices, counter-clockwise when seen
/// from outside, with their outward normals.
const FACES: [([usize; 4], Vec3); 6] = [
    ([0, 4, 6, 2], Vec3::NEG_X),
    ([1, 3, 7, 5], Vec3::X),
    ([0, 1, 5, 4], Vec3::NEG_Y),
    ([2, 6, 7, 3], Vec3::Y),
    ([0, 2, 3, 1], Vec3::NEG_Z),
    ([4, 5, 7, 6], Vec3::Z),
];

/// The point where the segment from `a` to `b` crosses the plane, given the
/// signed distances `fa` and `fb` of its ends (of opposite signs). Always
/// interpolates from the lower end, so the same edge gives exactly the same
/// point whichever direction it is walked in, which keeps the mesh watertight.
fn crossing(a: Vec3, fa: f32, b: Vec3, fb: f32) -> Vec3 {
    let (a, fa, b, fb) = if fa <= fb {
        (a, fa, b, fb)
    } else {
        (b, fb, a, fa)
    };
    a + (b - a) * (fa / (fa - fb))
}

impl WaterCube {
    /// Signed distance of each cube corner from the water plane at height
    /// `h`: negative below, positive above. The plane is
    /// g·(p - lowest corner) = h.
    fn corner_distances(&self, h: f32) -> [f32; 8] {
        let g = self.g();
        let lowest = 0.5 - g.signum() * 0.5;
        let offset = g.dot(lowest) + h;
        CORNERS.map(|c| g.dot(c) - offset)
    }

    /// The water surface at height `h`: the plane clipped to the cube, a
    /// convex polygon with 3 to 6 vertices, counter-clockwise when seen from
    /// above (its normal is `g`, normalized). Empty when `h` is outside
    /// (0, `max_height`), where the plane doesn't cut the cube.
    ///
    /// Its area equals [`WaterCube::surface_area`].
    pub fn surface_polygon(&self, h: f32) -> Polygon {
        let normal = self.g().normalize();
        let mut polygon = Polygon::empty(normal);
        if h.is_nan() || h <= 0.0 || h >= self.max_height() {
            return polygon;
        }

        // where the plane crosses the edges; a corner exactly on the plane is
        // found once per edge to a corner below it, so collect and dedupe
        let f = self.corner_distances(h);
        let mut points = [Vec3::ZERO; 12];
        let mut n = 0;
        for &(a, b) in &EDGES {
            if (f[a] < 0.0) != (f[b] < 0.0) {
                let p = crossing(CORNERS[a], f[a], CORNERS[b], f[b]);
                if points[..n].iter().all(|q| q.distance_squared(p) >= 1e-12) {
                    points[n] = p;
                    n += 1;
                }
            }
        }

        // order counter-clockwise around the normal: by angle in a basis
        // (u, v) of the plane with u × v = normal
        let center = points[..n].iter().sum::<Vec3>() / n as f32;
        let u = normal.any_orthonormal_vector();
        let v = normal.cross(u);
        let angle = |p: &Vec3| {
            let r = *p - center;
            r.dot(v).atan2(r.dot(u))
        };
        let points = &mut points[..n];
        points.sort_unstable_by(|p, q| angle(p).total_cmp(&angle(q)));

        for &p in points.iter() {
            polygon.push(p);
        }
        polygon.finish()
    }

    /// The faces of the body of water at height `h`: the water surface (see
    /// [`WaterCube::surface_polygon`]) and the parts of the 6 cube faces
    /// below it, each with an outward normal. Together they form a closed
    /// mesh whose volume equals [`WaterCube::volume`].
    ///
    /// Empty for `h <= 0`; the 6 full cube faces for `h >= max_height`.
    pub fn water_polygons(&self, h: f32) -> WaterPolygons {
        let mut body = WaterPolygons {
            faces: [Polygon::empty(Vec3::ZERO); 7],
            len: 0,
        };
        if h.is_nan() || h <= 0.0 {
            return body;
        }

        let surface = self.surface_polygon(h);
        if !surface.is_empty() {
            body.faces[0] = surface;
            body.len = 1;
        }

        // clip each cube face to the part below the plane (Sutherland-Hodgman
        // against a single plane); corners exactly on the plane count as below
        let f = self.corner_distances(h);
        for (corners, normal) in &FACES {
            let mut face = Polygon::empty(*normal);
            for i in 0..4 {
                let (a, b) = (corners[i], corners[(i + 1) % 4]);
                let (a_below, b_below) = (f[a] <= 0.0, f[b] <= 0.0);
                if a_below != b_below {
                    face.push(crossing(CORNERS[a], f[a], CORNERS[b], f[b]));
                }
                if b_below {
                    face.push(CORNERS[b]);
                }
            }
            let face = face.finish();
            if !face.is_empty() {
                body.faces[body.len] = face;
                body.len += 1;
            }
        }
        body
    }
}

/// A triangle mesh with flat-shaded faces, for collecting the polygons of
/// many cells into one mesh for a renderer: positions and normals per
/// vertex, and three vertex indices per triangle.
#[derive(Debug, Clone, Default)]
pub struct TriangleMesh {
    pub positions: Vec<Vec3>,
    pub normals: Vec<Vec3>,
    pub indices: Vec<u32>,
}

impl TriangleMesh {
    /// An empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Removes all vertices and triangles, keeping the allocated memory.
    pub fn clear(&mut self) {
        self.positions.clear();
        self.normals.clear();
        self.indices.clear();
    }

    /// Adds a polygon moved by `offset` (for example the cell's position in
    /// the grid), split into a fan of triangles. Its vertices get the
    /// polygon's normal, so each face is shaded flat. Triangles are
    /// counter-clockwise when seen from the side the normal points to.
    pub fn add_polygon(&mut self, polygon: &Polygon, offset: Vec3) {
        let vertices = polygon.vertices();
        if vertices.len() < 3 {
            return;
        }
        let first = self.positions.len() as u32;
        for &p in vertices {
            self.positions.push(p + offset);
            self.normals.push(polygon.normal());
        }
        for i in 1..vertices.len() as u32 - 1 {
            self.indices
                .extend_from_slice(&[first, first + i, first + i + 1]);
        }
    }

    /// Adds all faces of a body of water, moved by `offset`.
    pub fn add_water(&mut self, water: &WaterPolygons, offset: Vec3) {
        for face in water.faces() {
            self.add_polygon(face, offset);
        }
    }

    /// The number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::volume::{exact, normal_with_smallest, smallest_components, surface_area};
    use glam::vec3;
    use rand::RngExt;

    /// Volume enclosed by a closed mesh with outward-facing triangles
    /// (divergence theorem).
    fn enclosed_volume(mesh: &TriangleMesh) -> f32 {
        mesh.indices
            .chunks(3)
            .map(|t| {
                let (a, b, c) = (
                    mesh.positions[t[0] as usize],
                    mesh.positions[t[1] as usize],
                    mesh.positions[t[2] as usize],
                );
                a.dot(b.cross(c)) / 6.0
            })
            .sum()
    }

    fn directions(rng: &mut rand::rngs::ThreadRng) -> Vec<Vec3> {
        let mut gs = vec![
            vec3(0.0, 0.0, 1.0),
            vec3(0.0, -1.0, 0.0),
            vec3(0.0, 1.0, 1.0).normalize(),
            vec3(1.0, 1.0, 1.0).normalize(),
            vec3(-0.25, 0.5, -0.75).normalize(),
            vec3(0.2, 0.3, 0.9).normalize(),
        ];
        for _ in 0..300 {
            for smallest in smallest_components(rng) {
                gs.push(normal_with_smallest(rng, smallest));
            }
        }
        gs
    }

    #[test]
    fn surface_polygon_matches_surface_area() {
        let mut rng = rand::rng();
        for g in directions(&mut rng) {
            let cube = WaterCube::new(g);
            for _ in 0..20 {
                let h = rng.random::<f32>() * cube.max_height();
                let polygon = cube.surface_polygon(h);
                assert!(
                    (polygon.area() - surface_area(g, h)).abs() < 1e-4,
                    "g = {g}, h = {h}"
                );
                let lowest = 0.5 - g.signum() * 0.5;
                for &p in polygon.vertices() {
                    // inside the cube and on the plane
                    assert!(
                        p.cmpge(Vec3::splat(-1e-6)).all() && p.cmple(Vec3::splat(1.0 + 1e-6)).all()
                    );
                    assert!(
                        (g.dot(p - lowest) - h).abs() < 1e-5,
                        "g = {g}, h = {h}, p = {p}"
                    );
                }
                if !polygon.is_empty() {
                    assert!((3..=6).contains(&polygon.vertices().len()));
                }
            }
        }
    }

    #[test]
    fn water_polygons_enclose_the_volume() {
        let mut rng = rand::rng();
        let mut mesh = TriangleMesh::new();
        for g in directions(&mut rng) {
            let cube = WaterCube::new(g);
            for _ in 0..20 {
                let h = rng.random::<f32>() * cube.max_height();
                let water = cube.water_polygons(h);
                mesh.clear();
                mesh.add_water(&water, Vec3::ZERO);
                assert!(
                    (enclosed_volume(&mesh) - exact(g, h)).abs() < 1e-5,
                    "g = {g}, h = {h}: {} vs {}",
                    enclosed_volume(&mesh),
                    exact(g, h)
                );
                // closed and outward: the area vectors of the faces cancel
                let total: Vec3 = water.faces().iter().map(|f| f.normal() * f.area()).sum();
                assert!(total.length() < 1e-5, "g = {g}, h = {h}: {total}");
                // every face is counter-clockwise around its normal
                for face in water.faces() {
                    assert!(
                        face.area() > 0.0,
                        "g = {g}, h = {h}, face {:?} area {}",
                        face.vertices(),
                        face.area()
                    );
                }
            }
        }
    }

    #[test]
    fn empty_and_full_cubes() {
        let cube = WaterCube::new(vec3(0.3, -0.5, 0.8).normalize());
        let s = cube.max_height();
        for h in [-0.1, 0.0, f32::NAN] {
            assert!(cube.surface_polygon(h).is_empty());
            assert!(cube.water_polygons(h).faces().is_empty());
        }
        for h in [s, s + 0.1] {
            assert!(cube.surface_polygon(h).is_empty());
            let water = cube.water_polygons(h);
            assert_eq!(water.faces().len(), 6);
            assert!(water.faces().iter().all(|f| (f.area() - 1.0).abs() < 1e-6));
        }
    }

    #[test]
    fn axis_aligned_water_is_a_box() {
        let cube = WaterCube::new(vec3(0.0, 0.0, 1.0));
        let surface = cube.surface_polygon(0.25);
        assert_eq!(surface.vertices().len(), 4);
        assert!(surface.vertices().iter().all(|p| (p.z - 0.25).abs() < 1e-6));
        assert!((surface.area() - 1.0).abs() < 1e-6);
        // surface, bottom and four sides
        let water = cube.water_polygons(0.25);
        assert_eq!(water.faces().len(), 6);
    }

    #[test]
    fn triangle_mesh_layout() {
        let cube = WaterCube::new(vec3(0.3, 0.5, 0.8).normalize());
        let mut mesh = TriangleMesh::new();
        let surface = cube.surface_polygon(0.8);
        mesh.add_polygon(&surface, vec3(2.0, 0.0, 5.0));
        let n = surface.vertices().len();
        assert_eq!(mesh.positions.len(), n);
        assert_eq!(mesh.normals.len(), n);
        assert_eq!(mesh.triangle_count(), n - 2);
        assert!(mesh.indices.iter().all(|&i| (i as usize) < n));
        assert!(mesh.positions.iter().all(|p| p.x >= 2.0 && p.z >= 5.0));
    }
}
