use crate::spatial::quadtree::Point2D;

/// Returns the signed area of triangle ABC (×2).
///
/// Positive: C is left of A→B (counterclockwise turn)
/// Negative: C is right of A→B (clockwise turn)
/// Zero: A, B, C are collinear
///
/// Uses robust adaptive precision arithmetic for nearly-collinear points.
#[inline]
fn orientation(a: &(f64, f64), b: &(f64, f64), c: &(f64, f64)) -> f64 {
    robust::orient2d(
        robust::Coord { x: a.0, y: a.1 },
        robust::Coord { x: b.0, y: b.1 },
        robust::Coord { x: c.0, y: c.1 },
    )
}

// =============================================================================
// Orientation Predicate
// =============================================================================

/// Returns the sign of the orientation of three points in 2D.
///
/// - Returns `1` if C is to the left of line AB (counterclockwise)
/// - Returns `-1` if C is to the right of line AB (clockwise)
/// - Returns `0` if A, B, C are collinear
///
/// Uses robust adaptive precision arithmetic via the `robust` crate.
#[inline]
pub fn orient_2d(a: &Point2D, b: &Point2D, c: &Point2D) -> i32 {
    let det = robust::orient2d(
        robust::Coord { x: a.x, y: a.y },
        robust::Coord { x: b.x, y: b.y },
        robust::Coord { x: c.x, y: c.y },
    );

    if det > 0.0 {
        1
    } else if det < 0.0 {
        -1
    } else {
        0
    }
}

/// Returns the raw determinant value for orientation (not just sign).
///
/// Useful when you need the magnitude for error estimation or other purposes.
#[inline]
pub fn orient_2d_det(a: &Point2D, b: &Point2D, c: &Point2D) -> f64 {
    robust::orient2d(
        robust::Coord { x: a.x, y: a.y },
        robust::Coord { x: b.x, y: b.y },
        robust::Coord { x: c.x, y: c.y },
    )
}

// =============================================================================
// Edge Crossing
// =============================================================================

/// Tests whether segment AB crosses segment CD.
///
/// Returns:
/// - `1` if AB properly crosses CD (interiors intersect)
/// - `0` if any two vertices from different edges are the same (vertex touch)
/// - `-1` if no crossing occurs
///
/// A "proper crossing" means the segments intersect at a point that is
/// interior to both segments (not at an endpoint).
///
/// Properties:
/// - crossing_sign_2d(a, b, c, d) == crossing_sign_2d(b, a, c, d)
/// - crossing_sign_2d(a, b, c, d) == crossing_sign_2d(c, d, a, b)
/// - Returns 0 if a==c, a==d, b==c, or b==d
pub fn crossing_sign_2d(a: &Point2D, b: &Point2D, c: &Point2D, d: &Point2D) -> i32 {
    // Check for shared vertices
    if points_equal(a, c) || points_equal(a, d) || points_equal(b, c) || points_equal(b, d) {
        return 0;
    }

    // Check for degenerate edges
    if points_equal(a, b) || points_equal(c, d) {
        return -1;
    }

    // For a proper crossing:
    // 1. C and D must be on opposite sides of line AB
    // 2. A and B must be on opposite sides of line CD
    let abc = orient_2d(a, b, c);
    let abd = orient_2d(a, b, d);
    let cda = orient_2d(c, d, a);
    let cdb = orient_2d(c, d, b);

    // If both points of one segment are on the same side of the other segment,
    // there is no crossing
    if abc * abd > 0 || cda * cdb > 0 {
        return -1;
    }

    // If any orientation is zero, points are collinear - need special handling
    if abc == 0 || abd == 0 || cda == 0 || cdb == 0 {
        // Collinear case - check if segments overlap
        // For our purposes (polygon containment), we treat collinear overlap
        // as "no proper crossing"
        return -1;
    }

    // Opposite sides on both tests means proper crossing
    1
}

/// Checks if two points are equal within floating-point tolerance.
#[inline]
fn points_equal(a: &Point2D, b: &Point2D) -> bool {
    a.x == b.x && a.y == b.y
}

/// Tests whether segment AB crosses segment CD, considering vertex crossings.
///
/// Returns `true` if:
/// - AB properly crosses CD (crossing_sign_2d returns 1), OR
/// - AB and CD share a vertex and the edges "cross" at that vertex according to the vertex crossing
///   rules for point containment tests
///
/// The vertex crossing rules ensure that when counting edge crossings from
/// a point to determine polygon containment, the count is consistent even
/// when the ray passes through a polygon vertex.
pub fn edge_or_vertex_crossing_2d(a: &Point2D, b: &Point2D, c: &Point2D, d: &Point2D) -> bool {
    let sign = crossing_sign_2d(a, b, c, d);
    if sign > 0 {
        return true;
    }
    if sign < 0 {
        return false;
    }

    // sign == 0: edges share a vertex
    vertex_crossing_2d(a, b, c, d)
}

/// Returns true if point P lies strictly on the interior of segment AB.
/// Assumes P is already known to be collinear with A and B.
#[inline]
fn point_on_segment_interior(p: &Point2D, a: &Point2D, b: &Point2D) -> bool {
    // Use the coordinate with larger range for numerical stability
    let dx = b.x - a.x;
    let dy = b.y - a.y;

    let t = if dx.abs() > dy.abs() {
        if dx == 0.0 {
            return false;
        }
        (p.x - a.x) / dx
    } else {
        if dy == 0.0 {
            return false;
        }
        (p.y - a.y) / dy
    };

    // Strictly interior: t in (0, 1), excluding endpoints
    t > 0.0 && t < 1.0
}

/// Determines whether edges AB and CD "cross" when they share a vertex.
///
/// This is used for point-in-polygon tests where we need consistent
/// crossing counts even when the test ray passes through polygon vertices.
///
/// The rule is: a crossing occurs at a shared vertex V if and only if
/// the other endpoints are on opposite sides of the line through V.
fn vertex_crossing_2d(a: &Point2D, b: &Point2D, c: &Point2D, d: &Point2D) -> bool {
    // Case 1: Endpoints are equal (shared vertex)
    if points_equal(a, c) {
        return orient_2d(d, a, b) > 0;
    }
    if points_equal(a, d) {
        return orient_2d(c, a, b) > 0;
    }
    if points_equal(b, c) {
        return orient_2d(d, b, a) > 0;
    }
    if points_equal(b, d) {
        return orient_2d(c, b, a) > 0;
    }

    // Case 2: Vertex C lies on segment AB (not at endpoints)
    // This happens when the test ray passes through a polygon vertex
    if orient_2d(a, b, c) == 0 && point_on_segment_interior(c, a, b) {
        // Count as crossing if D is to the left of ray AB
        return orient_2d(a, b, d) > 0;
    }

    // Case 3: Vertex D lies on segment AB (not at endpoints)
    if orient_2d(a, b, d) == 0 && point_on_segment_interior(d, a, b) {
        // Count as crossing if C is to the left of ray AB
        return orient_2d(a, b, c) > 0;
    }

    // No vertex touch detected
    false
}

/// Returns the signed crossing of AB with CD.
///
/// Returns:
/// - `-1` if AB crosses CD from left to right (exiting a region)
/// - `+1` if AB crosses CD from right to left (entering a region)
/// - `0` if AB does not cross CD
///
/// This is useful for computing winding numbers and signed area.
pub fn signed_edge_crossing_2d(a: &Point2D, b: &Point2D, c: &Point2D, d: &Point2D) -> i32 {
    if !edge_or_vertex_crossing_2d(a, b, c, d) {
        return 0;
    }

    // Determine the sign based on orientation
    // If A is to the right of CD and B is to the left, we're crossing left-to-right (-1)
    // If A is to the left of CD and B is to the right, we're crossing right-to-left (+1)
    let cda = orient_2d(c, d, a);
    if cda > 0 {
        // A is left of CD, so we're entering from left (or B is right)
        -1
    } else {
        1
    }
}

// =============================================================================
// EdgeCrosser2D - Efficient Edge Chain Crossing Test
// =============================================================================

/// Efficiently tests whether a fixed edge AB crosses a chain of edges.
///
/// This is the 2D analog of S2EdgeCrosser. It caches computations for the
/// fixed edge AB and reuses them when testing against multiple edges CD1,
/// CD2, CD3, etc.
///
/// # Example
///
/// ```ignore
/// let mut crosser = EdgeCrosser2D::new(&a, &b);
/// for edge in polygon_edges {
///     if crosser.edge_or_vertex_crossing(&edge.0, &edge.1) {
///         crossing_count += 1;
///     }
/// }
/// ```
pub struct EdgeCrosser2D {
    // Fixed edge AB
    a: Point2D,
    b: Point2D,

    // Previous vertex in the edge chain
    c: Point2D,

    // Cached orientation of triangle ACB (negated for efficiency)
    // acb stores -orient_2d(a, c, b) so we can check quickly
    acb: i32,
}

impl EdgeCrosser2D {
    /// Creates a new edge crosser for testing edge AB against other edges.
    pub fn new(a: &Point2D, b: &Point2D) -> Self {
        Self {
            a: *a,
            b: *b,
            c: Point2D::ORIGIN,
            acb: 0,
        }
    }

    /// Creates a new edge crosser and initializes the chain at vertex c.
    pub fn new_with_start(a: &Point2D, b: &Point2D, c: &Point2D) -> Self {
        let mut crosser = Self::new(a, b);
        crosser.restart_at(c);
        crosser
    }

    /// Restarts the edge chain at vertex c.
    ///
    /// Call this when the chain jumps to a new location (not continuing
    /// from the previous edge's endpoint).
    pub fn restart_at(&mut self, c: &Point2D) {
        self.c = *c;
        // We store -orient(a, c, b) = orient(a, b, c) for our crossing test
        self.acb = orient_2d(&self.a, &self.b, c);
    }

    /// Tests whether edge AB crosses edge CD.
    ///
    /// Returns:
    /// - `1` if AB properly crosses CD
    /// - `0` if edges share a vertex
    /// - `-1` if no crossing
    pub fn crossing_sign(&mut self, c: &Point2D, d: &Point2D) -> i32 {
        if !points_equal(c, &self.c) {
            self.restart_at(c);
        }
        self.crossing_sign_chain(d)
    }

    /// Tests crossing assuming c is the previous vertex.
    ///
    /// This is more efficient when testing a chain of edges where each
    /// edge's first vertex is the previous edge's second vertex.
    pub fn crossing_sign_chain(&mut self, d: &Point2D) -> i32 {
        // Check for shared vertices with AB
        if points_equal(&self.a, &self.c)
            || points_equal(&self.a, d)
            || points_equal(&self.b, &self.c)
            || points_equal(&self.b, d)
        {
            // Save d as the new c
            let old_c = self.c;
            self.c = *d;
            self.acb = orient_2d(&self.a, &self.b, d);

            // Return 0 for shared vertex
            if points_equal(&self.a, &old_c)
                || points_equal(&self.a, d)
                || points_equal(&self.b, &old_c)
                || points_equal(&self.b, d)
            {
                return 0;
            }
        }

        // Check for degenerate edge CD
        if points_equal(&self.c, d) {
            return -1;
        }

        // Save original c for collinear check
        let original_c = self.c;

        // Orientation of D relative to AB
        let abd = orient_2d(&self.a, &self.b, d);

        // If C and D are on the same side of AB, no crossing
        if self.acb * abd > 0 {
            // Update for next edge in chain
            self.c = *d;
            self.acb = abd;
            return -1;
        }

        // Now check if A and B are on opposite sides of CD
        let cda = orient_2d(&self.c, d, &self.a);
        let cdb = orient_2d(&self.c, d, &self.b);

        // Update for next edge in chain
        let old_acb = self.acb;
        self.c = *d;
        self.acb = abd;

        if cda * cdb > 0 {
            return -1;
        }

        // Handle cases where a polygon vertex lies on the test ray
        // This is a vertex touch, not a collinear overlap
        if old_acb == 0 && point_on_segment_interior(&original_c, &self.a, &self.b) {
            return 0; // Vertex C touches ray AB
        }
        if abd == 0 && point_on_segment_interior(d, &self.a, &self.b) {
            return 0; // Vertex D touches ray AB
        }

        // Handle other collinear cases (segment overlap, point on extended line)
        if old_acb == 0 || abd == 0 || cda == 0 || cdb == 0 {
            return -1;
        }

        // Proper crossing
        1
    }

    /// Tests whether edge AB crosses edge CD, including vertex crossings.
    pub fn edge_or_vertex_crossing(&mut self, c: &Point2D, d: &Point2D) -> bool {
        if !points_equal(c, &self.c) {
            self.restart_at(c);
        }
        self.edge_or_vertex_crossing_chain(d)
    }

    /// Tests crossing (including vertex) assuming c is the previous vertex.
    pub fn edge_or_vertex_crossing_chain(&mut self, d: &Point2D) -> bool {
        // Save c before it's overwritten
        let old_c = self.c;
        let sign = self.crossing_sign_chain(d);

        if sign > 0 {
            return true;
        }
        if sign < 0 {
            return false;
        }

        // sign == 0: check vertex crossing
        vertex_crossing_2d(&self.a, &self.b, &old_c, d)
    }

    /// Returns the signed crossing value for winding number computation.
    pub fn signed_edge_crossing(&mut self, c: &Point2D, d: &Point2D) -> i32 {
        if !self.edge_or_vertex_crossing(c, d) {
            return 0;
        }

        let cda = orient_2d(c, d, &self.a);
        if cda > 0 {
            -1
        } else {
            1
        }
    }
}

// =============================================================================
// Point-in-Polygon
// =============================================================================

/// Determines if a point is inside a simple polygon using ray casting.
///
/// The polygon is defined by a slice of vertices in order (either CW or CCW).
/// The polygon is implicitly closed (last vertex connects to first).
///
/// This is the "brute force" method that tests all edges. For indexed
/// polygons with precomputed origin_inside, use the faster containment
/// method that only tests edges intersecting a query cell.
///
/// # Arguments
///
/// * `point` - The query point
/// * `vertices` - Polygon vertices in order
///
/// # Returns
///
/// `true` if the point is strictly inside the polygon (not on the boundary)
pub fn brute_force_contains_2d(point: &Point2D, vertices: &[Point2D]) -> bool {
    if vertices.len() < 3 {
        return false;
    }

    // Choose a reference direction for the ray
    // We use a point far to the right (positive x direction)
    let ray_end = Point2D::new(1e18, point.y);

    let mut crossings = 0;
    let n = vertices.len();

    for i in 0..n {
        let v0 = &vertices[i];
        let v1 = &vertices[(i + 1) % n];

        if edge_or_vertex_crossing_2d(point, &ray_end, v0, v1) {
            crossings += 1;
        }
    }

    // Odd number of crossings means inside
    crossings % 2 == 1
}

/// Computes whether a reference origin point is inside the polygon.
///
/// The origin is chosen as a point outside the polygon's bounding box
/// (specifically, at the minimum corner minus a small offset). This
/// ensures the origin is always outside convex polygons and most
/// practical polygons.
///
/// # Returns
///
/// A tuple of (origin, origin_inside) where:
/// - `origin` is the computed reference point
/// - `origin_inside` is true if the origin is inside the polygon
///
/// For typical polygons, origin_inside will be false since we place
/// the origin outside the bounding box. However, for polygons that
/// wrap around (e.g., the exterior of a hole), it could be true.
pub fn compute_origin_inside_2d(vertices: &[Point2D]) -> (Point2D, bool) {
    if vertices.is_empty() {
        return (Point2D::ORIGIN, false);
    }

    // Compute bounding box
    let mut min_x = vertices[0].x;
    let mut min_y = vertices[0].y;

    for v in vertices.iter().skip(1) {
        if v.x < min_x {
            min_x = v.x;
        }
        if v.y < min_y {
            min_y = v.y;
        }
    }

    // Place origin outside the bounding box
    // Use a small offset to avoid exact corner cases
    let origin = Point2D::new(min_x - 1.0, min_y - 1.0);

    // Check if origin is inside (should be false for most polygons)
    let inside = brute_force_contains_2d(&origin, vertices);

    (origin, inside)
}

/// Determines if a point is inside a polygon using precomputed origin information.
///
/// This method counts edge crossings between the origin and the query point,
/// then XORs with the origin_inside flag. This is more efficient when testing
/// many points against the same polygon, especially when only a subset of edges
/// need to be tested (e.g., edges that intersect a query cell).
///
/// # Arguments
///
/// * `point` - The query point
/// * `origin` - The precomputed reference origin
/// * `origin_inside` - Whether the origin is inside the polygon
/// * `vertices` - Polygon vertices in order
///
/// # Returns
///
/// `true` if the point is inside the polygon
pub fn contains_with_origin_2d(
    point: &Point2D,
    origin: &Point2D,
    origin_inside: bool,
    vertices: &[Point2D],
) -> bool {
    if vertices.len() < 3 {
        return false;
    }

    let mut inside = origin_inside;
    let n = vertices.len();

    let mut crosser = EdgeCrosser2D::new(origin, point);

    for i in 0..n {
        let v0 = &vertices[i];
        let v1 = &vertices[(i + 1) % n];

        if crosser.edge_or_vertex_crossing(v0, v1) {
            inside = !inside;
        }
    }

    inside
}

/// Determines if a point is inside a polygon, testing only specified edges.
///
/// This is the most efficient method when you have precomputed which edges
/// intersect the query region (e.g., from a quadtree cell's edge list).
///
/// # Arguments
///
/// * `point` - The query point
/// * `origin` - The precomputed reference origin
/// * `origin_inside` - Whether the origin is inside the polygon
/// * `vertices` - All polygon vertices
/// * `edge_indices` - Indices of edges to test (edge i connects vertex i to i+1)
///
/// # Returns
///
/// `true` if the point is inside the polygon
pub fn contains_with_edges_2d(
    point: &Point2D,
    origin: &Point2D,
    origin_inside: bool,
    vertices: &[Point2D],
    edge_indices: &[usize],
) -> bool {
    let mut inside = origin_inside;
    let n = vertices.len();

    let mut crosser = EdgeCrosser2D::new(origin, point);

    for &i in edge_indices {
        let v0 = &vertices[i];
        let v1 = &vertices[(i + 1) % n];

        if crosser.edge_or_vertex_crossing(v0, v1) {
            inside = !inside;
        }
    }

    inside
}
#[cfg(test)]
mod predicate_tests {
    use super::*;

    // -------------------------------------------------------------------------
    // Orientation Tests
    // -------------------------------------------------------------------------

    mod orientation_tests {
        use super::*;

        #[test]
        fn test_orient_2d_ccw() {
            // Counterclockwise triangle
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(1.0, 0.0);
            let c = Point2D::new(0.5, 1.0);

            assert_eq!(orient_2d(&a, &b, &c), 1);
        }

        #[test]
        fn test_orient_2d_cw() {
            // Clockwise triangle
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(0.5, 1.0);
            let c = Point2D::new(1.0, 0.0);

            assert_eq!(orient_2d(&a, &b, &c), -1);
        }

        #[test]
        fn test_orient_2d_collinear() {
            // Collinear points
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(1.0, 1.0);
            let c = Point2D::new(2.0, 2.0);

            assert_eq!(orient_2d(&a, &b, &c), 0);
        }

        #[test]
        fn test_orient_2d_rotation_invariant() {
            // Sign should be preserved under rotation of arguments
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(1.0, 0.0);
            let c = Point2D::new(0.5, 1.0);

            let sign_abc = orient_2d(&a, &b, &c);
            let sign_bca = orient_2d(&b, &c, &a);
            let sign_cab = orient_2d(&c, &a, &b);

            assert_eq!(sign_abc, sign_bca);
            assert_eq!(sign_bca, sign_cab);
        }

        #[test]
        fn test_orient_2d_swap_inverts() {
            // Swapping two arguments should negate the result
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(1.0, 0.0);
            let c = Point2D::new(0.5, 1.0);

            let sign_abc = orient_2d(&a, &b, &c);
            let sign_bac = orient_2d(&b, &a, &c);

            assert_eq!(sign_abc, -sign_bac);
        }

        #[test]
        fn test_orient_2d_nearly_collinear() {
            // Nearly collinear points should still give correct sign
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(1.0, 0.0);
            let c = Point2D::new(0.5, 1e-15); // Very slightly above the line

            // Should be CCW (positive)
            assert_eq!(orient_2d(&a, &b, &c), 1);

            let d = Point2D::new(0.5, -1e-15); // Very slightly below the line
            assert_eq!(orient_2d(&a, &b, &d), -1);
        }
    }

    // -------------------------------------------------------------------------
    // Edge Crossing Tests
    // -------------------------------------------------------------------------

    mod crossing_tests {
        use super::*;

        #[test]
        fn test_crossing_sign_proper_cross() {
            // Two segments that properly cross
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(2.0, 2.0);
            let c = Point2D::new(0.0, 2.0);
            let d = Point2D::new(2.0, 0.0);

            assert_eq!(crossing_sign_2d(&a, &b, &c, &d), 1);
        }

        #[test]
        fn test_crossing_sign_no_cross_parallel() {
            // Parallel segments
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(2.0, 0.0);
            let c = Point2D::new(0.0, 1.0);
            let d = Point2D::new(2.0, 1.0);

            assert_eq!(crossing_sign_2d(&a, &b, &c, &d), -1);
        }

        #[test]
        fn test_crossing_sign_no_cross_disjoint() {
            // Disjoint segments
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(1.0, 0.0);
            let c = Point2D::new(2.0, 0.0);
            let d = Point2D::new(3.0, 0.0);

            assert_eq!(crossing_sign_2d(&a, &b, &c, &d), -1);
        }

        #[test]
        fn test_crossing_sign_shared_vertex() {
            // Segments share a vertex
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(1.0, 1.0);
            let c = Point2D::new(1.0, 1.0);
            let d = Point2D::new(2.0, 0.0);

            assert_eq!(crossing_sign_2d(&a, &b, &c, &d), 0);
        }

        #[test]
        fn test_crossing_sign_t_intersection() {
            // T-intersection: one endpoint on the other segment
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(2.0, 0.0);
            let c = Point2D::new(1.0, 0.0);
            let d = Point2D::new(1.0, 1.0);

            // This is collinear (c is on AB), should return -1
            assert_eq!(crossing_sign_2d(&a, &b, &c, &d), -1);
        }

        #[test]
        fn test_crossing_sign_symmetric() {
            // crossing_sign(a,b,c,d) == crossing_sign(c,d,a,b)
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(2.0, 2.0);
            let c = Point2D::new(0.0, 2.0);
            let d = Point2D::new(2.0, 0.0);

            assert_eq!(
                crossing_sign_2d(&a, &b, &c, &d),
                crossing_sign_2d(&c, &d, &a, &b)
            );
        }

        #[test]
        fn test_edge_or_vertex_crossing_proper() {
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(2.0, 2.0);
            let c = Point2D::new(0.0, 2.0);
            let d = Point2D::new(2.0, 0.0);

            assert!(edge_or_vertex_crossing_2d(&a, &b, &c, &d));
        }

        #[test]
        fn test_edge_or_vertex_crossing_no_cross() {
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(1.0, 0.0);
            let c = Point2D::new(0.0, 1.0);
            let d = Point2D::new(1.0, 1.0);

            assert!(!edge_or_vertex_crossing_2d(&a, &b, &c, &d));
        }
    }

    // -------------------------------------------------------------------------
    // EdgeCrosser2D Tests
    // -------------------------------------------------------------------------

    mod edge_crosser_tests {
        use super::*;

        #[test]
        fn test_edge_crosser_basic() {
            let a = Point2D::new(0.0, 1.0);
            let b = Point2D::new(2.0, 1.0);

            let mut crosser = EdgeCrosser2D::new(&a, &b);

            // Edge that crosses AB
            let c = Point2D::new(1.0, 0.0);
            let d = Point2D::new(1.0, 2.0);
            assert_eq!(crosser.crossing_sign(&c, &d), 1);

            // Edge that doesn't cross AB
            let e = Point2D::new(3.0, 0.0);
            let f = Point2D::new(3.0, 2.0);
            assert_eq!(crosser.crossing_sign(&e, &f), -1);
        }

        #[test]
        fn test_edge_crosser_chain() {
            // Test crossing against a chain of edges
            let a = Point2D::new(0.5, 0.5);
            let b = Point2D::new(1.5, 0.5);

            // Square from (0,0) to (1,1)
            let square = [
                Point2D::new(0.0, 0.0),
                Point2D::new(1.0, 0.0),
                Point2D::new(1.0, 1.0),
                Point2D::new(0.0, 1.0),
            ];

            let mut crosser = EdgeCrosser2D::new(&a, &b);
            let mut crossings = 0;

            // Test each edge of the square
            for i in 0..4 {
                let v0 = &square[i];
                let v1 = &square[(i + 1) % 4];
                if crosser.edge_or_vertex_crossing(v0, v1) {
                    crossings += 1;
                }
            }

            // AB starts at (0.5, 0.5) which is inside the square,
            // and ends at (1.5, 0.5) which is outside to the right.
            // The ray only crosses the right edge (x=1), not the left edge.
            assert_eq!(crossings, 1);
        }

        #[test]
        fn test_edge_crosser_consistent_with_standalone() {
            // EdgeCrosser should give same results as standalone function
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(2.0, 2.0);
            let c = Point2D::new(0.0, 2.0);
            let d = Point2D::new(2.0, 0.0);

            let mut crosser = EdgeCrosser2D::new(&a, &b);
            let crosser_result = crosser.edge_or_vertex_crossing(&c, &d);
            let standalone_result = edge_or_vertex_crossing_2d(&a, &b, &c, &d);

            assert_eq!(crosser_result, standalone_result);
        }
    }

    // -------------------------------------------------------------------------
    // Point-in-Polygon Tests
    // -------------------------------------------------------------------------

    mod containment_tests {
        use super::*;

        fn unit_square() -> Vec<Point2D> {
            vec![
                Point2D::new(0.0, 0.0),
                Point2D::new(1.0, 0.0),
                Point2D::new(1.0, 1.0),
                Point2D::new(0.0, 1.0),
            ]
        }

        #[test]
        fn test_brute_force_inside() {
            let square = unit_square();
            let inside = Point2D::new(0.5, 0.5);

            assert!(brute_force_contains_2d(&inside, &square));
        }

        #[test]
        fn test_brute_force_outside() {
            let square = unit_square();
            let outside = Point2D::new(2.0, 0.5);

            assert!(!brute_force_contains_2d(&outside, &square));
        }

        #[test]
        fn test_brute_force_outside_all_directions() {
            let square = unit_square();

            // Test points outside in various directions
            assert!(!brute_force_contains_2d(&Point2D::new(-1.0, 0.5), &square)); // left
            assert!(!brute_force_contains_2d(&Point2D::new(2.0, 0.5), &square)); // right
            assert!(!brute_force_contains_2d(&Point2D::new(0.5, -1.0), &square)); // below
            assert!(!brute_force_contains_2d(&Point2D::new(0.5, 2.0), &square)); // above
            assert!(!brute_force_contains_2d(&Point2D::new(-1.0, -1.0), &square)); // corner
        }

        #[test]
        fn test_brute_force_triangle() {
            let triangle = vec![
                Point2D::new(0.0, 0.0),
                Point2D::new(2.0, 0.0),
                Point2D::new(1.0, 2.0),
            ];

            // Center of triangle
            assert!(brute_force_contains_2d(&Point2D::new(1.0, 0.5), &triangle));

            // Outside
            assert!(!brute_force_contains_2d(&Point2D::new(0.0, 1.0), &triangle));
            assert!(!brute_force_contains_2d(&Point2D::new(2.0, 1.0), &triangle));
        }

        #[test]
        fn test_compute_origin_inside() {
            let square = unit_square();
            let (origin, inside) = compute_origin_inside_2d(&square);

            // Origin should be outside the bounding box
            assert!(origin.x < 0.0);
            assert!(origin.y < 0.0);

            // Origin should be outside the polygon
            assert!(!inside);
        }

        #[test]
        fn test_contains_with_origin() {
            let square = unit_square();
            let (origin, origin_inside) = compute_origin_inside_2d(&square);

            let inside = Point2D::new(0.5, 0.5);
            let outside = Point2D::new(2.0, 0.5);

            assert!(contains_with_origin_2d(
                &inside,
                &origin,
                origin_inside,
                &square
            ));
            assert!(!contains_with_origin_2d(
                &outside,
                &origin,
                origin_inside,
                &square
            ));
        }

        #[test]
        fn test_contains_with_edges() {
            let square = unit_square();
            let (origin, origin_inside) = compute_origin_inside_2d(&square);

            let point = Point2D::new(0.5, 0.5);

            // Test with all edges
            let all_edges: Vec<usize> = (0..4).collect();
            assert!(contains_with_edges_2d(
                &point,
                &origin,
                origin_inside,
                &square,
                &all_edges
            ));

            // Test with subset of edges (those actually crossed by the test ray)
            // For a point at (0.5, 0.5) and origin at (-1, -1), the ray
            // crosses edges 0 (bottom) and 3 (left)
            let subset = vec![0, 3];
            assert!(contains_with_edges_2d(
                &point,
                &origin,
                origin_inside,
                &square,
                &subset
            ));
        }

        #[test]
        fn test_contains_complex_polygon() {
            // L-shaped polygon
            let l_shape = vec![
                Point2D::new(0.0, 0.0),
                Point2D::new(2.0, 0.0),
                Point2D::new(2.0, 1.0),
                Point2D::new(1.0, 1.0),
                Point2D::new(1.0, 2.0),
                Point2D::new(0.0, 2.0),
            ];

            // Points inside the L
            assert!(brute_force_contains_2d(&Point2D::new(0.5, 0.5), &l_shape));
            assert!(brute_force_contains_2d(&Point2D::new(0.5, 1.5), &l_shape));
            assert!(brute_force_contains_2d(&Point2D::new(1.5, 0.5), &l_shape));

            // Point in the "cutout" of the L
            assert!(!brute_force_contains_2d(&Point2D::new(1.5, 1.5), &l_shape));

            // Points outside
            assert!(!brute_force_contains_2d(&Point2D::new(3.0, 0.5), &l_shape));
            assert!(!brute_force_contains_2d(&Point2D::new(-1.0, 0.5), &l_shape));
        }

        #[test]
        fn test_degenerate_cases() {
            // Empty polygon
            let empty: Vec<Point2D> = vec![];
            assert!(!brute_force_contains_2d(&Point2D::new(0.0, 0.0), &empty));

            // Single point
            let single = vec![Point2D::new(1.0, 1.0)];
            assert!(!brute_force_contains_2d(&Point2D::new(1.0, 1.0), &single));

            // Line segment (2 points)
            let line = vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 1.0)];
            assert!(!brute_force_contains_2d(&Point2D::new(0.5, 0.5), &line));
        }
    }
}
