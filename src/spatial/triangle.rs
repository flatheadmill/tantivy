//! Triangle spatial indexing with BKD encoding.
//!
//! This module provides triangle encoding and decoding following Lucene's BKD approach
//! for spatial indexing. Triangles are encoded as 7-dimensional points with 8 canonical
//! orderings for deterministic storage and efficient spatial queries.

use std::cmp;

// Constants for triangle encoding following Lucene's BKD approach
const BYTES: usize = 4; // Integer bytes (i32 = 4 bytes)
const TOTAL_BYTES: usize = 7 * BYTES; // 7 dimensions * 4 bytes each = 28 bytes total

// Triangle encoding constants - these represent the 8 possible canonical forms
// Based on which vertex has minX/minY and maxX/maxY properties
// Each constant encodes a specific geometric configuration for deterministic encoding
const MINY_MINX_MAXY_MAXX_Y_X: i32 = 0; // Point A has both minX and minY
const MINY_MINX_Y_X_MAXY_MAXX: i32 = 1; // Point A has minX/minY, C has maxX/maxY
const MAXY_MINX_Y_X_MINY_MAXX: i32 = 2; // Point A has minX/maxY, C has maxX/minY
const MAXY_MINX_MINY_MAXX_Y_X: i32 = 3; // Point A has minX/maxY, B has maxX/minY
const Y_MINX_MINY_X_MAXY_MAXX: i32 = 4; // Point A has minX (middle Y), B has minY, C has maxX/maxY
const Y_MINX_MINY_MAXX_MAXY_X: i32 = 5; // Point A has minX (middle Y), B has minY/maxX
const MAXY_MINX_MINY_X_Y_MAXX: i32 = 6; // Point A has minX/maxY, B has minY
const MINY_MINX_Y_MAXX_MAXY_X: i32 = 7; // Point A has minX/minY, B has maxX

/// Triangle representation with edge flags for polygon reconstruction
///
/// This follows Lucene's approach where triangles from tessellated polygons
/// retain information about which edges belonged to the original polygon boundary.
/// This is crucial for accurate point-in-polygon queries.
/// Triangle representation with edge flags for polygon reconstruction.
///
/// This follows Lucene's approach where triangles from tessellated polygons
/// retain information about which edges belonged to the original polygon boundary.
/// This is crucial for accurate point-in-polygon queries.
#[derive(Debug, Clone, PartialEq)]
pub struct Triangle {
    pub a_x: i32,
    pub a_y: i32, // Vertex A coordinates
    pub b_x: i32,
    pub b_y: i32, // Vertex B coordinates
    pub c_x: i32,
    pub c_y: i32, // Vertex C coordinates
    pub ab: bool, // Edge AB belongs to original polygon boundary
    pub bc: bool, // Edge BC belongs to original polygon boundary
    pub ca: bool, // Edge CA belongs to original polygon boundary
}

impl Triangle {
    /// Create a new triangle with vertex coordinates and edge flags
    ///
    /// # Arguments
    /// * Vertex coordinates as i32 (following Lucene's integer coordinate system)
    /// * Edge flags indicating which edges are from the original polygon
    pub fn new(
        a_x: i32,
        a_y: i32,
        b_x: i32,
        b_y: i32,
        c_x: i32,
        c_y: i32,
        ab: bool,
        bc: bool,
        ca: bool,
    ) -> Self {
        Triangle {
            a_x,
            a_y,
            b_x,
            b_y,
            c_x,
            c_y,
            ab,
            bc,
            ca,
        }
    }

    /// Calculate the bounding box of the triangle
    ///
    /// Returns (min_x, min_y, max_x, max_y) following Lucene's approach where
    /// these values become the first 4 dimensions of the BKD encoding.
    pub fn bounding_box(&self) -> (i32, i32, i32, i32) {
        let min_x = cmp::min(self.a_x, cmp::min(self.b_x, self.c_x));
        let min_y = cmp::min(self.a_y, cmp::min(self.b_y, self.c_y));
        let max_x = cmp::max(self.a_x, cmp::max(self.b_x, self.c_x));
        let max_y = cmp::max(self.a_y, cmp::max(self.b_y, self.c_y));
        (min_x, min_y, max_x, max_y)
    }

    /// Test if a point is inside this triangle using Lucene's winding order method
    ///
    /// This follows Lucene's Component2D.pointInTriangle implementation using
    /// the orient() function to determine if a point is inside the triangle.
    /// The method first checks the bounding box, then uses orientation tests.
    pub fn contains_point(&self, x: i32, y: i32) -> bool {
        let (min_x, min_y, max_x, max_y) = self.bounding_box();

        // Check bounding box first - Lucene optimization for degenerate cases
        if x < min_x || x > max_x || y < min_y || y > max_y {
            return false;
        }

        // Use winding order method with orientation tests
        let a = orient(x, y, self.a_x, self.a_y, self.b_x, self.b_y);
        let b = orient(x, y, self.b_x, self.b_y, self.c_x, self.c_y);

        if a == 0 || b == 0 || (a < 0) == (b < 0) {
            let c = orient(x, y, self.c_x, self.c_y, self.a_x, self.a_y);
            c == 0 || (c < 0) == (b < 0 || a < 0)
        } else {
            false
        }
    }

    /// Test if this triangle intersects with another triangle
    ///
    /// Based on Lucene's intersectsTriangle approach: checks if any vertices
    /// of either triangle are contained in the other, or if any edges intersect.
    pub fn intersects_triangle(&self, other: &Triangle) -> bool {
        let (min_x, min_y, max_x, max_y) = self.bounding_box();
        let (other_min_x, other_min_y, other_max_x, other_max_y) = other.bounding_box();

        // Quick bounding box disjoint test
        if max_x < other_min_x || min_x > other_max_x || max_y < other_min_y || min_y > other_max_y
        {
            return false;
        }

        // Check if any vertex of other triangle is inside this triangle
        if self.contains_point(other.a_x, other.a_y)
            || self.contains_point(other.b_x, other.b_y)
            || self.contains_point(other.c_x, other.c_y)
        {
            return true;
        }

        // Check if any vertex of this triangle is inside other triangle
        if other.contains_point(self.a_x, self.a_y)
            || other.contains_point(self.b_x, self.b_y)
            || other.contains_point(self.c_x, self.c_y)
        {
            return true;
        }

        // Check if any edges intersect
        self.edges_intersect_triangle(other)
    }

    /// Helper method to check if any edges of this triangle intersect with another triangle
    fn edges_intersect_triangle(&self, other: &Triangle) -> bool {
        // Check if any edge of this triangle intersects any edge of the other triangle
        self.edge_intersects_edge(
            self.a_x, self.a_y, self.b_x, self.b_y, other.a_x, other.a_y, other.b_x, other.b_y,
        ) || self.edge_intersects_edge(
            self.a_x, self.a_y, self.b_x, self.b_y, other.b_x, other.b_y, other.c_x, other.c_y,
        ) || self.edge_intersects_edge(
            self.a_x, self.a_y, self.b_x, self.b_y, other.c_x, other.c_y, other.a_x, other.a_y,
        ) || self.edge_intersects_edge(
            self.b_x, self.b_y, self.c_x, self.c_y, other.a_x, other.a_y, other.b_x, other.b_y,
        ) || self.edge_intersects_edge(
            self.b_x, self.b_y, self.c_x, self.c_y, other.b_x, other.b_y, other.c_x, other.c_y,
        ) || self.edge_intersects_edge(
            self.b_x, self.b_y, self.c_x, self.c_y, other.c_x, other.c_y, other.a_x, other.a_y,
        ) || self.edge_intersects_edge(
            self.c_x, self.c_y, self.a_x, self.a_y, other.a_x, other.a_y, other.b_x, other.b_y,
        ) || self.edge_intersects_edge(
            self.c_x, self.c_y, self.a_x, self.a_y, other.b_x, other.b_y, other.c_x, other.c_y,
        ) || self.edge_intersects_edge(
            self.c_x, self.c_y, self.a_x, self.a_y, other.c_x, other.c_y, other.a_x, other.a_y,
        )
    }

    /// Test if two line segments intersect using orientation tests
    fn edge_intersects_edge(
        &self,
        ax: i32,
        ay: i32,
        bx: i32,
        by: i32,
        cx: i32,
        cy: i32,
        dx: i32,
        dy: i32,
    ) -> bool {
        let o1 = orient(ax, ay, bx, by, cx, cy);
        let o2 = orient(ax, ay, bx, by, dx, dy);
        let o3 = orient(cx, cy, dx, dy, ax, ay);
        let o4 = orient(cx, cy, dx, dy, bx, by);

        // General case - segments intersect if orientations are different
        if o1 != o2 && o3 != o4 {
            return true;
        }

        // Special cases - points are collinear and lie on the other segment
        (o1 == 0 && self.point_on_segment(ax, ay, cx, cy, bx, by))
            || (o2 == 0 && self.point_on_segment(ax, ay, dx, dy, bx, by))
            || (o3 == 0 && self.point_on_segment(cx, cy, ax, ay, dx, dy))
            || (o4 == 0 && self.point_on_segment(cx, cy, bx, by, dx, dy))
    }

    /// Check if point q lies on line segment pr (when they are collinear)
    fn point_on_segment(&self, px: i32, py: i32, qx: i32, qy: i32, rx: i32, ry: i32) -> bool {
        qx <= cmp::max(px, rx)
            && qx >= cmp::min(px, rx)
            && qy <= cmp::max(py, ry)
            && qy >= cmp::min(py, ry)
    }

    /// Test if this triangle intersects with a bounding box
    ///
    /// This follows Lucene's approach for spatial indexing where the first
    /// 4 dimensions of encoded triangles represent bounding boxes.
    pub fn intersects_bbox(
        &self,
        bbox_min_x: i32,
        bbox_min_y: i32,
        bbox_max_x: i32,
        bbox_max_y: i32,
    ) -> bool {
        let (min_x, min_y, max_x, max_y) = self.bounding_box();

        // Quick disjoint test
        if max_x < bbox_min_x || min_x > bbox_max_x || max_y < bbox_min_y || min_y > bbox_max_y {
            return false;
        }

        // Check if any triangle vertex is inside the bbox
        if (self.a_x >= bbox_min_x
            && self.a_x <= bbox_max_x
            && self.a_y >= bbox_min_y
            && self.a_y <= bbox_max_y)
            || (self.b_x >= bbox_min_x
                && self.b_x <= bbox_max_x
                && self.b_y >= bbox_min_y
                && self.b_y <= bbox_max_y)
            || (self.c_x >= bbox_min_x
                && self.c_x <= bbox_max_x
                && self.c_y >= bbox_min_y
                && self.c_y <= bbox_max_y)
        {
            return true;
        }

        // Check if any bbox corner is inside the triangle
        if self.contains_point(bbox_min_x, bbox_min_y)
            || self.contains_point(bbox_max_x, bbox_min_y)
            || self.contains_point(bbox_max_x, bbox_max_y)
            || self.contains_point(bbox_min_x, bbox_max_y)
        {
            return true;
        }

        // Check if triangle edges intersect bbox edges
        self.edge_intersects_bbox_edge(bbox_min_x, bbox_min_y, bbox_max_x, bbox_max_y)
    }

    /// Helper to check if triangle edges intersect bounding box edges
    fn edge_intersects_bbox_edge(
        &self,
        bbox_min_x: i32,
        bbox_min_y: i32,
        bbox_max_x: i32,
        bbox_max_y: i32,
    ) -> bool {
        // Check each triangle edge against each bbox edge
        self.edge_intersects_edge(
            self.a_x, self.a_y, self.b_x, self.b_y, bbox_min_x, bbox_min_y, bbox_max_x, bbox_min_y,
        ) || self.edge_intersects_edge(
            self.a_x, self.a_y, self.b_x, self.b_y, bbox_max_x, bbox_min_y, bbox_max_x, bbox_max_y,
        ) || self.edge_intersects_edge(
            self.a_x, self.a_y, self.b_x, self.b_y, bbox_max_x, bbox_max_y, bbox_min_x, bbox_max_y,
        ) || self.edge_intersects_edge(
            self.a_x, self.a_y, self.b_x, self.b_y, bbox_min_x, bbox_max_y, bbox_min_x, bbox_min_y,
        ) || self.edge_intersects_edge(
            self.b_x, self.b_y, self.c_x, self.c_y, bbox_min_x, bbox_min_y, bbox_max_x, bbox_min_y,
        ) || self.edge_intersects_edge(
            self.b_x, self.b_y, self.c_x, self.c_y, bbox_max_x, bbox_min_y, bbox_max_x, bbox_max_y,
        ) || self.edge_intersects_edge(
            self.b_x, self.b_y, self.c_x, self.c_y, bbox_max_x, bbox_max_y, bbox_min_x, bbox_max_y,
        ) || self.edge_intersects_edge(
            self.b_x, self.b_y, self.c_x, self.c_y, bbox_min_x, bbox_max_y, bbox_min_x, bbox_min_y,
        ) || self.edge_intersects_edge(
            self.c_x, self.c_y, self.a_x, self.a_y, bbox_min_x, bbox_min_y, bbox_max_x, bbox_min_y,
        ) || self.edge_intersects_edge(
            self.c_x, self.c_y, self.a_x, self.a_y, bbox_max_x, bbox_min_y, bbox_max_x, bbox_max_y,
        ) || self.edge_intersects_edge(
            self.c_x, self.c_y, self.a_x, self.a_y, bbox_max_x, bbox_max_y, bbox_min_x, bbox_max_y,
        ) || self.edge_intersects_edge(
            self.c_x, self.c_y, self.a_x, self.a_y, bbox_min_x, bbox_max_y, bbox_min_x, bbox_min_y,
        )
    }
}

/// Convert signed integer to sortable byte representation
///
/// This is equivalent to Lucene's NumericUtils.intToSortableBytes()
/// The key insight is flipping the sign bit so that:
/// - Negative numbers sort before positive numbers
/// - Within each group, natural ordering is preserved
/// This enables direct byte-wise comparison for spatial indexing
fn int_to_sortable_bytes(value: i32, bytes: &mut [u8], offset: usize) {
    let sortable = (value as u32) ^ 0x80000000; // Flip sign bit for proper sorting
    bytes[offset..offset + 4].copy_from_slice(&sortable.to_be_bytes());
}

/// Convert sortable bytes back to signed integer
///
/// Reverse operation of int_to_sortable_bytes() - flips the sign bit back
/// to recover the original signed integer value
fn sortable_bytes_to_int(bytes: &[u8], offset: usize) -> i32 {
    let mut array = [0u8; 4];
    array.copy_from_slice(&bytes[offset..offset + 4]);
    let sortable = u32::from_be_bytes(array);
    (sortable ^ 0x80000000) as i32 // Flip sign bit back to original
}

/// Geometric orientation test using cross product
///
/// Determines if three points form a counter-clockwise (CCW), clockwise (CW),
/// or collinear arrangement. This is fundamental for ensuring consistent
/// triangle orientation in the canonical form.
///
/// Returns: 1 for CCW, -1 for CW, 0 for collinear
/// Uses i64 arithmetic to prevent overflow during cross product calculation
fn orient(ax: i32, ay: i32, bx: i32, by: i32, cx: i32, cy: i32) -> i32 {
    // Cross product: (B-A) × (C-A) = (bx-ax)(cy-ay) - (by-ay)(cx-ax)
    let det = (bx as i64 - ax as i64) * (cy as i64 - ay as i64)
        - (by as i64 - ay as i64) * (cx as i64 - ax as i64);
    if det > 0 {
        1
    } else if det < 0 {
        -1
    } else {
        0
    }
}

/// Encode triangle into 7-dimensional point for BKD tree indexing
///
/// This is the core of Lucene's triangle encoding algorithm. The process:
/// 1. Rotate vertices to place the one with minimum X coordinate first
/// 2. Ensure counter-clockwise orientation for consistency
/// 3. Determine canonical encoding pattern based on bounding box properties
/// 4. Encode as 7D point: [minY, minX, maxY, maxX, y, x, metadata]
///
/// The first 4 dimensions form a 2D bounding box for efficient spatial indexing,
/// while dimensions 5-6 store the "floating" vertex, and dimension 7 stores
/// encoding metadata and edge flags.
pub fn encode_triangle(triangle: &Triangle) -> [u8; TOTAL_BYTES] {
    let mut bytes = [0u8; TOTAL_BYTES];

    // Work with mutable copies for canonical form transformation
    let mut a_x = triangle.a_x;
    let mut a_y = triangle.a_y;
    let mut b_x = triangle.b_x;
    let mut b_y = triangle.b_y;
    let mut c_x = triangle.c_x;
    let mut c_y = triangle.c_y;
    let mut ab = triangle.ab;
    let mut bc = triangle.bc;
    let mut ca = triangle.ca;

    // STEP 1: Rotate triangle to place vertex with minimum X coordinate first
    // This ensures deterministic encoding regardless of input vertex order
    if b_x < a_x || c_x < a_x {
        if b_x < c_x {
            // Rotate: a -> b -> c -> a (B becomes new A)
            let temp_x = a_x;
            let temp_y = a_y;
            let temp_bool = ab;
            a_x = b_x;
            a_y = b_y;
            ab = bc;
            b_x = c_x;
            b_y = c_y;
            bc = ca;
            c_x = temp_x;
            c_y = temp_y;
            ca = temp_bool;
        } else {
            // Rotate: a -> c -> b -> a (C becomes new A)
            let temp_x = a_x;
            let temp_y = a_y;
            let temp_bool = ab;
            a_x = c_x;
            a_y = c_y;
            ab = ca;
            c_x = b_x;
            c_y = b_y;
            ca = bc;
            b_x = temp_x;
            b_y = temp_y;
            bc = temp_bool;
        }
    } else if a_x == b_x && a_x == c_x {
        // Degenerate case: All points have same X coordinate
        // Use Y coordinate as tie-breaker to maintain deterministic ordering
        if b_y < a_y || c_y < a_y {
            if b_y < c_y {
                let temp_x = a_x;
                let temp_y = a_y;
                let temp_bool = ab;
                a_x = b_x;
                a_y = b_y;
                ab = bc;
                b_x = c_x;
                b_y = c_y;
                bc = ca;
                c_x = temp_x;
                c_y = temp_y;
                ca = temp_bool;
            } else {
                let temp_x = a_x;
                let temp_y = a_y;
                let temp_bool = ab;
                a_x = c_x;
                a_y = c_y;
                ab = ca;
                c_x = b_x;
                c_y = b_y;
                ca = bc;
                b_x = temp_x;
                b_y = temp_y;
                bc = temp_bool;
            }
        }
    }

    // STEP 2: Ensure counter-clockwise orientation for canonical form
    // This eliminates ambiguity between CW and CCW representations of the same triangle
    if orient(a_x, a_y, b_x, b_y, c_x, c_y) == -1 {
        // Triangle is clockwise - swap B and C to make it counter-clockwise
        let temp_x = b_x;
        let temp_y = b_y;
        let temp_bool = ab;
        ab = bc; // Edge AB becomes AC (since we're swapping B and C)
        b_x = c_x;
        b_y = c_y;
        c_x = temp_x;
        c_y = temp_y;
        ca = temp_bool;
    }

    // STEP 3: Calculate bounding box coordinates
    let min_x = a_x; // A has minimum X by construction
    let min_y = cmp::min(a_y, cmp::min(b_y, c_y));
    let max_x = cmp::max(a_x, cmp::max(b_x, c_x));
    let max_y = cmp::max(a_y, cmp::max(b_y, c_y));

    // STEP 4: Determine encoding pattern based on bounding box vertex properties
    // This is where Lucene's algorithm gets clever - it uses the geometric properties
    // of which vertex has min/max X/Y to determine how to encode the triangle
    let (bits, x, y) = if min_y == a_y {
        // Point A has minimum Y coordinate
        if max_y == b_y && max_x == b_x {
            // B has both maxY and maxX - C is the "floating" point
            (MINY_MINX_MAXY_MAXX_Y_X, c_x, c_y)
        } else if max_y == c_y && max_x == c_x {
            // C has both maxY and maxX - B is the "floating" point
            (MINY_MINX_Y_X_MAXY_MAXX, b_x, b_y)
        } else {
            // Mixed case - store C's X and B's Y as floating coordinates
            (MINY_MINX_Y_MAXX_MAXY_X, c_x, b_y)
        }
    } else if max_y == a_y {
        // Point A has maximum Y coordinate
        if min_y == b_y && max_x == b_x {
            (MAXY_MINX_MINY_MAXX_Y_X, c_x, c_y)
        } else if min_y == c_y && max_x == c_x {
            (MAXY_MINX_Y_X_MINY_MAXX, b_x, b_y)
        } else {
            (MAXY_MINX_MINY_X_Y_MAXX, b_x, c_y)
        }
    } else if max_x == b_x && min_y == b_y {
        // B has maxX and minY - A's Y coordinate becomes floating
        (Y_MINX_MINY_MAXX_MAXY_X, c_x, a_y)
    } else if max_x == c_x && max_y == c_y {
        // C has maxX and maxY - A's Y coordinate becomes floating
        (Y_MINX_MINY_X_MAXY_MAXX, b_x, a_y)
    } else {
        panic!("Could not encode the provided triangle - invalid geometric configuration");
    };

    // STEP 5: Pack edge flags into the encoding bits
    // Bits 0-2: encoding pattern, Bits 3-5: edge flags
    let mut bits = bits;
    bits |= if ab { 1 << 3 } else { 0 }; // Bit 3: edge AB
    bits |= if bc { 1 << 4 } else { 0 }; // Bit 4: edge BC
    bits |= if ca { 1 << 5 } else { 0 }; // Bit 5: edge CA

    // STEP 6: Encode as 7-dimensional point with sortable bytes
    // Dimensions 0-3: 2D bounding box (minY, minX, maxY, maxX)
    // Dimensions 4-5: floating coordinates (y, x)
    // Dimension 6: encoding metadata and edge flags
    int_to_sortable_bytes(min_y, &mut bytes, 0 * BYTES);
    int_to_sortable_bytes(min_x, &mut bytes, 1 * BYTES);
    int_to_sortable_bytes(max_y, &mut bytes, 2 * BYTES);
    int_to_sortable_bytes(max_x, &mut bytes, 3 * BYTES);
    int_to_sortable_bytes(y, &mut bytes, 4 * BYTES);
    int_to_sortable_bytes(x, &mut bytes, 5 * BYTES);
    int_to_sortable_bytes(bits, &mut bytes, 6 * BYTES);

    bytes
}

/// Decode 7-dimensional point back to triangle representation
///
/// Reverse process of encode_triangle() - reconstructs the original triangle
/// from its encoded 7D representation by:
/// 1. Extracting encoding pattern and edge flags from metadata
/// 2. Using pattern to determine vertex coordinate assignments
/// 3. Reconstructing triangle with proper orientation validation
pub fn decode_triangle(bytes: &[u8; TOTAL_BYTES]) -> Triangle {
    // Extract metadata from dimension 6
    let bits = sortable_bytes_to_int(bytes, 6 * BYTES);
    let t_code = bits & ((1 << 3) - 1); // Extract encoding pattern (bits 0-2)

    // Reconstruct triangle vertices based on encoding pattern
    // Each pattern corresponds to a specific geometric configuration
    let (a_x, a_y, b_x, b_y, c_x, c_y) = match t_code {
        MINY_MINX_MAXY_MAXX_Y_X => (
            sortable_bytes_to_int(bytes, 1 * BYTES), // aX = minX
            sortable_bytes_to_int(bytes, 0 * BYTES), // aY = minY
            sortable_bytes_to_int(bytes, 3 * BYTES), // bX = maxX
            sortable_bytes_to_int(bytes, 2 * BYTES), // bY = maxY
            sortable_bytes_to_int(bytes, 5 * BYTES), // cX = x (floating)
            sortable_bytes_to_int(bytes, 4 * BYTES), // cY = y (floating)
        ),
        MINY_MINX_Y_X_MAXY_MAXX => (
            sortable_bytes_to_int(bytes, 1 * BYTES), // aX = minX
            sortable_bytes_to_int(bytes, 0 * BYTES), // aY = minY
            sortable_bytes_to_int(bytes, 5 * BYTES), // bX = x (floating)
            sortable_bytes_to_int(bytes, 4 * BYTES), // bY = y (floating)
            sortable_bytes_to_int(bytes, 3 * BYTES), // cX = maxX
            sortable_bytes_to_int(bytes, 2 * BYTES), // cY = maxY
        ),
        MAXY_MINX_Y_X_MINY_MAXX => (
            sortable_bytes_to_int(bytes, 1 * BYTES), // aX = minX
            sortable_bytes_to_int(bytes, 2 * BYTES), // aY = maxY
            sortable_bytes_to_int(bytes, 5 * BYTES), // bX = x (floating)
            sortable_bytes_to_int(bytes, 4 * BYTES), // bY = y (floating)
            sortable_bytes_to_int(bytes, 3 * BYTES), // cX = maxX
            sortable_bytes_to_int(bytes, 0 * BYTES), // cY = minY
        ),
        MAXY_MINX_MINY_MAXX_Y_X => (
            sortable_bytes_to_int(bytes, 1 * BYTES), // aX = minX
            sortable_bytes_to_int(bytes, 2 * BYTES), // aY = maxY
            sortable_bytes_to_int(bytes, 3 * BYTES), // bX = maxX
            sortable_bytes_to_int(bytes, 0 * BYTES), // bY = minY
            sortable_bytes_to_int(bytes, 5 * BYTES), // cX = x (floating)
            sortable_bytes_to_int(bytes, 4 * BYTES), // cY = y (floating)
        ),
        Y_MINX_MINY_X_MAXY_MAXX => (
            sortable_bytes_to_int(bytes, 1 * BYTES), // aX = minX
            sortable_bytes_to_int(bytes, 4 * BYTES), // aY = y (floating)
            sortable_bytes_to_int(bytes, 5 * BYTES), // bX = x (floating)
            sortable_bytes_to_int(bytes, 0 * BYTES), // bY = minY
            sortable_bytes_to_int(bytes, 3 * BYTES), // cX = maxX
            sortable_bytes_to_int(bytes, 2 * BYTES), // cY = maxY
        ),
        Y_MINX_MINY_MAXX_MAXY_X => (
            sortable_bytes_to_int(bytes, 1 * BYTES), // aX = minX
            sortable_bytes_to_int(bytes, 4 * BYTES), // aY = y (floating)
            sortable_bytes_to_int(bytes, 3 * BYTES), // bX = maxX
            sortable_bytes_to_int(bytes, 0 * BYTES), // bY = minY
            sortable_bytes_to_int(bytes, 5 * BYTES), // cX = x (floating)
            sortable_bytes_to_int(bytes, 2 * BYTES), // cY = maxY
        ),
        MAXY_MINX_MINY_X_Y_MAXX => (
            sortable_bytes_to_int(bytes, 1 * BYTES), // aX = minX
            sortable_bytes_to_int(bytes, 2 * BYTES), // aY = maxY
            sortable_bytes_to_int(bytes, 5 * BYTES), // bX = x (floating)
            sortable_bytes_to_int(bytes, 0 * BYTES), // bY = minY
            sortable_bytes_to_int(bytes, 3 * BYTES), // cX = maxX
            sortable_bytes_to_int(bytes, 4 * BYTES), // cY = y (floating)
        ),
        MINY_MINX_Y_MAXX_MAXY_X => (
            sortable_bytes_to_int(bytes, 1 * BYTES), // aX = minX
            sortable_bytes_to_int(bytes, 0 * BYTES), // aY = minY
            sortable_bytes_to_int(bytes, 3 * BYTES), // bX = maxX
            sortable_bytes_to_int(bytes, 4 * BYTES), // bY = y (floating)
            sortable_bytes_to_int(bytes, 5 * BYTES), // cX = x (floating)
            sortable_bytes_to_int(bytes, 2 * BYTES), // cY = maxY
        ),
        _ => panic!(
            "Could not decode the provided triangle - invalid encoding pattern: {}",
            t_code
        ),
    };

    // Extract edge flags from metadata bits
    let ab = (bits & (1 << 3)) != 0; // Bit 3: edge AB flag
    let bc = (bits & (1 << 4)) != 0; // Bit 4: edge BC flag
    let ca = (bits & (1 << 5)) != 0; // Bit 5: edge CA flag

    // Verify that decoded triangle has proper orientation (CCW or collinear)
    // This serves as a sanity check for the decoding process
    assert!(
        orient(a_x, a_y, b_x, b_y, c_x, c_y) >= 0,
        "Decoded triangle has incorrect orientation - should be CCW or collinear"
    );

    Triangle::new(a_x, a_y, b_x, b_y, c_x, c_y, ab, bc, ca)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_triangle_bounding_box() {
        let triangle = Triangle::new(10, 5, 20, 15, 5, 20, true, false, true);
        let (min_x, min_y, max_x, max_y) = triangle.bounding_box();

        assert_eq!(min_x, 5); // min of 10, 20, 5
        assert_eq!(min_y, 5); // min of 5, 15, 20
        assert_eq!(max_x, 20); // max of 10, 20, 5
        assert_eq!(max_y, 20); // max of 5, 15, 20
    }

    #[test]
    fn test_point_in_triangle() {
        // Simple triangle: (0,0), (10,0), (5,10)
        let triangle = Triangle::new(0, 0, 10, 0, 5, 10, true, true, true);

        // Point inside triangle
        assert!(triangle.contains_point(5, 3));

        // Point outside triangle
        assert!(!triangle.contains_point(15, 5));
        assert!(!triangle.contains_point(5, 15));

        // Point on vertex
        assert!(triangle.contains_point(0, 0));
        assert!(triangle.contains_point(10, 0));
        assert!(triangle.contains_point(5, 10));

        // Point outside bounding box
        assert!(!triangle.contains_point(-5, 0));
        assert!(!triangle.contains_point(0, -5));
    }

    #[test]
    fn test_triangle_intersects_bbox() {
        let triangle = Triangle::new(5, 5, 15, 5, 10, 15, true, true, true);

        // Bbox completely contains triangle
        assert!(triangle.intersects_bbox(0, 0, 20, 20));

        // Bbox partially overlaps triangle
        assert!(triangle.intersects_bbox(12, 7, 18, 12));

        // Bbox disjoint from triangle
        assert!(!triangle.intersects_bbox(20, 20, 30, 30));
        assert!(!triangle.intersects_bbox(0, 0, 3, 3));
    }

    #[test]
    fn test_encode_decode_roundtrip() {
        let original = Triangle::new(100, -50, 200, 300, -100, 150, true, false, true);
        let encoded = encode_triangle(&original);
        let decoded = decode_triangle(&encoded);

        // The decoded triangle should have the same bounding box
        assert_eq!(original.bounding_box(), decoded.bounding_box());

        // The decoded triangle should contain the same vertices (in any order)
        let orig_vertices = [
            (original.a_x, original.a_y),
            (original.b_x, original.b_y),
            (original.c_x, original.c_y),
        ];
        let decoded_vertices = [
            (decoded.a_x, decoded.a_y),
            (decoded.b_x, decoded.b_y),
            (decoded.c_x, decoded.c_y),
        ];

        for vertex in orig_vertices.iter() {
            assert!(
                decoded_vertices.contains(vertex),
                "Decoded triangle missing vertex {:?}",
                vertex
            );
        }

        // Both triangles should have the same orientation (CCW or collinear)
        let orig_orient = orient(
            original.a_x,
            original.a_y,
            original.b_x,
            original.b_y,
            original.c_x,
            original.c_y,
        );
        let decoded_orient = orient(
            decoded.a_x,
            decoded.a_y,
            decoded.b_x,
            decoded.b_y,
            decoded.c_x,
            decoded.c_y,
        );
        assert!(
            orig_orient >= 0 && decoded_orient >= 0,
            "Both triangles should have CCW or collinear orientation"
        );
    }

    #[test]
    fn test_encode_decode_with_bounding_box() {
        let triangle = Triangle::new(10, 20, 30, 5, 15, 25, false, true, false);
        let encoded = encode_triangle(&triangle);

        // First 4 dimensions should be the bounding box
        let min_y = sortable_bytes_to_int(&encoded, 0 * BYTES);
        let min_x = sortable_bytes_to_int(&encoded, 1 * BYTES);
        let max_y = sortable_bytes_to_int(&encoded, 2 * BYTES);
        let max_x = sortable_bytes_to_int(&encoded, 3 * BYTES);

        let (expected_min_x, expected_min_y, expected_max_x, expected_max_y) =
            triangle.bounding_box();

        assert_eq!(min_x, expected_min_x);
        assert_eq!(min_y, expected_min_y);
        assert_eq!(max_x, expected_max_x);
        assert_eq!(max_y, expected_max_y);
    }

    #[test]
    fn test_triangles_intersect() {
        // Two overlapping triangles
        let triangle1 = Triangle::new(0, 0, 10, 0, 5, 10, true, true, true);
        let triangle2 = Triangle::new(3, 2, 8, 2, 6, 8, true, true, true);

        assert!(triangle1.intersects_triangle(&triangle2));
        assert!(triangle2.intersects_triangle(&triangle1));

        // Two disjoint triangles
        let triangle3 = Triangle::new(20, 20, 30, 20, 25, 30, true, true, true);
        assert!(!triangle1.intersects_triangle(&triangle3));
        assert!(!triangle3.intersects_triangle(&triangle1));
    }

    #[test]
    fn test_orientation_consistency() {
        // Test that our orient function matches Lucene's behavior
        assert_eq!(orient(0, 0, 1, 0, 0, 1), 1); // CCW
        assert_eq!(orient(0, 0, 0, 1, 1, 0), -1); // CW
        assert_eq!(orient(0, 0, 1, 1, 2, 2), 0); // Collinear
    }
}
