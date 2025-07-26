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
