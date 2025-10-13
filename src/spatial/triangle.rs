//! Canonical triangle encoding for spatial indexing.
//!
//! Encodes triangles as 7-dimensional points for BKD tree storage, using a deterministic canonical
//! form that preserves boundary edge information while enabling efficient bounding box queries.

use robust::{orient2d, Coord};

const MINY_MINX_MAXY_MAXX_Y_X: i32 = 0;
const MINY_MINX_Y_X_MAXY_MAXX: i32 = 1;
const MAXY_MINX_Y_X_MINY_MAXX: i32 = 2;
const MAXY_MINX_MINY_MAXX_Y_X: i32 = 3;
const Y_MINX_MINY_X_MAXY_MAXX: i32 = 4;
const Y_MINX_MINY_MAXX_MAXY_X: i32 = 5;
const MAXY_MINX_MINY_X_Y_MAXX: i32 = 6;
const MINY_MINX_Y_MAXX_MAXY_X: i32 = 7;

/// A triangle encoded in canonical form for spatial indexing.
///
/// Contains the bounding box, one vertex, boundary edge flags, and a reconstruction code that
/// together allow exact triangle recovery while optimizing for spatial query performance.
pub struct Triangle {
    words: [i32; 6],
    boundaries: [bool; 3],
    code: i32,
}

impl Triangle {
    /// Encodes a triangle into canonical form.
    ///
    /// Takes three vertices as [x0, y0, x1, y1, x2, y2] and edge boundary flags [ab, bc, ca]
    /// indicating which edges are polygon boundaries. Returns a canonically encoded triangle with
    /// consistent vertex ordering (minX first, CCW orientation) for deterministic representation.
    pub fn encode(triangle: [i32; 6], boundary: [bool; 3]) -> Self {
        let mut ax = triangle[0];
        let mut ay = triangle[1];
        let mut bx = triangle[2];
        let mut by = triangle[3];
        let mut cx = triangle[4];
        let mut cy = triangle[5];
        let mut ab = boundary[0];
        let mut bc = boundary[1];
        let mut ca = boundary[2];
        // rotate edges and place minX at the beginning
        if bx < ax || cx < ax {
            let temp_x = ax;
            let temp_y = ay;
            let temp_boundary = ab;
            if bx < cx {
                ax = bx;
                ay = by;
                ab = bc;
                bx = cx;
                by = cy;
                bc = ca;
                cx = temp_x;
                cy = temp_y;
                ca = temp_boundary;
            } else {
                ax = cx;
                ay = cy;
                ab = ca;
                cx = bx;
                cy = by;
                ca = bc;
                bx = temp_x;
                by = temp_y;
                bc = temp_boundary;
            }
        } else if ax == bx && ax == cx {
            // degenerated case, all points with same longitude
            // we need to prevent that ax is in the middle (not part of the MBS)
            if by < ay || cy < ay {
                let temp_x = ax;
                let temp_y = ay;
                let temp_boundary = ab;
                if by < cy {
                    ax = bx;
                    ay = by;
                    ab = bc;
                    bx = cx;
                    by = cy;
                    bc = ca;
                    cx = temp_x;
                    cy = temp_y;
                    ca = temp_boundary;
                } else {
                    ax = cx;
                    ay = cy;
                    ab = ca;
                    cx = bx;
                    cy = by;
                    ca = bc;
                    bx = temp_x;
                    by = temp_y;
                    bc = temp_boundary;
                }
            }
        }
        // change orientation if CW
        if orient2d(
            Coord { x: ax, y: ay },
            Coord { x: bx, y: by },
            Coord { x: cx, y: cy },
        ) < 0.0
        {
            let temp_x = bx;
            let temp_y = by;
            let temp_boundary = ab;
            // ax and ay do not change, ab becomes bc
            ab = bc;
            bx = cx;
            by = cy;
            // bc does not change, ca becomes ab
            cx = temp_x;
            cy = temp_y;
            ca = temp_boundary;
        }
        let min_x = ax;
        let min_y = ay.min(by).min(cy);
        let max_x = ax.max(bx).max(cx);
        let max_y = ay.max(by).max(cy);
        let (y, x, code) = if min_y == ay {
            if max_y == by && max_x == bx {
                (cy, cx, MINY_MINX_MAXY_MAXX_Y_X)
            } else if max_y == cy && max_x == cx {
                (by, bx, MINY_MINX_Y_X_MAXY_MAXX)
            } else {
                (by, cx, MINY_MINX_Y_MAXX_MAXY_X)
            }
        } else if max_y == ay {
            if min_y == by && max_x == bx {
                (cy, cx, MAXY_MINX_MINY_MAXX_Y_X)
            } else if min_y == cy && max_x == cx {
                (by, bx, MAXY_MINX_Y_X_MINY_MAXX)
            } else {
                (cy, bx, MAXY_MINX_MINY_X_Y_MAXX)
            }
        } else if max_x == bx && min_y == by {
            (ay, cx, Y_MINX_MINY_MAXX_MAXY_X)
        } else if max_x == cx && max_y == cy {
            (bx, ay, Y_MINX_MINY_X_MAXY_MAXX)
        } else {
            panic!("Could not encode the provided triangle");
        };
        Triangle {
            words: [min_y, min_x, max_y, max_x, y, x],
            boundaries: [ab, bc, ca],
            code: code,
        }
    }

    /// Decodes the triangle back to vertex coordinates and boundary flags.
    ///
    /// Returns vertices as [x0, y0, x1, y1, x2, y2] in canonical CCW order and boundary flags [ab,
    /// bc, ca]. The vertex order may differ from the original input to encode() due to canonical
    /// rotation.
    pub fn decode(&self) -> ([i32; 6], [bool; 3]) {
        let (ay, ax, by, bx, cy, cx) = match self.code {
            MINY_MINX_MAXY_MAXX_Y_X => (
                self.words[0],
                self.words[1],
                self.words[2],
                self.words[3],
                self.words[4],
                self.words[5],
            ),
            MINY_MINX_Y_X_MAXY_MAXX => (
                self.words[0],
                self.words[1],
                self.words[4],
                self.words[5],
                self.words[2],
                self.words[3],
            ),
            MAXY_MINX_Y_X_MINY_MAXX => (
                self.words[2],
                self.words[1],
                self.words[4],
                self.words[5],
                self.words[0],
                self.words[3],
            ),
            MAXY_MINX_MINY_MAXX_Y_X => (
                self.words[2],
                self.words[1],
                self.words[0],
                self.words[3],
                self.words[4],
                self.words[5],
            ),
            Y_MINX_MINY_X_MAXY_MAXX => (
                self.words[4],
                self.words[1],
                self.words[0],
                self.words[5],
                self.words[2],
                self.words[3],
            ),
            Y_MINX_MINY_MAXX_MAXY_X => (
                self.words[4],
                self.words[1],
                self.words[0],
                self.words[3],
                self.words[2],
                self.words[5],
            ),
            MAXY_MINX_MINY_X_Y_MAXX => (
                self.words[2],
                self.words[1],
                self.words[0],
                self.words[5],
                self.words[4],
                self.words[3],
            ),
            MINY_MINX_Y_MAXX_MAXY_X => (
                self.words[0],
                self.words[1],
                self.words[4],
                self.words[3],
                self.words[2],
                self.words[5],
            ),
            _ => panic!("Could not decode the provided triangle"),
        };
        ([ax, ay, bx, by, cx, cy], self.boundaries)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_triangle() {
        let test_cases = [
            ([1, 1, 2, 3, 4, 2], [true, false, false]),
            ([1, 1, 4, 2, 2, 3], [false, false, true]),
            ([4, 2, 1, 1, 2, 3], [false, true, false]),
            ([4, 2, 2, 3, 1, 1], [false, true, false]),
            ([2, 3, 1, 1, 4, 2], [true, false, false]),
            ([2, 3, 4, 2, 1, 1], [false, false, true]),
        ];
        let ccw_coords = [1, 1, 4, 2, 2, 3];
        let ccw_bounds = [false, false, true];
        for (coords, bounds) in test_cases {
            let triangle = Triangle::encode(coords, bounds);
            println!("words: {:?}, code: {:?}", triangle.words, triangle.code);
            let (decoded_coords, decoded_bounds) = triangle.decode();
            assert_eq!(decoded_coords, ccw_coords);
            assert_eq!(decoded_bounds, ccw_bounds);
        }
    }

    #[test]
    fn degenerate_triangle() {
        let test_cases = [
            (
                [1, 1, 1, 2, 1, 3],
                [true, false, false],
                [1, 1, 1, 2, 1, 3],
                [true, false, false],
            ),
            (
                [1, 2, 1, 1, 1, 3],
                [true, false, false],
                [1, 1, 1, 3, 1, 2],
                [false, false, true],
            ),
            (
                [1, 2, 1, 3, 1, 1],
                [false, false, true],
                [1, 1, 1, 2, 1, 3],
                [true, false, false],
            ),
        ];
        for (coords, bounds, ccw_coords, ccw_bounds) in test_cases {
            let triangle = Triangle::encode(coords, bounds);
            println!("words: {:?}, code: {:?}", triangle.words, triangle.code);
            let (decoded_coords, decoded_bounds) = triangle.decode();
            assert_eq!(decoded_coords, ccw_coords);
            assert_eq!(decoded_bounds, ccw_bounds);
        }
    }

    #[test]
    fn decode_triangle() {
        let test_cases = [
            [1, 1, 8, 6, 4, 4],
            [0, 0, 2, 1, 3, 3],
            [0, 3, 1, 1, 2, 0],
            [0, 3, 2, 0, 1, 2],
            [0, 2, 2, 0, 3, 3],
            [0, 2, 3, 0, 1, 3],
            [0, 3, 1, 0, 2, 1],
            [0, 0, 2, 1, 1, 2],
        ];
        for coords in test_cases {
            let triangle = Triangle::encode(coords, [true, true, true]);
            println!("words: {:?}, code: {:?}", triangle.words, triangle.code);
            let (decoded_coords, _) = triangle.decode();
            assert_eq!(decoded_coords, coords);
        }
    }
}
