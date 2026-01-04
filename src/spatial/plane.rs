//! HUSH

/// Robust orientation test (cross product sign)
#[inline]
fn orientation(a: &(f64, f64), b: &(f64, f64), c: &(f64, f64)) -> f64 {
    robust::orient2d(
        robust::Coord { x: a.0, y: a.1 },
        robust::Coord { x: b.0, y: b.1 },
        robust::Coord { x: c.0, y: c.1 },
    )
}

/// Cohen-Sutherland outcode for line clipping
const INSIDE: u8 = 0;
const LEFT: u8 = 1;
const RIGHT: u8 = 2;
const BOTTOM: u8 = 4;
const TOP: u8 = 8;

fn outcode(p: &(f64, f64), aabb: &[(f64, f64); 2]) -> u8 {
    let (x, y) = p;
    let [(min_x, min_y), (max_x, max_y)] = aabb;

    let mut code = INSIDE;
    if x < min_x {
        code |= LEFT;
    } else if x > max_x {
        code |= RIGHT;
    }
    if y < min_y {
        code |= BOTTOM;
    } else if y > max_y {
        code |= TOP;
    }
    code
}

#[inline]
fn on_segment(p: &(f64, f64), q: &(f64, f64), r: &(f64, f64)) -> bool {
    q.0 <= p.0.max(r.0) && q.0 >= p.0.min(r.0) && q.1 <= p.1.max(r.1) && q.1 >= p.1.min(r.1)
}

/// Check if line segment AB intersects line segment CD
fn segments_intersect(a: &(f64, f64), b: &(f64, f64), c: &(f64, f64), d: &(f64, f64)) -> bool {
    let o1 = orientation(a, b, c);
    let o2 = orientation(a, b, d);
    let o3 = orientation(c, d, a);
    let o4 = orientation(c, d, b);
    // Different signs means straddle
    if (o1 > 0.0) != (o2 > 0.0) && (o3 > 0.0) != (o4 > 0.0) {
        return true;
    }
    // Collinear cases
    if o1 == 0.0 && on_segment(a, c, b) {
        return true;
    }
    if o2 == 0.0 && on_segment(a, d, b) {
        return true;
    }
    if o3 == 0.0 && on_segment(c, a, d) {
        return true;
    }
    if o4 == 0.0 && on_segment(c, b, d) {
        return true;
    }
    false
}

/// Point-in-polygon with holes using winding number (non-zero rule for robustness)
fn point_in_polygon(p: &(f64, f64), polygon: &[Vec<(f64, f64)>]) -> bool {
    let mut winding = winding_number(p, &polygon[0]);
    for hole in &polygon[1..] {
        winding -= winding_number(p, hole);
    }
    winding != 0
}

/// Winding number for a single ring (positive for CCW)
fn winding_number(p: &(f64, f64), ring: &[(f64, f64)]) -> i32 {
    let n = ring.len();
    let mut winding = 0;
    let mut j = n - 1;
    for i in 0..n {
        let vi = &ring[i];
        let vj = &ring[j];
        if ((vi.1 > p.1) != (vj.1 > p.1))
            && (p.0 < vi.0 + (vj.0 - vi.0) * (p.1 - vi.1) / (vj.1 - vi.1))
        {
            winding += 1;
        }
        j = i;
    }
    winding
}

/// Polygon-AABB intersection predicate.
///
/// Use for the common case of searching for geometries within the bounding rectangle of a viewport
/// on a Web Mercator projection.
///
/// polygon: first ring is exterior, rest are holes.
/// aabb: [min_corner, max_corner]
pub fn polygon_aabb_intersect(polygon: &[Vec<(f64, f64)>], aabb: [(f64, f64); 2]) -> bool {
    let [(min_x, min_y), (max_x, max_y)] = aabb;

    // Step 1: Trivial MBR rejection (exterior only)
    let exterior = &polygon[0];
    let (mut poly_min_x, mut poly_min_y) = (f64::INFINITY, f64::INFINITY);
    let (mut poly_max_x, mut poly_max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &(x, y) in exterior {
        poly_min_x = poly_min_x.min(x);
        poly_min_y = poly_min_y.min(y);
        poly_max_x = poly_max_x.max(x);
        poly_max_y = poly_max_y.max(y);
    }
    if poly_max_x < min_x || poly_min_x > max_x || poly_max_y < min_y || poly_min_y > max_y {
        return false;
    }

    // Step 2: Any polygon vertex inside rectangle?
    for ring in polygon {
        for &(x, y) in ring {
            if x >= min_x && x <= max_x && y >= min_y && y <= max_y {
                return true;
            }
        }
    }

    // Step 3: Any rectangle corner inside polygon?
    let corners = [
        (min_x, min_y),
        (min_x, max_y),
        (max_x, min_y),
        (max_x, max_y),
    ];
    for corner in corners {
        if point_in_polygon(&corner, polygon) {
            return true;
        }
    }

    // Step 4: Any polygon edge intersects rectangle edge?
    let sides = [
        ((min_x, min_y), (min_x, max_y)), // Left
        ((min_x, max_y), (max_x, max_y)), // Top
        ((max_x, max_y), (max_x, min_y)), // Right
        ((max_x, min_y), (min_x, min_y)), // Bottom
    ];
    for ring in polygon {
        let n = ring.len();
        for i in 0..n {
            let a = ring[i];
            let b = ring[(i + 1) % n];
            let code_a = outcode(&a, &aabb);
            let code_b = outcode(&b, &aabb);
            if (code_a & code_b) != 0 {
                continue;
            }
            for &(c, d) in &sides {
                if segments_intersect(&a, &b, &c, &d) {
                    return true;
                }
            }
        }
    }

    false
}
