//! Exact arithmetic predicates for spherical geometry.
//!
//! This module provides robust geometric predicates that handle degenerate cases using exact
//! arithmetic and symbolic perturbation.
//!
//! s2predicates.cc

use super::math::cross;

/// Uses arbitrary-precision arithmetic and the "simulation of simplicity" technique in order to be
/// completely robust (i.e., to return consistent results for all possible inputs).
///
/// Below we define a floating-point type with enough precision so that it can represent the exact
/// determinant of any 3x3 matrix of floating-point numbers.  It uses ExactFloat, which is based on
/// the OpenSSL Bignum library and therefore has a permissive BSD-style license.  (At one time we
/// also supported an option based on MPFR, but that has an LGPL license and is
/// therefore not suited for some applications.)
pub fn expensive_sign(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3], perturb: bool) -> i32 {
    // Return zero if and only if two points are the same.  This ensures (1).
    if a == b || b == c || c == a {
        return 0;
    }

    // Next we try recomputing the determinant still using floating-point
    // arithmetic but in a more precise way.  This is more expensive than the
    // simple calculation done by TriageSign(), but it is still *much* cheaper
    // than using arbitrary-precision arithmetic.  This optimization is able to
    // compute the correct determinant sign in virtually all cases except when
    // the three points are truly collinear (e.g., three points on the equator).
    let det_sign = stable_sign(a, b, c);
    if det_sign != 0 {
        return det_sign;
    }

    // Otherwise fall back to exact arithmetic and symbolic permutations.
    exact_sign(a, b, c, perturb)
}

/// Compute the determinant using exact arithmetic and/or symbolic
/// permutations.  Requires that the three points are distinct.
fn exact_sign(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3], perturb: bool) -> i32 {
    let result = robust::orient3d(
        robust::Coord3D {
            x: a[0],
            y: a[1],
            z: a[2],
        },
        robust::Coord3D {
            x: b[0],
            y: b[1],
            z: b[2],
        },
        robust::Coord3D {
            x: c[0],
            y: c[1],
            z: c[2],
        },
        robust::Coord3D {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
    );

    if result > 0.0 {
        1
    } else if result < 0.0 {
        -1
    } else if perturb {
        symbolically_perturbed_sign(a, b, c)
    } else {
        0
    }
}

// Compute the determinant in a numerically stable way.  Unlike TriageSign(),
// this method can usually compute the correct determinant sign even when all
// three points are as collinear as possible.  For example if three points are
// spaced 1km apart along a random line on the Earth's surface using the
// nearest representable points, there is only a 0.4% chance that this method
// will not be able to find the determinant sign.  The probability of failure
// decreases as the points get closer together; if the collinear points are
// 1 meter apart, the failure rate drops to 0.0004%.
//
// This method could be extended to also handle nearly-antipodal points (and
// in fact an earlier version of this code did exactly that), but antipodal
// points are rare in practice so it seems better to simply fall back to
// exact arithmetic in that case.
fn stable_sign(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> i32 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let ca = [a[0] - c[0], a[1] - c[1], a[2] - c[2]];
    let ab2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    let bc2 = bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2];
    let ca2 = ca[0] * ca[0] + ca[1] * ca[1] + ca[2] * ca[2];

    // Now compute the determinant ((A-C)x(B-C)).C, where the vertices have been
    // cyclically permuted if necessary so that AB is the longest edge.  (This
    // minimizes the magnitude of cross product.)  At the same time we also
    // compute the maximum error in the determinant.  Using a similar technique
    // to the one used for kMaxDetError, the error is at most
    //
    //   |d| <= (3 + 6/sqrt(3)) * |A-C| * |B-C| * e
    //
    // where e = 0.5 * DBL_EPSILON.  If the determinant magnitude is larger than
    // this value then we know its sign with certainty.
    const DET_ERROR_MULTIPLIER: f64 = 3.2321 * f64::EPSILON; // see above
    let (det, max_error) = if ab2 >= bc2 && ab2 >= ca2 {
        // AB is the longest edge, so compute (A-C)x(B-C).C.
        let cross = cross(&ca, &bc);
        let det = -(cross[0] * c[0] + cross[1] * c[1] + cross[2] * c[2]);
        let max_error = DET_ERROR_MULTIPLIER * (ca2 * bc2).sqrt();
        (det, max_error)
    } else if bc2 >= ca2 {
        // BC is the longest edge, so compute (B-A)x(C-A).A.
        let cross = cross(&ab, &ca);
        let det = -(cross[0] * a[0] + cross[1] * a[1] + cross[2] * a[2]);
        let max_error = DET_ERROR_MULTIPLIER * (ab2 * ca2).sqrt();
        (det, max_error)
    } else {
        // CA is the longest edge, so compute (C-B)x(A-B).B.
        let cross = cross(&bc, &ab);
        let det = -(cross[0] * b[0] + cross[1] * b[1] + cross[2] * b[2]);
        let max_error = DET_ERROR_MULTIPLIER * (bc2 * ab2).sqrt();
        (det, max_error)
    };

    // Errors smaller than this value may not be accurate due to underflow.
    const MIN_NO_UNDERFLOW_ERROR: f64 = 3.2321 * f64::EPSILON * 1.4916681462400413e-154; // sqrt(DBL_MIN)
    if max_error < MIN_NO_UNDERFLOW_ERROR {
        return 0;
    }

    if det.abs() <= max_error {
        0
    } else if det > 0.0 {
        1
    } else {
        -1
    }
}

/// The following function returns the sign of the determinant of three points A, B, C under a
/// model where every possible S2Point is slightly perturbed by a unique infinitesmal amount such
/// that no three perturbed points are collinear and no four points are coplanar.  The
/// perturbations are so small that they do not change the sign of any determinant that was
/// non-zero before the perturbations, and therefore can be safely ignored unless the determinant
/// of three points is exactly zero (using multiple-precision arithmetic).
///
/// Since the symbolic perturbation of a given point is fixed (i.e., the perturbation is the same
/// for all calls to this method and does not depend on the other two arguments), the results of
/// this method are always self-consistent.  It will never return results that would correspond to
/// an "impossible" configuration of non-degenerate points.
///
/// Requirements:
/// - The 3x3 determinant of A, B, C must be exactly zero.
/// - The points must be distinct, with A < B < C in lexicographic order.
///
/// Returns:
/// - +1 or -1 according to the sign of the determinant after the symbolic perturbations are taken
///   into account.
///
/// Reference:
/// - "Simulation of Simplicity" (Edelsbrunner and Muecke, ACM Transactions on Graphics, 1990).
fn symbolically_perturbed_sign(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> i32 {
    // This method requires that the points are sorted in lexicographically increasing order.  This
    // is because every possible S2Point has its own symbolic perturbation such that if A < B then
    // the symbolic perturbation for A is much larger than the perturbation for B.
    //
    // Alternatively, we could sort the points in this method and keep track of the sign of the
    // permutation, but it is more efficient to do this before converting the inputs to the
    // multi-precision representation, and this also lets us re-use the result of the cross product
    // B x C.
    let mut points = [*a, *b, *c];
    let perm_sign = sort_with_parity(&mut points);

    let (pa, pb, pc) = (points[0], points[1], points[2]);

    // Compute b x c here instead of the caller as in C++.
    let bc = cross(&pb, &pc);

    // Check components in order of decreasing perturbation magnitude.
    let det_sign = sign_of(bc[2]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    let det_sign = sign_of(bc[1]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    let det_sign = sign_of(bc[0]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    let det_sign = sign_of(pc[0] * pa[1] - pc[1] * pa[0]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    let det_sign = sign_of(pc[0]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    let det_sign = sign_of(-pc[1]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    let det_sign = sign_of(pc[2] * pa[0] - pc[0] * pa[2]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    let det_sign = sign_of(pc[2]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    let det_sign = sign_of(pa[0] * pb[1] - pa[1] * pb[0]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    let det_sign = sign_of(-pb[0]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    let det_sign = sign_of(pb[1]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    let det_sign = sign_of(pa[0]);
    if det_sign != 0 {
        return perm_sign * det_sign;
    }

    perm_sign
}

// Sort the three points in lexicographic order, keeping track of the sign of the permutation.
// (Each exchange inverts the sign of the determinant.)
fn sort_with_parity(points: &mut [[f64; 3]; 3]) -> i32 {
    let mut parity = 1;

    if points[1] < points[0] {
        points.swap(0, 1);
        parity = -parity;
    }
    if points[2] < points[1] {
        points.swap(1, 2);
        parity = -parity;
    }
    if points[1] < points[0] {
        points.swap(0, 1);
        parity = -parity;
    }

    parity
}

/// Return the sign of a floating-point value.
#[inline]
fn sign_of(x: f64) -> i32 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}
