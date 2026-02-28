// =============================================================================
// Point2D - 2D Point and Vector Operations
// =============================================================================

use std::ops::{Add, Mul, Neg, Sub};

use crate::spatial::quadtree::bounds::Bounds;

/// A 2D point or vector with f64 coordinates.
///
/// Used for both positions and directions in the plane.
/// Implements basic vector operations: addition, subtraction, scaling,
/// dot product, cross product (2D determinant), and normalization.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(missing_docs)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

impl Point2D {
    /// Creates a new point from x and y coordinates.
    #[inline]
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// The origin point (0, 0).
    pub const ORIGIN: Self = Self { x: 0.0, y: 0.0 };

    /// Dot product with another point/vector.
    ///
    /// Returns a . b = a.x * b.x + a.y * b.y
    #[inline]
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// 2D cross product (perpendicular dot product).
    ///
    /// Returns the z-component of the 3D cross product if the points
    /// were embedded in 3D with z=0: a.x * b.y - a.y * b.x
    ///
    /// This is also the signed area of the parallelogram formed by the vectors.
    /// Positive if b is counterclockwise from a, negative if clockwise.
    #[inline]
    pub fn cross(&self, other: &Self) -> f64 {
        self.x * other.y - self.y * other.x
    }

    /// Squared Euclidean norm (length squared).
    #[inline]
    pub fn norm_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y
    }

    /// Euclidean norm (length).
    #[inline]
    pub fn norm(&self) -> f64 {
        self.norm_squared().sqrt()
    }

    /// Returns a normalized (unit length) vector in the same direction.
    ///
    /// Panics if the vector has zero length.
    #[inline]
    pub fn normalize(&self) -> Self {
        let n = self.norm();
        debug_assert!(n > 0.0, "cannot normalize zero vector");
        Self {
            x: self.x / n,
            y: self.y / n,
        }
    }

    /// Returns a vector perpendicular to this one (rotated 90 degrees CCW).
    ///
    /// For vector (x, y), returns (-y, x).
    #[inline]
    pub fn orthogonal(&self) -> Self {
        Self {
            x: -self.y,
            y: self.x,
        }
    }

    /// Scales the point by a scalar factor.
    #[inline]
    pub fn scale(&self, factor: f64) -> Self {
        Self {
            x: self.x * factor,
            y: self.y * factor,
        }
    }

    /// Linear interpolation between self and other.
    ///
    /// Returns self + t * (other - self) = (1-t) * self + t * other
    #[inline]
    pub fn lerp(&self, other: &Self, t: f64) -> Self {
        Self {
            x: self.x + t * (other.x - self.x),
            y: self.y + t * (other.y - self.y),
        }
    }

    /// Distance to another point.
    #[inline]
    pub fn distance_to(&self, other: &Self) -> f64 {
        (*self - *other).norm()
    }

    /// Squared distance to another point (avoids sqrt).
    #[inline]
    pub fn distance_squared_to(&self, other: &Self) -> f64 {
        (*self - *other).norm_squared()
    }
}

impl Add for Point2D {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub for Point2D {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Mul<f64> for Point2D {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: f64) -> Self {
        Self {
            x: self.x * scalar,
            y: self.y * scalar,
        }
    }
}

impl Mul<Point2D> for f64 {
    type Output = Point2D;

    #[inline]
    fn mul(self, point: Point2D) -> Point2D {
        Point2D {
            x: self * point.x,
            y: self * point.y,
        }
    }
}

impl Neg for Point2D {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}


// =============================================================================
// Interval - 1D Interval
// =============================================================================

/// A closed 1D interval [lo, hi].
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(missing_docs)]
pub struct Interval {
    pub lo: f64,
    pub hi: f64,
}

impl Interval {
    /// Creates a new interval. Automatically orders lo and hi.
    #[inline]
    pub fn new(a: f64, b: f64) -> Self {
        if a <= b {
            Self { lo: a, hi: b }
        } else {
            Self { lo: b, hi: a }
        }
    }

    /// Creates an empty interval.
    pub fn empty() -> Self {
        Self {
            lo: f64::INFINITY,
            hi: f64::NEG_INFINITY,
        }
    }

    /// Checks if the interval is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.lo > self.hi
    }

    /// Length of the interval.
    #[inline]
    pub fn length(&self) -> f64 {
        if self.is_empty() {
            0.0
        } else {
            self.hi - self.lo
        }
    }

    /// Center of the interval.
    #[inline]
    pub fn center(&self) -> f64 {
        (self.lo + self.hi) * 0.5
    }

    /// Checks if the interval contains a value.
    #[inline]
    pub fn contains(&self, value: f64) -> bool {
        value >= self.lo && value <= self.hi
    }

    /// Checks if this interval contains another interval.
    #[inline]
    pub fn contains_interval(&self, other: &Self) -> bool {
        if other.is_empty() {
            return true;
        }
        other.lo >= self.lo && other.hi <= self.hi
    }

    /// Checks if the intervals intersect.
    #[inline]
    pub fn intersects(&self, other: &Self) -> bool {
        self.lo <= other.hi && other.lo <= self.hi
    }

    /// Returns the intersection of two intervals.
    #[inline]
    pub fn intersection(&self, other: &Self) -> Self {
        Self {
            lo: self.lo.max(other.lo),
            hi: self.hi.min(other.hi),
        }
    }

    /// Returns the union of two intervals (smallest interval containing both).
    #[inline]
    pub fn union(&self, other: &Self) -> Self {
        if self.is_empty() {
            return *other;
        }
        if other.is_empty() {
            return *self;
        }
        Self {
            lo: self.lo.min(other.lo),
            hi: self.hi.max(other.hi),
        }
    }

    /// Expands the interval by the given amount on each side.
    #[inline]
    pub fn expanded(&self, margin: f64) -> Self {
        Self {
            lo: self.lo - margin,
            hi: self.hi + margin,
        }
    }
}

// =============================================================================
// Rect - Axis-Aligned Rectangle
// =============================================================================

/// An axis-aligned rectangle defined by x and y intervals.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(missing_docs)]
pub struct Rect {
    pub x: Interval,
    pub y: Interval,
}

impl Rect {
    /// Creates a new rectangle from x and y intervals.
    #[inline]
    pub fn new(x: Interval, y: Interval) -> Self {
        Self { x, y }
    }

    /// Creates a rectangle from corner coordinates.
    #[inline]
    pub fn from_coords(min_x: f64, min_y: f64, max_x: f64, max_y: f64) -> Self {
        Self {
            x: Interval::new(min_x, max_x),
            y: Interval::new(min_y, max_y),
        }
    }

    /// Creates a rectangle from a center point and half-widths.
    #[inline]
    pub fn from_center(center: &Point2D, half_width: f64, half_height: f64) -> Self {
        Self {
            x: Interval::new(center.x - half_width, center.x + half_width),
            y: Interval::new(center.y - half_height, center.y + half_height),
        }
    }

    /// Creates an empty rectangle.
    pub fn empty() -> Self {
        Self {
            x: Interval::empty(),
            y: Interval::empty(),
        }
    }

    /// Checks if the rectangle is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.x.is_empty() || self.y.is_empty()
    }

    /// Width of the rectangle.
    #[inline]
    pub fn width(&self) -> f64 {
        self.x.length()
    }

    /// Height of the rectangle.
    #[inline]
    pub fn height(&self) -> f64 {
        self.y.length()
    }

    /// Area of the rectangle.
    #[inline]
    pub fn area(&self) -> f64 {
        self.width() * self.height()
    }

    /// Center point of the rectangle.
    #[inline]
    pub fn center(&self) -> Point2D {
        Point2D {
            x: self.x.center(),
            y: self.y.center(),
        }
    }

    /// Low corner (min_x, min_y).
    #[inline]
    pub fn lo(&self) -> Point2D {
        Point2D {
            x: self.x.lo,
            y: self.y.lo,
        }
    }

    /// High corner (max_x, max_y).
    #[inline]
    pub fn hi(&self) -> Point2D {
        Point2D {
            x: self.x.hi,
            y: self.y.hi,
        }
    }

    /// Checks if the rectangle contains a point.
    #[inline]
    pub fn contains_point(&self, point: &Point2D) -> bool {
        self.x.contains(point.x) && self.y.contains(point.y)
    }

    /// Checks if this rectangle contains another rectangle.
    #[inline]
    pub fn contains_rect(&self, other: &Self) -> bool {
        self.x.contains_interval(&other.x) && self.y.contains_interval(&other.y)
    }

    /// Checks if the rectangles intersect.
    #[inline]
    pub fn intersects(&self, other: &Self) -> bool {
        self.x.intersects(&other.x) && self.y.intersects(&other.y)
    }

    /// Returns the intersection of two rectangles.
    #[inline]
    pub fn intersection(&self, other: &Self) -> Self {
        Self {
            x: self.x.intersection(&other.x),
            y: self.y.intersection(&other.y),
        }
    }

    /// Returns the union of two rectangles (smallest rectangle containing both).
    #[inline]
    pub fn union(&self, other: &Self) -> Self {
        Self {
            x: self.x.union(&other.x),
            y: self.y.union(&other.y),
        }
    }

    /// Expands the rectangle by the given margin on all sides.
    #[inline]
    pub fn expanded(&self, margin: f64) -> Self {
        Self {
            x: self.x.expanded(margin),
            y: self.y.expanded(margin),
        }
    }

    /// Converts to Bounds.
    #[inline]
    pub fn to_bounds(&self) -> Bounds {
        Bounds {
            min_x: self.x.lo,
            min_y: self.y.lo,
            max_x: self.x.hi,
            max_y: self.y.hi,
        }
    }

    /// Returns the corner point at the given index (0-3).
    ///
    /// Index layout (in Z-order):
    /// - 0: (lo.x, lo.y) - bottom-left
    /// - 1: (hi.x, lo.y) - bottom-right
    /// - 2: (lo.x, hi.y) - top-left
    /// - 3: (hi.x, hi.y) - top-right
    #[inline]
    pub fn corner(&self, index: usize) -> Point2D {
        debug_assert!(index < 4);
        Point2D {
            x: if index & 1 == 0 { self.x.lo } else { self.x.hi },
            y: if index & 2 == 0 { self.y.lo } else { self.y.hi },
        }
    }
}
