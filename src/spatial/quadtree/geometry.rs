// =============================================================================
// Point2D - 2D Point and Vector Operations
// =============================================================================

/// A 2D point or vector with f64 coordinates.
///
/// Used for both positions and directions in the plane.
/// Implements basic vector operations: addition, subtraction, scaling,
/// dot product, cross product (2D determinant), and normalization.
#[derive(Debug, Clone, Copy, PartialEq)]
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
// Bounds - Global Coordinate System Bounds
// =============================================================================

/// Axis-aligned bounding box representing the global coordinate space.
///
/// All points indexed by the quadtree must fall within these bounds.
/// The bounds are used to normalize coordinates to [0, 1] for cell ID encoding.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds {
    pub min_x: f64,
    pub min_y: f64,
    pub max_x: f64,
    pub max_y: f64,
}

impl Bounds {
    /// Creates new bounds from min and max coordinates.
    ///
    /// # Panics
    /// Panics if min >= max for either dimension.
    pub fn new(min_x: f64, min_y: f64, max_x: f64, max_y: f64) -> Self {
        assert!(min_x < max_x, "min_x must be less than max_x");
        assert!(min_y < max_y, "min_y must be less than max_y");
        Self {
            min_x,
            min_y,
            max_x,
            max_y,
        }
    }

    /// Creates bounds that encompass the entire f64 range (not recommended for indexing).
    pub fn infinite() -> Self {
        Self {
            min_x: f64::NEG_INFINITY,
            min_y: f64::NEG_INFINITY,
            max_x: f64::INFINITY,
            max_y: f64::INFINITY,
        }
    }

    /// Creates unit bounds [0, 1] x [0, 1].
    pub fn unit() -> Self {
        Self {
            min_x: 0.0,
            min_y: 0.0,
            max_x: 1.0,
            max_y: 1.0,
        }
    }

    /// Width of the bounds (max_x - min_x).
    #[inline]
    pub fn width(&self) -> f64 {
        self.max_x - self.min_x
    }

    /// Height of the bounds (max_y - min_y).
    #[inline]
    pub fn height(&self) -> f64 {
        self.max_y - self.min_y
    }

    /// Center point of the bounds.
    #[inline]
    pub fn center(&self) -> Point2D {
        Point2D {
            x: (self.min_x + self.max_x) * 0.5,
            y: (self.min_y + self.max_y) * 0.5,
        }
    }

    /// Normalizes a point from world coordinates to [0, 1] x [0, 1].
    ///
    /// Points exactly at min map to 0, points at max map to 1.
    #[inline]
    pub fn normalize(&self, point: &Point2D) -> Point2D {
        Point2D {
            x: (point.x - self.min_x) / self.width(),
            y: (point.y - self.min_y) / self.height(),
        }
    }

    /// Denormalizes a point from [0, 1] x [0, 1] back to world coordinates.
    #[inline]
    pub fn denormalize(&self, normalized: &Point2D) -> Point2D {
        Point2D {
            x: self.min_x + normalized.x * self.width(),
            y: self.min_y + normalized.y * self.height(),
        }
    }

    /// Checks if a point is contained within the bounds (inclusive).
    #[inline]
    pub fn contains_point(&self, point: &Point2D) -> bool {
        point.x >= self.min_x
            && point.x <= self.max_x
            && point.y >= self.min_y
            && point.y <= self.max_y
    }

    /// Converts to a Rect.
    #[inline]
    pub fn to_rect(&self) -> Rect {
        Rect {
            x: Interval::new(self.min_x, self.max_x),
            y: Interval::new(self.min_y, self.max_y),
        }
    }

    /// Returns the corner point at the given index (0-3).
    ///
    /// Index layout (in Z-order):
    /// - 0: (min_x, min_y) - bottom-left
    /// - 1: (max_x, min_y) - bottom-right
    /// - 2: (min_x, max_y) - top-left
    /// - 3: (max_x, max_y) - top-right
    #[inline]
    pub fn corner(&self, index: usize) -> Point2D {
        debug_assert!(index < 4);
        Point2D {
            x: if index & 1 == 0 {
                self.min_x
            } else {
                self.max_x
            },
            y: if index & 2 == 0 {
                self.min_y
            } else {
                self.max_y
            },
        }
    }
}

// =============================================================================
// Interval - 1D Interval
// =============================================================================

/// A closed 1D interval [lo, hi].
#[derive(Debug, Clone, Copy, PartialEq)]
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
#[cfg(test)]
mod tests {
    use super::*;
    use super::*;

    // -------------------------------------------------------------------------
    // Point2D Tests
    // -------------------------------------------------------------------------

    mod point2d_tests {
        use super::*;

        #[test]
        fn test_new_and_origin() {
            let p = Point2D::new(3.0, 4.0);
            assert_eq!(p.x, 3.0);
            assert_eq!(p.y, 4.0);

            assert_eq!(Point2D::ORIGIN.x, 0.0);
            assert_eq!(Point2D::ORIGIN.y, 0.0);
        }

        #[test]
        fn test_add_sub() {
            let a = Point2D::new(1.0, 2.0);
            let b = Point2D::new(3.0, 4.0);

            let sum = a + b;
            assert_eq!(sum.x, 4.0);
            assert_eq!(sum.y, 6.0);

            let diff = b - a;
            assert_eq!(diff.x, 2.0);
            assert_eq!(diff.y, 2.0);
        }

        #[test]
        fn test_scale_mul() {
            let p = Point2D::new(2.0, 3.0);

            let scaled = p.scale(2.0);
            assert_eq!(scaled.x, 4.0);
            assert_eq!(scaled.y, 6.0);

            let multiplied = p * 3.0;
            assert_eq!(multiplied.x, 6.0);
            assert_eq!(multiplied.y, 9.0);

            let multiplied2 = 4.0 * p;
            assert_eq!(multiplied2.x, 8.0);
            assert_eq!(multiplied2.y, 12.0);
        }

        #[test]
        fn test_neg() {
            let p = Point2D::new(3.0, -4.0);
            let neg = -p;
            assert_eq!(neg.x, -3.0);
            assert_eq!(neg.y, 4.0);
        }

        #[test]
        fn test_dot_product() {
            let a = Point2D::new(1.0, 2.0);
            let b = Point2D::new(3.0, 4.0);
            assert_eq!(a.dot(&b), 11.0); // 1*3 + 2*4 = 11
        }

        #[test]
        fn test_cross_product() {
            let a = Point2D::new(1.0, 0.0);
            let b = Point2D::new(0.0, 1.0);

            // a x b = 1*1 - 0*0 = 1 (b is CCW from a)
            assert_eq!(a.cross(&b), 1.0);

            // b x a = 0*0 - 1*1 = -1 (a is CW from b)
            assert_eq!(b.cross(&a), -1.0);
        }

        #[test]
        fn test_norm() {
            let p = Point2D::new(3.0, 4.0);
            assert_eq!(p.norm_squared(), 25.0);
            assert_eq!(p.norm(), 5.0);
        }

        #[test]
        fn test_normalize() {
            let p = Point2D::new(3.0, 4.0);
            let n = p.normalize();

            // Should be unit length
            let len = n.norm();
            assert!((len - 1.0).abs() < 1e-10);

            // Should preserve direction
            assert!((n.x - 0.6).abs() < 1e-10);
            assert!((n.y - 0.8).abs() < 1e-10);
        }

        #[test]
        fn test_orthogonal() {
            let p = Point2D::new(3.0, 4.0);
            let o = p.orthogonal();

            // Orthogonal should be perpendicular (dot product = 0)
            assert!((p.dot(&o)).abs() < 1e-10);

            // Orthogonal is (-y, x)
            assert_eq!(o.x, -4.0);
            assert_eq!(o.y, 3.0);
        }

        #[test]
        fn test_lerp() {
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(10.0, 20.0);

            let mid = a.lerp(&b, 0.5);
            assert_eq!(mid.x, 5.0);
            assert_eq!(mid.y, 10.0);

            let quarter = a.lerp(&b, 0.25);
            assert_eq!(quarter.x, 2.5);
            assert_eq!(quarter.y, 5.0);

            // t=0 should give a
            let at_a = a.lerp(&b, 0.0);
            assert_eq!(at_a.x, 0.0);
            assert_eq!(at_a.y, 0.0);

            // t=1 should give b
            let at_b = a.lerp(&b, 1.0);
            assert_eq!(at_b.x, 10.0);
            assert_eq!(at_b.y, 20.0);
        }

        #[test]
        fn test_distance() {
            let a = Point2D::new(0.0, 0.0);
            let b = Point2D::new(3.0, 4.0);

            assert_eq!(a.distance_to(&b), 5.0);
            assert_eq!(a.distance_squared_to(&b), 25.0);
        }
    }

    // -------------------------------------------------------------------------
    // Interval Tests
    // -------------------------------------------------------------------------

    mod interval_tests {
        use super::*;

        #[test]
        fn test_new_orders_correctly() {
            let i1 = Interval::new(1.0, 5.0);
            assert_eq!(i1.lo, 1.0);
            assert_eq!(i1.hi, 5.0);

            // Should auto-order
            let i2 = Interval::new(5.0, 1.0);
            assert_eq!(i2.lo, 1.0);
            assert_eq!(i2.hi, 5.0);
        }

        #[test]
        fn test_empty_interval() {
            let empty = Interval::empty();
            assert!(empty.is_empty());
            assert_eq!(empty.length(), 0.0);
        }

        #[test]
        fn test_length_and_center() {
            let i = Interval::new(2.0, 8.0);
            assert_eq!(i.length(), 6.0);
            assert_eq!(i.center(), 5.0);
        }

        #[test]
        fn test_contains() {
            let i = Interval::new(2.0, 8.0);

            assert!(i.contains(2.0)); // boundary
            assert!(i.contains(5.0)); // interior
            assert!(i.contains(8.0)); // boundary
            assert!(!i.contains(1.9)); // outside
            assert!(!i.contains(8.1)); // outside
        }

        #[test]
        fn test_contains_interval() {
            let outer = Interval::new(0.0, 10.0);
            let inner = Interval::new(2.0, 8.0);

            assert!(outer.contains_interval(&inner));
            assert!(!inner.contains_interval(&outer));
            assert!(outer.contains_interval(&outer)); // self-containment
        }

        #[test]
        fn test_intersects() {
            let a = Interval::new(0.0, 5.0);
            let b = Interval::new(3.0, 8.0);
            let c = Interval::new(6.0, 10.0);

            assert!(a.intersects(&b)); // overlap
            assert!(!a.intersects(&c)); // no overlap
            assert!(a.intersects(&a)); // self-intersection
        }

        #[test]
        fn test_intersection() {
            let a = Interval::new(0.0, 5.0);
            let b = Interval::new(3.0, 8.0);

            let inter = a.intersection(&b);
            assert_eq!(inter.lo, 3.0);
            assert_eq!(inter.hi, 5.0);

            let c = Interval::new(6.0, 10.0);
            let no_inter = a.intersection(&c);
            assert!(no_inter.is_empty());
        }

        #[test]
        fn test_union() {
            let a = Interval::new(0.0, 3.0);
            let b = Interval::new(5.0, 8.0);

            let u = a.union(&b);
            assert_eq!(u.lo, 0.0);
            assert_eq!(u.hi, 8.0);
        }

        #[test]
        fn test_expanded() {
            let i = Interval::new(2.0, 8.0);
            let expanded = i.expanded(1.0);

            assert_eq!(expanded.lo, 1.0);
            assert_eq!(expanded.hi, 9.0);
        }
    }

    // -------------------------------------------------------------------------
    // Bounds Tests
    // -------------------------------------------------------------------------

    mod bounds_tests {
        use super::*;

        #[test]
        fn test_new() {
            let b = Bounds::new(0.0, 0.0, 10.0, 20.0);
            assert_eq!(b.min_x, 0.0);
            assert_eq!(b.min_y, 0.0);
            assert_eq!(b.max_x, 10.0);
            assert_eq!(b.max_y, 20.0);
        }

        #[test]
        #[should_panic]
        fn test_new_invalid_x() {
            Bounds::new(10.0, 0.0, 0.0, 20.0); // min_x > max_x
        }

        #[test]
        fn test_width_height_center() {
            let b = Bounds::new(0.0, 0.0, 10.0, 20.0);
            assert_eq!(b.width(), 10.0);
            assert_eq!(b.height(), 20.0);

            let c = b.center();
            assert_eq!(c.x, 5.0);
            assert_eq!(c.y, 10.0);
        }

        #[test]
        fn test_normalize_denormalize() {
            let b = Bounds::new(100.0, 200.0, 200.0, 400.0);

            // Min corner should normalize to (0, 0)
            let p_min = Point2D::new(100.0, 200.0);
            let n_min = b.normalize(&p_min);
            assert!((n_min.x - 0.0).abs() < 1e-10);
            assert!((n_min.y - 0.0).abs() < 1e-10);

            // Max corner should normalize to (1, 1)
            let p_max = Point2D::new(200.0, 400.0);
            let n_max = b.normalize(&p_max);
            assert!((n_max.x - 1.0).abs() < 1e-10);
            assert!((n_max.y - 1.0).abs() < 1e-10);

            // Center should normalize to (0.5, 0.5)
            let p_center = Point2D::new(150.0, 300.0);
            let n_center = b.normalize(&p_center);
            assert!((n_center.x - 0.5).abs() < 1e-10);
            assert!((n_center.y - 0.5).abs() < 1e-10);

            // Round-trip
            let original = Point2D::new(175.0, 350.0);
            let normalized = b.normalize(&original);
            let recovered = b.denormalize(&normalized);
            assert!((recovered.x - original.x).abs() < 1e-10);
            assert!((recovered.y - original.y).abs() < 1e-10);
        }

        #[test]
        fn test_contains_point() {
            let b = Bounds::new(0.0, 0.0, 10.0, 10.0);

            assert!(b.contains_point(&Point2D::new(5.0, 5.0))); // interior
            assert!(b.contains_point(&Point2D::new(0.0, 0.0))); // corner
            assert!(b.contains_point(&Point2D::new(10.0, 10.0))); // corner
            assert!(!b.contains_point(&Point2D::new(-1.0, 5.0))); // outside
            assert!(!b.contains_point(&Point2D::new(11.0, 5.0))); // outside
        }

        #[test]
        fn test_corner() {
            let b = Bounds::new(0.0, 0.0, 10.0, 20.0);

            assert_eq!(b.corner(0), Point2D::new(0.0, 0.0)); // bottom-left
            assert_eq!(b.corner(1), Point2D::new(10.0, 0.0)); // bottom-right
            assert_eq!(b.corner(2), Point2D::new(0.0, 20.0)); // top-left
            assert_eq!(b.corner(3), Point2D::new(10.0, 20.0)); // top-right
        }
    }

    // -------------------------------------------------------------------------
    // Rect Tests
    // -------------------------------------------------------------------------

    mod rect_tests {
        use super::*;

        #[test]
        fn test_from_coords() {
            let r = Rect::from_coords(1.0, 2.0, 5.0, 8.0);
            assert_eq!(r.x.lo, 1.0);
            assert_eq!(r.x.hi, 5.0);
            assert_eq!(r.y.lo, 2.0);
            assert_eq!(r.y.hi, 8.0);
        }

        #[test]
        fn test_from_center() {
            let r = Rect::from_center(&Point2D::new(5.0, 5.0), 2.0, 3.0);
            assert_eq!(r.x.lo, 3.0);
            assert_eq!(r.x.hi, 7.0);
            assert_eq!(r.y.lo, 2.0);
            assert_eq!(r.y.hi, 8.0);
        }

        #[test]
        fn test_dimensions() {
            let r = Rect::from_coords(0.0, 0.0, 4.0, 3.0);
            assert_eq!(r.width(), 4.0);
            assert_eq!(r.height(), 3.0);
            assert_eq!(r.area(), 12.0);
        }

        #[test]
        fn test_center_lo_hi() {
            let r = Rect::from_coords(2.0, 4.0, 8.0, 10.0);
            assert_eq!(r.center(), Point2D::new(5.0, 7.0));
            assert_eq!(r.lo(), Point2D::new(2.0, 4.0));
            assert_eq!(r.hi(), Point2D::new(8.0, 10.0));
        }

        #[test]
        fn test_contains_point() {
            let r = Rect::from_coords(0.0, 0.0, 10.0, 10.0);

            assert!(r.contains_point(&Point2D::new(5.0, 5.0)));
            assert!(r.contains_point(&Point2D::new(0.0, 0.0)));
            assert!(!r.contains_point(&Point2D::new(-1.0, 5.0)));
        }

        #[test]
        fn test_contains_rect() {
            let outer = Rect::from_coords(0.0, 0.0, 10.0, 10.0);
            let inner = Rect::from_coords(2.0, 2.0, 8.0, 8.0);
            let partial = Rect::from_coords(5.0, 5.0, 15.0, 15.0);

            assert!(outer.contains_rect(&inner));
            assert!(!inner.contains_rect(&outer));
            assert!(!outer.contains_rect(&partial));
        }

        #[test]
        fn test_intersects() {
            let a = Rect::from_coords(0.0, 0.0, 5.0, 5.0);
            let b = Rect::from_coords(3.0, 3.0, 8.0, 8.0);
            let c = Rect::from_coords(6.0, 6.0, 10.0, 10.0);

            assert!(a.intersects(&b));
            assert!(!a.intersects(&c));
        }

        #[test]
        fn test_intersection() {
            let a = Rect::from_coords(0.0, 0.0, 5.0, 5.0);
            let b = Rect::from_coords(3.0, 3.0, 8.0, 8.0);

            let inter = a.intersection(&b);
            assert_eq!(inter.x.lo, 3.0);
            assert_eq!(inter.x.hi, 5.0);
            assert_eq!(inter.y.lo, 3.0);
            assert_eq!(inter.y.hi, 5.0);
        }

        #[test]
        fn test_union() {
            let a = Rect::from_coords(0.0, 0.0, 3.0, 3.0);
            let b = Rect::from_coords(5.0, 5.0, 8.0, 8.0);

            let u = a.union(&b);
            assert_eq!(u.x.lo, 0.0);
            assert_eq!(u.x.hi, 8.0);
            assert_eq!(u.y.lo, 0.0);
            assert_eq!(u.y.hi, 8.0);
        }

        #[test]
        fn test_corner() {
            let r = Rect::from_coords(0.0, 0.0, 10.0, 20.0);

            assert_eq!(r.corner(0), Point2D::new(0.0, 0.0)); // bottom-left
            assert_eq!(r.corner(1), Point2D::new(10.0, 0.0)); // bottom-right
            assert_eq!(r.corner(2), Point2D::new(0.0, 20.0)); // top-left
            assert_eq!(r.corner(3), Point2D::new(10.0, 20.0)); // top-right
        }
    }
}
