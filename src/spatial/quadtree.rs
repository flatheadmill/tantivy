//! Quadtree implementation for 2D spatial indexing.
//!
//! This module provides a multi-shape quadtree for Tantivy's spatial indexing,
//! based on the architecture of Google S2's MutableS2ShapeIndex but adapted
//! for planar (non-spherical) geometry.

use std::ops::{Add, Mul, Neg, Sub};

use smallvec::SmallVec;

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

// =============================================================================
// QuadtreeCellId - Z-Order Cell Encoding
// =============================================================================

/// Maximum level for quadtree cells (0-30).
/// At level 30 with typical geographic bounds, resolution is sub-centimeter.
pub const MAX_LEVEL: u8 = 30;

/// A cell ID using Z-order (Morton code) encoding.
///
/// The encoding uses 64 bits:
/// - Bits 63-1: Position (2 bits per level, up to 30 levels = 60 bits)
/// - Bit 0: Always 1 when shifted to mark the level (sentinel)
///
/// The cell ID format follows the pattern used by S2CellId:
/// - The position bits are stored at the high end
/// - A trailing 1-bit marks where the meaningful bits end
/// - The level can be computed by finding the position of the lowest set bit
///
/// Example for level 2 cell at position (1, 2) = binary 01, 10 -> Z-order 0110:
/// ```text
/// bit:  63      60 59                                    4  3  2  1  0
///       |-------|--|------------------------------------|----|----|----|
///       | 0110  |  |           ... zeros ...            |    | 1  |    |
///       position                                        sentinel
/// ```
/// The sentinel 1-bit is at position (62 - 2*level) = 62 - 4 = 58... actually
/// let me recalculate.
///
/// Actually, the encoding is:
/// - Position bits at the top (MSB)
/// - Sentinel 1-bit immediately after position bits
/// - All trailing bits are 0
///
/// For level L: the sentinel is at bit position 62 - 2*L
/// For level 0: position is in bits 63-62, sentinel at bit 61
/// For level 30: position uses all 60 bits (63-4), sentinel at bit 3
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QuadtreeCellId(pub u64);

impl QuadtreeCellId {
    /// The invalid/sentinel cell ID (value 0).
    pub const NONE: Self = Self(0);

    /// The root cell ID at level 0 covering the entire space.
    /// Position = 0, level = 0, so sentinel at bit 62.
    pub const ROOT: Self = Self(1 << 62);

    /// Creates a cell ID from a raw u64 value.
    #[inline]
    pub const fn from_raw(value: u64) -> Self {
        Self(value)
    }

    /// Returns the raw u64 value.
    #[inline]
    pub const fn to_raw(&self) -> u64 {
        self.0
    }

    /// Checks if this is a valid cell ID (non-zero with proper sentinel).
    #[inline]
    pub fn is_valid(&self) -> bool {
        self.0 != 0 && self.0.trailing_zeros() % 2 == 0
    }

    /// Returns the level of this cell (0 = root, MAX_LEVEL = finest).
    ///
    /// The level is determined by the position of the sentinel bit.
    #[inline]
    pub fn level(&self) -> u8 {
        debug_assert!(self.is_valid());
        // Sentinel is at bit (62 - 2*level), so level = (62 - sentinel_pos) / 2
        // trailing_zeros gives us the sentinel position
        let sentinel_pos = self.0.trailing_zeros();
        ((62 - sentinel_pos) / 2) as u8
    }

    /// Creates a cell ID from normalized coordinates at the given level.
    ///
    /// Coordinates should be in [0, 1) range. Values exactly at 1.0 are
    /// clamped to just below 1.0.
    pub fn from_normalized(x: f64, y: f64, level: u8) -> Self {
        debug_assert!(level <= MAX_LEVEL);

        // Clamp coordinates to valid range
        let x = x.clamp(0.0, 1.0 - f64::EPSILON);
        let y = y.clamp(0.0, 1.0 - f64::EPSILON);

        // Convert to integer coordinates at this level
        let scale = (1u64 << level) as f64;
        let ix = (x * scale) as u32;
        let iy = (y * scale) as u32;

        Self::from_ij(ix, iy, level)
    }

    /// Creates a cell ID from integer (i, j) coordinates at the given level.
    ///
    /// i is the x-coordinate (column), j is the y-coordinate (row).
    /// Both must be in range [0, 2^level).
    pub fn from_ij(i: u32, j: u32, level: u8) -> Self {
        debug_assert!(level <= MAX_LEVEL);
        let max_coord = 1u32 << level;
        debug_assert!(i < max_coord && j < max_coord);

        // Interleave i and j bits to create Z-order (Morton) code
        let position = interleave_bits(i, j);

        // Shift position to top and add sentinel
        // Position uses 2*level bits, starting at bit 62 (just below bit 63)
        let shift = 62 - 2 * level as u32;
        let id = (position << shift) | (1u64 << shift.saturating_sub(1).max(0));

        // Actually, let's reconsider the encoding:
        // For level L, we have 2*L position bits at the top, then a sentinel 1.
        // The sentinel is at bit (62 - 2*L).
        // But we need position bits in 63..(64-2*L), sentinel at (62-2*L).
        // Hmm, let me think again...

        // Simpler approach: position in high bits, sentinel marks end
        let sentinel_pos = 62 - 2 * level as u32;
        let id = (position << (sentinel_pos + 1)) | (1u64 << sentinel_pos);

        Self(id)
    }

    /// Extracts the (i, j) integer coordinates at this cell's level.
    pub fn to_ij(&self) -> (u32, u32) {
        let level = self.level();
        let sentinel_pos = 62 - 2 * level as u32;
        let position = self.0 >> (sentinel_pos + 1);
        deinterleave_bits(position, level)
    }

    /// Returns the center of this cell in normalized [0, 1] coordinates.
    pub fn to_center_normalized(&self) -> Point2D {
        let level = self.level();
        let (i, j) = self.to_ij();
        let scale = (1u64 << level) as f64;

        Point2D {
            x: (i as f64 + 0.5) / scale,
            y: (j as f64 + 0.5) / scale,
        }
    }

    /// Returns the center of this cell in world coordinates.
    pub fn to_center(&self, bounds: &Bounds) -> Point2D {
        bounds.denormalize(&self.to_center_normalized())
    }

    /// Returns the bounds of this cell in normalized [0, 1] coordinates.
    pub fn to_bounds_normalized(&self) -> Rect {
        let level = self.level();
        let (i, j) = self.to_ij();
        let scale = (1u64 << level) as f64;

        Rect {
            x: Interval::new(i as f64 / scale, (i + 1) as f64 / scale),
            y: Interval::new(j as f64 / scale, (j + 1) as f64 / scale),
        }
    }

    /// Returns the bounds of this cell in world coordinates.
    pub fn to_bounds(&self, bounds: &Bounds) -> Rect {
        let norm = self.to_bounds_normalized();
        Rect {
            x: Interval::new(
                bounds.min_x + norm.x.lo * bounds.width(),
                bounds.min_x + norm.x.hi * bounds.width(),
            ),
            y: Interval::new(
                bounds.min_y + norm.y.lo * bounds.height(),
                bounds.min_y + norm.y.hi * bounds.height(),
            ),
        }
    }

    /// Returns the parent cell (one level up).
    ///
    /// Returns None if this is already the root cell.
    pub fn parent(&self) -> Option<Self> {
        if self.level() == 0 {
            return None;
        }

        // Move sentinel two bits to the right and clear the old position bits
        let sentinel_pos = 62 - 2 * self.level() as u32;
        let new_sentinel_pos = sentinel_pos + 2;
        let mask = !((1u64 << (new_sentinel_pos + 1)) - 1);
        let id = (self.0 & mask) | (1u64 << new_sentinel_pos);

        Some(Self(id))
    }

    /// Returns the four child cells (one level down).
    ///
    /// Returns None if this is already at MAX_LEVEL.
    pub fn children(&self) -> Option<[Self; 4]> {
        let level = self.level();
        if level >= MAX_LEVEL {
            return None;
        }

        let (i, j) = self.to_ij();
        let new_level = level + 1;

        Some([
            Self::from_ij(2 * i, 2 * j, new_level),
            Self::from_ij(2 * i + 1, 2 * j, new_level),
            Self::from_ij(2 * i, 2 * j + 1, new_level),
            Self::from_ij(2 * i + 1, 2 * j + 1, new_level),
        ])
    }

    /// Returns the child at the given index (0-3).
    ///
    /// Index layout (Z-order):
    /// - 0: bottom-left
    /// - 1: bottom-right
    /// - 2: top-left
    /// - 3: top-right
    pub fn child(&self, index: usize) -> Option<Self> {
        debug_assert!(index < 4);
        let level = self.level();
        if level >= MAX_LEVEL {
            return None;
        }

        let (i, j) = self.to_ij();
        let di = (index & 1) as u32;
        let dj = ((index >> 1) & 1) as u32;

        Some(Self::from_ij(2 * i + di, 2 * j + dj, level + 1))
    }

    /// Returns the minimum cell ID in the range covered by this cell.
    ///
    /// This is the cell ID of the first (smallest) leaf cell contained
    /// within this cell. For leaf cells, returns self.
    pub fn range_min(&self) -> Self {
        // The minimum is this cell's position with the sentinel moved to MAX_LEVEL
        let level = self.level();
        let sentinel_pos = 62 - 2 * level as u32;
        let position = self.0 & !((1u64 << (sentinel_pos + 1)) - 1);

        // Set sentinel at MAX_LEVEL position
        let new_sentinel_pos = 62 - 2 * MAX_LEVEL as u32;
        Self(position | (1u64 << new_sentinel_pos))
    }

    /// Returns the maximum cell ID in the range covered by this cell.
    ///
    /// This is the cell ID of the last (largest) leaf cell contained
    /// within this cell. For leaf cells, returns self.
    pub fn range_max(&self) -> Self {
        // The maximum fills all descendant bits with 1s
        let level = self.level();
        let sentinel_pos = 62 - 2 * level as u32;

        // Keep position bits, fill lower bits with 1s up to sentinel
        let position = self.0 & !((1u64 << (sentinel_pos + 1)) - 1);
        let fill_mask = (1u64 << (sentinel_pos + 1)) - (1u64 << (62 - 2 * MAX_LEVEL as u32));

        // Set sentinel at MAX_LEVEL position
        let new_sentinel_pos = 62 - 2 * MAX_LEVEL as u32;
        Self(position | fill_mask | (1u64 << new_sentinel_pos))
    }

    /// Checks if this cell contains another cell (other is a descendant).
    pub fn contains(&self, other: &Self) -> bool {
        other.0 >= self.range_min().0 && other.0 <= self.range_max().0
    }

    /// Checks if this cell intersects another cell (one contains the other or they're equal).
    pub fn intersects(&self, other: &Self) -> bool {
        self.contains(other) || other.contains(self)
    }

    /// Returns the common ancestor of this cell and another.
    pub fn common_ancestor(&self, other: &Self) -> Self {
        // XOR to find differing bits, then mask to common prefix
        let diff = self.0 ^ other.0;
        if diff == 0 {
            return *self;
        }

        // Find the highest differing bit
        let highest_diff = 63 - diff.leading_zeros();

        // The common ancestor level is where the first difference occurs
        // Each level uses 2 bits starting from bit 62
        let diff_level = (62 - highest_diff) / 2;

        // Get parent at that level
        let mut result = *self;
        while result.level() > diff_level as u8 {
            result = result.parent().unwrap();
        }
        result
    }
}

/// Interleaves bits of two 32-bit values into a 64-bit Z-order code.
///
/// The result has bits of x in even positions and bits of y in odd positions.
fn interleave_bits(x: u32, y: u32) -> u64 {
    let mut x = x as u64;
    let mut y = y as u64;

    // Spread x bits: 0000abcd -> 0a0b0c0d
    x = (x | (x << 16)) & 0x0000_FFFF_0000_FFFF;
    x = (x | (x << 8)) & 0x00FF_00FF_00FF_00FF;
    x = (x | (x << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    x = (x | (x << 2)) & 0x3333_3333_3333_3333;
    x = (x | (x << 1)) & 0x5555_5555_5555_5555;

    // Spread y bits
    y = (y | (y << 16)) & 0x0000_FFFF_0000_FFFF;
    y = (y | (y << 8)) & 0x00FF_00FF_00FF_00FF;
    y = (y | (y << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    y = (y | (y << 2)) & 0x3333_3333_3333_3333;
    y = (y | (y << 1)) & 0x5555_5555_5555_5555;

    // Combine: x in even bits, y in odd bits
    x | (y << 1)
}

/// Deinterleaves a Z-order code back into (x, y) coordinates.
fn deinterleave_bits(z: u64, level: u8) -> (u32, u32) {
    let mask = (1u64 << (2 * level as u32)) - 1;
    let z = z & mask;

    // Extract even bits (x) and odd bits (y)
    let mut x = z & 0x5555_5555_5555_5555;
    let mut y = (z >> 1) & 0x5555_5555_5555_5555;

    // Compact x bits
    x = (x | (x >> 1)) & 0x3333_3333_3333_3333;
    x = (x | (x >> 2)) & 0x0F0F_0F0F_0F0F_0F0F;
    x = (x | (x >> 4)) & 0x00FF_00FF_00FF_00FF;
    x = (x | (x >> 8)) & 0x0000_FFFF_0000_FFFF;
    x = (x | (x >> 16)) & 0x0000_0000_FFFF_FFFF;

    // Compact y bits
    y = (y | (y >> 1)) & 0x3333_3333_3333_3333;
    y = (y | (y >> 2)) & 0x0F0F_0F0F_0F0F_0F0F;
    y = (y | (y >> 4)) & 0x00FF_00FF_00FF_00FF;
    y = (y | (y >> 8)) & 0x0000_FFFF_0000_FFFF;
    y = (y | (y >> 16)) & 0x0000_0000_FFFF_FFFF;

    (x as u32, y as u32)
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
// =============================================================================
// ClippedShape - Per-Shape Data Within a Cell
// =============================================================================

/// Per-shape data within a quadtree cell.
///
/// This structure tracks which edges of a shape intersect a particular cell,
/// and whether the shape's interior contains the cell center. This information
/// enables efficient point-in-polygon queries: start with `contains_center`,
/// then count edge crossings from the cell center to the query point.
///
/// # Memory Layout
///
/// The edge indices use SmallVec with inline storage for 2 edges. This optimizes
/// for the common case where a cell is crossed by few edges of each shape:
/// - 0-2 edges: no heap allocation (inline storage)
/// - 3+ edges: spills to heap
///
/// # Invariants
///
/// - Edge indices are stored in sorted order (ascending)
/// - Edge index `i` refers to the edge from vertex[i] to vertex[(i+1) % n]
#[derive(Debug, Clone)]
pub struct ClippedShape {
    /// Document ID (shape identifier within the segment).
    doc_id: u32,

    /// Whether the shape's interior contains the cell center.
    ///
    /// For shapes without an interior (e.g., polylines), this is always false.
    /// For polygons, this indicates whether the cell center is inside the polygon.
    contains_center: bool,

    /// Indices of edges that intersect this cell.
    ///
    /// Each index refers to an edge in the shape's vertex array:
    /// edge[i] connects vertex[i] to vertex[(i+1) % num_vertices].
    ///
    /// Edges are stored in sorted order for efficient lookup and merging.
    edges: SmallVec<[u16; 2]>,
}

impl ClippedShape {
    /// Creates a new ClippedShape with no edges.
    ///
    /// # Arguments
    ///
    /// * `doc_id` - The document ID of the shape
    /// * `contains_center` - Whether the shape interior contains the cell center
    #[inline]
    pub fn new(doc_id: u32, contains_center: bool) -> Self {
        Self {
            doc_id,
            contains_center,
            edges: SmallVec::new(),
        }
    }

    /// Creates a ClippedShape with preallocated edge capacity.
    #[inline]
    pub fn with_capacity(doc_id: u32, contains_center: bool, capacity: usize) -> Self {
        Self {
            doc_id,
            contains_center,
            edges: SmallVec::with_capacity(capacity),
        }
    }

    /// Returns the document ID.
    #[inline]
    pub fn doc_id(&self) -> u32 {
        self.doc_id
    }

    /// Returns whether the shape's interior contains the cell center.
    #[inline]
    pub fn contains_center(&self) -> bool {
        self.contains_center
    }

    /// Sets whether the shape's interior contains the cell center.
    #[inline]
    pub fn set_contains_center(&mut self, contains: bool) {
        self.contains_center = contains;
    }

    /// Returns the number of edges intersecting this cell.
    #[inline]
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    /// Returns the edge indices as a slice.
    ///
    /// The indices are in sorted order.
    #[inline]
    pub fn edges(&self) -> &[u16] {
        &self.edges
    }

    /// Returns the edge index at position `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i >= num_edges()`.
    #[inline]
    pub fn edge(&self, i: usize) -> u16 {
        self.edges[i]
    }

    /// Adds an edge index, maintaining sorted order.
    ///
    /// If the edge is already present, it is not added again.
    pub fn add_edge(&mut self, edge_idx: u16) {
        match self.edges.binary_search(&edge_idx) {
            Ok(_) => {} // Already present
            Err(pos) => self.edges.insert(pos, edge_idx),
        }
    }

    /// Adds multiple edge indices, maintaining sorted order.
    ///
    /// More efficient than calling `add_edge` repeatedly when adding
    /// multiple edges, as it sorts once at the end.
    pub fn add_edges(&mut self, edge_indices: &[u16]) {
        if edge_indices.is_empty() {
            return;
        }

        // For small additions to small vectors, individual inserts may be faster
        if edge_indices.len() <= 2 && self.edges.len() <= 4 {
            for &idx in edge_indices {
                self.add_edge(idx);
            }
            return;
        }

        // For larger additions, extend and sort
        let original_len = self.edges.len();
        self.edges.extend_from_slice(edge_indices);
        self.edges.sort_unstable();
        self.edges.dedup();

        // If no duplicates were possible, the sort was sufficient
        if self.edges.len() < original_len + edge_indices.len() {
            // Duplicates were removed by dedup
        }
    }

    /// Returns true if this shape contains the given edge index.
    #[inline]
    pub fn contains_edge(&self, edge_idx: u16) -> bool {
        self.edges.binary_search(&edge_idx).is_ok()
    }

    /// Returns true if this shape has no edges in the cell.
    ///
    /// Note: A shape with no edges can still have `contains_center = true`
    /// if the entire cell is contained within the shape's interior.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }

    /// Clears all edge indices, keeping doc_id and contains_center.
    #[inline]
    pub fn clear_edges(&mut self) {
        self.edges.clear();
    }
}

impl PartialEq for ClippedShape {
    fn eq(&self, other: &Self) -> bool {
        self.doc_id == other.doc_id
            && self.contains_center == other.contains_center
            && self.edges == other.edges
    }
}

impl Eq for ClippedShape {}

// =============================================================================
// QuadtreeCell - A Cell in the Index
// =============================================================================

/// A single cell in the quadtree index.
///
/// Each cell contains a set of ClippedShapes representing the shapes that
/// intersect this cell. Shapes are stored sorted by doc_id for efficient
/// lookup and merging.
///
/// # Cell Properties
///
/// - `cell_id`: Identifies the cell's position and level in the quadtree
/// - `shapes`: Shapes intersecting this cell, sorted by doc_id
///
/// # Query Usage
///
/// To test if a point is inside any shape:
/// 1. Find the cell containing the point
/// 2. For each ClippedShape in the cell: a. Start with `contains_center` b. Count edge crossings
///    from cell center to query point c. XOR the crossing parity with contains_center
#[derive(Debug, Clone)]
pub struct QuadtreeCell {
    /// The cell ID (encodes position and level).
    cell_id: QuadtreeCellId,

    /// Shapes intersecting this cell, sorted by doc_id.
    shapes: Vec<ClippedShape>,
}

impl QuadtreeCell {
    /// Creates a new empty cell.
    #[inline]
    pub fn new(cell_id: QuadtreeCellId) -> Self {
        Self {
            cell_id,
            shapes: Vec::new(),
        }
    }

    /// Creates a new cell with preallocated shape capacity.
    #[inline]
    pub fn with_capacity(cell_id: QuadtreeCellId, capacity: usize) -> Self {
        Self {
            cell_id,
            shapes: Vec::with_capacity(capacity),
        }
    }

    /// Returns the cell ID.
    #[inline]
    pub fn cell_id(&self) -> QuadtreeCellId {
        self.cell_id
    }

    /// Returns the number of shapes in this cell.
    #[inline]
    pub fn num_shapes(&self) -> usize {
        self.shapes.len()
    }

    /// Returns the shapes as a slice.
    #[inline]
    pub fn shapes(&self) -> &[ClippedShape] {
        &self.shapes
    }

    /// Returns a mutable reference to the shapes.
    #[inline]
    pub fn shapes_mut(&mut self) -> &mut Vec<ClippedShape> {
        &mut self.shapes
    }

    /// Returns the shape at the given index.
    ///
    /// # Panics
    ///
    /// Panics if `i >= num_shapes()`.
    #[inline]
    pub fn shape(&self, i: usize) -> &ClippedShape {
        &self.shapes[i]
    }

    /// Finds the ClippedShape for the given doc_id.
    ///
    /// Returns `None` if no shape with that doc_id is in this cell.
    /// Uses binary search for O(log n) lookup.
    pub fn find_shape(&self, doc_id: u32) -> Option<&ClippedShape> {
        self.shapes
            .binary_search_by_key(&doc_id, |s| s.doc_id)
            .ok()
            .map(|idx| &self.shapes[idx])
    }

    /// Finds the mutable ClippedShape for the given doc_id.
    pub fn find_shape_mut(&mut self, doc_id: u32) -> Option<&mut ClippedShape> {
        self.shapes
            .binary_search_by_key(&doc_id, |s| s.doc_id)
            .ok()
            .map(|idx| &mut self.shapes[idx])
    }

    /// Adds a shape to the cell, maintaining sorted order by doc_id.
    ///
    /// If a shape with the same doc_id already exists, it is replaced.
    pub fn add_shape(&mut self, shape: ClippedShape) {
        match self
            .shapes
            .binary_search_by_key(&shape.doc_id, |s| s.doc_id)
        {
            Ok(idx) => self.shapes[idx] = shape,        // Replace existing
            Err(idx) => self.shapes.insert(idx, shape), // Insert at sorted position
        }
    }

    /// Gets or creates a ClippedShape for the given doc_id.
    ///
    /// If no shape exists, creates one with `contains_center = false`
    /// and no edges.
    pub fn get_or_create_shape(&mut self, doc_id: u32) -> &mut ClippedShape {
        match self.shapes.binary_search_by_key(&doc_id, |s| s.doc_id) {
            Ok(idx) => &mut self.shapes[idx],
            Err(idx) => {
                self.shapes.insert(idx, ClippedShape::new(doc_id, false));
                &mut self.shapes[idx]
            }
        }
    }

    /// Removes the shape with the given doc_id.
    ///
    /// Returns the removed shape, or `None` if not found.
    pub fn remove_shape(&mut self, doc_id: u32) -> Option<ClippedShape> {
        self.shapes
            .binary_search_by_key(&doc_id, |s| s.doc_id)
            .ok()
            .map(|idx| self.shapes.remove(idx))
    }

    /// Returns the total number of edges across all shapes in this cell.
    pub fn total_edges(&self) -> usize {
        self.shapes.iter().map(|s| s.num_edges()).sum()
    }

    /// Returns true if this cell has no shapes.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.shapes.is_empty()
    }

    /// Clears all shapes from the cell.
    #[inline]
    pub fn clear(&mut self) {
        self.shapes.clear();
    }

    /// Returns an iterator over (doc_id, contains_center) pairs.
    pub fn doc_ids(&self) -> impl Iterator<Item = (u32, bool)> + '_ {
        self.shapes.iter().map(|s| (s.doc_id, s.contains_center))
    }

    /// Sorts shapes by doc_id (call after bulk insertion).
    ///
    /// This is useful when shapes were added without maintaining sort order.
    pub fn sort_shapes(&mut self) {
        self.shapes.sort_by_key(|s| s.doc_id);
    }
}

impl PartialEq for QuadtreeCell {
    fn eq(&self, other: &Self) -> bool {
        self.cell_id == other.cell_id && self.shapes == other.shapes
    }
}

impl Eq for QuadtreeCell {}

// =============================================================================
// PaddedCell - Cell with Padding for Edge Clipping
// =============================================================================

/// A quadtree cell with padding for edge clipping tolerance.
///
/// PaddedCell extends a cell's bounds slightly to handle numerical precision
/// issues when edges lie exactly on cell boundaries. It also provides
/// precomputed values for efficient recursive subdivision during index building.
///
/// # Usage in Index Building
///
/// During index construction, PaddedCell is used to:
/// 1. Clip edges against the padded bounds (ensures no edges are lost)
/// 2. Track entry/exit vertices for InteriorTracker state maintenance
/// 3. Efficiently subdivide into child cells
///
/// # Coordinate Systems
///
/// - `bounds`: The actual cell bounds in world coordinates
/// - `padded`: The padded bounds (expanded by padding amount)
/// - `middle`: The bounds shared by all four children (for subdivision)
///
/// # Differences from S2PaddedCell
///
/// Our 2D version is simpler:
/// - No face projection (single plane instead of cube faces)
/// - Z-order traversal (fixed child ordering, no Hilbert curve orientation)
/// - Simpler entry/exit vertex computation (corners, not sphere intersections)
#[derive(Debug, Clone)]
pub struct PaddedCell {
    /// The cell ID.
    id: QuadtreeCellId,

    /// The padding amount (absolute, in world coordinates).
    padding: f64,

    /// The actual cell bounds (without padding).
    bounds: Rect,

    /// The padded cell bounds (expanded by padding on all sides).
    padded: Rect,

    /// The center of the cell.
    center: Point2D,

    /// The "middle" rectangle shared by all four children.
    /// Computed lazily on first access.
    middle: Option<Rect>,
}

impl PaddedCell {
    /// Creates a PaddedCell from a cell ID and global bounds.
    ///
    /// # Arguments
    ///
    /// * `id` - The cell ID
    /// * `padding` - The padding amount (in world coordinates)
    /// * `global_bounds` - The global coordinate bounds for the quadtree
    pub fn new(id: QuadtreeCellId, padding: f64, global_bounds: &Bounds) -> Self {
        let bounds = id.to_bounds(global_bounds);
        let padded = bounds.expanded(padding);
        let center = bounds.center();

        Self {
            id,
            padding,
            bounds,
            padded,
            center,
            middle: None,
        }
    }

    /// Creates a PaddedCell for the root cell.
    pub fn root(padding: f64, global_bounds: &Bounds) -> Self {
        Self::new(QuadtreeCellId::ROOT, padding, global_bounds)
    }

    /// Creates a child PaddedCell from this parent.
    ///
    /// # Arguments
    ///
    /// * `child_index` - The child index (0-3) in Z-order:
    ///   - 0: bottom-left (min_x, min_y)
    ///   - 1: bottom-right (max_x, min_y)
    ///   - 2: top-left (min_x, max_y)
    ///   - 3: top-right (max_x, max_y)
    pub fn child(&self, child_index: usize) -> Option<Self> {
        debug_assert!(child_index < 4);

        let child_id = self.id.child(child_index)?;

        // Compute child bounds from parent (more efficient than from cell ID)
        let mid_x = self.bounds.x.center();
        let mid_y = self.bounds.y.center();

        let x_interval = if child_index & 1 == 0 {
            Interval::new(self.bounds.x.lo, mid_x)
        } else {
            Interval::new(mid_x, self.bounds.x.hi)
        };

        let y_interval = if child_index & 2 == 0 {
            Interval::new(self.bounds.y.lo, mid_y)
        } else {
            Interval::new(mid_y, self.bounds.y.hi)
        };

        let bounds = Rect::new(x_interval, y_interval);
        let padded = bounds.expanded(self.padding);
        let center = bounds.center();

        Some(Self {
            id: child_id,
            padding: self.padding,
            bounds,
            padded,
            center,
            middle: None,
        })
    }

    /// Returns all four children.
    pub fn children(&self) -> Option<[Self; 4]> {
        Some([
            self.child(0)?,
            self.child(1)?,
            self.child(2)?,
            self.child(3)?,
        ])
    }

    /// Returns the cell ID.
    #[inline]
    pub fn id(&self) -> QuadtreeCellId {
        self.id
    }

    /// Returns the cell level.
    #[inline]
    pub fn level(&self) -> u8 {
        self.id.level()
    }

    /// Returns the padding amount.
    #[inline]
    pub fn padding(&self) -> f64 {
        self.padding
    }

    /// Returns the actual cell bounds (without padding).
    #[inline]
    pub fn bounds(&self) -> &Rect {
        &self.bounds
    }

    /// Returns the padded cell bounds.
    #[inline]
    pub fn padded_bounds(&self) -> &Rect {
        &self.padded
    }

    /// Returns the center of the cell.
    #[inline]
    pub fn center(&self) -> Point2D {
        self.center
    }

    /// Returns the "middle" rectangle shared by all four children.
    ///
    /// This is the region where the parent's edges need to be tested
    /// against all four children (edges outside this region only need
    /// to be tested against a subset of children).
    ///
    /// The middle is computed lazily and cached.
    pub fn middle(&mut self) -> Rect {
        if let Some(middle) = self.middle {
            return middle;
        }

        let mid_x = self.bounds.x.center();
        let mid_y = self.bounds.y.center();

        // The middle is the region within padding distance of both split lines
        let middle = Rect::new(
            Interval::new(mid_x - self.padding, mid_x + self.padding),
            Interval::new(mid_y - self.padding, mid_y + self.padding),
        );

        self.middle = Some(middle);
        middle
    }

    /// Returns the entry vertex for Z-order traversal.
    ///
    /// This is the corner where the space-filling curve enters this cell.
    /// For Z-order (Morton code), this is always the bottom-left corner (0,0).
    #[inline]
    pub fn entry_vertex(&self) -> Point2D {
        // Z-order always enters at the minimum corner
        self.bounds.lo()
    }

    /// Returns the exit vertex for Z-order traversal.
    ///
    /// This is the corner where the space-filling curve exits this cell.
    /// For Z-order (Morton code), this is always the top-right corner (1,1).
    #[inline]
    pub fn exit_vertex(&self) -> Point2D {
        // Z-order always exits at the maximum corner
        self.bounds.hi()
    }

    /// Returns the (i, j) indices for the child at the given Z-order position.
    ///
    /// For Z-order traversal, the positions map directly to indices:
    /// - pos 0 -> (0, 0) -> child index 0
    /// - pos 1 -> (1, 0) -> child index 1
    /// - pos 2 -> (0, 1) -> child index 2
    /// - pos 3 -> (1, 1) -> child index 3
    #[inline]
    pub fn child_ij(&self, pos: usize) -> (usize, usize) {
        debug_assert!(pos < 4);
        (pos & 1, (pos >> 1) & 1)
    }

    /// Returns the child index for the given (i, j) indices.
    #[inline]
    pub fn ij_to_child_index(i: usize, j: usize) -> usize {
        debug_assert!(i < 2 && j < 2);
        i + (j << 1)
    }

    /// Finds the smallest cell that contains all descendants whose bounds
    /// intersect the given rectangle.
    ///
    /// This is useful for skipping initial subdivision steps when only
    /// one branch of the tree needs to be explored.
    ///
    /// # Preconditions
    ///
    /// The rectangle must intersect this cell's padded bounds.
    pub fn shrink_to_fit(&self, rect: &Rect, global_bounds: &Bounds) -> QuadtreeCellId {
        // Start with the current cell
        let mut result_id = self.id;
        let mut current_bounds = self.bounds;

        loop {
            let mid_x = current_bounds.x.center();
            let mid_y = current_bounds.y.center();

            // Determine which children the rect could intersect
            let x_lo = rect.x.lo <= mid_x + self.padding;
            let x_hi = rect.x.hi >= mid_x - self.padding;
            let y_lo = rect.y.lo <= mid_y + self.padding;
            let y_hi = rect.y.hi >= mid_y - self.padding;

            // Count how many quadrants are needed
            let quadrants = (x_lo as u8 + x_hi as u8) * (y_lo as u8 + y_hi as u8);

            if quadrants != 1 {
                // Rect spans multiple children or is in the middle padding zone
                break;
            }

            // Rect is entirely within one child
            let i = if x_hi && !x_lo { 1 } else { 0 };
            let j = if y_hi && !y_lo { 1 } else { 0 };

            let child_id = match result_id.child(Self::ij_to_child_index(i, j)) {
                Some(id) => id,
                None => break, // At max level
            };

            // Update bounds for next iteration
            let x_interval = if i == 0 {
                Interval::new(current_bounds.x.lo, mid_x)
            } else {
                Interval::new(mid_x, current_bounds.x.hi)
            };

            let y_interval = if j == 0 {
                Interval::new(current_bounds.y.lo, mid_y)
            } else {
                Interval::new(mid_y, current_bounds.y.hi)
            };

            current_bounds = Rect::new(x_interval, y_interval);
            result_id = child_id;
        }

        result_id
    }

    /// Tests if a point is contained in the padded bounds.
    #[inline]
    pub fn contains_point(&self, point: &Point2D) -> bool {
        self.padded.contains_point(point)
    }

    /// Tests if an edge (line segment) intersects the padded bounds.
    ///
    /// This is a conservative test - it may return true for edges that
    /// don't actually intersect, but will never return false for edges
    /// that do intersect.
    pub fn intersects_edge(&self, v0: &Point2D, v1: &Point2D) -> bool {
        // Fast rejection: check bounding box of edge against padded bounds
        let edge_bounds = Rect::new(
            Interval::new(v0.x.min(v1.x), v0.x.max(v1.x)),
            Interval::new(v0.y.min(v1.y), v0.y.max(v1.y)),
        );

        if !self.padded.intersects(&edge_bounds) {
            return false;
        }

        // If either endpoint is inside, the edge intersects
        if self.padded.contains_point(v0) || self.padded.contains_point(v1) {
            return true;
        }

        // Check if edge crosses any of the four sides of the padded rectangle
        let corners = [
            self.padded.corner(0), // bottom-left
            self.padded.corner(1), // bottom-right
            self.padded.corner(3), // top-right
            self.padded.corner(2), // top-left
        ];

        // Test edge against each side of the rectangle
        for i in 0..4 {
            let c0 = &corners[i];
            let c1 = &corners[(i + 1) % 4];
            if edge_or_vertex_crossing_2d(v0, v1, c0, c1) {
                return true;
            }
        }

        false
    }
}

// =============================================================================
// Unit Tests
// =============================================================================

#[cfg(test)]
mod tests {
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

    // -------------------------------------------------------------------------
    // QuadtreeCellId Tests
    // -------------------------------------------------------------------------

    mod cell_id_tests {
        use super::*;

        #[test]
        fn test_root_cell() {
            let root = QuadtreeCellId::ROOT;
            assert!(root.is_valid());
            assert_eq!(root.level(), 0);

            // Root should cover entire space
            let bounds = root.to_bounds_normalized();
            assert!((bounds.x.lo - 0.0).abs() < 1e-10);
            assert!((bounds.x.hi - 1.0).abs() < 1e-10);
            assert!((bounds.y.lo - 0.0).abs() < 1e-10);
            assert!((bounds.y.hi - 1.0).abs() < 1e-10);
        }

        #[test]
        fn test_none_cell() {
            let none = QuadtreeCellId::NONE;
            assert!(!none.is_valid());
        }

        #[test]
        fn test_from_ij_level() {
            // Level 1: 2x2 grid
            let cell_00 = QuadtreeCellId::from_ij(0, 0, 1);
            let cell_10 = QuadtreeCellId::from_ij(1, 0, 1);
            let cell_01 = QuadtreeCellId::from_ij(0, 1, 1);
            let cell_11 = QuadtreeCellId::from_ij(1, 1, 1);

            assert!(cell_00.is_valid());
            assert_eq!(cell_00.level(), 1);
            assert_eq!(cell_10.level(), 1);
            assert_eq!(cell_01.level(), 1);
            assert_eq!(cell_11.level(), 1);

            // Should be ordered by Z-order
            assert!(cell_00 < cell_10);
            assert!(cell_10 < cell_01);
            assert!(cell_01 < cell_11);
        }

        #[test]
        fn test_to_ij_roundtrip() {
            for level in 0..=10 {
                let max = 1u32 << level;
                for i in [0, max / 4, max / 2, max - 1].iter().filter(|&&x| x < max) {
                    for j in [0, max / 4, max / 2, max - 1].iter().filter(|&&x| x < max) {
                        let cell = QuadtreeCellId::from_ij(*i, *j, level);
                        let (ri, rj) = cell.to_ij();
                        assert_eq!(
                            (*i, *j),
                            (ri, rj),
                            "Roundtrip failed for level={}, i={}, j={}",
                            level,
                            i,
                            j
                        );
                    }
                }
            }
        }

        #[test]
        fn test_from_normalized() {
            // Corner points
            let bl = QuadtreeCellId::from_normalized(0.0, 0.0, 2);
            assert_eq!(bl.to_ij(), (0, 0));

            // Near top-right (but not exactly 1.0)
            let tr = QuadtreeCellId::from_normalized(0.999, 0.999, 2);
            assert_eq!(tr.to_ij(), (3, 3));

            // Center should be (2, 2) at level 2
            let center = QuadtreeCellId::from_normalized(0.5, 0.5, 2);
            assert_eq!(center.to_ij(), (2, 2));
        }

        #[test]
        fn test_to_center() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);

            // Root cell center should be (50, 50)
            let root_center = QuadtreeCellId::ROOT.to_center(&bounds);
            assert!((root_center.x - 50.0).abs() < 1e-10);
            assert!((root_center.y - 50.0).abs() < 1e-10);

            // Level 1 cell (0,0) should have center at (25, 25)
            let cell = QuadtreeCellId::from_ij(0, 0, 1);
            let center = cell.to_center(&bounds);
            assert!((center.x - 25.0).abs() < 1e-10);
            assert!((center.y - 25.0).abs() < 1e-10);
        }

        #[test]
        fn test_to_bounds() {
            let global = Bounds::new(0.0, 0.0, 100.0, 100.0);

            // Root cell should cover entire bounds
            let root_bounds = QuadtreeCellId::ROOT.to_bounds(&global);
            assert!((root_bounds.x.lo - 0.0).abs() < 1e-10);
            assert!((root_bounds.x.hi - 100.0).abs() < 1e-10);

            // Level 1 cell (0,0) should cover bottom-left quadrant
            let cell = QuadtreeCellId::from_ij(0, 0, 1);
            let cell_bounds = cell.to_bounds(&global);
            assert!((cell_bounds.x.lo - 0.0).abs() < 1e-10);
            assert!((cell_bounds.x.hi - 50.0).abs() < 1e-10);
            assert!((cell_bounds.y.lo - 0.0).abs() < 1e-10);
            assert!((cell_bounds.y.hi - 50.0).abs() < 1e-10);
        }

        #[test]
        fn test_parent() {
            let root = QuadtreeCellId::ROOT;
            assert!(root.parent().is_none());

            let level1 = QuadtreeCellId::from_ij(0, 0, 1);
            let parent = level1.parent().unwrap();
            assert_eq!(parent.level(), 0);

            let level2 = QuadtreeCellId::from_ij(2, 3, 2);
            let parent2 = level2.parent().unwrap();
            assert_eq!(parent2.level(), 1);
            assert_eq!(parent2.to_ij(), (1, 1)); // (2,3) -> (1,1) at parent level
        }

        #[test]
        fn test_children() {
            let root = QuadtreeCellId::ROOT;
            let children = root.children().unwrap();

            assert_eq!(children.len(), 4);
            for child in &children {
                assert_eq!(child.level(), 1);
            }

            // Children should be in Z-order
            assert_eq!(children[0].to_ij(), (0, 0));
            assert_eq!(children[1].to_ij(), (1, 0));
            assert_eq!(children[2].to_ij(), (0, 1));
            assert_eq!(children[3].to_ij(), (1, 1));

            // All children's parent should be root
            for child in &children {
                assert_eq!(child.parent().unwrap(), root);
            }
        }

        #[test]
        fn test_child() {
            let root = QuadtreeCellId::ROOT;

            let c0 = root.child(0).unwrap();
            let c1 = root.child(1).unwrap();
            let c2 = root.child(2).unwrap();
            let c3 = root.child(3).unwrap();

            assert_eq!(c0.to_ij(), (0, 0));
            assert_eq!(c1.to_ij(), (1, 0));
            assert_eq!(c2.to_ij(), (0, 1));
            assert_eq!(c3.to_ij(), (1, 1));
        }

        #[test]
        fn test_max_level_no_children() {
            let max_cell = QuadtreeCellId::from_ij(0, 0, MAX_LEVEL);
            assert!(max_cell.children().is_none());
            assert!(max_cell.child(0).is_none());
        }

        #[test]
        fn test_contains() {
            let root = QuadtreeCellId::ROOT;
            let level1 = QuadtreeCellId::from_ij(0, 0, 1);
            let level2 = QuadtreeCellId::from_ij(0, 0, 2);

            // Root contains everything
            assert!(root.contains(&root));
            assert!(root.contains(&level1));
            assert!(root.contains(&level2));

            // Parent contains child
            assert!(level1.contains(&level2));

            // Child doesn't contain parent
            assert!(!level2.contains(&level1));
            assert!(!level1.contains(&root));

            // Different branches don't contain each other
            let other = QuadtreeCellId::from_ij(1, 1, 1);
            assert!(!level1.contains(&other));
            assert!(!other.contains(&level1));
        }

        #[test]
        fn test_intersects() {
            let cell_a = QuadtreeCellId::from_ij(0, 0, 2);
            let cell_b = QuadtreeCellId::from_ij(0, 0, 3); // child of cell_a
            let cell_c = QuadtreeCellId::from_ij(3, 3, 2); // different cell

            assert!(cell_a.intersects(&cell_b)); // parent-child
            assert!(cell_b.intersects(&cell_a)); // child-parent
            assert!(!cell_a.intersects(&cell_c)); // disjoint
        }

        #[test]
        fn test_range_min_max() {
            let cell = QuadtreeCellId::from_ij(0, 0, 2);

            let min = cell.range_min();
            let max = cell.range_max();

            assert!(min.is_valid());
            assert!(max.is_valid());
            assert_eq!(min.level(), MAX_LEVEL);
            assert_eq!(max.level(), MAX_LEVEL);
            assert!(min < max || min == max);

            // All children should be in range
            if let Some(children) = cell.children() {
                for child in &children {
                    assert!(child.0 >= min.0);
                    assert!(child.0 <= max.0);
                }
            }
        }

        #[test]
        fn test_common_ancestor() {
            let a = QuadtreeCellId::from_ij(0, 0, 4);
            let b = QuadtreeCellId::from_ij(1, 0, 4);

            let ancestor = a.common_ancestor(&b);
            assert!(ancestor.contains(&a));
            assert!(ancestor.contains(&b));

            // Same cell should return itself
            assert_eq!(a.common_ancestor(&a), a);

            // Root should be common ancestor of any two cells in different quadrants
            let c = QuadtreeCellId::from_ij(0, 0, 1);
            let d = QuadtreeCellId::from_ij(1, 1, 1);
            let root_ancestor = c.common_ancestor(&d);
            assert_eq!(root_ancestor, QuadtreeCellId::ROOT);
        }

        #[test]
        fn test_ordering() {
            // Create all cells at level 2 (4x4 grid = 16 cells) in arbitrary order
            let mut cells: Vec<QuadtreeCellId> = (0..4)
                .flat_map(|j| (0..4).map(move |i| QuadtreeCellId::from_ij(i, j, 2)))
                .collect();

            // Sort by cell ID (uses Z-order/Morton code ordering)
            cells.sort();

            // After sorting, cells should be in strictly increasing order
            for i in 1..cells.len() {
                assert!(
                    cells[i - 1] < cells[i],
                    "After sorting: cells[{}] ({:?}, morton={}) should be < cells[{}] ({:?}, morton={})",
                    i - 1,
                    cells[i - 1].to_ij(),
                    cells[i - 1].0 >> 55, // Approximate morton for debugging
                    i,
                    cells[i].to_ij(),
                    cells[i].0 >> 55,
                );
            }

            // Verify the Z-order traversal sequence
            // Z-order interleaves x and y bits: morton = ...y1 x1 y0 x0
            // For level 2 (4x4 grid), the expected sequence is:
            let expected_z_order: [(u32, u32); 16] = [
                (0, 0),
                (1, 0),
                (0, 1),
                (1, 1), // Bottom-left quadrant
                (2, 0),
                (3, 0),
                (2, 1),
                (3, 1), // Bottom-right quadrant
                (0, 2),
                (1, 2),
                (0, 3),
                (1, 3), // Top-left quadrant
                (2, 2),
                (3, 2),
                (2, 3),
                (3, 3), // Top-right quadrant
            ];

            for (idx, cell) in cells.iter().enumerate() {
                assert_eq!(
                    cell.to_ij(),
                    expected_z_order[idx],
                    "Cell at sorted index {} should be {:?}, got {:?}",
                    idx,
                    expected_z_order[idx],
                    cell.to_ij()
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // Bit Interleaving Tests
    // -------------------------------------------------------------------------

    mod interleave_tests {
        use super::*;

        #[test]
        fn test_interleave_basic() {
            // (0, 0) -> 0
            assert_eq!(interleave_bits(0, 0), 0);

            // (1, 0) -> 1 (x in bit 0)
            assert_eq!(interleave_bits(1, 0), 0b01);

            // (0, 1) -> 2 (y in bit 1)
            assert_eq!(interleave_bits(0, 1), 0b10);

            // (1, 1) -> 3
            assert_eq!(interleave_bits(1, 1), 0b11);

            // (2, 0) -> 4 (x bit 1 goes to bit 2)
            assert_eq!(interleave_bits(2, 0), 0b0100);

            // (0, 2) -> 8 (y bit 1 goes to bit 3)
            assert_eq!(interleave_bits(0, 2), 0b1000);
        }

        #[test]
        fn test_deinterleave_roundtrip() {
            for level in 1..=10 {
                let max = 1u32 << level;
                for x in [0, 1, max / 2, max - 1] {
                    for y in [0, 1, max / 2, max - 1] {
                        let z = interleave_bits(x, y);
                        let (dx, dy) = deinterleave_bits(z, level);
                        assert_eq!(
                            (x, y),
                            (dx, dy),
                            "Roundtrip failed for level={}, x={}, y={}",
                            level,
                            x,
                            y
                        );
                    }
                }
            }
        }
    }
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
#[cfg(test)]
mod index_structure_tests {
    use super::*;

    // -------------------------------------------------------------------------
    // ClippedShape Tests
    // -------------------------------------------------------------------------

    mod clipped_shape_tests {
        use super::*;

        #[test]
        fn test_new() {
            let shape = ClippedShape::new(42, true);
            assert_eq!(shape.doc_id(), 42);
            assert!(shape.contains_center());
            assert_eq!(shape.num_edges(), 0);
            assert!(shape.is_empty());
        }

        #[test]
        fn test_with_capacity() {
            let shape = ClippedShape::with_capacity(10, false, 8);
            assert_eq!(shape.doc_id(), 10);
            assert!(!shape.contains_center());
            assert_eq!(shape.num_edges(), 0);
        }

        #[test]
        fn test_add_edge_maintains_order() {
            let mut shape = ClippedShape::new(1, false);

            shape.add_edge(5);
            shape.add_edge(2);
            shape.add_edge(8);
            shape.add_edge(1);

            assert_eq!(shape.num_edges(), 4);
            assert_eq!(shape.edges(), &[1, 2, 5, 8]);
        }

        #[test]
        fn test_add_edge_no_duplicates() {
            let mut shape = ClippedShape::new(1, false);

            shape.add_edge(3);
            shape.add_edge(3); // Duplicate
            shape.add_edge(5);
            shape.add_edge(3); // Duplicate

            assert_eq!(shape.num_edges(), 2);
            assert_eq!(shape.edges(), &[3, 5]);
        }

        #[test]
        fn test_add_edges_bulk() {
            let mut shape = ClippedShape::new(1, false);

            shape.add_edges(&[5, 2, 8, 1, 5]); // 5 appears twice

            assert_eq!(shape.num_edges(), 4);
            assert_eq!(shape.edges(), &[1, 2, 5, 8]);
        }

        #[test]
        fn test_contains_edge() {
            let mut shape = ClippedShape::new(1, false);
            shape.add_edges(&[2, 4, 6, 8]);

            assert!(shape.contains_edge(2));
            assert!(shape.contains_edge(4));
            assert!(shape.contains_edge(6));
            assert!(shape.contains_edge(8));
            assert!(!shape.contains_edge(1));
            assert!(!shape.contains_edge(3));
            assert!(!shape.contains_edge(5));
        }

        #[test]
        fn test_set_contains_center() {
            let mut shape = ClippedShape::new(1, false);
            assert!(!shape.contains_center());

            shape.set_contains_center(true);
            assert!(shape.contains_center());

            shape.set_contains_center(false);
            assert!(!shape.contains_center());
        }

        #[test]
        fn test_clear_edges() {
            let mut shape = ClippedShape::new(1, true);
            shape.add_edges(&[1, 2, 3]);

            assert_eq!(shape.num_edges(), 3);
            assert!(shape.contains_center());

            shape.clear_edges();

            assert_eq!(shape.num_edges(), 0);
            assert!(shape.contains_center()); // Preserved
        }

        #[test]
        fn test_equality() {
            let mut shape1 = ClippedShape::new(1, true);
            shape1.add_edges(&[1, 2, 3]);

            let mut shape2 = ClippedShape::new(1, true);
            shape2.add_edges(&[1, 2, 3]);

            let mut shape3 = ClippedShape::new(1, true);
            shape3.add_edges(&[1, 2, 4]);

            assert_eq!(shape1, shape2);
            assert_ne!(shape1, shape3);
        }

        #[test]
        fn test_inline_vs_heap() {
            // SmallVec<[u16; 2]> stores up to 2 elements inline
            let mut shape = ClippedShape::new(1, false);

            // These should be inline
            shape.add_edge(1);
            shape.add_edge(2);
            assert_eq!(shape.num_edges(), 2);

            // This forces heap allocation
            shape.add_edge(3);
            assert_eq!(shape.num_edges(), 3);
            assert_eq!(shape.edges(), &[1, 2, 3]);
        }
    }

    // -------------------------------------------------------------------------
    // QuadtreeCell Tests
    // -------------------------------------------------------------------------

    mod quadtree_cell_tests {
        use super::*;

        #[test]
        fn test_new() {
            let cell = QuadtreeCell::new(QuadtreeCellId::ROOT);
            assert_eq!(cell.cell_id(), QuadtreeCellId::ROOT);
            assert_eq!(cell.num_shapes(), 0);
            assert!(cell.is_empty());
        }

        #[test]
        fn test_add_shape_maintains_order() {
            let mut cell = QuadtreeCell::new(QuadtreeCellId::ROOT);

            cell.add_shape(ClippedShape::new(30, false));
            cell.add_shape(ClippedShape::new(10, true));
            cell.add_shape(ClippedShape::new(50, false));
            cell.add_shape(ClippedShape::new(20, true));

            assert_eq!(cell.num_shapes(), 4);

            let doc_ids: Vec<u32> = cell.shapes().iter().map(|s| s.doc_id()).collect();
            assert_eq!(doc_ids, vec![10, 20, 30, 50]);
        }

        #[test]
        fn test_add_shape_replaces_existing() {
            let mut cell = QuadtreeCell::new(QuadtreeCellId::ROOT);

            let mut shape1 = ClippedShape::new(10, false);
            shape1.add_edge(1);
            cell.add_shape(shape1);

            assert_eq!(cell.num_shapes(), 1);
            assert!(!cell.shapes()[0].contains_center());

            // Replace with new shape having same doc_id
            let mut shape2 = ClippedShape::new(10, true);
            shape2.add_edge(2);
            cell.add_shape(shape2);

            assert_eq!(cell.num_shapes(), 1);
            assert!(cell.shapes()[0].contains_center());
            assert_eq!(cell.shapes()[0].edges(), &[2]);
        }

        #[test]
        fn test_find_shape() {
            let mut cell = QuadtreeCell::new(QuadtreeCellId::ROOT);

            cell.add_shape(ClippedShape::new(10, true));
            cell.add_shape(ClippedShape::new(20, false));
            cell.add_shape(ClippedShape::new(30, true));

            let found = cell.find_shape(20);
            assert!(found.is_some());
            assert_eq!(found.unwrap().doc_id(), 20);
            assert!(!found.unwrap().contains_center());

            assert!(cell.find_shape(15).is_none());
            assert!(cell.find_shape(25).is_none());
        }

        #[test]
        fn test_get_or_create_shape() {
            let mut cell = QuadtreeCell::new(QuadtreeCellId::ROOT);

            // Create new
            let shape = cell.get_or_create_shape(10);
            assert_eq!(shape.doc_id(), 10);
            assert!(!shape.contains_center());
            shape.set_contains_center(true);

            // Get existing
            let shape_again = cell.get_or_create_shape(10);
            assert!(shape_again.contains_center());

            assert_eq!(cell.num_shapes(), 1);
        }

        #[test]
        fn test_remove_shape() {
            let mut cell = QuadtreeCell::new(QuadtreeCellId::ROOT);

            cell.add_shape(ClippedShape::new(10, true));
            cell.add_shape(ClippedShape::new(20, false));
            cell.add_shape(ClippedShape::new(30, true));

            let removed = cell.remove_shape(20);
            assert!(removed.is_some());
            assert_eq!(removed.unwrap().doc_id(), 20);

            assert_eq!(cell.num_shapes(), 2);
            assert!(cell.find_shape(20).is_none());

            // Remove non-existent
            let not_found = cell.remove_shape(25);
            assert!(not_found.is_none());
        }

        #[test]
        fn test_total_edges() {
            let mut cell = QuadtreeCell::new(QuadtreeCellId::ROOT);

            let mut shape1 = ClippedShape::new(10, false);
            shape1.add_edges(&[1, 2, 3]);
            cell.add_shape(shape1);

            let mut shape2 = ClippedShape::new(20, false);
            shape2.add_edges(&[4, 5]);
            cell.add_shape(shape2);

            let shape3 = ClippedShape::new(30, true); // No edges
            cell.add_shape(shape3);

            assert_eq!(cell.total_edges(), 5);
        }

        #[test]
        fn test_doc_ids_iterator() {
            let mut cell = QuadtreeCell::new(QuadtreeCellId::ROOT);

            cell.add_shape(ClippedShape::new(10, true));
            cell.add_shape(ClippedShape::new(20, false));
            cell.add_shape(ClippedShape::new(30, true));

            let doc_ids: Vec<(u32, bool)> = cell.doc_ids().collect();
            assert_eq!(doc_ids, vec![(10, true), (20, false), (30, true)]);
        }

        #[test]
        fn test_clear() {
            let mut cell = QuadtreeCell::new(QuadtreeCellId::ROOT);

            cell.add_shape(ClippedShape::new(10, true));
            cell.add_shape(ClippedShape::new(20, false));

            assert_eq!(cell.num_shapes(), 2);

            cell.clear();

            assert_eq!(cell.num_shapes(), 0);
            assert!(cell.is_empty());
        }
    }

    // -------------------------------------------------------------------------
    // PaddedCell Tests
    // -------------------------------------------------------------------------

    mod padded_cell_tests {
        use super::*;

        fn test_bounds() -> Bounds {
            Bounds::new(0.0, 0.0, 100.0, 100.0)
        }

        #[test]
        fn test_new_root() {
            let bounds = test_bounds();
            let cell = PaddedCell::new(QuadtreeCellId::ROOT, 0.5, &bounds);

            assert_eq!(cell.id(), QuadtreeCellId::ROOT);
            assert_eq!(cell.level(), 0);
            assert_eq!(cell.padding(), 0.5);

            // Root bounds should cover entire space
            assert!((cell.bounds().x.lo - 0.0).abs() < 1e-10);
            assert!((cell.bounds().x.hi - 100.0).abs() < 1e-10);
            assert!((cell.bounds().y.lo - 0.0).abs() < 1e-10);
            assert!((cell.bounds().y.hi - 100.0).abs() < 1e-10);

            // Padded bounds should be expanded
            assert!((cell.padded_bounds().x.lo - (-0.5)).abs() < 1e-10);
            assert!((cell.padded_bounds().x.hi - 100.5).abs() < 1e-10);

            // Center
            assert!((cell.center().x - 50.0).abs() < 1e-10);
            assert!((cell.center().y - 50.0).abs() < 1e-10);
        }

        #[test]
        fn test_child() {
            let bounds = test_bounds();
            let root = PaddedCell::new(QuadtreeCellId::ROOT, 0.5, &bounds);

            // Child 0: bottom-left
            let child0 = root.child(0).unwrap();
            assert_eq!(child0.level(), 1);
            assert!((child0.bounds().x.lo - 0.0).abs() < 1e-10);
            assert!((child0.bounds().x.hi - 50.0).abs() < 1e-10);
            assert!((child0.bounds().y.lo - 0.0).abs() < 1e-10);
            assert!((child0.bounds().y.hi - 50.0).abs() < 1e-10);

            // Child 1: bottom-right
            let child1 = root.child(1).unwrap();
            assert!((child1.bounds().x.lo - 50.0).abs() < 1e-10);
            assert!((child1.bounds().x.hi - 100.0).abs() < 1e-10);
            assert!((child1.bounds().y.lo - 0.0).abs() < 1e-10);
            assert!((child1.bounds().y.hi - 50.0).abs() < 1e-10);

            // Child 2: top-left
            let child2 = root.child(2).unwrap();
            assert!((child2.bounds().x.lo - 0.0).abs() < 1e-10);
            assert!((child2.bounds().x.hi - 50.0).abs() < 1e-10);
            assert!((child2.bounds().y.lo - 50.0).abs() < 1e-10);
            assert!((child2.bounds().y.hi - 100.0).abs() < 1e-10);

            // Child 3: top-right
            let child3 = root.child(3).unwrap();
            assert!((child3.bounds().x.lo - 50.0).abs() < 1e-10);
            assert!((child3.bounds().x.hi - 100.0).abs() < 1e-10);
            assert!((child3.bounds().y.lo - 50.0).abs() < 1e-10);
            assert!((child3.bounds().y.hi - 100.0).abs() < 1e-10);
        }

        #[test]
        fn test_children() {
            let bounds = test_bounds();
            let root = PaddedCell::new(QuadtreeCellId::ROOT, 0.5, &bounds);

            let children = root.children().unwrap();
            assert_eq!(children.len(), 4);

            for (i, child) in children.iter().enumerate() {
                assert_eq!(child.level(), 1);
                assert_eq!(child.id(), root.id().child(i).unwrap());
            }
        }

        #[test]
        fn test_entry_exit_vertices() {
            let bounds = test_bounds();
            let root = PaddedCell::new(QuadtreeCellId::ROOT, 0.5, &bounds);

            let entry = root.entry_vertex();
            let exit = root.exit_vertex();

            // Z-order: entry at (0,0), exit at (1,1)
            assert_eq!(entry, Point2D::new(0.0, 0.0));
            assert_eq!(exit, Point2D::new(100.0, 100.0));
        }

        #[test]
        fn test_middle() {
            let bounds = test_bounds();
            let mut root = PaddedCell::new(QuadtreeCellId::ROOT, 0.5, &bounds);

            let middle = root.middle();

            // Middle should be the region within padding of the center
            assert!((middle.x.lo - (50.0 - 0.5)).abs() < 1e-10);
            assert!((middle.x.hi - (50.0 + 0.5)).abs() < 1e-10);
            assert!((middle.y.lo - (50.0 - 0.5)).abs() < 1e-10);
            assert!((middle.y.hi - (50.0 + 0.5)).abs() < 1e-10);

            // Should be cached (same result)
            let middle2 = root.middle();
            assert_eq!(middle, middle2);
        }

        #[test]
        fn test_child_ij() {
            let bounds = test_bounds();
            let cell = PaddedCell::new(QuadtreeCellId::ROOT, 0.5, &bounds);

            assert_eq!(cell.child_ij(0), (0, 0));
            assert_eq!(cell.child_ij(1), (1, 0));
            assert_eq!(cell.child_ij(2), (0, 1));
            assert_eq!(cell.child_ij(3), (1, 1));
        }

        #[test]
        fn test_ij_to_child_index() {
            assert_eq!(PaddedCell::ij_to_child_index(0, 0), 0);
            assert_eq!(PaddedCell::ij_to_child_index(1, 0), 1);
            assert_eq!(PaddedCell::ij_to_child_index(0, 1), 2);
            assert_eq!(PaddedCell::ij_to_child_index(1, 1), 3);
        }

        #[test]
        fn test_contains_point() {
            let bounds = test_bounds();
            let cell = PaddedCell::new(QuadtreeCellId::ROOT, 0.5, &bounds);

            // Inside bounds
            assert!(cell.contains_point(&Point2D::new(50.0, 50.0)));

            // Inside padding but outside bounds
            assert!(cell.contains_point(&Point2D::new(-0.25, 50.0)));

            // Outside padding
            assert!(!cell.contains_point(&Point2D::new(-1.0, 50.0)));
        }

        #[test]
        fn test_intersects_edge() {
            let bounds = test_bounds();
            let cell = PaddedCell::new(QuadtreeCellId::ROOT, 0.5, &bounds);

            // Edge entirely inside
            let v0_inside = Point2D::new(10.0, 10.0);
            let v1_inside = Point2D::new(90.0, 90.0);
            assert!(cell.intersects_edge(&v0_inside, &v1_inside));

            // Edge crossing boundary
            let v0_cross = Point2D::new(-10.0, 50.0);
            let v1_cross = Point2D::new(10.0, 50.0);
            assert!(cell.intersects_edge(&v0_cross, &v1_cross));

            // Edge entirely outside
            let v0_outside = Point2D::new(-10.0, -10.0);
            let v1_outside = Point2D::new(-5.0, -5.0);
            assert!(!cell.intersects_edge(&v0_outside, &v1_outside));
        }

        #[test]
        fn test_shrink_to_fit_single_child() {
            let bounds = test_bounds();
            let root = PaddedCell::new(QuadtreeCellId::ROOT, 0.5, &bounds);

            // Rect entirely in bottom-left quadrant
            let rect = Rect::from_coords(5.0, 5.0, 20.0, 20.0);
            let shrunk = root.shrink_to_fit(&rect, &bounds);

            // Should descend to the bottom-left child
            assert!(shrunk.level() >= 1);

            // The shrunk cell should contain the rect
            let shrunk_bounds = shrunk.to_bounds(&bounds);
            assert!(shrunk_bounds.contains_rect(&rect));
        }

        #[test]
        fn test_shrink_to_fit_spanning_multiple() {
            let bounds = test_bounds();
            let root = PaddedCell::new(QuadtreeCellId::ROOT, 0.5, &bounds);

            // Rect spanning multiple quadrants
            let rect = Rect::from_coords(40.0, 40.0, 60.0, 60.0);
            let shrunk = root.shrink_to_fit(&rect, &bounds);

            // Should stay at root because rect spans all quadrants
            assert_eq!(shrunk, QuadtreeCellId::ROOT);
        }

        #[test]
        fn test_padding_consistency() {
            let bounds = test_bounds();
            let root = PaddedCell::new(QuadtreeCellId::ROOT, 1.0, &bounds);

            // All children should have the same padding
            let children = root.children().unwrap();
            for child in &children {
                assert_eq!(child.padding(), 1.0);

                let grandchildren = child.children().unwrap();
                for gc in &grandchildren {
                    assert_eq!(gc.padding(), 1.0);
                }
            }
        }
    }
}
