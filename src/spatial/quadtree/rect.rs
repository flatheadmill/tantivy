use crate::spatial::quadtree::{interval::Interval, point::Point2D};

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
