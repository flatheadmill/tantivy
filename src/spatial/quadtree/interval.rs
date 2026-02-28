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
