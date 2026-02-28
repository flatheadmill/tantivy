// =============================================================================
// Bounds - Global Coordinate System Bounds
// =============================================================================

use crate::spatial::quadtree::point::Point2D;

/// Axis-aligned bounding box representing the global coordinate space.
///
/// All points indexed by the quadtree must fall within these bounds.
/// The bounds are used to normalize coordinates to [0, 1] for cell ID encoding.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(missing_docs)]
pub struct Bounds {
    pub lo: Point2D,
    pub hi: Point2D,
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
              lo: Point2D { x: min_x, y: min_y },
              hi: Point2D { x: max_x, y: max_y },
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
          x: if index & 1 == 0 { self.lo.x } else { self.hi.x },
          y: if index & 2 == 0 { self.lo.y } else { self.hi.y },
      }
  }
}
