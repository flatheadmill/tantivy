// =============================================================================
// PaddedQuadtreeCell - Cell with Padding for Edge Clipping
// =============================================================================

use crate::spatial::{quadtree::{edge_or_vertex_crossing_2d, interval::Interval, point::Point2D, rect::Rect, QuadtreeCellId}, PaddedCell};

/// A quadtree cell with padding for edge clipping tolerance.
///
/// PaddedQuadtreeCell extends a cell's bounds slightly to handle numerical precision
/// issues when edges lie exactly on cell boundaries. It also provides
/// precomputed values for efficient recursive subdivision during index building.
///
/// # Usage in Index Building
///
/// During index construction, PaddedQuadtreeCell is used to:
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
/// # Differences from S2PaddedQuadtreeCell
///
/// Our 2D version is simpler:
/// - No face projection (single plane instead of cube faces)
/// - Z-order traversal (fixed child ordering, no Hilbert curve orientation)
/// - Simpler entry/exit vertex computation (corners, not sphere intersections)
#[derive(Debug, Clone)]
pub struct PaddedQuadtreeCell {
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

impl PaddedCell for PaddedQuadtreeCell {
    type Point = Point2D;
    type CellId = QuadtreeCellId;
    type Rect = Rect;

    /// Returns the cell ID.
    #[inline]
    fn id(&self) -> QuadtreeCellId {
        self.id
    }


    /// Returns the cell level.
    #[inline]
    fn level(&self) -> u8 {
        self.id.level()
    }

    /// Returns the actual cell bounds (without padding).
    #[inline]
    fn bounds(&self) -> Rect {
        self.bounds
    }

    /// Returns the padded cell bounds.
    #[inline]
    fn padded_bounds(&self) -> Rect {
        self.padded
    }

    /// Returns the center of the cell.
    #[inline]
    fn center(&self) -> Point2D {
        self.center
    }

    /// Creates a child PaddedQuadtreeCell from this parent.
    fn child(&self, child_index: usize) -> Option<Self> {
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
    fn children(&self) -> Option<[Self; 4]> {
        Some([
            self.child(0)?,
            self.child(1)?,
            self.child(2)?,
            self.child(3)?,
        ])
    }
    /// Returns the "middle" rectangle shared by all four children.
    ///
    /// This is the region where the parent's edges need to be tested
    /// against all four children (edges outside this region only need
    /// to be tested against a subset of children).
    ///
    /// The middle is computed lazily and cached.
    fn middle(&mut self) -> Option<Rect> {
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
        Some(middle)
    }

    /// Tests if a point is contained in the padded bounds.
    #[inline]
    fn contains_point(&self, point: &Point2D) -> bool {
        self.padded.contains_point(point)
    }

    /// Returns the entry vertex for Z-order traversal.
    ///
    /// This is the corner where the space-filling curve enters this cell.
    /// For Z-order (Morton code), this is always the bottom-left corner (0,0).
    #[inline]
    fn entry_vertex(&self) -> Point2D {
        // Z-order always enters at the minimum corner
        self.bounds.lo()
    }

    /// Returns the exit vertex for Z-order traversal.
    ///
    /// This is the corner where the space-filling curve exits this cell.
    /// For Z-order (Morton code), this is always the top-right corner (1,1).
    #[inline]
    fn exit_vertex(&self) -> Point2D {
        // Z-order always exits at the maximum corner
        self.bounds.hi()
    }
    #[inline]
  fn curve_pos_to_child(&self, pos: usize) -> usize {
      pos
  }

    #[inline]
  fn child_to_curve_pos(&self, child: usize) -> usize {
      child
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
    fn shrink_to_fit(&self, bounds: &Rect) -> QuadtreeCellId {
        // Start with the current cell
        let mut result_id = self.id;
        let mut current_bounds = self.bounds;

        loop {
            let mid_x = current_bounds.x.center();
            let mid_y = current_bounds.y.center();

            // Determine which children the rect could intersect
            let x_lo = bounds.x.lo <= mid_x + self.padding;
            let x_hi = bounds.x.hi >= mid_x - self.padding;
            let y_lo = bounds.y.lo <= mid_y + self.padding;
            let y_hi = bounds.y.hi >= mid_y - self.padding;

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
}

impl PaddedQuadtreeCell {
    /// Creates a PaddedQuadtreeCell from a cell ID and global bounds.
    ///
    /// # Arguments
    ///
    /// * `id` - The cell ID
    /// * `padding` - The padding amount (in world coordinates)
    /// * `global_bounds` - The global coordinate bounds for the quadtree
    pub fn new(id: QuadtreeCellId, padding: f64, global_bounds: &Rect) -> Self {
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

    /// Creates a PaddedQuadtreeCell for the root cell.
    pub fn root(padding: f64, global_bounds: &Rect) -> Self {
        Self::new(QuadtreeCellId::ROOT, padding, global_bounds)
    }


    /// Returns the padding amount.
    #[inline]
    fn padding(&self) -> f64 {
        self.padding
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
