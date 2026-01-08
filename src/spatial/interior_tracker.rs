//! HUSH

use crate::spatial::surface::Surface;

#[derive(Debug)]
pub struct InteriorTracker<S: Surface> {
    origin: S::Point,

    /// Previous focus point (start of current DrawTo segment).
    a: S::Point,

    /// Current focus point.
    b: S::Point,

    /// The next expected cell ID (for optimization).
    /// When the caller moves to this cell's entry vertex, we can skip MoveTo.
    next_cellid: S::CellId,

    /// Edge crosser for the current A->B segment.
    crosser: S::Crosser,

    /// Set of shape IDs whose interiors contain the current focus.
    shape_ids: Vec<u32>,

    /// Whether any shapes are being tracked.
    is_active: bool,
}

impl<S: Surface> InteriorTracker<S> {
    pub fn new(surface: &S) -> Self {
        let origin = surface.origin();
        let next_cellid = S::cell_begin(S::MAX_LEVEL);

        Self {
            origin,
            a: origin,
            b: origin,
            next_cellid,
            crosser: S::new_crosser(&origin, &origin),
            shape_ids: Vec::new(),
            is_active: false,
        }
    }

    /// Returns the origin point (start of Z-order curve).
    #[inline]
    pub fn origin(&self) -> S::Point {
        self.origin
    }

    /// Returns the current focus point.
    #[inline]
    pub fn focus(&self) -> S::Point {
        self.b
    }

    /// Returns true if any shapes are being tracked.
    #[inline]
    pub fn is_active(&self) -> bool {
        self.is_active
    }

    /// Adds a shape whose interior should be tracked.
    ///
    /// # Arguments
    ///
    /// * `shape_id` - The document/shape ID
    /// * `contains_focus` - Whether the current focus point is inside this shape
    ///
    /// If the focus is being moved via DrawTo, you can also call this with
    /// `contains_focus` at the OLD focus position, then call TestEdge for
    /// edges that might cross the DrawTo line.
    pub fn add_shape(&mut self, shape_id: u32, contains_focus: bool) {
        self.is_active = true;
        if contains_focus {
            self.toggle_shape(shape_id);
        }
    }

    #[inline]
    pub fn shape_ids(&self) -> &[u32] {
        &self.shape_ids
    }

    /// Moves the focus to the given point without testing edge crossings.
    ///
    /// Use this when you know there are no edges between the old and new focus.
    /// For example, moving within a cell that has no edges.
    #[inline]
    pub fn move_to(&mut self, point: &S::Point) {
        self.b = *point;
    }

    /// Prepares to move the focus to the given point.
    ///
    /// After calling this, call `test_edge()` for every edge that might cross
    /// the line segment from the old focus to the new focus.
    pub fn draw_to(&mut self, point: &S::Point) {
        self.a = self.b;
        self.b = *point;
        self.crosser = S::new_crosser(&self.a, &self.b);
    }

    /// Tests whether an edge crosses the current DrawTo segment.
    ///
    /// If the edge crosses, the containment state for `shape_id` is toggled.
    ///
    /// # Arguments
    ///
    /// * `shape_id` - The shape that owns this edge
    /// * `v0` - First vertex of the edge
    /// * `v1` - Second vertex of the edge
    #[inline]
    pub fn test_edge(&mut self, shape_id: u32, v0: &S::Point, v1: &S::Point) {
      if S::test_crossing(&mut self.crosser, v0, v1) {
          self.toggle_shape(shape_id);
      }
    }

    /// Toggles the containment state for a shape.
    ///
    /// If the shape is currently in the set, it is removed.
    /// If the shape is not in the set, it is added.
    #[inline]
     fn toggle_shape(&mut self, shape_id: u32) {
      if let Some(pos) = self.shape_ids.iter().position(|&id| id == shape_id) {
          self.shape_ids.remove(pos);
      } else {
          let pos = self.shape_ids.iter().position(|&id| id > shape_id).unwrap_or(self.shape_ids.len());
          self.shape_ids.insert(pos, shape_id);
      }
  }

    /// Indicates that the focus is positioned at the entry vertex of the given cell.
    ///
    /// Call this after processing each cell. When traversing in Z-order, the
    /// exit vertex of one cell is often the entry vertex of the next cell,
    /// allowing us to skip the MoveTo call.
    pub fn set_next_cellid(&mut self, next_cellid: S::CellId) {
        self.next_cellid = S::cell_range_min(next_cellid);
    }

    /// Returns true if the focus is already at the entry vertex of the given cell.
    ///
    /// This optimization avoids redundant DrawTo/TestEdge calls when traversing
    /// cells in Z-order, since adjacent cells often share vertices.
    #[inline]
    pub fn at_cellid(&self, cellid: S::CellId) -> bool {
        S::cell_range_min(cellid) == self.next_cellid
    }
}
