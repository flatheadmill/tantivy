//! HUSH

use crate::spatial::quadtree::{
    brute_force_contains_2d, Rect, EdgeCrosser2D, Point2D, QuadtreeCellId,
};

#[derive(Debug)]
pub struct InteriorTracker {
    /// The origin point (start of Z-order curve, outside all shapes).
    origin: Point2D,

    /// Previous focus point (start of current DrawTo segment).
    a: Point2D,

    /// Current focus point.
    b: Point2D,

    /// The next expected cell ID (for optimization).
    /// When the caller moves to this cell's entry vertex, we can skip MoveTo.
    next_cellid: QuadtreeCellId,

    /// Edge crosser for the current A->B segment.
    crosser: EdgeCrosser2D,

    /// Set of shape IDs whose interiors contain the current focus.
    shape_ids: Vec<u32>,

    /// Whether any shapes are being tracked.
    is_active: bool,

    /// Saved state for incremental updates.
    saved_ids: Vec<u32>,
}

impl InteriorTracker {
    /// Creates a new InteriorTracker with the origin at the bottom-left corner.
    ///
    /// The origin is placed at (bounds.min_x - 1, bounds.min_y - 1) to ensure
    /// it is outside all shapes indexed within the bounds.
    pub fn new(bounds: &Rect) -> Self {
        let origin = Point2D::new(bounds.min_x - 1.0, bounds.min_y - 1.0);
        let next_cellid = QuadtreeCellId::begin(QuadtreeCellId::MAX_LEVEL);

        Self {
            origin,
            a: origin,
            b: origin,
            next_cellid,
            crosser: EdgeCrosser2D::new(&origin, &origin),
            shape_ids: Vec::new(),
            is_active: false,
            saved_ids: Vec::new(),
        }
    }

    /// Returns the origin point (start of Z-order curve).
    #[inline]
    pub fn origin(&self) -> Point2D {
        self.origin
    }

    /// Returns the current focus point.
    #[inline]
    pub fn focus(&self) -> Point2D {
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
    pub fn move_to(&mut self, point: &Point2D) {
        self.b = *point;
    }

    /// Prepares to move the focus to the given point.
    ///
    /// After calling this, call `test_edge()` for every edge that might cross
    /// the line segment from the old focus to the new focus.
    pub fn draw_to(&mut self, point: &Point2D) {
        self.a = self.b;
        self.b = *point;
        self.crosser = EdgeCrosser2D::new(&self.a, &self.b);
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
    pub fn test_edge(&mut self, shape_id: u32, v0: &Point2D, v1: &Point2D) {
        if self.crosser.edge_or_vertex_crossing(v0, v1) {
            self.toggle_shape(shape_id);
        }
    }

    /// Toggles the containment state for a shape.
    ///
    /// If the shape is currently in the set, it is removed.
    /// If the shape is not in the set, it is added.
    #[inline]
    pub fn toggle_shape(&mut self, shape_id: u32) {
        self.shape_ids.toggle(shape_id);
    }

    /// Indicates that the focus is positioned at the entry vertex of the given cell.
    ///
    /// Call this after processing each cell. When traversing in Z-order, the
    /// exit vertex of one cell is often the entry vertex of the next cell,
    /// allowing us to skip the MoveTo call.
    pub fn set_next_cellid(&mut self, next_cellid: QuadtreeCellId) {
        self.next_cellid = next_cellid.range_min();
    }

    /// Returns true if the focus is already at the entry vertex of the given cell.
    ///
    /// This optimization avoids redundant DrawTo/TestEdge calls when traversing
    /// cells in Z-order, since adjacent cells often share vertices.
    #[inline]
    pub fn at_cellid(&self, cellid: QuadtreeCellId) -> bool {
        cellid.range_min() == self.next_cellid
    }
    /// Sets the partial shape ID.
    #[inline]
    pub fn set_partial_shape_id(&mut self, shape_id: i32) {
        self.partial_shape_id = shape_id;
    }
}

// =============================================================================
// Helper: Determine Initial Containment at Origin
// =============================================================================

/// Computes whether a shape's interior contains the tracker's origin point.
///
/// This is used when adding a shape to the InteriorTracker. For most shapes,
/// the origin (bottom-left corner of bounds) will be outside, so this typically
/// returns false.
///
/// # Arguments
///
/// * `origin` - The InteriorTracker's origin point
/// * `vertices` - The shape's vertices in order
///
/// # Returns
///
/// `true` if the origin is inside the shape's interior
pub fn contains_tracker_origin(origin: &Point2D, vertices: &[Point2D]) -> bool {
    if vertices.len() < 3 {
        return false;
    }

    // Use brute force contains check
    brute_force_contains_2d(origin, vertices)
}
