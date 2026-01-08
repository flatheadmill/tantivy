//! HUSH

use crate::spatial::quadtree::{
    brute_force_contains_2d, Bounds, EdgeCrosser2D, Point2D, QuadtreeCellId,
};

/// An efficient set of shape IDs optimized for very small cardinality.
///
/// Following S2's observation, the number of shapes containing a given point
/// is typically very small (0, 1, or 2). Using a sorted Vec and linear search
/// is significantly faster than a HashSet for this case.
///
/// # Performance
///
/// - Lookup: O(n) linear search, but n is typically <= 2
/// - Insert: O(n) to maintain sorted order
/// - Toggle: O(n) for both lookup and potential insert/remove
///
/// For larger sets, consider upgrading to BTreeSet, but in practice the
/// number of overlapping polygons at any point is bounded.
#[derive(Debug, Clone, Default)]
pub struct ShapeIdSet {
    /// Shape IDs in sorted order.
    ids: Vec<u32>,
}

impl ShapeIdSet {
    /// Creates an empty set.
    #[inline]
    pub fn new() -> Self {
        Self { ids: Vec::new() }
    }

    /// Creates a set with preallocated capacity.
    #[inline]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            ids: Vec::with_capacity(capacity),
        }
    }

    /// Returns the number of shape IDs in the set.
    #[inline]
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// Returns true if the set is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    /// Returns the shape IDs as a slice (in sorted order).
    #[inline]
    pub fn as_slice(&self) -> &[u32] {
        &self.ids
    }

    /// Returns true if the set contains the given shape ID.
    #[inline]
    pub fn contains(&self, shape_id: u32) -> bool {
        // Linear search is faster than binary search for small n
        self.ids.iter().any(|&id| id == shape_id)
    }

    /// Inserts a shape ID into the set, maintaining sorted order.
    ///
    /// Returns true if the shape was newly inserted, false if already present.
    pub fn insert(&mut self, shape_id: u32) -> bool {
        // Fast path for empty set
        if self.ids.is_empty() {
            self.ids.push(shape_id);
            return true;
        }

        // Fast path for first element
        if self.ids[0] == shape_id {
            return false;
        }

        // Linear search for insertion point
        let mut pos = 0;
        while pos < self.ids.len() && self.ids[pos] < shape_id {
            pos += 1;
        }

        if pos < self.ids.len() && self.ids[pos] == shape_id {
            return false; // Already present
        }

        self.ids.insert(pos, shape_id);
        true
    }

    /// Removes a shape ID from the set.
    ///
    /// Returns true if the shape was present, false otherwise.
    pub fn remove(&mut self, shape_id: u32) -> bool {
        // Fast path for empty set
        if self.ids.is_empty() {
            return false;
        }

        // Fast path for first element
        if self.ids[0] == shape_id {
            self.ids.remove(0);
            return true;
        }

        // Linear search
        for i in 1..self.ids.len() {
            if self.ids[i] == shape_id {
                self.ids.remove(i);
                return true;
            }
            if self.ids[i] > shape_id {
                return false; // Not found (sorted order)
            }
        }

        false
    }

    /// Toggles a shape ID: removes it if present, inserts if absent.
    ///
    /// Returns true if the shape is now in the set, false if removed.
    pub fn toggle(&mut self, shape_id: u32) -> bool {
        // Fast path for empty set
        if self.ids.is_empty() {
            self.ids.push(shape_id);
            return true;
        }

        // Fast path for first element
        if self.ids[0] == shape_id {
            self.ids.remove(0);
            return false;
        }

        // Linear search
        let mut pos = 0;
        while pos < self.ids.len() && self.ids[pos] < shape_id {
            pos += 1;
        }

        if pos < self.ids.len() && self.ids[pos] == shape_id {
            self.ids.remove(pos);
            false
        } else {
            self.ids.insert(pos, shape_id);
            true
        }
    }

    /// Clears all shape IDs from the set.
    #[inline]
    pub fn clear(&mut self) {
        self.ids.clear();
    }

    /// Returns an iterator over the shape IDs.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = u32> + '_ {
        self.ids.iter().copied()
    }

    /// Returns the position of the first element >= shape_id.
    fn lower_bound(&self, shape_id: u32) -> usize {
        let mut pos = 0;
        while pos < self.ids.len() && self.ids[pos] < shape_id {
            pos += 1;
        }
        pos
    }

    /// Saves and clears state for shape IDs < limit.
    ///
    /// Returns the saved state for later restoration.
    pub fn save_and_clear_before(&mut self, limit_shape_id: u32) -> Vec<u32> {
        let pos = self.lower_bound(limit_shape_id);
        let saved: Vec<u32> = self.ids[..pos].to_vec();
        self.ids.drain(..pos);
        saved
    }

    /// Restores previously saved state for shape IDs < limit.
    pub fn restore_before(&mut self, limit_shape_id: u32, saved: &[u32]) {
        // Remove any IDs < limit that were added since save
        let pos = self.lower_bound(limit_shape_id);
        self.ids.drain(..pos);

        // Insert saved IDs at the beginning
        self.ids.splice(0..0, saved.iter().copied());
    }
}

impl PartialEq for ShapeIdSet {
    fn eq(&self, other: &Self) -> bool {
        self.ids == other.ids
    }
}

impl Eq for ShapeIdSet {}

// =============================================================================
// InteriorTracker - Tracks Shape Containment During Cell Traversal
// =============================================================================

/// Tracks which shapes contain the current focus point during Z-order traversal.
///
/// The InteriorTracker maintains the set of shape IDs whose interiors contain
/// the current "focus" point. As the focus moves through cells in Z-order,
/// edge crossings are detected and the containment state is toggled accordingly.
///
/// # Usage Pattern
///
/// ```ignore
/// let mut tracker = InteriorTracker::new(&bounds);
///
/// // Register shapes with their initial containment at origin
/// for (shape_id, contains_origin) in shapes_with_interiors {
///     tracker.add_shape(shape_id, contains_origin);
/// }
///
/// // Traverse cells in Z-order
/// for cell in cells_in_z_order {
///     // Move to entry vertex
///     if !tracker.at_cellid(cell.id()) {
///         tracker.draw_to(&cell.entry_vertex());
///         for edge in edges_crossing_focus_line {
///             tracker.test_edge(edge.shape_id, &edge.v0, &edge.v1);
///         }
///     }
///
///     // Now shape_ids() contains shapes with focus at entry vertex
///     tracker.move_to(&cell.center());
///
///     // Record which shapes contain the cell center
///     for shape_id in tracker.shape_ids() {
///         // This shape contains the cell center
///     }
///
///     // Prepare for next cell
///     tracker.move_to(&cell.exit_vertex());
///     tracker.set_next_cellid(cell.id().next());
/// }
/// ```
///
/// # Differences from S2 InteriorTracker
///
/// | S2 (Spherical) | Quadtree (Planar) |
/// |----------------|-------------------|
/// | Origin at face 0 corner on sphere | Origin at bottom-left of bounds |
/// | S2EdgeCrosser (great circles) | EdgeCrosser2D (line segments) |
/// | Hilbert curve traversal | Z-order (Morton) traversal |
/// | S2CellId with 6 faces | QuadtreeCellId single plane |
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
    shape_ids: ShapeIdSet,

    /// Whether any shapes are being tracked.
    is_active: bool,

    /// Saved state for incremental updates.
    saved_ids: Vec<u32>,

    /// Saved is_active state.
    saved_is_active: bool,

    /// If non-negative, indicates only some edges of this shape are being
    /// added, so its interior should not be processed yet.
    partial_shape_id: i32,
}

impl InteriorTracker {
    /// Creates a new InteriorTracker with the origin at the bottom-left corner.
    ///
    /// The origin is placed at (bounds.min_x - 1, bounds.min_y - 1) to ensure
    /// it is outside all shapes indexed within the bounds.
    pub fn new(bounds: &Bounds) -> Self {
        let origin = Point2D::new(bounds.min_x - 1.0, bounds.min_y - 1.0);
        let next_cellid = QuadtreeCellId::begin(QuadtreeCellId::MAX_LEVEL);

        Self {
            origin,
            a: origin,
            b: origin,
            next_cellid,
            crosser: EdgeCrosser2D::new(&origin, &origin),
            shape_ids: ShapeIdSet::new(),
            is_active: false,
            saved_ids: Vec::new(),
            saved_is_active: false,
            partial_shape_id: -1,
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

    /// Removes a shape from tracking.
    ///
    /// This is used during incremental updates when shapes are removed.
    pub fn remove_shape(&mut self, shape_id: u32) {
        self.shape_ids.remove(shape_id);
    }

    /// Returns true if the given shape contains the current focus.
    #[inline]
    pub fn shape_contains_focus(&self, shape_id: u32) -> bool {
        self.shape_ids.contains(shape_id)
    }

    /// Returns the set of shape IDs that contain the current focus.
    #[inline]
    pub fn shape_ids(&self) -> &ShapeIdSet {
        &self.shape_ids
    }

    /// Returns the shape IDs as a slice.
    #[inline]
    pub fn shape_ids_slice(&self) -> &[u32] {
        self.shape_ids.as_slice()
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

    /// Returns the partial shape ID, or -1 if none.
    ///
    /// When only some edges of a shape are being added (during batched updates),
    /// this indicates which shape's interior should not be processed yet.
    #[inline]
    pub fn partial_shape_id(&self) -> i32 {
        self.partial_shape_id
    }

    /// Sets the partial shape ID.
    #[inline]
    pub fn set_partial_shape_id(&mut self, shape_id: i32) {
        self.partial_shape_id = shape_id;
    }

    /// Clears the partial shape ID.
    #[inline]
    pub fn clear_partial_shape_id(&mut self) {
        self.partial_shape_id = -1;
    }

    /// Resets the tracker to its initial state at the origin.
    pub fn reset(&mut self) {
        self.a = self.origin;
        self.b = self.origin;
        self.next_cellid = QuadtreeCellId::begin(QuadtreeCellId::MAX_LEVEL);
        self.crosser = EdgeCrosser2D::new(&self.origin, &self.origin);
        self.shape_ids.clear();
        self.is_active = false;
        self.saved_ids.clear();
        self.partial_shape_id = -1;
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

/// Computes whether a shape's interior contains the origin, given an existing
/// index cell that contains the origin.
///
/// This is more efficient than brute_force when we have an existing index,
/// as we only need to test the edges in the cell containing the origin.
///
/// # Arguments
///
/// * `origin` - The InteriorTracker's origin point
/// * `cell_center` - The center of the cell containing the origin
/// * `contains_center` - Whether the shape contains the cell center
/// * `vertices` - The shape's vertices
/// * `edge_indices` - The edges that intersect the cell
///
/// # Returns
///
/// `true` if the origin is inside the shape's interior
pub fn contains_origin_via_cell(
    origin: &Point2D,
    cell_center: &Point2D,
    contains_center: bool,
    vertices: &[Point2D],
    edge_indices: &[u16],
) -> bool {
    if vertices.len() < 3 {
        return false;
    }

    let n = vertices.len();
    let mut inside = contains_center;

    let mut crosser = EdgeCrosser2D::new(cell_center, origin);

    for &edge_idx in edge_indices {
        let i = edge_idx as usize;
        let v0 = &vertices[i];
        let v1 = &vertices[(i + 1) % n];

        if crosser.edge_or_vertex_crossing(v0, v1) {
            inside = !inside;
        }
    }

    inside
}
