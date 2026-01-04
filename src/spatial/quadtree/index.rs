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
