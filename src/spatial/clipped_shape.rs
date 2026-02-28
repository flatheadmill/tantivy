// =============================================================================
// ClippedShape - Per-Shape Data Within a Cell
// =============================================================================
use smallvec::SmallVec;

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
    pub geometry_id: u32,

    /// Whether the shape's interior contains the cell center.
    ///
    /// For shapes without an interior (e.g., polylines), this is always false.
    /// For polygons, this indicates whether the cell center is inside the polygon.
    pub contains_center: bool,

    /// Indices of edges that intersect this cell.
    ///
    /// Each index refers to an edge in the shape's vertex array:
    /// edge[i] connects vertex[i] to vertex[(i+1) % num_vertices].
    ///
    /// Edges are stored in sorted order for efficient lookup and merging.
    pub edges: SmallVec<[u16; 2]>,
}

impl ClippedShape {
    /// Creates a new ClippedShape with no edges.
    ///
    /// # Arguments
    ///
    /// * `geometry_id` - The document ID of the shape
    /// * `contains_center` - Whether the shape interior contains the cell center
    #[inline]
    pub fn new(geometry_id: u32, contains_center: bool) -> Self {
        Self {
            geometry_id,
            contains_center,
            edges: SmallVec::new(),
        }
    }

    /// Creates a ClippedShape with preallocated edge capacity.
    #[inline]
    pub fn with_capacity(geometry_id: u32, contains_center: bool, capacity: usize) -> Self {
        Self {
            geometry_id,
            contains_center,
            edges: SmallVec::with_capacity(capacity),
        }
    }

    /// Returns the document ID.
    #[inline]
    pub fn geometry_id(&self) -> u32 {
        self.geometry_id
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

    /// Clears all edge indices, keeping geometry_id and contains_center.
    #[inline]
    pub fn clear_edges(&mut self) {
        self.edges.clear();
    }
}

impl PartialEq for ClippedShape {
    fn eq(&self, other: &Self) -> bool {
        self.geometry_id == other.geometry_id
            && self.contains_center == other.contains_center
            && self.edges == other.edges
    }
}

impl Eq for ClippedShape {}
