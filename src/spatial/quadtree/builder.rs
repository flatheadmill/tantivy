//! Multi-Shape Index Builder
//!
//! This module implements the quadtree index construction algorithm, following
//! the architecture of S2's MutableS2ShapeIndex but adapted for planar 2D geometry.

use std::collections::BTreeMap;
use std::io::{self, Write};

use common::CountingWriter;

use crate::directory::WritePtr;
use crate::spatial::quadtree::{
    contains_tracker_origin, Bounds, ClippedShape, InteriorTracker, PaddedCell, Point2D,
    QuadtreeCell, QuadtreeCellId, Rect, ShapeIdSet,
};

// =============================================================================
// Configuration Constants
// =============================================================================

/// Default maximum edges per cell before subdivision.
/// S2 uses 10 by default; we follow the same.
pub const DEFAULT_MAX_EDGES_PER_CELL: usize = 10;

/// Default padding for numerical tolerance.
/// Edges within this distance of cell boundaries are included in both cells.
pub const DEFAULT_CELL_PADDING: f64 = 1e-14;

/// Ratio of cell size to edge length that determines when an edge is "long".
/// Long edges don't count toward the max_edges_per_cell limit.
pub const CELL_SIZE_TO_LONG_EDGE_RATIO: f64 = 1.0;

/// Minimum fraction of short edges required to trigger subdivision.
/// This bounds the total index size to O(n) edges.
pub const MIN_SHORT_EDGE_FRACTION: f64 = 0.1;

// =============================================================================
// FaceEdge - Original Edge Data with Metadata
// =============================================================================

/// Edge data for index construction.
///
/// A FaceEdge represents a single edge from a shape, with metadata needed
/// for index building. Unlike S2 which projects to cube faces, we work
/// directly in the 2D coordinate space.
#[derive(Debug, Clone)]
struct FaceEdge {
    /// Shape/document ID.
    shape_id: u32,

    /// Edge index within the shape (0 to num_edges-1).
    edge_id: u16,

    /// Whether this shape has an interior (polygons do, polylines don't).
    has_interior: bool,

    /// First level at which this edge is considered "long".
    /// Long edges don't count toward max_edges_per_cell.
    max_level: u8,

    /// First endpoint of the edge.
    v0: Point2D,

    /// Second endpoint of the edge.
    v1: Point2D,
}

impl FaceEdge {
    /// Creates a new FaceEdge.
    fn new(
        shape_id: u32,
        edge_id: u16,
        has_interior: bool,
        v0: Point2D,
        v1: Point2D,
        global_bounds: &Bounds,
    ) -> Self {
        let max_level = Self::compute_max_level(&v0, &v1, global_bounds);
        Self {
            shape_id,
            edge_id,
            has_interior,
            max_level,
            v0,
            v1,
        }
    }

    /// Computes the first level at which this edge is considered "long".
    ///
    /// An edge is "long" at a given level if the cell edge length at that level
    /// is less than the edge length multiplied by CELL_SIZE_TO_LONG_EDGE_RATIO.
    fn compute_max_level(v0: &Point2D, v1: &Point2D, global_bounds: &Bounds) -> u8 {
        let dx = v1.x - v0.x;
        let dy = v1.y - v0.y;
        let edge_length = (dx * dx + dy * dy).sqrt();

        if edge_length == 0.0 {
            return QuadtreeCellId::MAX_LEVEL;
        }

        // The cell edge length at which this edge becomes "long"
        let max_cell_edge = edge_length * CELL_SIZE_TO_LONG_EDGE_RATIO;

        // Find the level where cell edge <= max_cell_edge
        let span = (global_bounds.max_x - global_bounds.min_x)
            .max(global_bounds.max_y - global_bounds.min_y);

        let mut level = 0u8;
        let mut cell_edge = span;

        while cell_edge > max_cell_edge && level < QuadtreeCellId::MAX_LEVEL {
            cell_edge /= 2.0;
            level += 1;
        }

        level
    }
}

// =============================================================================
// ClippedEdge - Edge with Clipped Bounds During Subdivision
// =============================================================================

/// A clipped edge during recursive subdivision.
///
/// During index construction, edges are recursively clipped to child cells.
/// The ClippedEdge stores a reference to the original FaceEdge plus the
/// clipped bounding box.
#[derive(Debug, Clone)]
struct ClippedEdge<'a> {
    /// Reference to the original edge data.
    face_edge: &'a FaceEdge,

    /// Clipped bounding rectangle in world coordinates.
    /// This is progressively shrunk as we descend the tree.
    bound: Rect,
}

impl<'a> ClippedEdge<'a> {
    /// Creates a ClippedEdge from a FaceEdge with initial bounds.
    fn from_face_edge(face_edge: &'a FaceEdge) -> Self {
        // FIX: Use from_coords with min/max instead of non-existent from_point_pair
        let bound = Rect::from_coords(
            face_edge.v0.x.min(face_edge.v1.x),
            face_edge.v0.y.min(face_edge.v1.y),
            face_edge.v0.x.max(face_edge.v1.x),
            face_edge.v0.y.max(face_edge.v1.y),
        );
        Self { face_edge, bound }
    }

    /// Creates a ClippedEdge with a specified bound.
    fn with_bound(face_edge: &'a FaceEdge, bound: Rect) -> Self {
        Self { face_edge, bound }
    }
}

// =============================================================================
// ShapeData - Per-Shape Information for Building
// =============================================================================

/// Information about a shape being indexed.
#[derive(Debug, Clone)]
struct ShapeData {
    /// Document/shape ID.
    geometry_id: u32,

    /// Vertices of the shape (closed polygon: first == last).
    vertices: Vec<Point2D>,

    /// Whether this shape has an interior (dimension 2).
    has_interior: bool,
}

// =============================================================================
// QuadtreeIndexBuilder - Main Builder
// =============================================================================

/// Configuration for index building.
#[derive(Debug, Clone)]
pub struct BuilderOptions {
    /// Maximum edges per cell before subdivision.
    pub max_edges_per_cell: usize,

    /// Padding for numerical tolerance.
    pub cell_padding: f64,
}

impl Default for BuilderOptions {
    fn default() -> Self {
        Self {
            max_edges_per_cell: DEFAULT_MAX_EDGES_PER_CELL,
            cell_padding: DEFAULT_CELL_PADDING,
        }
    }
}

/// Builds a quadtree index from a collection of shapes.
///
/// # Usage
///
/// ```ignore
/// let mut builder = QuadtreeIndexBuilder::new(bounds);
///
/// // Add shapes (polygons with geometry_id)
/// builder.add_shape(0, &polygon1_vertices);
/// builder.add_shape(1, &polygon2_vertices);
///
/// // Build the index
/// let index = builder.build();
///
/// // Query the index
/// for cell in index.cells() {
///     // Process cells...
/// }
/// ```
pub struct QuadtreeIndexBuilder {
    /// Global coordinate bounds.
    bounds: Bounds,

    /// Builder options.
    options: BuilderOptions,

    /// Shapes to be indexed.
    shapes: Vec<ShapeData>,

    /// Edges from all shapes.
    edges: Vec<FaceEdge>,
}

impl QuadtreeIndexBuilder {
    /// Creates a new builder with the given global bounds.
    pub fn new(bounds: Bounds) -> Self {
        Self {
            bounds,
            options: BuilderOptions::default(),
            shapes: Vec::new(),
            edges: Vec::new(),
        }
    }

    /// Creates a new builder with custom options.
    pub fn with_options(bounds: Bounds, options: BuilderOptions) -> Self {
        Self {
            bounds,
            options,
            shapes: Vec::new(),
            edges: Vec::new(),
        }
    }

    /// Returns the global bounds.
    pub fn bounds(&self) -> &Bounds {
        &self.bounds
    }

    /// Returns the builder options.
    pub fn options(&self) -> &BuilderOptions {
        &self.options
    }

    /// Adds a polygon shape to the index.
    ///
    /// # Arguments
    ///
    /// * `geometry_id` - The document ID for this shape
    /// * `vertices` - The polygon vertices in order (will be closed automatically)
    ///
    /// # Panics
    ///
    /// Panics if vertices has fewer than 3 points.
    pub fn add_shape(&mut self, geometry_id: u32, vertices: &[Point2D]) {
        assert!(vertices.len() >= 3, "Polygon must have at least 3 vertices");

        let has_interior = true; // Polygons have interiors

        // Store shape data
        self.shapes.push(ShapeData {
            geometry_id,
            vertices: vertices.to_vec(),
            has_interior,
        });

        // Create edges
        let n = vertices.len();
        for i in 0..n {
            let v0 = vertices[i];
            let v1 = vertices[(i + 1) % n];

            self.edges.push(FaceEdge::new(
                geometry_id,
                i as u16,
                has_interior,
                v0,
                v1,
                &self.bounds,
            ));
        }
    }

    /// Adds a polyline shape (no interior) to the index.
    ///
    /// # Arguments
    ///
    /// * `geometry_id` - The document ID for this shape
    /// * `vertices` - The polyline vertices in order
    pub fn add_polyline(&mut self, geometry_id: u32, vertices: &[Point2D]) {
        if vertices.len() < 2 {
            return; // Polylines need at least 2 vertices
        }

        let has_interior = false; // Polylines don't have interiors

        self.shapes.push(ShapeData {
            geometry_id,
            vertices: vertices.to_vec(),
            has_interior,
        });

        // Create edges (n-1 edges for n vertices)
        for i in 0..vertices.len() - 1 {
            self.edges.push(FaceEdge::new(
                geometry_id,
                i as u16,
                has_interior,
                vertices[i],
                vertices[i + 1],
                &self.bounds,
            ));
        }
    }

    /// Returns the number of shapes added.
    pub fn num_shapes(&self) -> usize {
        self.shapes.len()
    }

    /// Returns the number of edges across all shapes.
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    /// Builds the quadtree index.
    ///
    /// This consumes the builder and returns the final index.
    pub fn build(self) -> QuadtreeIndex {
        if self.edges.is_empty() {
            // Handle empty index (shapes with interiors but no edges)
            return self.build_empty_index();
        }

        // Initialize the interior tracker
        let mut tracker = InteriorTracker::new(&self.bounds);

        // Register shapes with interiors
        for shape in &self.shapes {
            if shape.has_interior {
                let contains_origin = contains_tracker_origin(&tracker.origin(), &shape.vertices);
                tracker.add_shape(shape.geometry_id, contains_origin);
            }
        }

        // Create initial clipped edges
        let clipped_edges: Vec<ClippedEdge> =
            self.edges.iter().map(ClippedEdge::from_face_edge).collect();

        // Compute bounding box of all edges
        // FIX: Use union instead of union_rect
        let mut all_bounds = Rect::empty();
        for edge in &clipped_edges {
            all_bounds = all_bounds.union(&edge.bound);
        }

        // Create root padded cell
        let root_pcell = PaddedCell::root(self.options.cell_padding, &self.bounds);

        // Try to shrink to a tighter starting cell
        let start_id = root_pcell.shrink_to_fit(&all_bounds, &self.bounds);
        let start_pcell = PaddedCell::new(start_id, self.options.cell_padding, &self.bounds);

        // Build the cell map
        let mut cells: BTreeMap<QuadtreeCellId, QuadtreeCell> = BTreeMap::new();

        // Recursively build the index
        self.update_edges(&start_pcell, &clipped_edges, &mut tracker, &mut cells);

        QuadtreeIndex {
            bounds: self.bounds,
            cells,
        }
    }

    /// Builds an empty index (for shapes with no edges but possibly with interiors).
    fn build_empty_index(self) -> QuadtreeIndex {
        QuadtreeIndex {
            bounds: self.bounds,
            cells: BTreeMap::new(),
        }
    }

    /// Recursively processes edges and creates index cells.
    fn update_edges(
        &self,
        pcell: &PaddedCell,
        edges: &[ClippedEdge],
        tracker: &mut InteriorTracker,
        cells: &mut BTreeMap<QuadtreeCellId, QuadtreeCell>,
    ) {
        // Try to create an index cell at this level
        if self.make_index_cell(pcell, edges, tracker, cells) {
            return;
        }

        // Need to subdivide - get the middle region
        let mut pcell_mut = pcell.clone();
        let middle = pcell_mut.middle();

        // Distribute edges to children
        // child_edges[i][j] where i=x_index, j=y_index
        let mut child_edges: [[Vec<ClippedEdge>; 2]; 2] = Default::default();

        // Reserve space for efficiency
        let num_edges = edges.len();
        for i in 0..2 {
            for j in 0..2 {
                child_edges[i][j].reserve(num_edges);
            }
        }

        // Clip each edge to the appropriate children
        for edge in edges {
            self.clip_edge_to_children(edge, &middle, &mut child_edges);
        }

        // Process children in Z-order (Morton code order)
        for pos in 0..4 {
            let (i, j) = pcell.child_ij(pos);
            let child_idx = PaddedCell::ij_to_child_index(i, j);

            if child_edges[i][j].is_empty() && tracker.shape_ids().is_empty() {
                continue;
            }

            if let Some(child_pcell) = pcell.child(child_idx) {
                self.update_edges(&child_pcell, &child_edges[i][j], tracker, cells);
            }
        }
    }

    /// Clips an edge to the four child quadrants.
    fn clip_edge_to_children<'a>(
        &self,
        edge: &ClippedEdge<'a>,
        middle: &Rect,
        child_edges: &mut [[Vec<ClippedEdge<'a>>; 2]; 2],
    ) {
        // Determine which children the edge's bound intersects
        let bound = &edge.bound;

        // X-axis clipping (unused variables prefixed with underscore)
        let _x_in_lo = bound.x.lo <= middle.x.hi; // Intersects left children
        let _x_in_hi = bound.x.hi >= middle.x.lo; // Intersects right children

        // Y-axis clipping
        let _y_in_lo = bound.y.lo <= middle.y.hi; // Intersects bottom children
        let _y_in_hi = bound.y.hi >= middle.y.lo; // Intersects top children

        // Fast path: edge is entirely in one child
        if bound.x.hi <= middle.x.lo {
            // Entirely in left children
            self.clip_v_axis(edge, middle, 0, child_edges);
        } else if bound.x.lo >= middle.x.hi {
            // Entirely in right children
            self.clip_v_axis(edge, middle, 1, child_edges);
        } else if bound.y.hi <= middle.y.lo {
            // Entirely in bottom children
            child_edges[0][0].push(self.clip_u_bound(edge, true, middle.x.hi));
            child_edges[1][0].push(self.clip_u_bound(edge, false, middle.x.lo));
        } else if bound.y.lo >= middle.y.hi {
            // Entirely in top children
            child_edges[0][1].push(self.clip_u_bound(edge, true, middle.x.hi));
            child_edges[1][1].push(self.clip_u_bound(edge, false, middle.x.lo));
        } else {
            // Edge spans multiple quadrants - clip to each
            let left = self.clip_u_bound(edge, true, middle.x.hi);
            self.clip_v_axis(&left, middle, 0, child_edges);

            let right = self.clip_u_bound(edge, false, middle.x.lo);
            self.clip_v_axis(&right, middle, 1, child_edges);
        }
    }

    /// Clips an edge along the V (Y) axis and adds to children.
    fn clip_v_axis<'a>(
        &self,
        edge: &ClippedEdge<'a>,
        middle: &Rect,
        i: usize,
        child_edges: &mut [[Vec<ClippedEdge<'a>>; 2]; 2],
    ) {
        if edge.bound.y.hi <= middle.y.lo {
            // Entirely in bottom child
            child_edges[i][0].push(edge.clone());
        } else if edge.bound.y.lo >= middle.y.hi {
            // Entirely in top child
            child_edges[i][1].push(edge.clone());
        } else {
            // Spans both children
            child_edges[i][0].push(self.clip_v_bound(edge, true, middle.y.hi));
            child_edges[i][1].push(self.clip_v_bound(edge, false, middle.y.lo));
        }
    }

    /// Clips the U (X) bound of an edge.
    ///
    /// If `clip_hi` is true, clips the high end to `u`.
    /// If `clip_hi` is false, clips the low end to `u`.
    fn clip_u_bound<'a>(&self, edge: &ClippedEdge<'a>, clip_hi: bool, u: f64) -> ClippedEdge<'a> {
        let mut new_bound = edge.bound;

        if clip_hi {
            if edge.bound.x.hi <= u {
                return edge.clone(); // No clipping needed
            }
            new_bound.x.hi = u.min(edge.bound.x.hi);
        } else {
            if edge.bound.x.lo >= u {
                return edge.clone(); // No clipping needed
            }
            new_bound.x.lo = u.max(edge.bound.x.lo);
        }

        // Interpolate V bound based on where the edge crosses U = u
        let fe = edge.face_edge;
        if (fe.v1.x - fe.v0.x).abs() > f64::EPSILON {
            let t = (u - fe.v0.x) / (fe.v1.x - fe.v0.x);
            if t > 0.0 && t < 1.0 {
                let v_at_u = fe.v0.y + t * (fe.v1.y - fe.v0.y);
                // Clamp to the existing bound
                let v_clamped = v_at_u.clamp(edge.bound.y.lo, edge.bound.y.hi);

                // Update the v bound based on which endpoint we're keeping
                let slope_positive = (fe.v1.x > fe.v0.x) == (fe.v1.y > fe.v0.y);
                if clip_hi {
                    // Keeping the low-u side
                    if slope_positive {
                        new_bound.y.hi = new_bound.y.hi.min(v_clamped + f64::EPSILON);
                    } else {
                        new_bound.y.lo = new_bound.y.lo.max(v_clamped - f64::EPSILON);
                    }
                } else {
                    // Keeping the high-u side
                    if slope_positive {
                        new_bound.y.lo = new_bound.y.lo.max(v_clamped - f64::EPSILON);
                    } else {
                        new_bound.y.hi = new_bound.y.hi.min(v_clamped + f64::EPSILON);
                    }
                }
            }
        }

        ClippedEdge::with_bound(edge.face_edge, new_bound)
    }

    /// Clips the V (Y) bound of an edge.
    fn clip_v_bound<'a>(&self, edge: &ClippedEdge<'a>, clip_hi: bool, v: f64) -> ClippedEdge<'a> {
        let mut new_bound = edge.bound;

        if clip_hi {
            if edge.bound.y.hi <= v {
                return edge.clone();
            }
            new_bound.y.hi = v.min(edge.bound.y.hi);
        } else {
            if edge.bound.y.lo >= v {
                return edge.clone();
            }
            new_bound.y.lo = v.max(edge.bound.y.lo);
        }

        // Interpolate U bound based on where the edge crosses V = v
        let fe = edge.face_edge;
        if (fe.v1.y - fe.v0.y).abs() > f64::EPSILON {
            let t = (v - fe.v0.y) / (fe.v1.y - fe.v0.y);
            if t > 0.0 && t < 1.0 {
                let u_at_v = fe.v0.x + t * (fe.v1.x - fe.v0.x);
                let u_clamped = u_at_v.clamp(edge.bound.x.lo, edge.bound.x.hi);

                let slope_positive = (fe.v1.x > fe.v0.x) == (fe.v1.y > fe.v0.y);
                if clip_hi {
                    if slope_positive {
                        new_bound.x.hi = new_bound.x.hi.min(u_clamped + f64::EPSILON);
                    } else {
                        new_bound.x.lo = new_bound.x.lo.max(u_clamped - f64::EPSILON);
                    }
                } else {
                    if slope_positive {
                        new_bound.x.lo = new_bound.x.lo.max(u_clamped - f64::EPSILON);
                    } else {
                        new_bound.x.hi = new_bound.x.hi.min(u_clamped + f64::EPSILON);
                    }
                }
            }
        }

        ClippedEdge::with_bound(edge.face_edge, new_bound)
    }

    /// Attempts to create an index cell at this level.
    ///
    /// Returns true if a cell was created, false if subdivision is needed.
    fn make_index_cell(
        &self,
        pcell: &PaddedCell,
        edges: &[ClippedEdge],
        tracker: &mut InteriorTracker,
        cells: &mut BTreeMap<QuadtreeCellId, QuadtreeCell>,
    ) -> bool {
        // Early exit if nothing to do
        if edges.is_empty() && tracker.shape_ids().is_empty() {
            return true;
        }

        // Check if we should subdivide based on edge count
        if edges.len() > self.options.max_edges_per_cell {
            // Count "short" edges (edges that are short relative to cell size)
            let max_short = std::cmp::max(
                self.options.max_edges_per_cell,
                (MIN_SHORT_EDGE_FRACTION * (edges.len() + tracker.shape_ids().len()) as f64)
                    as usize,
            );

            let mut short_count = 0;
            for edge in edges {
                if pcell.level() < edge.face_edge.max_level {
                    short_count += 1;
                    if short_count > max_short {
                        return false; // Need to subdivide
                    }
                }
            }
        }

        // Create the index cell

        // Move the tracker to the cell center and test edges
        if tracker.is_active() && !edges.is_empty() {
            if !tracker.at_cellid(pcell.id()) {
                tracker.move_to(&pcell.entry_vertex());
            }
            tracker.draw_to(&pcell.center());
            self.test_all_edges(edges, tracker);
        }

        // FIX: Collect shape IDs into a Vec to end the immutable borrow of tracker
        let containing_shape_ids: Vec<u32> = tracker.shape_ids().iter().collect();
        let num_shapes = self.count_shapes_vec(edges, &containing_shape_ids);

        if num_shapes == 0 {
            // No shapes to record
            if tracker.is_active() && !edges.is_empty() {
                tracker.draw_to(&pcell.exit_vertex());
                self.test_all_edges(edges, tracker);
                tracker.set_next_cellid(pcell.id().next());
            }
            return true;
        }

        // Build the cell
        let mut cell = QuadtreeCell::with_capacity(pcell.id(), num_shapes);

        // Merge edge shapes and containing shapes
        let mut edge_idx = 0;
        let mut containing_iter = containing_shape_ids.iter().peekable();

        while edge_idx < edges.len() || containing_iter.peek().is_some() {
            let edge_shape_id = if edge_idx < edges.len() {
                edges[edge_idx].face_edge.shape_id
            } else {
                u32::MAX
            };

            let containing_shape_id = containing_iter.peek().copied().copied().unwrap_or(u32::MAX);

            if containing_shape_id < edge_shape_id {
                // Shape contains center but has no edges in this cell
                let clipped = ClippedShape::new(containing_shape_id, true);
                cell.add_shape(clipped);
                containing_iter.next();
            } else {
                // Shape has edges in this cell
                let shape_id = edge_shape_id;
                let mut clipped = ClippedShape::new(shape_id, false);

                // Collect all edges for this shape
                while edge_idx < edges.len() && edges[edge_idx].face_edge.shape_id == shape_id {
                    clipped.add_edge(edges[edge_idx].face_edge.edge_id);
                    edge_idx += 1;
                }

                // Check if this shape also contains the center
                if containing_shape_id == shape_id {
                    clipped.set_contains_center(true);
                    containing_iter.next();
                }

                cell.add_shape(clipped);
            }
        }

        cells.insert(pcell.id(), cell);

        // Move tracker to exit vertex - now safe because containing_iter is done
        if tracker.is_active() && !edges.is_empty() {
            tracker.draw_to(&pcell.exit_vertex());
            self.test_all_edges(edges, tracker);
            tracker.set_next_cellid(pcell.id().next());
        }

        true
    }

    /// Tests all edges with interiors against the tracker.
    fn test_all_edges(&self, edges: &[ClippedEdge], tracker: &mut InteriorTracker) {
        for edge in edges {
            if edge.face_edge.has_interior {
                tracker.test_edge(
                    edge.face_edge.shape_id,
                    &edge.face_edge.v0,
                    &edge.face_edge.v1,
                );
            }
        }
    }

    /// Counts the number of distinct shapes in edges plus containing shapes.
    fn count_shapes(&self, edges: &[ClippedEdge], containing: &ShapeIdSet) -> usize {
        let mut count = 0;
        let mut last_shape_id: Option<u32> = None;
        let mut containing_iter = containing.iter().peekable();

        for edge in edges {
            let shape_id = edge.face_edge.shape_id;

            if last_shape_id != Some(shape_id) {
                count += 1;
                last_shape_id = Some(shape_id);

                // Skip containing shapes up to and including this one
                while let Some(c_id) = containing_iter.peek().copied() {
                    if c_id > shape_id {
                        break;
                    }
                    if c_id < shape_id {
                        count += 1;
                    }
                    containing_iter.next();
                }
            }
        }

        // Count remaining containing shapes
        count += containing_iter.count();

        count
    }

    /// Counts the number of distinct shapes in edges plus containing shapes (Vec version).
    fn count_shapes_vec(&self, edges: &[ClippedEdge], containing: &[u32]) -> usize {
        let mut count = 0;
        let mut last_shape_id: Option<u32> = None;
        let mut containing_iter = containing.iter().peekable();

        for edge in edges {
            let shape_id = edge.face_edge.shape_id;

            if last_shape_id != Some(shape_id) {
                count += 1;
                last_shape_id = Some(shape_id);

                // Skip containing shapes up to and including this one
                while let Some(&&c_id) = containing_iter.peek() {
                    if c_id > shape_id {
                        break;
                    }
                    if c_id < shape_id {
                        count += 1;
                    }
                    containing_iter.next();
                }
            }
        }

        // Count remaining containing shapes
        count += containing_iter.count();

        count
    }
}

// =============================================================================
// QuadtreeIndex - The Built Index
// =============================================================================

/// A built quadtree spatial index.
///
/// The index is immutable once built. For updates, create a new builder
/// and rebuild the index.
#[derive(Debug)]
pub struct QuadtreeIndex {
    /// Global coordinate bounds.
    bounds: Bounds,

    /// Index cells, keyed by cell ID (in Z-order).
    cells: BTreeMap<QuadtreeCellId, QuadtreeCell>,
}

impl QuadtreeIndex {
    /// Returns the global bounds.
    pub fn bounds(&self) -> &Bounds {
        &self.bounds
    }

    /// Returns the number of cells in the index.
    pub fn num_cells(&self) -> usize {
        self.cells.len()
    }

    /// Returns true if the index is empty (no cells).
    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    /// Returns the cell for the given cell ID, if it exists.
    pub fn get_cell(&self, id: QuadtreeCellId) -> Option<&QuadtreeCell> {
        self.cells.get(&id)
    }

    /// Finds the cell containing the given point.
    ///
    /// Returns the most specific (deepest) cell that contains the point.
    pub fn locate_point(&self, point: &Point2D) -> Option<&QuadtreeCell> {
        // FIX: Use from_xy instead of non-existent from_point
        let target =
            QuadtreeCellId::from_xy(point.x, point.y, QuadtreeCellId::MAX_LEVEL, &self.bounds);

        // Search for the cell containing this point
        // The cell ID at max level is a descendant of any ancestor cell
        for (cell_id, cell) in self.cells.range(..=target).rev() {
            if cell_id.contains(&target) {
                return Some(cell);
            }
            // If we've gone past any possible ancestor, stop
            if cell_id.range_max() < target {
                break;
            }
        }

        None
    }

    /// Returns an iterator over all cells in Z-order.
    pub fn cells(&self) -> impl Iterator<Item = (&QuadtreeCellId, &QuadtreeCell)> {
        self.cells.iter()
    }

    /// Returns an iterator over cells that may intersect the given bounds.
    // FIX: Add '_ to capture the lifetime properly
    pub fn cells_intersecting(
        &self,
        query_bounds: &Rect,
    ) -> impl Iterator<Item = &QuadtreeCell> + '_ {
        // Copy the Rect (it implements Copy) to avoid lifetime issues
        let query = *query_bounds;
        // For now, simple implementation that filters all cells
        // A more efficient implementation would use the tree structure
        self.cells.values().filter(move |cell| {
            let cell_bounds = cell.cell_id().to_bounds(&self.bounds);
            cell_bounds.intersects(&query)
        })
    }

    /// Serializes the quadtree index.
    ///
    /// Format:
    /// - Cells in cell_id order (variable length)
    /// - Cell index: (cell_id: u64, offset: u64) pairs
    /// - Footer (48 bytes): cell_count, cell_index_offset, bounds, version, magic
    pub fn serialize(&self, write: &mut CountingWriter<WritePtr>) -> io::Result<()> {
        const VERSION: u16 = 1;
        const MAGIC: u16 = 0x5154; // "QT"

        let mut index_entries: Vec<(u64, u64)> = Vec::with_capacity(self.cells.len());

        // Write cells in cell_id order
        for (cell_id, cell) in &self.cells {
            index_entries.push((cell_id.to_raw(), write.written_bytes()));

            // num_shapes
            write.write_all(&(cell.num_shapes() as u32).to_le_bytes())?;

            for shape in cell.shapes() {
                // geometry_id
                write.write_all(&shape.geometry_id().to_le_bytes())?;
                // flags
                let flags: u8 = if shape.contains_center() { 1 } else { 0 };
                write.write_all(&[flags])?;
                // num_edges
                write.write_all(&(shape.num_edges() as u16).to_le_bytes())?;
                // edges
                for &edge in shape.edges() {
                    write.write_all(&edge.to_le_bytes())?;
                }
            }
        }

        // Write cell index
        let cell_index_offset = write.written_bytes();
        for (cell_id, offset) in &index_entries {
            write.write_all(&cell_id.to_le_bytes())?;
            write.write_all(&offset.to_le_bytes())?;
        }

        // Write footer (48 bytes)
        write.write_all(&(self.cells.len() as u32).to_le_bytes())?;
        write.write_all(&cell_index_offset.to_le_bytes())?;
        write.write_all(&self.bounds.min_x.to_le_bytes())?;
        write.write_all(&self.bounds.min_y.to_le_bytes())?;
        write.write_all(&self.bounds.max_x.to_le_bytes())?;
        write.write_all(&self.bounds.max_y.to_le_bytes())?;
        write.write_all(&VERSION.to_le_bytes())?;
        write.write_all(&MAGIC.to_le_bytes())?;

        Ok(())
    }
}
