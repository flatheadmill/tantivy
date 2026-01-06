//! K-way streaming merge for quadtree segments.
//!
//! This module implements the merge operation that combines multiple quadtree
//! segments into a single merged segment. The merge handles:
//!
//! - K-way merging using a min-heap on cell_id
//! - Cell splitting when combined edges exceed threshold
//! - Delete filtering and doc_id remapping
//! - Advisory collapse candidate detection

use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};
use std::io;

use super::collapse_detector::CollapseDetector;
use crate::spatial::quadtree::{
    Bounds, CacheStats, ClippedShape, GeometryCache, InputIterator, Point2D, QuadtreeCell,
    QuadtreeCellId, Rect,
};

// =============================================================================
// MergeOptions
// =============================================================================

/// Options for merge operation.
#[derive(Debug, Clone)]
pub struct MergeOptions {
    /// Maximum edges per cell before splitting.
    pub max_edges_per_cell: usize,

    /// Number of geometries to cache during split operations.
    pub geometry_cache_size: usize,

    /// Threshold for collapse detection (combined edges of siblings).
    /// Set to 0 to disable collapse detection.
    pub collapse_threshold: usize,

    /// Global bounds for the coordinate system.
    pub bounds: Bounds,
}

impl MergeOptions {
    /// Creates options with the given bounds and default thresholds.
    pub fn new(bounds: Bounds) -> Self {
        Self {
            max_edges_per_cell: 10,
            geometry_cache_size: 64,
            collapse_threshold: 25,
            bounds,
        }
    }

    /// Sets the maximum edges per cell.
    pub fn with_max_edges(mut self, max_edges: usize) -> Self {
        self.max_edges_per_cell = max_edges;
        self
    }

    /// Sets the geometry cache size.
    pub fn with_cache_size(mut self, size: usize) -> Self {
        self.geometry_cache_size = size;
        self
    }

    /// Sets the collapse threshold.
    pub fn with_collapse_threshold(mut self, threshold: usize) -> Self {
        self.collapse_threshold = threshold;
        self
    }
}

// =============================================================================
// MergeStats
// =============================================================================

/// Statistics from a merge operation.
#[derive(Debug, Default, Clone)]
pub struct MergeStats {
    /// Total cells read from all inputs.
    pub cells_input: u64,

    /// Total cells written to output.
    pub cells_output: u64,

    /// Number of cells merged (combined with same cell_id).
    pub cells_merged: u64,

    /// Number of cells that were split.
    pub cells_split: u64,

    /// Number of geometries loaded for splits.
    pub geometries_loaded: u64,

    /// Cache statistics from geometry access.
    pub cache_stats: CacheStats,
}

// =============================================================================
// MergeHeapEntry
// =============================================================================

/// Entry in the merge heap.
struct MergeHeapEntry {
    /// Cell ID for ordering.
    cell_id: QuadtreeCellId,

    /// Index of the input iterator.
    input_idx: usize,
}

impl Eq for MergeHeapEntry {}

impl PartialEq for MergeHeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.cell_id == other.cell_id
    }
}

impl Ord for MergeHeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        // Min-heap: reverse ordering
        other.cell_id.cmp(&self.cell_id)
    }
}

impl PartialOrd for MergeHeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// =============================================================================
// merge_cells - Combine Two Cells with Same ID
// =============================================================================

/// Merges two cells with the same cell_id.
///
/// Combines shapes from both cells:
/// - If the same doc_id appears in both, edges are merged
/// - contains_center is OR of both
fn merge_cells(mut a: QuadtreeCell, b: QuadtreeCell) -> QuadtreeCell {
    debug_assert_eq!(a.cell_id(), b.cell_id());

    for shape in b.shapes() {
        if let Some(existing) = a.find_shape_mut(shape.doc_id()) {
            // Merge edges into existing shape
            for &edge in shape.edges() {
                existing.add_edge(edge);
            }
            // contains_center is OR
            if shape.contains_center() {
                existing.set_contains_center(true);
            }
        } else {
            // Add new shape
            let mut new_shape = ClippedShape::with_capacity(
                shape.doc_id(),
                shape.contains_center(),
                shape.num_edges(),
            );
            for &edge in shape.edges() {
                new_shape.add_edge(edge);
            }
            a.add_shape(new_shape);
        }
    }

    a
}

// =============================================================================
// split_cell - Split an Over-Full Cell into Children
// =============================================================================

/// Splits a cell into four children when it exceeds the edge threshold.
///
/// For each shape in the cell:
/// 1. Load the geometry from the cache
/// 2. For each edge, determine which child cells it intersects
/// 3. Recompute contains_center for each child
fn split_cell<G: GeometryCache>(
    cell: QuadtreeCell,
    geometry: &mut G,
    bounds: &Bounds,
    stats: &mut MergeStats,
) -> io::Result<[QuadtreeCell; 4]> {
    let parent_id = cell.cell_id();
    let child_ids = match parent_id.children() {
        Some(ids) => ids,
        None => {
            // At max level, cannot split - return as single child
            return Ok([
                cell,
                QuadtreeCell::new(QuadtreeCellId::ROOT), // Empty placeholder
                QuadtreeCell::new(QuadtreeCellId::ROOT),
                QuadtreeCell::new(QuadtreeCellId::ROOT),
            ]);
        }
    };

    let mut children = [
        QuadtreeCell::new(child_ids[0]),
        QuadtreeCell::new(child_ids[1]),
        QuadtreeCell::new(child_ids[2]),
        QuadtreeCell::new(child_ids[3]),
    ];

    let child_bounds: [Rect; 4] = [
        child_ids[0].to_bounds(bounds),
        child_ids[1].to_bounds(bounds),
        child_ids[2].to_bounds(bounds),
        child_ids[3].to_bounds(bounds),
    ];

    let parent_center = parent_id.to_bounds(bounds).center();

    for shape in cell.shapes() {
        let doc_id = shape.doc_id();

        // Load geometry if shape has edges
        let vertices = if !shape.edges().is_empty() {
            stats.geometries_loaded += 1;
            geometry.get(doc_id)?
        } else {
            &[] as &[Point2D]
        };

        let num_vertices = vertices.len();

        // Distribute edges to children
        for &edge_idx in shape.edges() {
            let idx = edge_idx as usize;
            if idx >= num_vertices {
                continue; // Invalid edge index
            }

            let v0 = &vertices[idx];
            let v1 = &vertices[(idx + 1) % num_vertices];

            // Check which children this edge intersects
            for (i, child_rect) in child_bounds.iter().enumerate() {
                if edge_intersects_rect(v0, v1, child_rect) {
                    children[i].get_or_create_shape(doc_id).add_edge(edge_idx);
                }
            }
        }

        // Propagate contains_center to children
        if shape.contains_center() {
            // For each child, compute contains_center by counting crossings
            // from parent center to child center
            for (i, child_rect) in child_bounds.iter().enumerate() {
                let child_center = child_rect.center();

                // Start with parent's contains_center, count crossings
                let mut inside = true; // Parent contains center

                for &edge_idx in shape.edges() {
                    let idx = edge_idx as usize;
                    if idx >= num_vertices {
                        continue;
                    }

                    let v0 = &vertices[idx];
                    let v1 = &vertices[(idx + 1) % num_vertices];

                    // Does this edge cross the ray from parent_center to child_center?
                    if edge_crosses_ray(&parent_center, &child_center, v0, v1) {
                        inside = !inside;
                    }
                }

                if inside {
                    children[i]
                        .get_or_create_shape(doc_id)
                        .set_contains_center(true);
                }
            }
        }
    }

    Ok(children)
}

/// Tests if an edge intersects a rectangle.
///
/// Uses a conservative bounding box test followed by segment-rectangle intersection.
fn edge_intersects_rect(v0: &Point2D, v1: &Point2D, rect: &Rect) -> bool {
    // Quick bounding box rejection
    let edge_min_x = v0.x.min(v1.x);
    let edge_max_x = v0.x.max(v1.x);
    let edge_min_y = v0.y.min(v1.y);
    let edge_max_y = v0.y.max(v1.y);

    if edge_max_x < rect.x.lo || edge_min_x > rect.x.hi {
        return false;
    }
    if edge_max_y < rect.y.lo || edge_min_y > rect.y.hi {
        return false;
    }

    // If either endpoint is inside, intersection
    if rect.contains_point(v0) || rect.contains_point(v1) {
        return true;
    }

    // Check if edge crosses any side of rectangle
    let corners = [
        Point2D::new(rect.x.lo, rect.y.lo), // bottom-left
        Point2D::new(rect.x.hi, rect.y.lo), // bottom-right
        Point2D::new(rect.x.hi, rect.y.hi), // top-right
        Point2D::new(rect.x.lo, rect.y.hi), // top-left
    ];

    for i in 0..4 {
        let c0 = &corners[i];
        let c1 = &corners[(i + 1) % 4];
        if segments_intersect(v0, v1, c0, c1) {
            return true;
        }
    }

    false
}

/// Tests if two line segments intersect.
fn segments_intersect(a0: &Point2D, a1: &Point2D, b0: &Point2D, b1: &Point2D) -> bool {
    let d1 = cross_product_sign(b0, b1, a0);
    let d2 = cross_product_sign(b0, b1, a1);
    let d3 = cross_product_sign(a0, a1, b0);
    let d4 = cross_product_sign(a0, a1, b1);

    // Standard segment intersection test
    if ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0)) {
        return true;
    }

    // Collinear cases
    if d1 == 0 && on_segment(b0, b1, a0) {
        return true;
    }
    if d2 == 0 && on_segment(b0, b1, a1) {
        return true;
    }
    if d3 == 0 && on_segment(a0, a1, b0) {
        return true;
    }
    if d4 == 0 && on_segment(a0, a1, b1) {
        return true;
    }

    false
}

/// Returns the sign of the cross product (p1-p0) x (p2-p0).
fn cross_product_sign(p0: &Point2D, p1: &Point2D, p2: &Point2D) -> i32 {
    let cross = (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x);
    if cross > 1e-14 {
        1
    } else if cross < -1e-14 {
        -1
    } else {
        0
    }
}

/// Tests if point p is on segment (a, b) when collinear.
fn on_segment(a: &Point2D, b: &Point2D, p: &Point2D) -> bool {
    p.x >= a.x.min(b.x) - 1e-14
        && p.x <= a.x.max(b.x) + 1e-14
        && p.y >= a.y.min(b.y) - 1e-14
        && p.y <= a.y.max(b.y) + 1e-14
}

/// Tests if an edge crosses the ray from `from` to `to`.
///
/// Used for tracking contains_center during split.
fn edge_crosses_ray(from: &Point2D, to: &Point2D, v0: &Point2D, v1: &Point2D) -> bool {
    // Use the same crossing logic as point-in-polygon
    let d0 = cross_product_sign(from, to, v0);
    let d1 = cross_product_sign(from, to, v1);

    // Edge crosses the line through from->to
    if (d0 > 0 && d1 < 0) || (d0 < 0 && d1 > 0) {
        // Check if crossing point is between from and to
        let d2 = cross_product_sign(v0, v1, from);
        let d3 = cross_product_sign(v0, v1, to);

        // The crossing is on the ray if from and to are on opposite sides of edge
        if (d2 <= 0 && d3 >= 0) || (d2 >= 0 && d3 <= 0) {
            return true;
        }
    }

    false
}

// =============================================================================
// merge - Main K-Way Merge Function
// =============================================================================

/// Performs a K-way merge of cell iterators.
///
/// This is the main merge function that combines multiple input segments
/// into a single output stream of cells.
///
/// # Arguments
///
/// * `inputs` - Vector of input iterators (wrapped CellReaders)
/// * `geometry` - Geometry cache for split operations
/// * `options` - Merge configuration
/// * `emit` - Callback for emitting merged cells
///
/// # Returns
///
/// Statistics about the merge operation.
///
/// # Errors
///
/// Returns an error if:
/// - An input iterator encounters an I/O error
/// - The emit callback returns an error
/// - Geometry lookup fails during split
pub fn merge<G, F>(
    mut inputs: Vec<InputIterator>,
    geometry: &mut G,
    options: &MergeOptions,
    mut emit: F,
) -> io::Result<MergeStats>
where
    G: GeometryCache,
    F: FnMut(QuadtreeCell) -> io::Result<()>,
{
    let mut stats = MergeStats::default();
    let mut collapse_detector = CollapseDetector::new(options.collapse_threshold);

    // Pending cells from splits, kept in cell_id order
    let mut pending: VecDeque<QuadtreeCell> = VecDeque::new();

    // Initialize heap with first cell from each input
    let mut heap = BinaryHeap::new();
    for (idx, input) in inputs.iter().enumerate() {
        if let Some(cell_id) = input.current_cell_id() {
            stats.cells_input += 1;
            heap.push(MergeHeapEntry {
                cell_id,
                input_idx: idx,
            });
        }
    }

    while let Some(entry) = heap.pop() {
        let cell_id = entry.cell_id;

        // Emit pending cells that precede this cell_id
        while let Some(pending_cell) = pending.front() {
            if pending_cell.cell_id() >= cell_id {
                break;
            }
            let cell = pending.pop_front().unwrap();
            collapse_detector.observe_cell(&cell);
            emit(cell)?;
            stats.cells_output += 1;
        }

        // Take the cell from this input
        let mut merged = inputs[entry.input_idx].take().unwrap();

        // Advance the input and re-add to heap if it has more
        if let Some(next_id) = inputs[entry.input_idx].current_cell_id() {
            stats.cells_input += 1;
            heap.push(MergeHeapEntry {
                cell_id: next_id,
                input_idx: entry.input_idx,
            });
        }

        // Collect other inputs with the same cell_id
        while let Some(top) = heap.peek() {
            if top.cell_id != cell_id {
                break;
            }
            let other_entry = heap.pop().unwrap();
            let other_cell = inputs[other_entry.input_idx].take().unwrap();

            merged = merge_cells(merged, other_cell);
            stats.cells_merged += 1;

            // Advance and re-add
            if let Some(next_id) = inputs[other_entry.input_idx].current_cell_id() {
                stats.cells_input += 1;
                heap.push(MergeHeapEntry {
                    cell_id: next_id,
                    input_idx: other_entry.input_idx,
                });
            }
        }

        // Check if split is needed
        if merged.total_edges() > options.max_edges_per_cell {
            let children = split_cell(merged, geometry, &options.bounds, &mut stats)?;
            stats.cells_split += 1;

            // Add non-empty children to pending queue
            for child in children {
                if !child.is_empty() && child.cell_id() != QuadtreeCellId::ROOT {
                    // Insert in sorted order
                    let pos = pending
                        .iter()
                        .position(|c| c.cell_id() > child.cell_id())
                        .unwrap_or(pending.len());
                    pending.insert(pos, child);
                }
            }
        } else {
            // Emit directly if it precedes all pending
            if pending
                .front()
                .map_or(true, |p| merged.cell_id() < p.cell_id())
            {
                collapse_detector.observe_cell(&merged);
                emit(merged)?;
                stats.cells_output += 1;
            } else {
                // Insert in sorted order
                let pos = pending
                    .iter()
                    .position(|c| c.cell_id() > merged.cell_id())
                    .unwrap_or(pending.len());
                pending.insert(pos, merged);
            }
        }
    }

    // Emit remaining pending cells
    for cell in pending {
        collapse_detector.observe_cell(&cell);
        emit(cell)?;
        stats.cells_output += 1;
    }

    stats.cache_stats = geometry.stats();
    Ok(stats)
}

// =============================================================================
// Unit Tests
// =============================================================================

#[cfg(test)]
mod merge_tests {
    use std::io::Cursor;

    use super::*;

    fn test_bounds() -> Bounds {
        Bounds::new(0.0, 0.0, 100.0, 100.0)
    }

    // -------------------------------------------------------------------------
    // merge_cells Tests
    // -------------------------------------------------------------------------

    mod merge_cells_tests {
        use super::*;

        #[test]
        fn test_merge_disjoint_shapes() {
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 3, &test_bounds());

            let mut cell1 = QuadtreeCell::new(cell_id);
            cell1.add_shape(ClippedShape::new(1, true));

            let mut cell2 = QuadtreeCell::new(cell_id);
            cell2.add_shape(ClippedShape::new(2, false));

            let merged = merge_cells(cell1, cell2);
            assert_eq!(merged.num_shapes(), 2);
            assert!(merged.find_shape(1).is_some());
            assert!(merged.find_shape(2).is_some());
        }

        #[test]
        fn test_merge_same_shape() {
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 3, &test_bounds());

            let mut cell1 = QuadtreeCell::new(cell_id);
            let mut shape1 = ClippedShape::new(1, false);
            shape1.add_edge(0);
            shape1.add_edge(1);
            cell1.add_shape(shape1);

            let mut cell2 = QuadtreeCell::new(cell_id);
            let mut shape2 = ClippedShape::new(1, true); // Same doc_id, different edges
            shape2.add_edge(2);
            shape2.add_edge(3);
            cell2.add_shape(shape2);

            let merged = merge_cells(cell1, cell2);
            assert_eq!(merged.num_shapes(), 1);

            let shape = merged.find_shape(1).unwrap();
            assert!(shape.contains_center()); // OR of false and true
            assert_eq!(shape.num_edges(), 4);
            assert_eq!(shape.edges(), &[0, 1, 2, 3]);
        }
    }

    // -------------------------------------------------------------------------
    // edge_intersects_rect Tests
    // -------------------------------------------------------------------------

    mod edge_intersects_tests {
        use super::*;

        fn make_rect() -> Rect {
            Rect::from_coords(10.0, 10.0, 20.0, 20.0)
        }

        #[test]
        fn test_edge_inside() {
            let rect = make_rect();
            let v0 = Point2D::new(12.0, 12.0);
            let v1 = Point2D::new(18.0, 18.0);

            assert!(edge_intersects_rect(&v0, &v1, &rect));
        }

        #[test]
        fn test_edge_outside() {
            let rect = make_rect();
            let v0 = Point2D::new(0.0, 0.0);
            let v1 = Point2D::new(5.0, 5.0);

            assert!(!edge_intersects_rect(&v0, &v1, &rect));
        }

        #[test]
        fn test_edge_crossing() {
            let rect = make_rect();
            let v0 = Point2D::new(5.0, 15.0);
            let v1 = Point2D::new(25.0, 15.0);

            assert!(edge_intersects_rect(&v0, &v1, &rect));
        }

        #[test]
        fn test_endpoint_inside() {
            let rect = make_rect();
            let v0 = Point2D::new(15.0, 15.0); // Inside
            let v1 = Point2D::new(50.0, 50.0); // Outside

            assert!(edge_intersects_rect(&v0, &v1, &rect));
        }
    }

    // -------------------------------------------------------------------------
    // split_cell Tests
    // -------------------------------------------------------------------------

    mod split_cell_tests {
        use super::*;
        use crate::spatial::quadtree::{InMemoryGeometryReader, LruGeometryCache};

        fn make_square_geometry() -> InMemoryGeometryReader {
            let mut reader = InMemoryGeometryReader::new();

            // Square from (20, 20) to (80, 80)
            reader.add(
                1,
                vec![
                    Point2D::new(20.0, 20.0),
                    Point2D::new(80.0, 20.0),
                    Point2D::new(80.0, 80.0),
                    Point2D::new(20.0, 80.0),
                ],
            );

            reader
        }

        #[test]
        fn test_split_distributes_edges() {
            let bounds = test_bounds();
            let mut stats = MergeStats::default();
            let reader = make_square_geometry();
            let mut cache = LruGeometryCache::new(reader, 10);

            // Create a cell at root with a square shape
            let mut cell = QuadtreeCell::new(QuadtreeCellId::ROOT);
            let mut shape = ClippedShape::new(1, true);
            shape.add_edge(0); // Bottom edge: (20,20)-(80,20)
            shape.add_edge(1); // Right edge: (80,20)-(80,80)
            shape.add_edge(2); // Top edge: (80,80)-(20,80)
            shape.add_edge(3); // Left edge: (20,80)-(20,20)
            cell.add_shape(shape);

            let children = split_cell(cell, &mut cache, &bounds, &mut stats).unwrap();

            // All children should have some edges from the square
            // The square spans all four quadrants
            let total_edges: usize = children.iter().map(|c| c.total_edges()).sum();
            assert!(total_edges > 0);

            // Each child should have some portion of the edges
            for child in &children {
                if !child.is_empty() {
                    let shape = child.find_shape(1);
                    assert!(shape.is_some());
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Full Merge Integration Tests
    // -------------------------------------------------------------------------

    mod merge_integration_tests {
        use super::*;
        use crate::spatial::quadtree::segment::QuadtreeSegment;
        use crate::spatial::quadtree::{
            DeleteBitSet, DocIdMap, InMemoryGeometryReader, LruGeometryCache,
        };

        fn create_test_segment_bytes(
            bounds: &Bounds,
            cells: Vec<(QuadtreeCellId, Vec<(u32, bool, Vec<u16>)>)>,
        ) -> Vec<u8> {
            let mut buf = Vec::new();
            let mut cell_index = Vec::new();

            // Sort cells by cell_id
            let mut cells = cells;
            cells.sort_by_key(|(id, _)| id.to_raw());

            // Write cell data, track offsets
            for (cell_id, shapes) in &cells {
                let offset = buf.len() as u64;
                cell_index.push((cell_id.to_raw(), offset));

                buf.extend_from_slice(&(shapes.len() as u32).to_le_bytes());
                for (doc_id, contains_center, edges) in shapes {
                    buf.extend_from_slice(&doc_id.to_le_bytes());
                    buf.push(if *contains_center { 1 } else { 0 });
                    buf.extend_from_slice(&(edges.len() as u16).to_le_bytes());
                    for edge in edges {
                        buf.extend_from_slice(&edge.to_le_bytes());
                    }
                }
            }

            // Write cell index
            let cell_index_offset = buf.len() as u64;
            for (cell_id, data_offset) in &cell_index {
                buf.extend_from_slice(&cell_id.to_le_bytes());
                buf.extend_from_slice(&data_offset.to_le_bytes());
            }

            // Write footer
            buf.extend_from_slice(&(cells.len() as u32).to_le_bytes());
            buf.extend_from_slice(&cell_index_offset.to_le_bytes());
            buf.extend_from_slice(&bounds.min_x.to_le_bytes());
            buf.extend_from_slice(&bounds.min_y.to_le_bytes());
            buf.extend_from_slice(&bounds.max_x.to_le_bytes());
            buf.extend_from_slice(&bounds.max_y.to_le_bytes());
            buf.extend_from_slice(&1u16.to_le_bytes()); // version
            buf.extend_from_slice(&0x5154u16.to_le_bytes()); // magic "QT"

            buf
        }

        #[test]
        fn test_merge_empty_inputs() {
            let bounds = test_bounds();
            let options = MergeOptions::new(bounds.clone());

            let inputs: Vec<InputIterator> = Vec::new();
            let mut reader = InMemoryGeometryReader::new();
            let mut cache = LruGeometryCache::new(reader, 10);

            let mut output = Vec::new();
            let stats = merge(inputs, &mut cache, &options, |cell| {
                output.push(cell);
                Ok(())
            })
            .unwrap();

            assert_eq!(stats.cells_output, 0);
            assert!(output.is_empty());
        }

        #[test]
        fn test_merge_single_input() {
            let bounds = test_bounds();
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 3, &bounds);

            let bytes =
                create_test_segment_bytes(&bounds, vec![(cell_id, vec![(1, true, vec![0, 1, 2])])]);
            let segment = QuadtreeSegment::new(&bytes).unwrap();
            let input = InputIterator::unfiltered(&segment, 10);

            let mut geom_reader = InMemoryGeometryReader::new();
            let mut cache = LruGeometryCache::new(geom_reader, 10);
            let options = MergeOptions::new(bounds);

            let mut output = Vec::new();
            let stats = merge(vec![input], &mut cache, &options, |cell| {
                output.push(cell);
                Ok(())
            })
            .unwrap();

            assert_eq!(stats.cells_input, 1);
            assert_eq!(stats.cells_output, 1);
            assert_eq!(output.len(), 1);
            assert_eq!(output[0].cell_id(), cell_id);
        }

        #[test]
        fn test_merge_two_inputs_different_cells() {
            let bounds = test_bounds();
            let cell_id1 = QuadtreeCellId::from_xy(25.0, 25.0, 3, &bounds);
            let cell_id2 = QuadtreeCellId::from_xy(75.0, 75.0, 3, &bounds);

            let bytes1 =
                create_test_segment_bytes(&bounds, vec![(cell_id1, vec![(1, true, vec![0])])]);
            let bytes2 =
                create_test_segment_bytes(&bounds, vec![(cell_id2, vec![(2, false, vec![1])])]);

            let segment1 = QuadtreeSegment::new(&bytes1).unwrap();
            let segment2 = QuadtreeSegment::new(&bytes2).unwrap();
            let input1 = InputIterator::unfiltered(&segment1, 10);
            let input2 = InputIterator::unfiltered(&segment2, 10);

            let mut geom_reader = InMemoryGeometryReader::new();
            let mut cache = LruGeometryCache::new(geom_reader, 10);
            let options = MergeOptions::new(bounds);

            let mut output = Vec::new();
            let stats = merge(vec![input1, input2], &mut cache, &options, |cell| {
                output.push(cell);
                Ok(())
            })
            .unwrap();

            assert_eq!(stats.cells_output, 2);
            assert_eq!(output.len(), 2);

            // Should be in cell_id order
            assert!(output[0].cell_id() < output[1].cell_id());
        }

        #[test]
        fn test_merge_same_cell_combines() {
            let bounds = test_bounds();
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 3, &bounds);

            let bytes1 =
                create_test_segment_bytes(&bounds, vec![(cell_id, vec![(1, false, vec![0, 1])])]);
            let bytes2 =
                create_test_segment_bytes(&bounds, vec![(cell_id, vec![(1, true, vec![2, 3])])]);

            let segment1 = QuadtreeSegment::new(&bytes1).unwrap();
            let segment2 = QuadtreeSegment::new(&bytes2).unwrap();
            let input1 = InputIterator::unfiltered(&segment1, 10);
            let input2 = InputIterator::unfiltered(&segment2, 10);

            let mut geom_reader = InMemoryGeometryReader::new();
            let mut cache = LruGeometryCache::new(geom_reader, 10);
            let options = MergeOptions::new(bounds);

            let mut output = Vec::new();
            let stats = merge(vec![input1, input2], &mut cache, &options, |cell| {
                output.push(cell);
                Ok(())
            })
            .unwrap();

            assert_eq!(stats.cells_merged, 1);
            assert_eq!(stats.cells_output, 1);
            assert_eq!(output.len(), 1);

            let shape = output[0].find_shape(1).unwrap();
            assert!(shape.contains_center()); // OR of false and true
            assert_eq!(shape.num_edges(), 4);
        }

        #[test]
        fn test_delete_filtering() {
            let bounds = test_bounds();
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 3, &bounds);

            let bytes = create_test_segment_bytes(
                &bounds,
                vec![(
                    cell_id,
                    vec![(1, true, vec![0]), (2, true, vec![1]), (3, true, vec![2])],
                )],
            );

            let segment = QuadtreeSegment::new(&bytes).unwrap();

            // Delete doc_id 2
            let mut deletes = DeleteBitSet::new();
            deletes.delete(2);

            let doc_id_map = DocIdMap::identity(10);
            let input = InputIterator::new(&segment, deletes, doc_id_map);

            let mut geom_reader = InMemoryGeometryReader::new();
            let mut cache = LruGeometryCache::new(geom_reader, 10);
            let options = MergeOptions::new(bounds);

            let mut output = Vec::new();
            merge(vec![input], &mut cache, &options, |cell| {
                output.push(cell);
                Ok(())
            })
            .unwrap();

            assert_eq!(output.len(), 1);
            assert_eq!(output[0].num_shapes(), 2); // Only 1 and 3
            assert!(output[0].find_shape(1).is_some());
            assert!(output[0].find_shape(2).is_none()); // Deleted
            assert!(output[0].find_shape(3).is_some());
        }

        #[test]
        fn test_doc_id_remapping() {
            let bounds = test_bounds();
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 3, &bounds);

            let bytes = create_test_segment_bytes(
                &bounds,
                vec![(cell_id, vec![(0, true, vec![0]), (1, false, vec![1])])],
            );

            let segment = QuadtreeSegment::new(&bytes).unwrap();

            // Remap: 0 -> 100, 1 -> 101
            let doc_id_map = DocIdMap::with_offset(2, 100);
            let input = InputIterator::new(&segment, DeleteBitSet::new(), doc_id_map);

            let mut geom_reader = InMemoryGeometryReader::new();
            let mut cache = LruGeometryCache::new(geom_reader, 10);
            let options = MergeOptions::new(bounds);

            let mut output = Vec::new();
            merge(vec![input], &mut cache, &options, |cell| {
                output.push(cell);
                Ok(())
            })
            .unwrap();

            assert_eq!(output.len(), 1);
            assert!(output[0].find_shape(100).is_some());
            assert!(output[0].find_shape(101).is_some());
            assert!(output[0].find_shape(0).is_none()); // Old ID not present
        }
    }
}
