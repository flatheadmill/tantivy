//! Point-in-polygon query for quadtree spatial index.

use std::io;

use common::BitSet;

use super::merge::GeometryCache;
use super::segment::{CellView, QuadtreeSegment, ShapeView};
use crate::spatial::quadtree::{edge_or_vertex_crossing_2d, Point2D, QuadtreeCellId};

// =============================================================================
// QuadtreeSegment Query Extensions
// =============================================================================

impl<'a> QuadtreeSegment<'a> {
    /// Finds the cell containing the given point.
    ///
    /// Computes cell_id at max level, then walks up parents until
    /// a cell is found in the index.
    pub fn locate_point(&self, point: &Point2D) -> Option<CellView<'a>> {
        if !self.global_bounds().contains_point(point) {
            return None;
        }

        let mut id = QuadtreeCellId::from_xy(
            point.x,
            point.y,
            QuadtreeCellId::MAX_LEVEL,
            self.global_bounds(),
        );

        loop {
            if let Some(cell) = self.find_cell(id) {
                return Some(cell);
            }
            id = id.parent()?;
        }
    }
}

// =============================================================================
// Contains Point Query
// =============================================================================

/// Returns document IDs of shapes containing the query point.
pub fn contains_point<G: GeometryCache>(
    segment: &QuadtreeSegment,
    point: &Point2D,
    geometry: &mut G,
    results: &mut BitSet,
) -> io::Result<()> {
    let cell = match segment.locate_point(point) {
        Some(c) => c,
        None => return Ok(()),
    };

    let cell_bounds = cell.cell_id().to_bounds(segment.global_bounds());
    let cell_center = cell_bounds.center();

    for shape in cell.shapes() {
        if test_shape_contains_point(&shape, &cell_center, point, geometry)? {
            results.insert(shape.doc_id());
        }
    }

    Ok(())
}

/// Tests if a shape contains the query point.
///
/// Uses ray-casting from cell center (known containment) to query point.
/// Each edge crossing flips inside/outside state.
fn test_shape_contains_point<G: GeometryCache>(
    shape: &ShapeView,
    cell_center: &Point2D,
    point: &Point2D,
    geometry: &mut G,
) -> io::Result<bool> {
    let mut inside = shape.contains_center();

    if shape.num_edges() == 0 {
        return Ok(inside);
    }

    let vertices = geometry.get(shape.doc_id())?;
    let n = vertices.len();

    for edge_idx in shape.edges() {
        let i = edge_idx as usize;
        let v0 = &vertices[i];
        let v1 = &vertices[(i + 1) % n];

        if edge_or_vertex_crossing_2d(cell_center, point, v0, v1) {
            inside = !inside;
        }
    }

    Ok(inside)
}
