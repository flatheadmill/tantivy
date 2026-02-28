// =============================================================================
// Cell - A Cell in the Index
// =============================================================================

use crate::spatial::ClippedShape;

/// A single cell in the quadtree index.
///
/// Each cell contains a set of ClippedShapes representing the shapes that
/// intersect this cell. Shapes are stored sorted by geometry_id for efficient
/// lookup and merging.
///
/// # Cell Properties
///
/// - `cell_id`: Identifies the cell's position and level in the quadtree
/// - `shapes`: Shapes intersecting this cell, sorted by geometry_id
///
/// # Query Usage
///
/// To test if a point is inside any shape:
/// 1. Find the cell containing the point
/// 2. For each ClippedShape in the cell: a. Start with `contains_center` b. Count edge crossings
///    from cell center to query point c. XOR the crossing parity with contains_center
#[derive(Debug, Clone)]
pub struct Cell<C: Copy> {
    /// The cell ID (encodes position and level).
    cell_id: C,

    /// Shapes intersecting this cell, sorted by geometry_id.
    shapes: Vec<ClippedShape>,
}

impl<C: Copy> Cell<C> {
    /// Creates a new empty cell.
    #[inline]
    pub fn new(cell_id: C) -> Self {
        Self {
            cell_id,
            shapes: Vec::new(),
        }
    }

    /// Creates a new cell with preallocated shape capacity.
    #[inline]
    pub fn with_capacity(cell_id: C, capacity: usize) -> Self {
        Self {
            cell_id,
            shapes: Vec::with_capacity(capacity),
        }
    }

    /// Returns the cell ID.
    #[inline]
    pub fn cell_id(&self) -> C {
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

    /// Finds the ClippedShape for the given geometry_id.
    ///
    /// Returns `None` if no shape with that geometry_id is in this cell.
    /// Uses binary search for O(log n) lookup.
    pub fn find_shape(&self, geometry_id: u32) -> Option<&ClippedShape> {
        self.shapes
            .binary_search_by_key(&geometry_id, |s| s.geometry_id)
            .ok()
            .map(|idx| &self.shapes[idx])
    }

    /// Finds the mutable ClippedShape for the given geometry_id.
    pub fn find_shape_mut(&mut self, geometry_id: u32) -> Option<&mut ClippedShape> {
        self.shapes
            .binary_search_by_key(&geometry_id, |s| s.geometry_id)
            .ok()
            .map(|idx| &mut self.shapes[idx])
    }

    /// Adds a shape to the cell, maintaining sorted order by geometry_id.
    ///
    /// If a shape with the same geometry_id already exists, it is replaced.
    pub fn add_shape(&mut self, shape: ClippedShape) {
        match self
            .shapes
            .binary_search_by_key(&shape.geometry_id, |s| s.geometry_id)
        {
            Ok(idx) => self.shapes[idx] = shape,        // Replace existing
            Err(idx) => self.shapes.insert(idx, shape), // Insert at sorted position
        }
    }

    /// Gets or creates a ClippedShape for the given geometry_id.
    ///
    /// If no shape exists, creates one with `contains_center = false`
    /// and no edges.
    pub fn get_or_create_shape(&mut self, geometry_id: u32) -> &mut ClippedShape {
        match self.shapes.binary_search_by_key(&geometry_id, |s| s.geometry_id) {
            Ok(idx) => &mut self.shapes[idx],
            Err(idx) => {
                self.shapes.insert(idx, ClippedShape::new(geometry_id, false));
                &mut self.shapes[idx]
            }
        }
    }

    /// Removes the shape with the given geometry_id.
    ///
    /// Returns the removed shape, or `None` if not found.
    pub fn remove_shape(&mut self, geometry_id: u32) -> Option<ClippedShape> {
        self.shapes
            .binary_search_by_key(&geometry_id, |s| s.geometry_id)
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

    /// Returns an iterator over (geometry_id, contains_center) pairs.
    pub fn doc_ids(&self) -> impl Iterator<Item = (u32, bool)> + '_ {
        self.shapes.iter().map(|s| (s.geometry_id, s.contains_center))
    }

    /// Sorts shapes by geometry_id (call after bulk insertion).
    ///
    /// This is useful when shapes were added without maintaining sort order.
    pub fn sort_shapes(&mut self) {
        self.shapes.sort_by_key(|s| s.geometry_id);
    }
}

impl<C: Copy> PartialEq for Cell<C> {
    fn eq(&self, other: &Self) -> bool {
        self.cell_id == other.cell_id && self.shapes == other.shapes
    }
}

impl<C: Copy> Eq for Cell<C> {}
