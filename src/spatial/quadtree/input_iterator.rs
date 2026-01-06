//! HUSH

// =============================================================================
// DocIdMap - Document ID Remapping
// =============================================================================

use std::io::{Read, Seek};
use crate::spatial::quadtree::{serialization::CellReader, Bounds, ClippedShape, QuadtreeCell, QuadtreeCellId};

/// Maps document IDs from a source segment to a merged segment.
///
/// During segment merge, each input segment's doc_ids are remapped
/// to the merged segment's doc_id space.
#[derive(Debug, Clone)]
pub struct DocIdMap {
    /// Maps old doc_id to new doc_id.
    /// If an entry is None, the document was deleted.
    mapping: Vec<Option<u32>>,
}

impl DocIdMap {
    /// Creates a new identity mapping for the given number of documents.
    pub fn identity(num_docs: u32) -> Self {
        Self {
            mapping: (0..num_docs).map(Some).collect(),
        }
    }

    /// Creates a mapping with an offset.
    ///
    /// Each doc_id is mapped to doc_id + offset.
    pub fn with_offset(num_docs: u32, offset: u32) -> Self {
        Self {
            mapping: (0..num_docs).map(|id| Some(id + offset)).collect(),
        }
    }

    /// Creates an empty mapping.
    pub fn empty() -> Self {
        Self {
            mapping: Vec::new(),
        }
    }

    /// Sets the mapping for a document.
    pub fn set(&mut self, old_id: u32, new_id: Option<u32>) {
        let idx = old_id as usize;
        if idx >= self.mapping.len() {
            self.mapping.resize(idx + 1, None);
        }
        self.mapping[idx] = new_id;
    }

    /// Gets the new doc_id for an old doc_id.
    ///
    /// Returns None if the document was deleted or not mapped.
    pub fn get(&self, old_id: u32) -> Option<u32> {
        self.mapping.get(old_id as usize).copied().flatten()
    }

    /// Returns true if the document is deleted (mapped to None).
    pub fn is_deleted(&self, old_id: u32) -> bool {
        self.get(old_id).is_none()
    }

    /// Returns the number of mappings.
    pub fn len(&self) -> usize {
        self.mapping.len()
    }

    /// Returns true if empty.
    pub fn is_empty(&self) -> bool {
        self.mapping.is_empty()
    }
}

// =============================================================================
// DeleteBitSet - Deleted Document Tracking
// =============================================================================

/// A simple bitset for tracking deleted documents.
#[derive(Debug, Clone)]
pub struct DeleteBitSet {
    /// Bits, one per document. True = deleted.
    bits: Vec<bool>,
}

impl DeleteBitSet {
    /// Creates an empty bitset.
    pub fn new() -> Self {
        Self { bits: Vec::new() }
    }

    /// Creates a bitset with the given capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            bits: vec![false; capacity],
        }
    }

    /// Marks a document as deleted.
    pub fn delete(&mut self, doc_id: u32) {
        let idx = doc_id as usize;
        if idx >= self.bits.len() {
            self.bits.resize(idx + 1, false);
        }
        self.bits[idx] = true;
    }

    /// Returns true if the document is deleted.
    pub fn is_deleted(&self, doc_id: u32) -> bool {
        self.bits.get(doc_id as usize).copied().unwrap_or(false)
    }

    /// Returns the number of deleted documents.
    pub fn num_deleted(&self) -> usize {
        self.bits.iter().filter(|&&b| b).count()
    }
}

impl Default for DeleteBitSet {
    fn default() -> Self {
        Self::new()
    }
}

// =============================================================================
// InputIterator - Wrapper for Delete Filtering and Doc ID Remapping
// =============================================================================

/// Wraps a CellReader with delete filtering and doc_id remapping.
///
/// This iterator:
/// 1. Filters out shapes from deleted documents
/// 2. Remaps doc_ids to the merged segment's doc_id space
/// 3. Skips cells that become empty after filtering
pub struct InputIterator<R: Read + Seek> {
    /// Underlying cell reader.
    reader: CellReader<R>,

    /// Delete bitmap.
    deletes: DeleteBitSet,

    /// Document ID mapping.
    doc_id_map: DocIdMap,

    /// Current cell (after filtering).
    current: Option<QuadtreeCell>,

    /// Whether an error has occurred.
    errored: bool,
}

impl<R: Read + Seek> InputIterator<R> {
    /// Creates a new InputIterator.
    pub fn new(reader: CellReader<R>, deletes: DeleteBitSet, doc_id_map: DocIdMap) -> Self {
        let mut iter = Self {
            reader,
            deletes,
            doc_id_map,
            current: None,
            errored: false,
        };
        iter.advance();
        iter
    }

    /// Creates an InputIterator with no filtering (identity mapping).
    pub fn unfiltered(reader: CellReader<R>, num_docs: u32) -> Self {
        Self::new(reader, DeleteBitSet::new(), DocIdMap::identity(num_docs))
    }

    /// Returns the current cell, if any.
    pub fn current(&self) -> Option<&QuadtreeCell> {
        self.current.as_ref()
    }

    /// Returns the cell_id of the current cell, if any.
    pub fn current_cell_id(&self) -> Option<QuadtreeCellId> {
        self.current.as_ref().map(|c| c.cell_id())
    }

    /// Returns true if there are more cells.
    pub fn has_more(&self) -> bool {
        self.current.is_some()
    }

    /// Advances to the next non-empty cell.
    pub fn advance(&mut self) {
        if self.errored {
            self.current = None;
            return;
        }

        loop {
            match self.reader.next() {
                Some(Ok(cell)) => {
                    let filtered = self.filter_cell(cell);
                    if !filtered.is_empty() {
                        self.current = Some(filtered);
                        return;
                    }
                    // Cell became empty - skip it
                }
                Some(Err(_)) => {
                    self.errored = true;
                    self.current = None;
                    return;
                }
                None => {
                    self.current = None;
                    return;
                }
            }
        }
    }

    /// Takes the current cell and advances to the next.
    pub fn take(&mut self) -> Option<QuadtreeCell> {
        let cell = self.current.take();
        if cell.is_some() {
            self.advance();
        }
        cell
    }

    /// Filters a cell by removing deleted shapes and remapping doc_ids.
    fn filter_cell(&self, cell: QuadtreeCell) -> QuadtreeCell {
        let mut filtered = QuadtreeCell::new(cell.cell_id());

        for shape in cell.shapes() {
            let old_id = shape.doc_id();

            // Skip deleted documents
            if self.deletes.is_deleted(old_id) {
                continue;
            }

            // Remap doc_id
            if let Some(new_id) = self.doc_id_map.get(old_id) {
                let mut new_shape =
                    ClippedShape::with_capacity(new_id, shape.contains_center(), shape.num_edges());
                for &edge in shape.edges() {
                    new_shape.add_edge(edge);
                }
                filtered.add_shape(new_shape);
            }
            // If doc_id not in map, skip the shape
        }

        filtered
    }

    /// Returns the bounds from the underlying reader.
    pub fn bounds(&self) -> Bounds {
        self.reader.bounds()
    }
}

#[cfg(test)]
mod tests {
    // -------------------------------------------------------------------------
    // DocIdMap Tests
    // -------------------------------------------------------------------------
        use super::*;

    mod doc_id_map_tests {
        use super::*;

        #[test]
        fn test_identity() {
            let map = DocIdMap::identity(10);
            for i in 0..10 {
                assert_eq!(map.get(i), Some(i));
            }
            assert_eq!(map.get(10), None);
        }

        #[test]
        fn test_with_offset() {
            let map = DocIdMap::with_offset(5, 100);
            assert_eq!(map.get(0), Some(100));
            assert_eq!(map.get(4), Some(104));
        }

        #[test]
        fn test_set() {
            let mut map = DocIdMap::empty();
            map.set(5, Some(10));
            map.set(10, None); // Deleted

            assert_eq!(map.get(5), Some(10));
            assert!(map.is_deleted(10));
            assert!(map.is_deleted(100)); // Not mapped
        }
    }

    // -------------------------------------------------------------------------
    // DeleteBitSet Tests
    // -------------------------------------------------------------------------

    mod delete_bitset_tests {
        use super::*;

        #[test]
        fn test_delete_and_check() {
            let mut bits = DeleteBitSet::new();
            assert!(!bits.is_deleted(5));

            bits.delete(5);
            assert!(bits.is_deleted(5));
            assert!(!bits.is_deleted(4));
            assert!(!bits.is_deleted(6));
        }

        #[test]
        fn test_num_deleted() {
            let mut bits = DeleteBitSet::new();
            bits.delete(1);
            bits.delete(5);
            bits.delete(10);

            assert_eq!(bits.num_deleted(), 3);
        }
    }
}
