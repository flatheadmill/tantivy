//! HUSH

use crate::spatial::quadtree::segment::{CellView, QuadtreeSegment};
use crate::spatial::quadtree::{Bounds, ClippedShape, QuadtreeCell, QuadtreeCellId};

// =============================================================================
// DocIdMap - Document ID Remapping
// =============================================================================

#[derive(Debug, Clone)]
pub struct DocIdMap {
    mapping: Vec<Option<u32>>,
}

impl DocIdMap {
    pub fn identity(num_docs: u32) -> Self {
        Self {
            mapping: (0..num_docs).map(Some).collect(),
        }
    }

    pub fn with_offset(num_docs: u32, offset: u32) -> Self {
        Self {
            mapping: (0..num_docs).map(|id| Some(id + offset)).collect(),
        }
    }

    pub fn empty() -> Self {
        Self {
            mapping: Vec::new(),
        }
    }

    pub fn set(&mut self, old_id: u32, new_id: Option<u32>) {
        let idx = old_id as usize;
        if idx >= self.mapping.len() {
            self.mapping.resize(idx + 1, None);
        }
        self.mapping[idx] = new_id;
    }

    pub fn get(&self, old_id: u32) -> Option<u32> {
        self.mapping.get(old_id as usize).copied().flatten()
    }

    pub fn is_deleted(&self, old_id: u32) -> bool {
        self.get(old_id).is_none()
    }

    pub fn len(&self) -> usize {
        self.mapping.len()
    }

    pub fn is_empty(&self) -> bool {
        self.mapping.is_empty()
    }
}

// =============================================================================
// DeleteBitSet - Deleted Document Tracking
// =============================================================================

#[derive(Debug, Clone)]
pub struct DeleteBitSet {
    bits: Vec<bool>,
}

impl DeleteBitSet {
    pub fn new() -> Self {
        Self { bits: Vec::new() }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            bits: vec![false; capacity],
        }
    }

    pub fn delete(&mut self, doc_id: u32) {
        let idx = doc_id as usize;
        if idx >= self.bits.len() {
            self.bits.resize(idx + 1, false);
        }
        self.bits[idx] = true;
    }

    pub fn is_deleted(&self, doc_id: u32) -> bool {
        self.bits.get(doc_id as usize).copied().unwrap_or(false)
    }

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

pub struct InputIterator<'a> {
    segment: &'a QuadtreeSegment<'a>,
    index: u32,
    deletes: DeleteBitSet,
    doc_id_map: DocIdMap,
    current: Option<QuadtreeCell>,
}

impl<'a> InputIterator<'a> {
    pub fn new(
        segment: &'a QuadtreeSegment<'a>,
        deletes: DeleteBitSet,
        doc_id_map: DocIdMap,
    ) -> Self {
        let mut iter = Self {
            segment,
            index: 0,
            deletes,
            doc_id_map,
            current: None,
        };
        iter.advance();
        iter
    }

    pub fn unfiltered(segment: &'a QuadtreeSegment<'a>, num_docs: u32) -> Self {
        Self::new(segment, DeleteBitSet::new(), DocIdMap::identity(num_docs))
    }

    pub fn current(&self) -> Option<&QuadtreeCell> {
        self.current.as_ref()
    }

    pub fn current_cell_id(&self) -> Option<QuadtreeCellId> {
        self.current.as_ref().map(|c| c.cell_id())
    }

    pub fn has_more(&self) -> bool {
        self.current.is_some()
    }

    pub fn advance(&mut self) {
        loop {
            if self.index >= self.segment.cell_count() {
                self.current = None;
                return;
            }

            if let Some(view) = self.segment.get_cell_by_index(self.index) {
                self.index += 1;
                let filtered = self.filter_cell(view);
                if !filtered.is_empty() {
                    self.current = Some(filtered);
                    return;
                }
            } else {
                self.current = None;
                return;
            }
        }
    }

    pub fn take(&mut self) -> Option<QuadtreeCell> {
        let cell = self.current.take();
        if cell.is_some() {
            self.advance();
        }
        cell
    }

    fn filter_cell(&self, view: CellView) -> QuadtreeCell {
        let mut filtered = QuadtreeCell::new(view.cell_id());

        for shape in view.shapes() {
            let old_id = shape.doc_id();

            if self.deletes.is_deleted(old_id) {
                continue;
            }

            if let Some(new_id) = self.doc_id_map.get(old_id) {
                let mut new_shape = ClippedShape::with_capacity(
                    new_id,
                    shape.contains_center(),
                    shape.num_edges() as usize,
                );
                for edge in shape.edges() {
                    new_shape.add_edge(edge);
                }
                filtered.add_shape(new_shape);
            }
        }

        filtered
    }

    pub fn bounds(&self) -> Bounds {
        self.segment.global_bounds().clone()
    }
}

#[cfg(test)]
mod tests {
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
            map.set(10, None);

            assert_eq!(map.get(5), Some(10));
            assert!(map.is_deleted(10));
            assert!(map.is_deleted(100));
        }
    }

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
