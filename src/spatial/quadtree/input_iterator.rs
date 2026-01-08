//! HUSH

use crate::spatial::quadtree::geometry_server::GeometryServer;
use crate::spatial::quadtree::segment::{CellView, QuadtreeSegment};
use crate::spatial::quadtree::{Bounds, ClippedShape, QuadtreeCell, QuadtreeCellId};
use crate::DocId;

// =============================================================================
// InputIterator - Wrapper for Delete Filtering and Doc ID Remapping
// =============================================================================

pub struct InputIterator<'a> {
    segment: &'a QuadtreeSegment<'a>,
    index: u32,
    doc_id_map: &'a Vec<DocId>,
    current: Option<QuadtreeCell>,
    geometry_server: GeometryServer,
}

impl<'a> InputIterator<'a> {
    pub fn new(
    segment: &'a QuadtreeSegment<'a>,
    doc_id_map: &'a Vec<DocId>,
    geometry_server: GeometryServer,
    ) -> Self {
        let mut iter = Self {
            segment,
            index: 0,
            doc_id_map,
            current: None,
            geometry_server,
        };
        iter.advance();
        iter
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
            let old_geometry_id = shape.geometry_id();
            let old_doc_id = self.geometry_server.get_doc_id(old_geometry_id);

            if self.doc_id_map[old_doc_id as usize] != u32::MAX {
                let mut new_shape = ClippedShape::with_capacity(
                    old_geometry_id,
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
