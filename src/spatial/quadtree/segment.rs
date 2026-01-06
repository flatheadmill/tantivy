//! HUSH
use std::io;

use crate::spatial::quadtree::{Bounds, QuadtreeCellId};

const FOOTER_SIZE: usize = 48;
const MAGIC: u16 = 0x5154; // "QT"
const VERSION: u16 = 1;
const CELL_INDEX_ENTRY_SIZE: usize = 16;

pub struct QuadtreeSegment<'a> {
    data: &'a [u8],
    cell_index_offset: u64,
    cell_count: u32,
    global_bounds: Bounds,
}

impl<'a> QuadtreeSegment<'a> {
    pub fn new(data: &'a [u8]) -> io::Result<Self> {
        if data.len() < FOOTER_SIZE {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "data too short"));
        }

        let footer_start = data.len() - FOOTER_SIZE;
        let footer = &data[footer_start..];

        let magic = u16::from_le_bytes(footer[46..48].try_into().unwrap());
        if magic != MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "bad magic"));
        }

        let version = u16::from_le_bytes(footer[44..46].try_into().unwrap());
        if version != VERSION {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "bad version"));
        }

        let cell_count = u32::from_le_bytes(footer[0..4].try_into().unwrap());
        let cell_index_offset = u64::from_le_bytes(footer[4..12].try_into().unwrap());

        let min_x = f64::from_le_bytes(footer[12..20].try_into().unwrap());
        let min_y = f64::from_le_bytes(footer[20..28].try_into().unwrap());
        let max_x = f64::from_le_bytes(footer[28..36].try_into().unwrap());
        let max_y = f64::from_le_bytes(footer[36..44].try_into().unwrap());
        let global_bounds = Bounds::new(min_x, min_y, max_x, max_y);

        Ok(Self {
            data,
            cell_index_offset,
            cell_count,
            global_bounds,
        })
    }

    pub fn global_bounds(&self) -> &Bounds {
        &self.global_bounds
    }

    pub fn cell_count(&self) -> u32 {
        self.cell_count
    }

    fn cell_index_entry(&self, index: u32) -> (u64, u64) {
        let offset = self.cell_index_offset as usize + (index as usize) * CELL_INDEX_ENTRY_SIZE;
        let entry = &self.data[offset..offset + CELL_INDEX_ENTRY_SIZE];
        let cell_id = u64::from_le_bytes(entry[0..8].try_into().unwrap());
        let data_offset = u64::from_le_bytes(entry[8..16].try_into().unwrap());
        (cell_id, data_offset)
    }

    pub fn find_cell(&self, target: QuadtreeCellId) -> Option<CellView<'a>> {
        if self.cell_count == 0 {
            return None;
        }

        let target_raw = target.to_raw();
        let mut lo: u32 = 0;
        let mut hi: u32 = self.cell_count;

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let (cell_id, offset) = self.cell_index_entry(mid);

            match cell_id.cmp(&target_raw) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    let next_offset = if mid + 1 < self.cell_count {
                        self.cell_index_entry(mid + 1).1
                    } else {
                        self.cell_index_offset
                    };
                    let len = (next_offset - offset) as usize;
                    return Some(CellView::new(
                        QuadtreeCellId::from_raw(cell_id),
                        &self.data[offset as usize..offset as usize + len],
                    ));
                }
            }
        }
        None
    }
}

pub struct CellView<'a> {
    cell_id: QuadtreeCellId,
    data: &'a [u8],
}

impl<'a> CellView<'a> {
    fn new(cell_id: QuadtreeCellId, data: &'a [u8]) -> Self {
        Self { cell_id, data }
    }

    pub fn cell_id(&self) -> QuadtreeCellId {
        self.cell_id
    }

    pub fn num_shapes(&self) -> u32 {
        u32::from_le_bytes(self.data[0..4].try_into().unwrap())
    }

    pub fn shapes(&self) -> ShapeIterator<'a> {
        ShapeIterator {
            data: &self.data[4..],
            remaining: self.num_shapes(),
        }
    }
}

pub struct ShapeIterator<'a> {
    data: &'a [u8],
    remaining: u32,
}

impl<'a> Iterator for ShapeIterator<'a> {
    type Item = ShapeView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }

        let doc_id = u32::from_le_bytes(self.data[0..4].try_into().unwrap());
        let flags = self.data[4];
        let contains_center = (flags & 1) != 0;
        let num_edges = u16::from_le_bytes(self.data[5..7].try_into().unwrap());

        let edges_size = (num_edges as usize) * 2;
        let shape_size = 7 + edges_size;

        let edges_data = &self.data[7..shape_size];
        self.data = &self.data[shape_size..];
        self.remaining -= 1;

        Some(ShapeView {
            doc_id,
            contains_center,
            num_edges,
            edges_data,
        })
    }
}

pub struct ShapeView<'a> {
    doc_id: u32,
    contains_center: bool,
    num_edges: u16,
    edges_data: &'a [u8],
}

impl<'a> ShapeView<'a> {
    pub fn doc_id(&self) -> u32 {
        self.doc_id
    }
    pub fn contains_center(&self) -> bool {
        self.contains_center
    }
    pub fn num_edges(&self) -> u16 {
        self.num_edges
    }

    pub fn edge(&self, i: usize) -> u16 {
        u16::from_le_bytes(self.edges_data[i * 2..i * 2 + 2].try_into().unwrap())
    }

    pub fn edges(&self) -> impl Iterator<Item = u16> + '_ {
        (0..self.num_edges as usize).map(|i| self.edge(i))
    }
}
