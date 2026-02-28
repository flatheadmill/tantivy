use std::io::{self, Read, Write};

use std::collections::HashMap;

use crate::spatial::quadtree::{QuadtreeCell, QuadtreeCellId};



// =============================================================================
// Collapse Candidate
// =============================================================================

/// A collapse candidate is a parent cell whose children could be merged.
///
/// During merge operations, if all four children of a cell have few edges
/// combined, they could be collapsed into the parent. This is recorded
/// as advisory information for the next merge cycle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CollapseCandidate {
    /// Parent cell ID (the cell that would exist after collapse).
    pub parent_id: u64,

    /// Combined edge count of all four children.
    pub combined_edges: u16,

    /// Flags (reserved for future use).
    pub flags: u16,
}

impl CollapseCandidate {
    /// Size of a serialized collapse candidate in bytes.
    pub const SIZE: usize = 12;

    /// Creates a new collapse candidate.
    pub fn new(parent_id: QuadtreeCellId, combined_edges: u16) -> Self {
        Self {
            parent_id: parent_id.to_raw(),
            combined_edges,
            flags: 0,
        }
    }

    /// Returns the parent cell ID.
    pub fn parent_cell_id(&self) -> QuadtreeCellId {
        QuadtreeCellId::from_raw(self.parent_id)
    }

    /// Writes to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.parent_id.to_le_bytes())?;
        writer.write_all(&self.combined_edges.to_le_bytes())?;
        writer.write_all(&self.flags.to_le_bytes())?;
        Ok(())
    }

    /// Reads from a reader.
    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut buf8 = [0u8; 8];
        let mut buf2 = [0u8; 2];

        reader.read_exact(&mut buf8)?;
        let parent_id = u64::from_le_bytes(buf8);

        reader.read_exact(&mut buf2)?;
        let combined_edges = u16::from_le_bytes(buf2);

        reader.read_exact(&mut buf2)?;
        let flags = u16::from_le_bytes(buf2);

        Ok(Self {
            parent_id,
            combined_edges,
            flags,
        })
    }
}

// =============================================================================
// Collapse Detector
// =============================================================================

/// Detects cells that are candidates for collapse during merge.
///
/// As cells are observed during a merge operation, the CollapseDetector
/// tracks sibling groups. When all four children of a parent are seen
/// and their combined edge count is below the threshold, the parent
/// is recorded as a collapse candidate.
///
/// # Usage
///
/// ```ignore
/// let mut detector = CollapseDetector::new(threshold);
///
/// for cell in cells_in_order {
///     detector.observe_cell(&cell);
/// }
///
/// let candidates = detector.finish();
/// ```
#[derive(Debug)]
pub struct CollapseDetector {
    /// Tracks sibling groups as we see cells.
    /// Key: parent cell ID raw value
    /// Value: (children_seen_mask, combined_edges)
    sibling_groups: HashMap<u64, (u8, usize)>,

    /// Completed candidates.
    candidates: Vec<CollapseCandidate>,

    /// Threshold for collapse (combined edges must be <= this).
    collapse_threshold: usize,
}

impl CollapseDetector {
    /// Creates a new CollapseDetector.
    ///
    /// # Arguments
    ///
    /// * `collapse_threshold` - Maximum combined edges for collapse eligibility
    pub fn new(collapse_threshold: usize) -> Self {
        Self {
            sibling_groups: HashMap::new(),
            candidates: Vec::new(),
            collapse_threshold,
        }
    }

    /// Creates a detector that never detects collapse (threshold = 0).
    pub fn disabled() -> Self {
        Self::new(0)
    }

    /// Returns the collapse threshold.
    pub fn threshold(&self) -> usize {
        self.collapse_threshold
    }

    /// Observes a cell during merge traversal.
    ///
    /// Call this for each cell in cell_id order. When all four siblings
    /// of a parent are seen and qualify for collapse, the parent is
    /// recorded as a candidate.
    pub fn observe_cell(&mut self, cell: &QuadtreeCell) {
        if self.collapse_threshold == 0 {
            return; // Disabled
        }

        let cell_id = cell.cell_id();
        let parent = match cell_id.parent() {
            Some(p) => p,
            None => return, // Root cell has no parent
        };

        let child_index = cell_id.child_position() as u8;
        let edge_count = cell.total_edges();

        let parent_raw = parent.to_raw();
        let entry = self.sibling_groups.entry(parent_raw).or_insert((0, 0));

        entry.0 |= 1 << child_index;
        entry.1 += edge_count;

        // All four siblings seen?
        if entry.0 == 0b1111 {
            let combined = entry.1;
            self.sibling_groups.remove(&parent_raw);

            if combined <= self.collapse_threshold {
                self.candidates
                    .push(CollapseCandidate::new(parent, combined as u16));
            }
        }
    }

    /// Observes a cell by its components.
    ///
    /// This is useful when you have cell_id and edge count separately.
    pub fn observe(&mut self, cell_id: QuadtreeCellId, edge_count: usize) {
        if self.collapse_threshold == 0 {
            return;
        }

        let parent = match cell_id.parent() {
            Some(p) => p,
            None => return,
        };

        let child_index = cell_id.child_position() as u8;
        let parent_raw = parent.to_raw();
        let entry = self.sibling_groups.entry(parent_raw).or_insert((0, 0));

        entry.0 |= 1 << child_index;
        entry.1 += edge_count;

        if entry.0 == 0b1111 {
            let combined = entry.1;
            self.sibling_groups.remove(&parent_raw);

            if combined <= self.collapse_threshold {
                self.candidates
                    .push(CollapseCandidate::new(parent, combined as u16));
            }
        }
    }

    /// Returns the number of candidates detected so far.
    pub fn num_candidates(&self) -> usize {
        self.candidates.len()
    }

    /// Returns the number of pending sibling groups.
    pub fn num_pending(&self) -> usize {
        self.sibling_groups.len()
    }

    /// Finishes detection and returns the candidates.
    ///
    /// Pending sibling groups (incomplete) are discarded.
    pub fn finish(self) -> Vec<CollapseCandidate> {
        self.candidates
    }

    /// Clears all state, keeping the threshold.
    pub fn reset(&mut self) {
        self.sibling_groups.clear();
        self.candidates.clear();
    }
}


#[cfg(test)]
mod serialization_tests {
    use std::io::Cursor;

    use super::*;


    // -------------------------------------------------------------------------
    // Collapse Candidate Tests
    // -------------------------------------------------------------------------

    mod collapse_candidate_tests {
        use super::*;

        #[test]
        fn test_collapse_candidate_roundtrip() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 5, &bounds);
            let parent_id = cell_id.parent().unwrap();

            let candidate = CollapseCandidate::new(parent_id, 15);

            let mut buf = Vec::new();
            candidate.write_to(&mut buf).unwrap();

            assert_eq!(buf.len(), CollapseCandidate::SIZE);

            let mut cursor = Cursor::new(&buf);
            let decoded = CollapseCandidate::read_from(&mut cursor).unwrap();

            assert_eq!(decoded, candidate);
            assert_eq!(decoded.parent_cell_id(), parent_id);
            assert_eq!(decoded.combined_edges, 15);
        }
    }

    // -------------------------------------------------------------------------
    // Collapse Detector Tests
    // -------------------------------------------------------------------------

    mod collapse_detector_tests {
        use super::*;

        fn make_cell(bounds: &Bounds, x: f64, y: f64, level: u8, num_edges: usize) -> QuadtreeCell {
            let cell_id = QuadtreeCellId::from_xy(x, y, level, bounds);
            let mut cell = QuadtreeCell::new(cell_id);

            if num_edges > 0 {
                let mut shape = ClippedShape::new(1, false);
                for i in 0..num_edges {
                    shape.add_edge(i as u16);
                }
                cell.add_shape(shape);
            }

            cell
        }

        #[test]
        fn test_detector_no_collapse() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let mut detector = CollapseDetector::new(10);

            // Create four siblings with too many edges
            let cells = [
                make_cell(&bounds, 25.0, 25.0, 2, 5), // child 0
                make_cell(&bounds, 75.0, 25.0, 2, 5), // child 1
                make_cell(&bounds, 25.0, 75.0, 2, 5), // child 2
                make_cell(&bounds, 75.0, 75.0, 2, 5), // child 3
            ];

            for cell in &cells {
                detector.observe_cell(cell);
            }

            let candidates = detector.finish();
            assert!(candidates.is_empty()); // 20 edges > threshold of 10
        }

        #[test]
        fn test_detector_collapse_detected() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let mut detector = CollapseDetector::new(10);

            // Create four siblings with few edges
            // These coordinates map to the four children of the level-1 cell at (0,0):
            // - (12.5, 12.5) -> level-2 cell (0, 0)
            // - (37.5, 12.5) -> level-2 cell (1, 0)
            // - (12.5, 37.5) -> level-2 cell (0, 1)
            // - (37.5, 37.5) -> level-2 cell (1, 1)
            // All share parent (0, 0) at level 1
            let cells = [
                make_cell(&bounds, 12.5, 12.5, 2, 2),
                make_cell(&bounds, 37.5, 12.5, 2, 2),
                make_cell(&bounds, 12.5, 37.5, 2, 2),
                make_cell(&bounds, 37.5, 37.5, 2, 2),
            ];

            for cell in &cells {
                detector.observe_cell(cell);
            }

            let candidates = detector.finish();
            assert_eq!(candidates.len(), 1); // 8 edges <= threshold of 10
            assert_eq!(candidates[0].combined_edges, 8);
        }

        #[test]
        fn test_detector_disabled() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let mut detector = CollapseDetector::disabled();

            let cells = [
                make_cell(&bounds, 25.0, 25.0, 2, 1),
                make_cell(&bounds, 75.0, 25.0, 2, 1),
                make_cell(&bounds, 25.0, 75.0, 2, 1),
                make_cell(&bounds, 75.0, 75.0, 2, 1),
            ];

            for cell in &cells {
                detector.observe_cell(cell);
            }

            let candidates = detector.finish();
            assert!(candidates.is_empty()); // Disabled, never detects
        }

        #[test]
        fn test_detector_incomplete_siblings() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let mut detector = CollapseDetector::new(10);

            // Only three siblings - incomplete
            // These are three of the four children of level-1 cell (0,0)
            let cells = [
                make_cell(&bounds, 12.5, 12.5, 2, 1), // level-2 cell (0, 0)
                make_cell(&bounds, 37.5, 12.5, 2, 1), // level-2 cell (1, 0)
                make_cell(&bounds, 12.5, 37.5, 2, 1), /* level-2 cell (0, 1)
                                                       * Missing: (37.5, 37.5) which would be
                                                       * level-2 cell (1, 1) */
            ];

            for cell in &cells {
                detector.observe_cell(cell);
            }

            assert_eq!(detector.num_pending(), 1); // One incomplete group

            let candidates = detector.finish();
            assert!(candidates.is_empty()); // Incomplete, no collapse
        }
    }
}
