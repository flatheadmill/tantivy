//! Block kd-tree spatial indexing for triangulated polygons.
//!
//! Implements an immutable bulk-loaded spatial index using recursive median partitioning on
//! bounding box dimensions. Each leaf stores up to 512 triangles with delta-compressed coordinates
//! and doc IDs. The tree provides three query types (intersects, within, contains) that use exact
//! integer arithmetic for geometric predicates and accumulate results in bit sets for efficient
//! deduplication across leaves.
//!
//! The serialized format stores compressed leaf pages followed by the tree structure (leaf and
//! branch nodes), enabling zero-copy access through memory-mapped segments without upfront
//! decompression.
use std::io;
use std::io::Write;
use std::marker::PhantomData;

use common::{BitSet, CountingWriter};

use crate::directory::WritePtr;
use crate::spatial::envelope::{Bounds, Envelope, LeafCompression, Spatial};

pub(crate) struct SpreadSurvey {
    min: f64,
    max: f64,
}

impl SpreadSurvey {
    pub(crate) fn survey(&mut self, value: f64) {
        self.min = self.min.min(value);
        self.max = self.max.max(value);
    }
    pub(crate) fn spread(&self) -> f64 {
        self.max - self.min
    }
}

impl Default for SpreadSurvey {
    fn default() -> Self {
        SpreadSurvey {
            min: f64::MAX,
            max: f64::MIN,
        }
    }
}

struct EnvelopeSurvey<S> {
    bounds: Vec<f64>,
    _marker: PhantomData<S>,
}

impl<S: Spatial> EnvelopeSurvey<S> {
    fn survey(&mut self, envelope: &Envelope<S::Bounds>) {
        let dims = S::DIMENSIONS;
        for i in 0..dims {
            self.bounds[i] = self.bounds[i].min(envelope.bounds.get(i));
            self.bounds[i + dims] = self.bounds[i + dims].max(envelope.bounds.get(i + dims));
        }
    }
    fn bounds(&self) -> Vec<f64> {
        self.bounds.clone()
    }
}

impl<S: Spatial> Default for EnvelopeSurvey<S> {
    fn default() -> Self {
        let mut bounds = vec![f64::MAX; S::COORDINATES];
        for i in S::DIMENSIONS..S::COORDINATES {
            bounds[i] = f64::MIN;
        }
        EnvelopeSurvey {
            bounds,
            _marker: PhantomData,
        }
    }
}

enum BuildNode {
    Branch {
        bounds: Vec<f64>,
        left: Box<BuildNode>,
        right: Box<BuildNode>,
    },
    Leaf {
        bounds: Vec<f64>,
        pos: u64,
        len: u16,
    },
}

// Leaf pages are first the count of triangles, followed by delta encoded doc_ids, followed by
// the delta encoded words in order. We will then have the length of the page. We build a tree
// after the pages with leaf nodes and branch nodes. Leaf nodes will contain the bounding box
// of the leaf followed position and length of the page. The leaf node is a level of direction
// to store the position and length of the page in a format that is easy to read directly from
// the mapping.

// We do not compress the tree nodes. We read them directly from the mapping.

//
fn write_leaf_pages<S: Spatial>(
    envelopes: &mut [Envelope<S::Bounds>],
    write: &mut CountingWriter<WritePtr>,
) -> io::Result<BuildNode> {
    // If less than 512 triangles we are at a leaf, otherwise we still in the inner nodes.
    if envelopes.len() <= S::Compression::PAGE_SIZE {
        let mut survey = EnvelopeSurvey::<S>::default();
        for envelope in envelopes.iter() {
            survey.survey(envelope);
        }
        let pos = write.written_bytes();
        S::Compression::compress(envelopes, write)?;
        let len = write.written_bytes() - pos;
        Ok(BuildNode::Leaf {
            bounds: survey.bounds(),
            pos,
            len: len as u16,
        })
    } else {
        let mut spreads: Vec<SpreadSurvey> = (0..S::COORDINATES)
            .map(|_| SpreadSurvey::default())
            .collect();
        let mut survey = EnvelopeSurvey::<S>::default();
        for envelope in envelopes.iter() {
            for i in 0..S::COORDINATES {
                spreads[i].survey(envelope.bounds.get(i));
            }
            survey.survey(envelope);
        }
        let (dimension, _) = spreads
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.spread().total_cmp(&b.spread()))
            .unwrap();

        // Partition the envelopes.
        let mid = envelopes.len() / 2;
        envelopes.select_nth_unstable_by(mid, |a, b| {
            a.bounds.get(dimension).total_cmp(&b.bounds.get(dimension))
        });
        let partition = envelopes[mid].bounds.get(dimension);
        let mut split_point = mid + 1;
        while split_point < envelopes.len()
            && envelopes[split_point].bounds.get(dimension) == partition
        {
            split_point += 1;
        }
        // If we reached the end of triangles then all of the triangles share the partition value
        // for the dimension. We handle this degeneracy by splitting at the midpoint so that we
        // won't have a leaf with zero triangles.
        if split_point == envelopes.len() {
            split_point = mid; // Force split at midpoint index
        } else {
            // Our partition does not sort the triangles, it only partitions. We have scan our right
            // partition to find all the midpoint values and move them to the left partition.
            let mut reverse = envelopes.len() - 1;
            loop {
                // Scan backwards looking for the partition value.
                while envelopes[reverse].bounds.get(dimension) != partition {
                    reverse -= 1;
                }
                // If we have reached the split point then we are done.
                if reverse <= split_point {
                    break;
                }
                // Swap the midpoint value with our current split point.
                envelopes.swap(split_point, reverse);
                // Move the split point up one.
                split_point += 1;
                // We know that what was at the split point was not the midpoint value.
                reverse -= 1;
            }
        }
        // Split into left and write partitions and create child nodes.
        let (left, right) = envelopes.split_at_mut(split_point);
        let left_node = write_leaf_pages::<S>(left, write)?;
        let right_node = write_leaf_pages::<S>(right, write)?;
        // Return an inner node.
        Ok(BuildNode::Branch {
            bounds: survey.bounds(),
            left: Box::new(left_node),
            right: Box::new(right_node),
        })
    }
}

fn write_leaf_nodes<S: Spatial>(
    node: &BuildNode,
    write: &mut CountingWriter<WritePtr>,
) -> io::Result<()> {
    match node {
        BuildNode::Branch {
            bounds: _,
            left,
            right,
        } => {
            write_leaf_nodes::<S>(right, write)?;
            write_leaf_nodes::<S>(left, write)?;
        }
        BuildNode::Leaf { bounds, pos, len } => {
            for &dimension in bounds.iter() {
                write.write_all(dimension.to_le_bytes().as_ref())?;
            }
            write.write_all(&pos.to_le_bytes())?;
            write.write_all(&len.to_le_bytes())?;
        }
    }
    Ok(())
}

fn leaf_node_size<S: Spatial>() -> usize {
    S::COORDINATES * 8 + 10 // bounds + pos(8) + len(2)
}

fn branch_node_size<S: Spatial>() -> usize {
    S::COORDINATES * 8 + 8 // bounds + left(4) + right(4)
}

fn write_branch_nodes<S: Spatial>(
    node: &BuildNode,
    branch_offset: &mut i32,
    leaf_offset: &mut i32,
    write: &mut CountingWriter<WritePtr>,
) -> io::Result<i32> {
    match node {
        BuildNode::Leaf { .. } => {
            let pos = *leaf_offset;
            *leaf_offset -= 1;
            Ok(pos * leaf_node_size::<S>() as i32)
        }
        BuildNode::Branch {
            bounds,
            left,
            right,
        } => {
            let left = write_branch_nodes::<S>(left, branch_offset, leaf_offset, write)?;
            let right = write_branch_nodes::<S>(right, branch_offset, leaf_offset, write)?;
            for &val in bounds {
                write.write_all(val.to_le_bytes().as_ref())?;
            }
            write.write_all(&left.to_le_bytes())?;
            write.write_all(&right.to_le_bytes())?;
            let pos = *branch_offset;
            *branch_offset += 1;
            Ok(pos * branch_node_size::<S>() as i32)
        }
    }
}

const VERSION: u16 = 1u16;

/// Builds and serializes a block kd-tree for spatial indexing of triangles.
///
/// Takes a collection of envelopes and constructs a complete block kd-tree, writing both the
/// compressed leaf pages and tree structure to the output. The tree uses recursive median
/// partitioning on the dimension with maximum spread, storing up to 512 envelopes per leaf.
///
/// The output format consists of:
/// - Version header (u16)
/// - Compressed leaf pages (delta-encoded doc_ids and envelopes coordinates)
/// - 32-byte aligned tree structure (leaf nodes, then branch nodes)
/// - Footer with envelopes count, root offset, and branch position
///
/// The `envelopes` slice will be reordered during tree construction as partitioning sorts by the
/// selected dimension at each level.
pub fn write_tree<S: Spatial>(
    envelopes: &mut [Envelope<S::Bounds>],
    write: &mut CountingWriter<WritePtr>,
) -> io::Result<()> {
    write.write_all(&VERSION.to_le_bytes())?;

    let tree = write_leaf_pages::<S>(envelopes, write)?;
    write_leaf_nodes::<S>(&tree, write)?;
    let branch_position = write.written_bytes();
    let mut branch_offset: i32 = 0;
    let mut leaf_offset: i32 = -1;
    let root = write_branch_nodes::<S>(&tree, &mut branch_offset, &mut leaf_offset, write)?;
    write.write_all(&envelopes.len().to_le_bytes())?;
    write.write_all(&root.to_le_bytes())?;
    write.write_all(&branch_position.to_le_bytes())?;
    Ok(())
}

struct BranchNode {
    left: i32,
    right: i32,
}

struct LeafNode {
    pos: u64,
    len: u16,
}

/// A read-only view into a serialized block kd-tree segment.
///
/// Provides access to the tree structure and compressed leaf data through memory-mapped or
/// buffered byte slices. The segment contains compressed leaf pages followed by the tree structure
/// (leaf nodes and branch nodes), with a footer containing metadata for locating the root and
/// interpreting offsets.
pub struct Segment<'a, S: Spatial> {
    data: &'a [u8],
    branch_position: u64,
    /// Offset to the root of the tree, used as the starting point for traversal.
    pub root_offset: i32,
    _marker: PhantomData<S>,
}

impl<'a, S: Spatial> Segment<'a, S> {
    /// Creates a new segment from serialized block kd-tree data.
    ///
    /// Reads the footer metadata from the last 12 bytes to locate the tree structure and root
    /// node.
    pub fn new(data: &'a [u8]) -> Self {
        Segment {
            data,
            branch_position: u64::from_le_bytes(data[data.len() - 8..].try_into().unwrap()),
            root_offset: i32::from_le_bytes(
                data[data.len() - 12..data.len() - 8].try_into().unwrap(),
            ),
            _marker: PhantomData,
        }
    }
    #[inline(always)]
    fn bounds(&self, offset: i32) -> Vec<f64> {
        let sizeof_bounds = S::COORDINATES * 8;
        let byte_offset = (self.branch_position as i64 + offset as i64) as usize;
        let bytes = &self.data[byte_offset..byte_offset + sizeof_bounds];
        let mut bounds = Vec::with_capacity(S::COORDINATES);
        for i in 0..S::COORDINATES {
            let start = i * 8;
            bounds.push(f64::from_le_bytes(
                bytes[start..start + 8].try_into().unwrap(),
            ));
        }
        bounds
    }
    #[inline(always)]
    fn branch_node(&self, offset: i32) -> BranchNode {
        let sizeof_bounds = S::COORDINATES * 8;
        let byte_offset = (self.branch_position as i64 + offset as i64) as usize;
        let child_offset = byte_offset + sizeof_bounds;
        let bytes = &self.data[child_offset..child_offset + 8];
        BranchNode {
            left: i32::from_le_bytes(bytes[0..4].try_into().unwrap()),
            right: i32::from_le_bytes(bytes[4..8].try_into().unwrap()),
        }
    }
    #[inline(always)]
    fn leaf_node(&self, offset: i32) -> LeafNode {
        let sizeof_bounds = S::COORDINATES * 8;
        let byte_offset = (self.branch_position as i64 + offset as i64) as usize;
        let data_offset = byte_offset + sizeof_bounds;
        let bytes = &self.data[data_offset..data_offset + 10];
        LeafNode {
            pos: u64::from_le_bytes(bytes[0..8].try_into().unwrap()),
            len: u16::from_le_bytes(bytes[8..10].try_into().unwrap()),
        }
    }
    fn leaf_page(&self, leaf_node: &LeafNode) -> &[u8] {
        &self.data[(leaf_node.pos as usize)..(leaf_node.pos as usize + leaf_node.len as usize)]
    }
}

fn bounds_within<S: Spatial>(bounds: &[f64], query: &[f64]) -> bool {
    for i in 0..S::DIMENSIONS {
        if bounds[i] < query[i] {
            return false;
        }
        if bounds[i + S::DIMENSIONS] > query[i + S::DIMENSIONS] {
            return false;
        }
    }
    true
}

fn bounds_intersects<S: Spatial>(bounds: &[f64], query: &[f64]) -> bool {
    for i in 0..S::DIMENSIONS {
        if bounds[i + S::DIMENSIONS] < query[i] || bounds[i] > query[i + S::DIMENSIONS] {
            return false; // disjoint in this dimension
        }
    }
    true
}

fn collect_all_docs<S: Spatial>(
    segment: &Segment<S>,
    offset: i32,
    docs: &mut BitSet,
) -> io::Result<()> {
    if offset < 0 {
        let leaf_node = segment.leaf_node(offset);
        let envelopes = S::Compression::decompress(segment.leaf_page(&leaf_node))?;
        for envelope in &envelopes {
            docs.insert(envelope.doc_id);
        }
    } else {
        let branch_node = segment.branch_node(offset);
        collect_all_docs::<S>(segment, branch_node.left, docs)?;
        collect_all_docs::<S>(segment, branch_node.right, docs)?;
    }
    Ok(())
}

/// Finds documents with triangles that intersect the query bounding box.
///
/// Traverses the tree starting at `offset` (typically `segment.root_offset`), pruning subtrees
/// whose bounding boxes don't intersect the query. When a node's bbox is entirely within the
/// query, all its documents are bulk-collected. Otherwise, individual triangles are tested using
/// exact geometric predicates.
///
/// The query is `[min_y, min_x, max_y, max_x]` in integer coordinates. Documents are inserted into
/// the `result` BitSet, which automatically deduplicates when the same document appears in
/// multiple leaves.
pub fn search_intersects<S: Spatial>(
    segment: &Segment<S>,
    offset: i32,
    query: &[f64],
    docs: &mut BitSet,
) -> io::Result<()> {
    let bounds = segment.bounds(offset);
    if !bounds_intersects::<S>(&bounds, query) {
        // skip
    } else if bounds_within::<S>(&bounds, query) {
        collect_all_docs::<S>(segment, offset, docs)?;
    } else if offset < 0 {
        let leaf_node = segment.leaf_node(offset);
        let envelopes = S::Compression::decompress(segment.leaf_page(&leaf_node))?;
        for envelope in &envelopes {
            let env_bounds: Vec<_> = (0..S::COORDINATES)
                .map(|i| envelope.bounds.get(i))
                .collect();
            if bounds_intersects::<S>(&env_bounds, query) {
                docs.insert(envelope.doc_id);
            }
        }
    } else {
        let branch_node = segment.branch_node(offset);
        search_intersects(segment, branch_node.left, query, docs)?;
        search_intersects(segment, branch_node.right, query, docs)?;
    }
    Ok(())
}

/// Finds documents where all triangles are within the query bounding box.
///
/// Traverses the tree starting at `offset` (typically `segment.root_offset`), testing each
/// triangle to determine if it lies entirely within the query bounds. Uses two `BitSet` instances
/// to track state: `include` accumulates candidate documents, while `exclude` marks documents that
/// have at least one triangle extending outside the query.
///
/// The query is `[min_y, min_x, max_y, max_x]` in integer coordinates. The final result is
/// documents in `include` that are NOT in `exclude` - the caller must compute this difference.
pub fn search_within<S: Spatial>(
    segment: &Segment<S>,
    offset: i32,
    query: &[f64],
    docs: &mut BitSet,
) -> io::Result<()> {
    let bounds = segment.bounds(offset);
    if !bounds_intersects::<S>(&bounds, query) {
        // skip
    } else if bounds_within::<S>(&bounds, query) {
        collect_all_docs::<S>(segment, offset, docs)?;
    } else if offset < 0 {
        let leaf_node = segment.leaf_node(offset);
        let envelopes = S::Compression::decompress(segment.leaf_page(&leaf_node))?;
        for envelope in &envelopes {
            let env_bounds: Vec<_> = (0..S::COORDINATES)
                .map(|i| envelope.bounds.get(i))
                .collect();
            if bounds_within::<S>(&env_bounds, query) {
                docs.insert(envelope.doc_id);
            }
        }
    } else {
        let branch_node = segment.branch_node(offset);
        search_within(segment, branch_node.left, query, docs)?;
        search_within(segment, branch_node.right, query, docs)?;
    }
    Ok(())
}

/// Finds documents whose polygons contain the query bounding box.
///
/// Traverses the tree starting at `offset` (typically `segment.root_offset`), testing each
/// triangle using three-state logic: `Candidate` (query might be contained), `NotWithin` (boundary
/// edge crosses query), or `Disjoint` (no overlap). Only boundary edges are tested for crossing -
/// internal tessellation edges are ignored.
///
/// The query is `[min_y, min_x, max_y, max_x]` in integer coordinates. Uses two `BitSet`
/// instances: `include` accumulates candidates, `excluded` marks documents with disqualifying
/// boundary crossings. The final result is documents in `include` that are NOT in `exclude`.
pub fn search_contains<S: Spatial>(
    segment: &Segment<S>,
    offset: i32,
    query: &[f64],
    docs: &mut BitSet,
) -> io::Result<()> {
    let bounds = segment.bounds(offset);
    if !bounds_within::<S>(&bounds, query) {
        // skip
    } else if offset < 0 {
        let leaf_node = segment.leaf_node(offset);
        let envelopes = S::Compression::decompress(segment.leaf_page(&leaf_node))?;
        for envelope in &envelopes {
            let env_bounds: Vec<_> = (0..S::COORDINATES)
                .map(|i| envelope.bounds.get(i))
                .collect();
            if bounds_within::<S>(query, &env_bounds) {
                docs.insert(envelope.doc_id);
            }
        }
    } else {
        let branch_node = segment.branch_node(offset);
        search_contains::<S>(segment, branch_node.left, query, docs)?;
        search_contains::<S>(segment, branch_node.right, query, docs)?;
    }
    Ok(())
}

/// HUSH
pub struct LeafPageIterator<'a, S: Spatial> {
    segment: &'a Segment<'a, S>,
    descent_stack: Vec<i32>,
}

impl<'a, S: Spatial> LeafPageIterator<'a, S> {
    /// HUSH
    pub fn new(segment: &'a Segment<'a, S>) -> Self {
        Self {
            segment,
            descent_stack: vec![segment.root_offset],
        }
    }
}

impl<'a, S: Spatial> Iterator for LeafPageIterator<'a, S> {
    type Item = io::Result<Vec<Envelope<S::Bounds>>>;

    fn next(&mut self) -> Option<Self::Item> {
        let offset = self.descent_stack.pop()?;
        if offset < 0 {
            let leaf_node = self.segment.leaf_node(offset);
            let leaf_page = self.segment.leaf_page(&leaf_node);
            match S::Compression::decompress(leaf_page) {
                Ok(envelopes) => Some(Ok(envelopes)),
                Err(e) => Some(Err(e)),
            }
        } else {
            let branch_node = self.segment.branch_node(offset);
            self.descent_stack.push(branch_node.right);
            self.descent_stack.push(branch_node.left);
            self.next()
        }
    }
}
