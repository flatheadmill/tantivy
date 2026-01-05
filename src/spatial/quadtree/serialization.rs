use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};

// =============================================================================
// Varint Encoding Utilities
// =============================================================================

/// Maximum bytes needed to encode a u64 as varint.
pub const MAX_VARINT_LEN: usize = 10;

/// Encodes a u64 as a variable-length integer.
///
/// Uses LEB128 encoding: each byte contains 7 bits of data and 1 continuation bit.
/// Values 0-127 use 1 byte, 128-16383 use 2 bytes, etc.
///
/// # Arguments
///
/// * `value` - The value to encode
/// * `buf` - Buffer to write to (must have at least MAX_VARINT_LEN bytes)
///
/// # Returns
///
/// The number of bytes written.
#[inline]
pub fn encode_varint(mut value: u64, buf: &mut [u8]) -> usize {
    let mut i = 0;
    while value >= 0x80 {
        buf[i] = (value as u8) | 0x80;
        value >>= 7;
        i += 1;
    }
    buf[i] = value as u8;
    i + 1
}

/// Decodes a variable-length integer from bytes.
///
/// # Arguments
///
/// * `buf` - Buffer to read from
///
/// # Returns
///
/// A tuple of (decoded value, bytes consumed), or None if buffer is too short.
#[inline]
pub fn decode_varint(buf: &[u8]) -> Option<(u64, usize)> {
    let mut value: u64 = 0;
    let mut shift = 0;

    for (i, &byte) in buf.iter().enumerate() {
        if i >= MAX_VARINT_LEN {
            return None; // Overflow protection
        }

        value |= ((byte & 0x7F) as u64) << shift;

        if byte & 0x80 == 0 {
            return Some((value, i + 1));
        }

        shift += 7;
    }

    None // Buffer too short
}

/// Writes a varint to a writer.
#[inline]
pub fn write_varint<W: Write>(writer: &mut W, value: u64) -> io::Result<usize> {
    let mut buf = [0u8; MAX_VARINT_LEN];
    let len = encode_varint(value, &mut buf);
    writer.write_all(&buf[..len])?;
    Ok(len)
}

/// Reads a varint from a reader.
#[inline]
pub fn read_varint<R: Read>(reader: &mut R) -> io::Result<u64> {
    let mut value: u64 = 0;
    let mut shift = 0;
    let mut buf = [0u8; 1];

    loop {
        reader.read_exact(&mut buf)?;
        let byte = buf[0];

        value |= ((byte & 0x7F) as u64) << shift;

        if byte & 0x80 == 0 {
            return Ok(value);
        }

        shift += 7;
        if shift >= 64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "varint too long",
            ));
        }
    }
}

// =============================================================================
// Header Structure
// =============================================================================

/// Magic number for quadtree index files: "QTRE" in little-endian.
pub const QUADTREE_MAGIC: u32 = 0x45525451; // "QTRE" as little-endian

/// Current format version.
pub const FORMAT_VERSION: u16 = 1;

/// Header at the start of a quadtree index file.
///
/// The header contains metadata needed to interpret the file contents.
/// It is written first and contains the global bounds required for
/// cell ID interpretation.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadtreeHeader {
    /// Magic number for file type identification.
    pub magic: u32,

    /// Format version.
    pub version: u16,

    /// Flags (reserved for future use).
    pub flags: u16,

    /// Global bounds - minimum X coordinate.
    pub bounds_min_x: f64,

    /// Global bounds - minimum Y coordinate.
    pub bounds_min_y: f64,

    /// Global bounds - maximum X coordinate.
    pub bounds_max_x: f64,

    /// Global bounds - maximum Y coordinate.
    pub bounds_max_y: f64,
}

impl QuadtreeHeader {
    /// Header size in bytes.
    pub const SIZE: usize = 40;

    /// Creates a new header with the given bounds.
    pub fn new(bounds: &Bounds) -> Self {
        Self {
            magic: QUADTREE_MAGIC,
            version: FORMAT_VERSION,
            flags: 0,
            bounds_min_x: bounds.min_x,
            bounds_min_y: bounds.min_y,
            bounds_max_x: bounds.max_x,
            bounds_max_y: bounds.max_y,
        }
    }

    /// Returns the bounds from this header.
    pub fn bounds(&self) -> Bounds {
        Bounds::new(
            self.bounds_min_x,
            self.bounds_min_y,
            self.bounds_max_x,
            self.bounds_max_y,
        )
    }

    /// Writes the header to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.magic.to_le_bytes())?;
        writer.write_all(&self.version.to_le_bytes())?;
        writer.write_all(&self.flags.to_le_bytes())?;
        writer.write_all(&self.bounds_min_x.to_le_bytes())?;
        writer.write_all(&self.bounds_min_y.to_le_bytes())?;
        writer.write_all(&self.bounds_max_x.to_le_bytes())?;
        writer.write_all(&self.bounds_max_y.to_le_bytes())?;
        Ok(())
    }

    /// Reads a header from a reader.
    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut buf4 = [0u8; 4];
        let mut buf2 = [0u8; 2];
        let mut buf8 = [0u8; 8];

        reader.read_exact(&mut buf4)?;
        let magic = u32::from_le_bytes(buf4);

        if magic != QUADTREE_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid magic number: expected 0x{:08X}, got 0x{:08X}",
                    QUADTREE_MAGIC, magic
                ),
            ));
        }

        reader.read_exact(&mut buf2)?;
        let version = u16::from_le_bytes(buf2);

        reader.read_exact(&mut buf2)?;
        let flags = u16::from_le_bytes(buf2);

        reader.read_exact(&mut buf8)?;
        let bounds_min_x = f64::from_le_bytes(buf8);

        reader.read_exact(&mut buf8)?;
        let bounds_min_y = f64::from_le_bytes(buf8);

        reader.read_exact(&mut buf8)?;
        let bounds_max_x = f64::from_le_bytes(buf8);

        reader.read_exact(&mut buf8)?;
        let bounds_max_y = f64::from_le_bytes(buf8);

        Ok(Self {
            magic,
            version,
            flags,
            bounds_min_x,
            bounds_min_y,
            bounds_max_x,
            bounds_max_y,
        })
    }
}

// =============================================================================
// Footer Structure
// =============================================================================

/// Footer at the end of a quadtree index file.
///
/// The footer contains statistics and offsets for efficient file access.
/// It is written last after all cells are written.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadtreeFooter {
    /// Number of cells in the index.
    pub num_cells: u32,

    /// Total number of shapes across all cells.
    pub num_shapes: u32,

    /// Total number of edges across all cells.
    pub num_edges: u32,

    /// Number of collapse candidates.
    pub num_collapse_candidates: u32,

    /// Offset to the collapse candidates array (from file start).
    pub collapse_candidates_offset: u64,

    /// Offset to the cell data (from file start).
    pub cell_data_offset: u64,

    /// Size of cell data in bytes.
    pub cell_data_size: u64,

    /// CRC32 checksum of cell data (optional, 0 if not computed).
    pub checksum: u32,

    /// Reserved for future use.
    pub reserved: u32,
}

impl QuadtreeFooter {
    /// Footer size in bytes.
    pub const SIZE: usize = 48;

    /// Creates an empty footer.
    pub fn new() -> Self {
        Self {
            num_cells: 0,
            num_shapes: 0,
            num_edges: 0,
            num_collapse_candidates: 0,
            collapse_candidates_offset: 0,
            cell_data_offset: QuadtreeHeader::SIZE as u64,
            cell_data_size: 0,
            checksum: 0,
            reserved: 0,
        }
    }

    /// Writes the footer to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.num_cells.to_le_bytes())?;
        writer.write_all(&self.num_shapes.to_le_bytes())?;
        writer.write_all(&self.num_edges.to_le_bytes())?;
        writer.write_all(&self.num_collapse_candidates.to_le_bytes())?;
        writer.write_all(&self.collapse_candidates_offset.to_le_bytes())?;
        writer.write_all(&self.cell_data_offset.to_le_bytes())?;
        writer.write_all(&self.cell_data_size.to_le_bytes())?;
        writer.write_all(&self.checksum.to_le_bytes())?;
        writer.write_all(&self.reserved.to_le_bytes())?;
        Ok(())
    }

    /// Reads a footer from a reader positioned at the footer start.
    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut buf4 = [0u8; 4];
        let mut buf8 = [0u8; 8];

        reader.read_exact(&mut buf4)?;
        let num_cells = u32::from_le_bytes(buf4);

        reader.read_exact(&mut buf4)?;
        let num_shapes = u32::from_le_bytes(buf4);

        reader.read_exact(&mut buf4)?;
        let num_edges = u32::from_le_bytes(buf4);

        reader.read_exact(&mut buf4)?;
        let num_collapse_candidates = u32::from_le_bytes(buf4);

        reader.read_exact(&mut buf8)?;
        let collapse_candidates_offset = u64::from_le_bytes(buf8);

        reader.read_exact(&mut buf8)?;
        let cell_data_offset = u64::from_le_bytes(buf8);

        reader.read_exact(&mut buf8)?;
        let cell_data_size = u64::from_le_bytes(buf8);

        reader.read_exact(&mut buf4)?;
        let checksum = u32::from_le_bytes(buf4);

        reader.read_exact(&mut buf4)?;
        let reserved = u32::from_le_bytes(buf4);

        Ok(Self {
            num_cells,
            num_shapes,
            num_edges,
            num_collapse_candidates,
            collapse_candidates_offset,
            cell_data_offset,
            cell_data_size,
            checksum,
            reserved,
        })
    }
}

impl Default for QuadtreeFooter {
    fn default() -> Self {
        Self::new()
    }
}

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

use std::collections::HashMap;

use crate::spatial::quadtree::{Bounds, ClippedShape, QuadtreeCell, QuadtreeCellId};

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

// =============================================================================
// CellWriter - Streaming Cell Output
// =============================================================================

/// Writes quadtree cells to a streaming output with delta encoding.
///
/// CellWriter handles the serialization of cells in cell_id order,
/// applying delta encoding for doc_ids and edge indices to achieve
/// compression.
///
/// # Usage
///
/// ```ignore
/// let file = File::create("index.qt")?;
/// let mut writer = CellWriter::new(file, bounds, options)?;
///
/// for cell in cells {
///     writer.write_cell(&cell)?;
/// }
///
/// let footer = writer.finish()?;
/// ```
///
/// # Format
///
/// Each cell is encoded as:
/// - cell_id: u64 (8 bytes, little-endian)
/// - num_shapes: varint
/// - For each shape:
///   - doc_id_delta: varint (delta from previous doc_id in this cell)
///   - flags: u8 (bit 0 = contains_center)
///   - num_edges: varint
///   - For each edge:
///     - edge_id_delta: varint (delta from previous edge in this shape)
pub struct CellWriter<W: Write + Seek> {
    /// Buffered writer for output.
    writer: BufWriter<W>,

    /// Collapse detector for advisory collapse candidates.
    collapse_detector: CollapseDetector,

    /// Header (written at start).
    header: QuadtreeHeader,

    /// Footer statistics (accumulated during writes).
    footer: QuadtreeFooter,

    /// Position where cell data started.
    cell_data_start: u64,

    /// Previous cell_id for validation (cells must be in order).
    prev_cell_id: Option<u64>,
}

/// Options for CellWriter.
#[derive(Debug, Clone)]
pub struct CellWriterOptions {
    /// Threshold for collapse detection (0 to disable).
    pub collapse_threshold: usize,

    /// Buffer size for writes.
    pub buffer_size: usize,
}

impl Default for CellWriterOptions {
    fn default() -> Self {
        Self {
            collapse_threshold: 25, // ~2.5x max_edges_per_cell / 4
            buffer_size: 64 * 1024, // 64KB
        }
    }
}

impl<W: Write + Seek> CellWriter<W> {
    /// Creates a new CellWriter.
    ///
    /// Writes the header immediately.
    pub fn new(writer: W, bounds: &Bounds, options: CellWriterOptions) -> io::Result<Self> {
        let mut buffered = BufWriter::with_capacity(options.buffer_size, writer);

        // Write header
        let header = QuadtreeHeader::new(bounds);
        header.write_to(&mut buffered)?;

        let cell_data_start = buffered.stream_position()?;

        Ok(Self {
            writer: buffered,
            collapse_detector: CollapseDetector::new(options.collapse_threshold),
            header,
            footer: QuadtreeFooter::new(),
            cell_data_start,
            prev_cell_id: None,
        })
    }

    /// Creates a CellWriter with default options.
    pub fn with_defaults(writer: W, bounds: &Bounds) -> io::Result<Self> {
        Self::new(writer, bounds, CellWriterOptions::default())
    }

    /// Returns the bounds from the header.
    pub fn bounds(&self) -> Bounds {
        self.header.bounds()
    }

    /// Returns the number of cells written so far.
    pub fn num_cells(&self) -> u32 {
        self.footer.num_cells
    }

    /// Returns the number of shapes written so far.
    pub fn num_shapes(&self) -> u32 {
        self.footer.num_shapes
    }

    /// Returns the number of edges written so far.
    pub fn num_edges(&self) -> u32 {
        self.footer.num_edges
    }

    /// Writes a single cell.
    ///
    /// Cells must be written in cell_id order (ascending).
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The cell is out of order
    /// - An I/O error occurs
    pub fn write_cell(&mut self, cell: &QuadtreeCell) -> io::Result<()> {
        let cell_id = cell.cell_id().to_raw();

        // Validate ordering
        if let Some(prev) = self.prev_cell_id {
            if cell_id <= prev {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("cells must be in ascending order: {} <= {}", cell_id, prev),
                ));
            }
        }
        self.prev_cell_id = Some(cell_id);

        // Observe for collapse detection
        self.collapse_detector.observe_cell(cell);

        // Write cell_id
        self.writer.write_all(&cell_id.to_le_bytes())?;

        // Write num_shapes
        write_varint(&mut self.writer, cell.num_shapes() as u64)?;

        // Write shapes with delta encoding
        let mut prev_doc_id: u32 = 0;

        for shape in cell.shapes() {
            // Delta-encoded doc_id
            let doc_id = shape.doc_id();
            let delta = doc_id.saturating_sub(prev_doc_id);
            write_varint(&mut self.writer, delta as u64)?;
            prev_doc_id = doc_id;

            // Flags
            let flags: u8 = if shape.contains_center() { 1 } else { 0 };
            self.writer.write_all(&[flags])?;

            // Num edges
            let num_edges = shape.num_edges();
            write_varint(&mut self.writer, num_edges as u64)?;

            // Delta-encoded edge indices
            let mut prev_edge: u16 = 0;
            for &edge in shape.edges() {
                let edge_delta = edge.saturating_sub(prev_edge);
                write_varint(&mut self.writer, edge_delta as u64)?;
                prev_edge = edge;
            }

            self.footer.num_shapes += 1;
            self.footer.num_edges += num_edges as u32;
        }

        self.footer.num_cells += 1;

        Ok(())
    }

    /// Finishes writing and returns the footer.
    ///
    /// This writes the collapse candidates and footer, and flushes the buffer.
    pub fn finish(mut self) -> io::Result<QuadtreeFooter> {
        // Record cell data size
        let cell_data_end = self.writer.stream_position()?;
        self.footer.cell_data_offset = self.cell_data_start;
        self.footer.cell_data_size = cell_data_end - self.cell_data_start;

        // Write collapse candidates
        let candidates = self.collapse_detector.finish();
        self.footer.collapse_candidates_offset = cell_data_end;
        self.footer.num_collapse_candidates = candidates.len() as u32;

        for candidate in &candidates {
            candidate.write_to(&mut self.writer)?;
        }

        // Write footer
        self.footer.write_to(&mut self.writer)?;

        // Flush
        self.writer.flush()?;

        Ok(self.footer)
    }

    /// Aborts writing and returns the inner writer.
    ///
    /// The file will be in an incomplete state.
    pub fn abort(self) -> W {
        match self.writer.into_inner() {
            Ok(w) => w,
            Err(e) => e
                .into_inner()
                .into_inner()
                .unwrap_or_else(|_| panic!("Failed to recover inner writer")),
        }
    }
}

// =============================================================================
// CellReader - Streaming Cell Iterator
// =============================================================================

/// Reads quadtree cells from a file as a streaming iterator.
///
/// CellReader handles the deserialization of cells, decoding the
/// delta-encoded doc_ids and edge indices.
///
/// # Usage
///
/// ```ignore
/// let file = File::open("index.qt")?;
/// let reader = CellReader::new(file)?;
///
/// println!("Bounds: {:?}", reader.bounds());
/// println!("Num cells: {}", reader.num_cells());
///
/// for cell_result in reader {
///     let cell = cell_result?;
///     // process cell
/// }
/// ```
pub struct CellReader<R: Read + Seek> {
    /// Buffered reader.
    reader: BufReader<R>,

    /// Header from the file.
    header: QuadtreeHeader,

    /// Footer from the file.
    footer: QuadtreeFooter,

    /// Number of cells remaining to read.
    cells_remaining: u32,

    /// Whether an error has occurred.
    errored: bool,
}

impl<R: Read + Seek> CellReader<R> {
    /// Opens a CellReader from a seekable reader.
    ///
    /// Reads the header and footer to initialize state.
    pub fn new(mut reader: R) -> io::Result<Self> {
        // Read header from start
        reader.seek(SeekFrom::Start(0))?;
        let header = QuadtreeHeader::read_from(&mut reader)?;

        // Read footer from end
        reader.seek(SeekFrom::End(-(QuadtreeFooter::SIZE as i64)))?;
        let footer = QuadtreeFooter::read_from(&mut reader)?;

        // Seek to start of cell data
        reader.seek(SeekFrom::Start(footer.cell_data_offset))?;

        // Extract values before moving footer
        let cells_remaining = footer.num_cells;

        Ok(Self {
            reader: BufReader::new(reader),
            header,
            footer,
            cells_remaining,
            errored: false,
        })
    }

    /// Returns the header.
    pub fn header(&self) -> &QuadtreeHeader {
        &self.header
    }

    /// Returns the footer.
    pub fn footer(&self) -> &QuadtreeFooter {
        &self.footer
    }

    /// Returns the bounds from the header.
    pub fn bounds(&self) -> Bounds {
        self.header.bounds()
    }

    /// Returns the total number of cells in the file.
    pub fn num_cells(&self) -> u32 {
        self.footer.num_cells
    }

    /// Returns the number of cells remaining to read.
    pub fn cells_remaining(&self) -> u32 {
        self.cells_remaining
    }

    /// Returns the total number of shapes in the file.
    pub fn num_shapes(&self) -> u32 {
        self.footer.num_shapes
    }

    /// Returns the total number of edges in the file.
    pub fn num_edges(&self) -> u32 {
        self.footer.num_edges
    }

    /// Reads the collapse candidates.
    ///
    /// This seeks to the collapse candidates section and reads all candidates.
    pub fn read_collapse_candidates(&mut self) -> io::Result<Vec<CollapseCandidate>> {
        let num = self.footer.num_collapse_candidates as usize;
        if num == 0 {
            return Ok(Vec::new());
        }

        // Save current position
        let current_pos = self.reader.stream_position()?;

        // Seek to collapse candidates
        self.reader
            .seek(SeekFrom::Start(self.footer.collapse_candidates_offset))?;

        let mut candidates = Vec::with_capacity(num);
        for _ in 0..num {
            candidates.push(CollapseCandidate::read_from(&mut self.reader)?);
        }

        // Restore position
        self.reader.seek(SeekFrom::Start(current_pos))?;

        Ok(candidates)
    }

    /// Reads the next cell.
    fn read_cell(&mut self) -> io::Result<QuadtreeCell> {
        // Read cell_id
        let mut buf8 = [0u8; 8];
        self.reader.read_exact(&mut buf8)?;
        let cell_id = QuadtreeCellId::from_raw(u64::from_le_bytes(buf8));

        // Read num_shapes
        let num_shapes = read_varint(&mut self.reader)? as usize;

        let mut cell = QuadtreeCell::with_capacity(cell_id, num_shapes);
        let mut prev_doc_id: u32 = 0;

        for _ in 0..num_shapes {
            // Read doc_id delta
            let delta = read_varint(&mut self.reader)? as u32;
            let doc_id = prev_doc_id.wrapping_add(delta);
            prev_doc_id = doc_id;

            // Read flags
            let mut flags_buf = [0u8; 1];
            self.reader.read_exact(&mut flags_buf)?;
            let contains_center = (flags_buf[0] & 1) != 0;

            // Read num_edges
            let num_edges = read_varint(&mut self.reader)? as usize;

            let mut shape = ClippedShape::with_capacity(doc_id, contains_center, num_edges);

            // Read edge deltas
            let mut prev_edge: u16 = 0;
            for _ in 0..num_edges {
                let edge_delta = read_varint(&mut self.reader)? as u16;
                let edge = prev_edge.wrapping_add(edge_delta);
                // Note: We can't use add_edge here as it maintains sorted order,
                // but the edges are already sorted from the writer.
                // We'll add directly to maintain the order.
                shape.add_edge(edge);
                prev_edge = edge;
            }

            cell.add_shape(shape);
        }

        Ok(cell)
    }

    /// Resets the reader to the beginning of cell data.
    pub fn reset(&mut self) -> io::Result<()> {
        self.reader
            .seek(SeekFrom::Start(self.footer.cell_data_offset))?;
        self.cells_remaining = self.footer.num_cells;
        self.errored = false;
        Ok(())
    }
}

impl<R: Read + Seek> Iterator for CellReader<R> {
    type Item = io::Result<QuadtreeCell>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.errored || self.cells_remaining == 0 {
            return None;
        }

        match self.read_cell() {
            Ok(cell) => {
                self.cells_remaining -= 1;
                Some(Ok(cell))
            }
            Err(e) => {
                self.errored = true;
                Some(Err(e))
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.cells_remaining as usize;
        (remaining, Some(remaining))
    }
}

impl<R: Read + Seek> ExactSizeIterator for CellReader<R> {}

// =============================================================================
// Convenience Functions
// =============================================================================

/// Writes an iterator of cells to a writer.
///
/// This is a convenience function for writing all cells at once.
pub fn write_cells<W, I>(writer: W, bounds: &Bounds, cells: I) -> io::Result<QuadtreeFooter>
where
    W: Write + Seek,
    I: IntoIterator<Item = QuadtreeCell>,
{
    let mut cell_writer = CellWriter::with_defaults(writer, bounds)?;

    for cell in cells {
        cell_writer.write_cell(&cell)?;
    }

    cell_writer.finish()
}

/// Reads all cells from a reader into a vector.
///
/// This is a convenience function for reading all cells at once.
pub fn read_all_cells<R: Read + Seek>(reader: R) -> io::Result<(Bounds, Vec<QuadtreeCell>)> {
    let cell_reader = CellReader::new(reader)?;
    let bounds = cell_reader.bounds();
    let cells: Result<Vec<_>, _> = cell_reader.collect();
    Ok((bounds, cells?))
}

// =============================================================================
// Unit Tests for Phase 6
// =============================================================================

#[cfg(test)]
mod serialization_tests {
    use std::io::Cursor;

    use super::*;

    // -------------------------------------------------------------------------
    // Varint Tests
    // -------------------------------------------------------------------------

    mod varint_tests {
        use super::*;

        #[test]
        fn test_encode_decode_small() {
            let mut buf = [0u8; MAX_VARINT_LEN];

            // Single byte values (0-127)
            for value in 0..128u64 {
                let len = encode_varint(value, &mut buf);
                assert_eq!(len, 1);
                let (decoded, consumed) = decode_varint(&buf).unwrap();
                assert_eq!(decoded, value);
                assert_eq!(consumed, 1);
            }
        }

        #[test]
        fn test_encode_decode_two_bytes() {
            let mut buf = [0u8; MAX_VARINT_LEN];

            // Two byte values (128-16383)
            for value in [128u64, 255, 1000, 10000, 16383] {
                let len = encode_varint(value, &mut buf);
                assert_eq!(len, 2);
                let (decoded, consumed) = decode_varint(&buf).unwrap();
                assert_eq!(decoded, value);
                assert_eq!(consumed, 2);
            }
        }

        #[test]
        fn test_encode_decode_large() {
            let mut buf = [0u8; MAX_VARINT_LEN];

            let test_values = [
                16384u64,    // 3 bytes
                1_000_000,   // 3 bytes
                268_435_456, // 4 bytes
                u32::MAX as u64,
                u64::MAX,
            ];

            for value in test_values {
                let len = encode_varint(value, &mut buf);
                assert!(len <= MAX_VARINT_LEN);
                let (decoded, consumed) = decode_varint(&buf).unwrap();
                assert_eq!(decoded, value);
                assert_eq!(consumed, len);
            }
        }

        #[test]
        fn test_write_read_varint() {
            let mut buf = Vec::new();

            let values = [0u64, 1, 127, 128, 255, 1000, u32::MAX as u64, u64::MAX];

            for &value in &values {
                write_varint(&mut buf, value).unwrap();
            }

            let mut cursor = Cursor::new(&buf);
            for &expected in &values {
                let decoded = read_varint(&mut cursor).unwrap();
                assert_eq!(decoded, expected);
            }
        }
    }

    // -------------------------------------------------------------------------
    // Header Tests
    // -------------------------------------------------------------------------

    mod header_tests {
        use super::*;

        #[test]
        fn test_header_roundtrip() {
            let bounds = Bounds::new(-180.0, -90.0, 180.0, 90.0);
            let header = QuadtreeHeader::new(&bounds);

            let mut buf = Vec::new();
            header.write_to(&mut buf).unwrap();

            assert_eq!(buf.len(), QuadtreeHeader::SIZE);

            let mut cursor = Cursor::new(&buf);
            let decoded = QuadtreeHeader::read_from(&mut cursor).unwrap();

            assert_eq!(decoded, header);
            assert_eq!(decoded.bounds(), bounds);
        }

        #[test]
        fn test_header_invalid_magic() {
            let mut buf = vec![0u8; QuadtreeHeader::SIZE];
            buf[0..4].copy_from_slice(&0xDEADBEEFu32.to_le_bytes());

            let mut cursor = Cursor::new(&buf);
            let result = QuadtreeHeader::read_from(&mut cursor);

            assert!(result.is_err());
            assert!(result.unwrap_err().to_string().contains("invalid magic"));
        }
    }

    // -------------------------------------------------------------------------
    // Footer Tests
    // -------------------------------------------------------------------------

    mod footer_tests {
        use super::*;

        #[test]
        fn test_footer_roundtrip() {
            let mut footer = QuadtreeFooter::new();
            footer.num_cells = 1234;
            footer.num_shapes = 5678;
            footer.num_edges = 9012;
            footer.num_collapse_candidates = 10;
            footer.collapse_candidates_offset = 12345;
            footer.cell_data_offset = QuadtreeHeader::SIZE as u64;
            footer.cell_data_size = 67890;
            footer.checksum = 0xABCD1234;

            let mut buf = Vec::new();
            footer.write_to(&mut buf).unwrap();

            assert_eq!(buf.len(), QuadtreeFooter::SIZE);

            let mut cursor = Cursor::new(&buf);
            let decoded = QuadtreeFooter::read_from(&mut cursor).unwrap();

            assert_eq!(decoded, footer);
        }
    }

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

    // -------------------------------------------------------------------------
    // CellWriter/CellReader Integration Tests
    // -------------------------------------------------------------------------

    mod writer_reader_tests {
        use super::*;

        fn create_test_cells(bounds: &Bounds) -> Vec<QuadtreeCell> {
            let mut cells = Vec::new();

            // Create cells at various positions
            let positions = [
                (10.0, 10.0, 3),
                (30.0, 10.0, 3),
                (50.0, 50.0, 2),
                (90.0, 90.0, 3),
            ];

            for (x, y, level) in positions {
                let cell_id = QuadtreeCellId::from_xy(x, y, level, bounds);
                let mut cell = QuadtreeCell::new(cell_id);

                // Add a shape with some edges
                let mut shape = ClippedShape::new(1, true);
                shape.add_edge(0);
                shape.add_edge(1);
                shape.add_edge(2);
                cell.add_shape(shape);

                // Add another shape
                let mut shape2 = ClippedShape::new(5, false);
                shape2.add_edge(10);
                shape2.add_edge(15);
                cell.add_shape(shape2);

                cells.push(cell);
            }

            // Sort by cell_id as required
            cells.sort_by_key(|c| c.cell_id().to_raw());

            cells
        }

        #[test]
        fn test_write_read_empty() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);

            let mut buf = Cursor::new(Vec::new());
            let writer = CellWriter::with_defaults(&mut buf, &bounds).unwrap();
            let footer = writer.finish().unwrap();

            assert_eq!(footer.num_cells, 0);
            assert_eq!(footer.num_shapes, 0);
            assert_eq!(footer.num_edges, 0);

            buf.set_position(0);
            let reader = CellReader::new(&mut buf).unwrap();

            assert_eq!(reader.num_cells(), 0);
            assert_eq!(reader.bounds(), bounds);

            let cells: Vec<_> = reader.collect::<Result<_, _>>().unwrap();
            assert!(cells.is_empty());
        }

        #[test]
        fn test_write_read_single_cell() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 3, &bounds);

            let mut cell = QuadtreeCell::new(cell_id);
            let mut shape = ClippedShape::new(42, true);
            shape.add_edge(5);
            shape.add_edge(10);
            shape.add_edge(15);
            cell.add_shape(shape);

            let mut buf = Cursor::new(Vec::new());
            let mut writer = CellWriter::with_defaults(&mut buf, &bounds).unwrap();
            writer.write_cell(&cell).unwrap();
            let footer = writer.finish().unwrap();

            assert_eq!(footer.num_cells, 1);
            assert_eq!(footer.num_shapes, 1);
            assert_eq!(footer.num_edges, 3);

            buf.set_position(0);
            let reader = CellReader::new(&mut buf).unwrap();

            let cells: Vec<_> = reader.collect::<Result<_, _>>().unwrap();
            assert_eq!(cells.len(), 1);

            let read_cell = &cells[0];
            assert_eq!(read_cell.cell_id(), cell_id);
            assert_eq!(read_cell.num_shapes(), 1);

            let read_shape = &read_cell.shapes()[0];
            assert_eq!(read_shape.doc_id(), 42);
            assert!(read_shape.contains_center());
            assert_eq!(read_shape.edges(), &[5, 10, 15]);
        }

        #[test]
        fn test_write_read_multiple_cells() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let cells = create_test_cells(&bounds);

            let mut buf = Cursor::new(Vec::new());
            let mut writer = CellWriter::with_defaults(&mut buf, &bounds).unwrap();

            for cell in &cells {
                writer.write_cell(cell).unwrap();
            }

            let footer = writer.finish().unwrap();

            assert_eq!(footer.num_cells, cells.len() as u32);

            buf.set_position(0);
            let reader = CellReader::new(&mut buf).unwrap();

            let read_cells: Vec<_> = reader.collect::<Result<_, _>>().unwrap();
            assert_eq!(read_cells.len(), cells.len());

            for (original, read) in cells.iter().zip(read_cells.iter()) {
                assert_eq!(read.cell_id(), original.cell_id());
                assert_eq!(read.num_shapes(), original.num_shapes());

                for (orig_shape, read_shape) in original.shapes().iter().zip(read.shapes().iter()) {
                    assert_eq!(read_shape.doc_id(), orig_shape.doc_id());
                    assert_eq!(read_shape.contains_center(), orig_shape.contains_center());
                    assert_eq!(read_shape.edges(), orig_shape.edges());
                }
            }
        }

        #[test]
        fn test_write_read_with_collapse_candidates() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);

            // Create four sibling cells with few edges
            let parent = QuadtreeCellId::from_xy(50.0, 50.0, 1, &bounds);
            let children = parent.children().unwrap();

            let mut cells: Vec<QuadtreeCell> = children
                .iter()
                .map(|&id| {
                    let mut cell = QuadtreeCell::new(id);
                    let shape = ClippedShape::new(1, false);
                    cell.add_shape(shape);
                    cell
                })
                .collect();

            cells.sort_by_key(|c| c.cell_id().to_raw());

            let options = CellWriterOptions {
                collapse_threshold: 10, // Low threshold
                buffer_size: 1024,
            };

            let mut buf = Cursor::new(Vec::new());
            let mut writer = CellWriter::new(&mut buf, &bounds, options).unwrap();

            for cell in &cells {
                writer.write_cell(cell).unwrap();
            }

            let footer = writer.finish().unwrap();

            assert!(footer.num_collapse_candidates > 0);

            buf.set_position(0);
            let mut reader = CellReader::new(&mut buf).unwrap();

            let candidates = reader.read_collapse_candidates().unwrap();
            assert_eq!(candidates.len(), footer.num_collapse_candidates as usize);
        }

        #[test]
        fn test_write_out_of_order_error() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);

            let cell1 = QuadtreeCell::new(QuadtreeCellId::from_xy(50.0, 50.0, 3, &bounds));
            let cell2 = QuadtreeCell::new(QuadtreeCellId::from_xy(10.0, 10.0, 3, &bounds));

            // cell2 might have a smaller cell_id than cell1

            let mut buf = Cursor::new(Vec::new());
            let mut writer = CellWriter::with_defaults(&mut buf, &bounds).unwrap();

            // Write in the correct order based on cell_id
            let (first, second) = if cell1.cell_id().to_raw() < cell2.cell_id().to_raw() {
                (&cell1, &cell2)
            } else {
                (&cell2, &cell1)
            };

            writer.write_cell(first).unwrap();
            writer.write_cell(second).unwrap();

            // Now try to write first again (out of order)
            let result = writer.write_cell(first);
            assert!(result.is_err());
        }

        #[test]
        fn test_reader_reset() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let cells = create_test_cells(&bounds);

            let mut buf = Cursor::new(Vec::new());
            let mut writer = CellWriter::with_defaults(&mut buf, &bounds).unwrap();
            for cell in &cells {
                writer.write_cell(cell).unwrap();
            }
            writer.finish().unwrap();

            buf.set_position(0);
            let mut reader = CellReader::new(&mut buf).unwrap();

            // Read all cells
            let first_read: Vec<_> = reader.by_ref().collect::<Result<_, _>>().unwrap();
            assert_eq!(first_read.len(), cells.len());

            // Reset and read again
            reader.reset().unwrap();

            let second_read: Vec<_> = reader.collect::<Result<_, _>>().unwrap();
            assert_eq!(second_read.len(), cells.len());

            for (a, b) in first_read.iter().zip(second_read.iter()) {
                assert_eq!(a.cell_id(), b.cell_id());
            }
        }

        #[test]
        fn test_convenience_functions() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let cells = create_test_cells(&bounds);

            let mut buf = Cursor::new(Vec::new());
            let footer = write_cells(&mut buf, &bounds, cells.clone()).unwrap();

            assert_eq!(footer.num_cells, cells.len() as u32);

            buf.set_position(0);
            let (read_bounds, read_cells) = read_all_cells(&mut buf).unwrap();

            assert_eq!(read_bounds, bounds);
            assert_eq!(read_cells.len(), cells.len());
        }

        #[test]
        fn test_delta_encoding_efficiency() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 3, &bounds);

            let mut cell = QuadtreeCell::new(cell_id);

            // Add shapes with consecutive doc_ids
            for doc_id in [100u32, 101, 102, 103, 104] {
                let mut shape = ClippedShape::new(doc_id, doc_id % 2 == 0);
                // Add consecutive edge indices
                for edge in 0..10u16 {
                    shape.add_edge(edge);
                }
                cell.add_shape(shape);
            }

            let mut buf = Cursor::new(Vec::new());
            let mut writer = CellWriter::with_defaults(&mut buf, &bounds).unwrap();
            writer.write_cell(&cell).unwrap();
            writer.finish().unwrap();

            // Verify file is reasonably compact (delta encoding should help)
            let file_size = buf.get_ref().len();
            // Header (48) + Footer (48) + cell_id (8) + data
            // With delta encoding, data should be much smaller than raw
            assert!(
                file_size < 200,
                "File size {} is larger than expected for delta-encoded data",
                file_size
            );

            buf.set_position(0);
            let (_, read_cells) = read_all_cells(&mut buf).unwrap();
            assert_eq!(read_cells.len(), 1);
            assert_eq!(read_cells[0].num_shapes(), 5);
        }
    }

    // -------------------------------------------------------------------------
    // Edge Cases
    // -------------------------------------------------------------------------

    mod edge_case_tests {
        use super::*;

        #[test]
        fn test_empty_shapes() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 3, &bounds);

            let mut cell = QuadtreeCell::new(cell_id);

            // Shape with no edges but contains_center = true
            let shape = ClippedShape::new(1, true);
            cell.add_shape(shape);

            let mut buf = Cursor::new(Vec::new());
            let mut writer = CellWriter::with_defaults(&mut buf, &bounds).unwrap();
            writer.write_cell(&cell).unwrap();
            let footer = writer.finish().unwrap();

            assert_eq!(footer.num_shapes, 1);
            assert_eq!(footer.num_edges, 0);

            buf.set_position(0);
            let (_, read_cells) = read_all_cells(&mut buf).unwrap();

            let read_shape = &read_cells[0].shapes()[0];
            assert_eq!(read_shape.doc_id(), 1);
            assert!(read_shape.contains_center());
            assert_eq!(read_shape.num_edges(), 0);
        }

        #[test]
        fn test_large_doc_ids() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 3, &bounds);

            let mut cell = QuadtreeCell::new(cell_id);

            // Large doc_id values
            for doc_id in [0u32, 1000, 100_000, 1_000_000, u32::MAX - 1] {
                let shape = ClippedShape::new(doc_id, false);
                cell.add_shape(shape);
            }

            let mut buf = Cursor::new(Vec::new());
            let mut writer = CellWriter::with_defaults(&mut buf, &bounds).unwrap();
            writer.write_cell(&cell).unwrap();
            writer.finish().unwrap();

            buf.set_position(0);
            let (_, read_cells) = read_all_cells(&mut buf).unwrap();

            let doc_ids: Vec<u32> = read_cells[0].shapes().iter().map(|s| s.doc_id()).collect();
            assert_eq!(doc_ids, vec![0, 1000, 100_000, 1_000_000, u32::MAX - 1]);
        }

        #[test]
        fn test_many_edges() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
            let cell_id = QuadtreeCellId::from_xy(50.0, 50.0, 3, &bounds);

            let mut cell = QuadtreeCell::new(cell_id);
            let mut shape = ClippedShape::new(1, false);

            // Add many edges
            for i in 0..1000u16 {
                shape.add_edge(i);
            }

            cell.add_shape(shape);

            let mut buf = Cursor::new(Vec::new());
            let mut writer = CellWriter::with_defaults(&mut buf, &bounds).unwrap();
            writer.write_cell(&cell).unwrap();
            let footer = writer.finish().unwrap();

            assert_eq!(footer.num_edges, 1000);

            buf.set_position(0);
            let (_, read_cells) = read_all_cells(&mut buf).unwrap();

            let read_shape = &read_cells[0].shapes()[0];
            assert_eq!(read_shape.num_edges(), 1000);

            for (i, &edge) in read_shape.edges().iter().enumerate() {
                assert_eq!(edge, i as u16);
            }
        }
    }
}
