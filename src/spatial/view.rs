//! Iterative access to geometry from multiple sources.
//!
//! This module provides a uniform interface for reading geometry whether the source is document
//! storage, index leaf pages, or streaming conversion from external types like GeoR ust polygons
//! or WKT. The PointListServer trait abstracts coordinate access, and View types pr ovide
//! iterator-based traversal of rings, linestrings, and geometry collections.
//!
//! The same View types work for deserializing user output, running predicates durin g query, and
//! converting between geometry representations. Implementations of PointListServer adapt to the
//! source be it raw byte slices, XOR-compressed streams, or foreign geometry types, while views
//! remain oblivious to storage format.

use std::cell::RefCell;
use std::collections::HashMap;
use std::io;
use std::rc::Rc;

// const POINT: u32 = 0;
// const MULTIPOINT: u32 = 1;
// const LINESTRING: u32 = 2;
// const MULTILINESTRING: u32 = 3;
const POLYGON: u32 = 4;
const MULTIPOLYGON: u32 = 5;
const GEOMETRYCOLLECTION: u32 = 6;

/// Indicates the type of geometry associated with a bounding area/volume and can potentially
/// create a tree of types for multi geometries and geometry collections.
#[derive(Clone)]
pub enum GeometrySnippet {
    /// Point simply uses one of values in the bounding area/volume.
    Point,
    /// An array of points.
    MultiPoint(u32),
    /// A single index into the point arrays array.
    LineString(u32),
    /// First index into the point arrays array is the exterior, subsequent indexes are the holes.
    Polygon(Vec<u32>),
    /// An array of line strings.
    MultiLineString(Vec<GeometrySnippet>),
    /// An array of polygons.
    MultiPolygon(Vec<GeometrySnippet>),
    /// An array of heterogenous geometries.
    GeometryCollection(Vec<GeometrySnippet>),
}

/// Parses a single geometry from the descriptor stream into a [`GeometrySnippet`].
///
/// The `kind` parameter specifies the geometry type—either read from the stream by the caller
/// or passed implicitly for Multi* children. The function consumes descriptors from `iter` and
/// appends point array lengths to `lengths`, returning a snippet with indices into that array.
///
/// Will recurse as necessary for multi-points, multi-linestrings, multi-polygons and geometry
/// collections.
fn descriptor_to_geometry_snippet<'a>(
    kind: u32,
    iter: &mut impl Iterator<Item = &'a u32>,
    lengths: &mut Vec<u32>,
) -> io::Result<GeometrySnippet> {
    match kind {
        POLYGON => {
            let ring_count = *iter
                .next()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "expected ring count"))?;
            let mut indices = Vec::with_capacity(ring_count as usize);
            for _ in 0..ring_count {
                let length = *iter.next().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "expected ring length
",
                    )
                })?;
                indices.push(lengths.len() as u32);
                lengths.push(length);
            }
            Ok(GeometrySnippet::Polygon(indices))
        }
        MULTIPOLYGON => {
            let count = *iter.next().ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "expected polygon count")
            })?;
            let mut children = Vec::with_capacity(count as usize);
            for _ in 0..count {
                children.push(descriptor_to_geometry_snippet(POLYGON, iter, lengths)?);
            }
            Ok(GeometrySnippet::MultiPolygon(children))
        }
        GEOMETRYCOLLECTION => {
            let count = *iter.next().ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "expected geometry count")
            })?;
            let mut children = Vec::with_capacity(count as usize);
            for _ in 0..count {
                let child_kind = *iter.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "expected geometry type")
                })?;
                children.push(descriptor_to_geometry_snippet(child_kind, iter, lengths)?);
            }
            Ok(GeometrySnippet::GeometryCollection(children))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported geometry type",
        )),
    }
}

/// Converts a flat descriptor stream into an indexed length tables while converting the
/// descriptors to reference the indexed table instead of recording the length.
///
/// Descriptors arrive as a flat sequence, below is polygon wiht a hold followed by a line string:
///
/// ```text
/// [type, ring_count, len, len, type, len, ...]
/// ```
/// The separation allows random access to any geometry's rings without scanning the entire
/// descriptor stream. XOR compression adds a third `offsets` vector for byte positi ons; raw and
/// bitshuffle formats derive offsets from `lengths` directly. For both binary encodings and XOR
/// encodings we subsequently zip the lengths array into an array of tuples of length, byte offset
/// and byte length.
pub fn descriptors_to_geometry_snippets(
    descriptors: &[u32],
) -> io::Result<(Vec<GeometrySnippet>, Vec<u32>)> {
    let mut iter = descriptors.iter();
    let mut lengths = Vec::new();
    let mut snippets = Vec::new();

    while let Some(&kind) = iter.next() {
        snippets.push(descriptor_to_geometry_snippet(
            kind,
            &mut iter,
            &mut lengths,
        )?);
    }

    Ok((snippets, lengths))
}

/// Provides indexed access to point coordinate arrays.
///
/// Views use this trait to retrieve point data for rings, linestrings, and multipoi nts without
/// binding to a specific storage format. Implementations range from direct byte sli ce access in
/// GeometryBuffer to caching decoders that hold decoded rings in an LRU cache with reference
/// counting.
///
/// The trait takes `&self` rather than `&mut self` so that multiple views can share a reference
/// to the same server. Implementations requiring internal mutation for caching or decode buffers
/// use interior mutability through RefCell.
pub trait PointListServer {
    /// Returns an iterator over the coordinates in the point list at the given index.
    ///
    /// Each point list represents a sequence of coordinates for a ring, linestring, or multipoint
    /// geometry. The index corresponds to entries in the point arrays table built during
    /// descriptor parsing. Exterior rings, interior rings, and linestrings each occupy their own
    /// index.
    ///
    /// Callers should not assume the iterator is reusable or clonable. For formats requiring
    /// stateful decode like XOR compression, the iterator consumes from an internal buffer and
    /// cannot be rewound.
    ///
    /// Returns `io::Result` because compressed formats like XOR may encounter decode errors. For
    /// trusted internal data where errors indicate bugs rather than runtime conditions,
    /// implementations may panic on malformed input instead.
    fn get_point_list(&self, index: u32) -> io::Result<impl Iterator<Item = (f64, f64)>>;
}

/// HUSH
pub struct GeometryBuffer {
    /// Geometry type information.
    pub snippet: GeometrySnippet,
    /// A tuple of point count, byte offset and byte length.
    pub point_arrays: Vec<(u32, u32, u32)>,
    /// Serialized points or vectors.
    pub bytes: Vec<u8>,
}

impl PointListServer for GeometryBuffer {
    fn get_point_list(&self, index: u32) -> io::Result<impl Iterator<Item = (f64, f64)>> {
        let (_count, offset, len) = self.point_arrays[index as usize];
        let slice = &self.bytes[offset as usize..(offset + len) as usize];
        Ok(slice.chunks_exact(16).map(|chunk| {
            let x = f64::from_le_bytes(chunk[0..8].try_into().unwrap());
            let y = f64::from_le_bytes(chunk[8..16].try_into().unwrap());
            (x, y)
        }))
    }
}

/// Interface for decoding compressed ring data.
///
/// CachingPointListServer delegates to a DecodeSource on cache misses. Implementations handle
/// format-specific decoding such as XOR delta decompression or bitshuffle reconstruction. The
/// returned Vec is wrapped in Rc and cached for subsequent access to the same ring.
pub trait DecodeSource {
    /// Decodes the ring at the given index into coordinate pairs.
    ///
    /// Called by CachingPointListServer when the requested ring is not in cache. The
    /// implementation reads from its underlying storage, decompresses as needed, and returns
    /// owned coordinate data. Errors indicate malformed or truncated data in the source.
    fn decode_ring(&self, index: u32) -> io::Result<Vec<[f64; 2]>>;
}

/// Point list server with caching for decoded rings.
///
/// Wraps a DecodeSource and caches decoded coordinate arrays in memory. Uses interi or mutability
/// through RefCell so that multiple views can share a reference while the cache updates
/// transparently. Decoded rings are reference-counted, allowing multiple outstanding iterators
/// without lifetime conflicts.
pub struct CachingPointListServer<S> {
    cache: RefCell<HashMap<u32, Rc<Vec<[f64; 2]>>>>,
    source: S,
}

impl<S: DecodeSource> CachingPointListServer<S> {
    /// Returns the coordinates for the ring at the given index.
    ///
    /// On cache hit, clones the Rc handle and returns immediately. On cache miss, decodes
    /// through the underlying source, caches the result, and returns an Rc to the new entry.
    /// RefCell borrows are kept brief to avoid runtime panics from overlapping mutable access.
    pub fn get_point_list(&self, index: u32) -> io::Result<Rc<Vec<[f64; 2]>>> {
        // Brief borrow - just lookup
        if let Some(data) = self.cache.borrow().get(&index) {
            return Ok(Rc::clone(data));
        }

        // Decode outside of borrow
        let decoded = self.source.decode_ring(index)?;
        let rc = Rc::new(decoded);

        // Brief borrow - just insert
        self.cache.borrow_mut().insert(index, Rc::clone(&rc));

        Ok(rc)
    }
}

/// HUSH
pub enum GeometryView<'a, P: PointListServer> {
    /// HUSH
    Point, // coords from envelope, not from PointListServer
    /// HUSH
    LineString(LineStringView<'a, P>),
    /// HUSH
    Polygon(PolygonView<'a, P>),
    /// HUSH
    MultiPoint(MultiPointView<'a, P>),
    /// HUSH
    MultiLineString(MultiLineStringView<'a, P>),
    /// HUSH
    MultiPolygon(MultiPolygonView<'a, P>),
    /// HUSH
    GeometryCollection(GeometryCollectionView<'a, P>),
}

/// HUSH
pub struct MultiPointView<'a, P: PointListServer> {
    index: u32,
    point_server: &'a P,
}

/// HUSH
impl<'a, P: PointListServer> MultiPointView<'a, P> {
    /// HUSH
    pub fn points(&mut self) -> io::Result<impl Iterator<Item = (f64, f64)> + use<'_, P>> {
        self.point_server.get_point_list(self.index)
    }
}

/// HUSH
pub struct LineStringView<'a, P: PointListServer> {
    index: u32,
    point_server: &'a P,
}

/// HUSH
impl<'a, P: PointListServer> LineStringView<'a, P> {
    /// HUSH
    pub fn line_string(&mut self) -> io::Result<impl Iterator<Item = (f64, f64)> + use<'_, P>> {
        self.point_server.get_point_list(self.index)
    }
}

/// HUSH
pub struct MultiLineStringView<'a, P: PointListServer> {
    indices: Vec<u32>,
    point_server: &'a P,
}

/// HUSH
impl<'a, P: PointListServer> MultiLineStringView<'a, P> {
    /// HUSH
    pub fn line_strings(&mut self) -> impl Iterator<Item = LineStringView<'a, P>> + '_ {
        self.indices.iter().map(|&index| LineStringView::<'a, P> {
            index,
            point_server: self.point_server,
        })
    }
}

/// HUSH
pub struct RingView<'a, P: PointListServer> {
    index: u32,
    point_server: &'a P,
}

/// HUSH
impl<'a, P: PointListServer> RingView<'a, P> {
    /// HUSH
    pub fn line_string(&mut self) -> io::Result<impl Iterator<Item = (f64, f64)> + use<'_, P>> {
        self.point_server.get_point_list(self.index)
    }
}

/// HUSH
pub struct PolygonView<'a, P: PointListServer> {
    rings: Vec<u32>,
    point_server: &'a P,
}

impl<'a, P: PointListServer> PolygonView<'a, P> {
    /// HUSH
    pub fn exterior(&mut self) -> io::Result<impl Iterator<Item = (f64, f64)> + use<'_, P>> {
        self.point_server.get_point_list(0)
    }
    /// HUSH
    pub fn interor(&mut self) -> impl Iterator<Item = RingView<'a, P>> + '_ {
        self.rings[1..].iter().map(|&index| RingView::<'a, P> {
            index,
            point_server: self.point_server,
        })
    }
}

/// HUSH
pub struct MultiPolygonView<'a, P: PointListServer> {
    polygons: Vec<Vec<u32>>,
    point_server: &'a P,
}

impl<'a, P: PointListServer> MultiPolygonView<'a, P> {
    /// HUSH
    pub fn line_strings(&mut self) -> impl Iterator<Item = PolygonView<'a, P>> + '_ {
        self.polygons.iter().map(|rings| PolygonView::<'a, P> {
            rings: rings.clone(),
            point_server: self.point_server,
        })
    }
}

/// HUSH
pub struct GeometryCollectionView<'a, P: PointListServer> {
    children: Vec<GeometrySnippet>,
    point_server: &'a P,
}

impl<'a, P: PointListServer> GeometryCollectionView<'a, P> {
    /// HUSH
    pub fn geometries(&self) -> impl Iterator<Item = GeometryView<'a, P>> + '_ {
        let point_server = self.point_server;
        self.children
            .iter()
            .map(move |snippet| snippet_to_view(snippet, point_server))
    }
}

fn snippet_to_view<'a, P: PointListServer>(
    snippet: &GeometrySnippet,
    point_server: &'a P,
) -> GeometryView<'a, P> {
    match snippet {
        GeometrySnippet::Point => GeometryView::Point,
        GeometrySnippet::LineString(index) => GeometryView::LineString(LineStringView {
            index: *index,
            point_server,
        }),
        GeometrySnippet::Polygon(rings) => GeometryView::Polygon(PolygonView {
            rings: rings.clone(),
            point_server,
        }),
        GeometrySnippet::MultiPoint(index) => {
            // Extract indices from children
            GeometryView::MultiPoint(MultiPointView {
                index: *index,
                point_server,
            })
        }
        GeometrySnippet::MultiLineString(children) => {
            let indices = children
                .iter()
                .filter_map(|s| {
                    if let GeometrySnippet::LineString(i) = s {
                        Some(*i)
                    } else {
                        None
                    }
                })
                .collect();
            GeometryView::MultiLineString(MultiLineStringView {
                indices,
                point_server,
            })
        }
        GeometrySnippet::MultiPolygon(children) => {
            let polygons = children
                .iter()
                .filter_map(|s| {
                    if let GeometrySnippet::Polygon(rings) = s {
                        Some(rings.clone())
                    } else {
                        None
                    }
                })
                .collect();
            GeometryView::MultiPolygon(MultiPolygonView {
                polygons,
                point_server,
            })
        }
        GeometrySnippet::GeometryCollection(children) => {
            GeometryView::GeometryCollection(GeometryCollectionView {
                children: children.clone(),
                point_server,
            })
        }
    }
}
