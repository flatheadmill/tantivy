//! Bounding areas and volumes envelopes for spatial indexing.
//!
//! Envelopes store the axis-aligned bounding area or volume of a document's geometry. The BVH tree
//! indexes envelopes for fast filtering; the scorer retrieves stored geometry and r uns polygon
//! predicates for verification. One envelope per document. Predictable memory for use
//! with bitsuffle+zstd compression.

use std::io;
use std::marker::PhantomData;

use crate::spatial::bvh::SpreadSurvey;
use crate::spatial::delta::{compress, decompress, Compressible};
use crate::spatial::point::GeoPoint;
use crate::DocId;

type XYSurfaceI32 = Plane<QuantizedPoint, XYBounds<QuantizedPoint>>;
/// HUSH
pub type SpatialI32 = SpatialDefinition<
    i32,
    QuantizedPoint,
    XYBounds<QuantizedPoint>,
    XYSurfaceI32,
    DeltaCompression,
>;

/// HUSH
pub trait Spatial {
/// HUSH
    type Coord: Coordinate;
/// HUSH
    type Point: Point<Coord = Self::Coord>;
/// HUSH
    type Bounds: Bounds<Coord = Self::Coord, Point = Self::Point>;
/// HUSH
    type Surface: Surface<Coord = Self::Coord, Point = Self::Point, Bounds = Self::Bounds>;
/// HUSH
    type Compression: LeafCompression<Self::Bounds>;
/// HUSH
    const COORDINATES: usize;
/// HUSH
    const DIMENSIONS: usize;
}

/// HUSH
pub struct SpatialDefinition<C, P, B, S, LC>
where
    C: Coordinate,
    P: Point<Coord = C>,
    B: Bounds<Coord = C, Point = P>,
    S: Surface<Coord = C, Bounds = B, Point = P>,
    LC: LeafCompression<B>,
{
    _marker: PhantomData<(C, P, B, S, LC)>,
}

impl<C, P, B, S, LC> Spatial for SpatialDefinition<C, P, B, S, LC>
where
    C: Coordinate,
    P: Point<Coord = C>,
    B: Bounds<Coord = C, Point = P>,
    S: Surface<Coord = C, Bounds = B, Point = P>,
    LC: LeafCompression<B>,
{
    type Coord = C;
    type Bounds = B;
    type Point = P;
    type Surface = S;
    type Compression = LC;
    const COORDINATES: usize = Self::Bounds::COORDINATES;
    const DIMENSIONS: usize = Self::COORDINATES / 2;
}

/// HUSH
pub trait Coordinate: Copy + PartialOrd + PartialEq {
    /// HUSH
    type LeBytes: AsRef<[u8]>;
    /// HUSH
    const BYTE_SIZE: usize;
    /// HUSH
    const MIN: Self;
    /// HUSH
    const MAX: Self;
    /// HUSH
    fn minimum(self, other: Self) -> Self;
    /// HUSH
    fn maximum(self, other: Self) -> Self;
    /// HUSH
    fn subtract(self, other: Self) -> Self;
    /// HUSH
    fn compare(&self, other: &Self) -> std::cmp::Ordering;
    /// HUSH
    fn to_le_bytes(self) -> Self::LeBytes;
    /// HUSH
    fn from_le_bytes(bytes: &[u8]) -> Self;
}

/// HUSH
pub trait Point: Copy {
/// HUSH
    type Coord: Coordinate;
/// HUSH
    fn from_geo(geo: GeoPoint) -> Self;
/// HUSH
    fn x(&self) -> Self::Coord;
/// HUSH
    fn y(&self) -> Self::Coord;
}

/// HUSH
#[derive(Copy, Clone)]
pub struct QuantizedPoint(pub i32, pub i32);

impl Point for QuantizedPoint {
    type Coord = i32;
    fn from_geo(geo: GeoPoint) -> Self {
        QuantizedPoint(
            (geo.lon / (360.0 / (1i64 << 32) as f64)).floor() as i32,
            (geo.lat / (180.0 / (1i64 << 32) as f64)).floor() as i32,
        )
    }
    fn x(&self) -> Self::Coord {
        self.0
    }
    fn y(&self) -> Self::Coord {
        self.1
    }
}

/// HUSH
pub trait Bounds: Copy {
/// HUSH
    type Coord: Coordinate;
/// HUSH
    type Point: Point<Coord = Self::Coord>;
/// HUSH
    const COORDINATES: usize;
/// HUSH
    fn empty() -> Self;
/// HUSH
    fn get(&self, index: usize) -> Self::Coord;
/// HUSH
    fn set(&mut self, index: usize, value: Self::Coord);
/// HUSH
    fn extend_by_point(&mut self, point: Self::Point);
/// HUSH
    fn extend_by_line(&mut self, from: Self::Point, to: Self::Point);
}

/// HUSH
#[derive(Copy, Clone)]
pub struct XYBounds<P: Point> {
    bounds: [P::Coord; 4],
}

impl<P: Point> Bounds for XYBounds<P> {
    type Coord = P::Coord;
    type Point = P;
    const COORDINATES: usize = 4;
    fn empty() -> Self {
        XYBounds {
            bounds: [Self::Coord::MIN; 4],
        }
    }
    fn get(&self, index: usize) -> Self::Coord {
        self.bounds[index]
    }
    fn set(&mut self, index: usize, value: Self::Coord) {
        self.bounds[index] = value
    }
    fn extend_by_point(&mut self, point: Self::Point) {
        self.bounds[0] = self.bounds[0].minimum(point.y()); // min_y
        self.bounds[1] = self.bounds[1].minimum(point.x()); // min_x
        self.bounds[2] = self.bounds[2].maximum(point.y()); // max_y
        self.bounds[3] = self.bounds[3].maximum(point.x()); // max_x
    }
    fn extend_by_line(&mut self, from: Self::Point, to: Self::Point) {
        // No arc bulge in 2D - just extend for both endpoints
        self.extend_by_point(from);
        self.extend_by_point(to);
    }
}

impl Coordinate for i32 {
    type LeBytes = [u8; 4];
    const BYTE_SIZE: usize = 4;
    const MIN: Self = i32::MIN;
    const MAX: Self = i32::MAX;
    fn minimum(self, other: Self) -> Self {
        std::cmp::Ord::min(self, other)
    }
    fn maximum(self, other: Self) -> Self {
        std::cmp::Ord::max(self, other)
    }
    fn subtract(self, other: Self) -> Self {
        self - other
    }
    fn compare(&self, other: &Self) -> std::cmp::Ordering {
        std::cmp::Ord::cmp(self, other)
    }
    fn to_le_bytes(self) -> Self::LeBytes {
        i32::to_le_bytes(self)
    }
    fn from_le_bytes(data: &[u8]) -> Self {
        i32::from_le_bytes(data.try_into().unwrap())
    }
}

impl Coordinate for f64 {
    type LeBytes = [u8; 8];
    const BYTE_SIZE: usize = 8;
    const MIN: Self = f64::NEG_INFINITY;
    const MAX: Self = f64::INFINITY;
    fn minimum(self, other: Self) -> Self {
        self.min(other)
    }
    fn maximum(self, other: Self) -> Self {
        self.max(other)
    }
    fn subtract(self, other: Self) -> Self {
        self - other
    }
    fn compare(&self, other: &Self) -> std::cmp::Ordering {
        self.total_cmp(other)
    }
    fn to_le_bytes(self) -> Self::LeBytes {
        f64::to_le_bytes(self)
    }
    fn from_le_bytes(data: &[u8]) -> Self {
        f64::from_le_bytes(data.try_into().unwrap())
    }
}

/// Marker for delta compression (i32).
pub struct DeltaCompression;

/// Compression strategy for leaf pages.
pub trait LeafCompression<B: Bounds> {
/// HUSH
    const PAGE_SIZE: usize;
/// HUSH
    fn compress<W: io::Write>(envelopes: &mut [Envelope<B>], write: &mut W) -> io::Result<()>;
/// HUSH
    fn decompress(data: &[u8]) -> io::Result<Vec<Envelope<B>>>;
}

/// HUSH
pub trait Surface {
/// HUSH
    type Coord: Coordinate;
/// HUSH
    type Point: Point<Coord = Self::Coord>;
/// HUSH
    type Bounds: Bounds<Coord = Self::Coord>;
/// HUSH
    const DIMENSIONS: usize;
/// HUSH
    fn create_bounds(point: Self::Point) -> Vec<Self::Coord>;
}

/// HUSH
pub struct Plane<P: Point, B: Bounds>(PhantomData<(P, B)>);

impl<P: Point, B: Bounds<Coord = P::Coord>> Surface for Plane<P, B> {
    type Coord = P::Coord;
    type Point = P;
    type Bounds = B;
    const DIMENSIONS: usize = 2;
    fn create_bounds(point: Self::Point) -> Vec<Self::Coord> {
        // [min_y, min_x, max_y, max_x]
        vec![point.y(), point.x(), point.y(), point.x()]
    }
}

/// Envelope for a 2D bounding area with integer coordinates.
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct Envelope<B: Bounds> {
    /// Bounding area as min x, min y, max x, max y.
    pub bounds: B,
    /// The document id associated with this bounding area.
    pub doc_id: DocId,
}

impl<B: Bounds> Envelope<B> {
/// HUSH
    pub fn from_bounds(doc_id: DocId, bounds: B) -> Self {
        Envelope { doc_id, bounds }
    }

/// HUSH
    pub fn skeleton(doc_id: u32) -> Self {
        Envelope {
            bounds: B::empty(),
            doc_id,
        }
    }
}

struct CompressibleCoord<'a, B: Bounds<Coord = i32>> {
    envelopes: &'a [Envelope<B>],
    dimension: usize,
}

impl<'a, B: Bounds<Coord = i32>> CompressibleCoord<'a, B> {
    fn new(envelopes: &'a [Envelope<B>], dimension: usize) -> Self {
        CompressibleCoord {
            envelopes,
            dimension,
        }
    }
}

impl<'a, B: Bounds<Coord = i32>> Compressible for CompressibleCoord<'a, B> {
    type Value = i32;
    fn len(&self) -> usize {
        self.envelopes.len()
    }
    fn get(&self, i: usize) -> i32 {
        self.envelopes[i].bounds.get(self.dimension)
    }
}

struct CompressibleDocId<'a, B: Bounds> {
    envelopes: &'a [Envelope<B>],
}

impl<'a, B: Bounds> Compressible for CompressibleDocId<'a, B> {
    type Value = u32;
    fn len(&self) -> usize {
        self.envelopes.len()
    }
    fn get(&self, i: usize) -> u32 {
        self.envelopes[i].doc_id
    }
}

/// Compression strategy for leaf pages.
impl<B: Bounds<Coord = i32>> LeafCompression<B> for DeltaCompression {
    const PAGE_SIZE: usize = 512;
    fn compress<W: io::Write>(envelopes: &mut [Envelope<B>], write: &mut W) -> io::Result<()> {
        let mut spreads: Vec<SpreadSurvey<i32>> = (0..B::COORDINATES)
            .map(|_| SpreadSurvey::default())
            .collect();
        for envelope in envelopes.iter() {
            for i in 0..B::COORDINATES {
                spreads[i].survey(envelope.bounds.get(i));
            }
        }
        let (dimension, _) = spreads
            .iter()
            .enumerate()
            .max_by_key(|(_, s)| s.spread())
            .unwrap();
        write.write_all(&(envelopes.len() as u16).to_le_bytes())?;
        envelopes.sort_by(|a, b| a.bounds.get(dimension).compare(&b.bounds.get(dimension)));
        compress(&CompressibleDocId { envelopes }, write)?;
        for i in 0..B::COORDINATES {
            compress(&CompressibleCoord::new(envelopes, i), write)?;
        }
        Ok(())
    }
    fn decompress(data: &[u8]) -> io::Result<Vec<Envelope<B>>> {
        use common::BinarySerializable;
        let mut cursor = data;
        let count: usize = u16::deserialize(&mut cursor)? as usize;
        let mut offset: usize = 0;
        let mut envelopes: Vec<Envelope<B>> = Vec::with_capacity(count);
        offset += decompress::<u32, _>(&cursor[offset..], count, |_, doc_id| {
            envelopes.push(Envelope::skeleton(doc_id))
        })?;
        for i in 0..B::COORDINATES {
            offset += decompress::<i32, _>(&cursor[offset..], count, |j, word| {
                envelopes[j].bounds.set(i, word);
            })?;
        }
        Ok(envelopes)
    }
}
