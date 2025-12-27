//! Bounding areas and volumes envelopes for spatial indexing.
//!
//! Envelopes store the axis-aligned bounding area or volume of a document's geometry. The BVH tree
//! indexes envelopes for fast filtering; the scorer retrieves stored geometry and r uns polygon
//! predicates for verification. One envelope per document. Predictable memory for use
//! with bitsuffle+zstd compression.

use std::io;
use std::marker::PhantomData;

#[cfg(feature = "zstd-compression")]
use zstd;

#[cfg(feature = "zstd-compression")]
use crate::spatial::bitshuffle::{bitshuffle, bitunshuffle};
use crate::spatial::bvh::SpreadSurvey;
use crate::spatial::delta::{compress, decompress, Compressible};
use crate::spatial::point::GeoPoint;
use crate::spatial::xor::{compress_f64, decompress_f64};
use crate::DocId;

/// HUSH
pub type SpatialF64 = SpatialDefinition<XY, XYBounds<XY>, RawCompression>;

/// HUSH
pub trait Spatial {
    /// HUSH
    type Point: Point;
    /// HUSH
    type Bounds: Bounds<Point = Self::Point>;
    /// HUSH
    type Compression: LeafCompression<Self::Bounds>;
    /// HUSH
    const COORDINATES: usize;
    /// HUSH
    const DIMENSIONS: usize;
}

/// HUSH
pub struct SpatialDefinition<P, B, LC>
where
    P: Point,
    B: Bounds<Point = P>,
    LC: LeafCompression<B>,
{
    _marker: PhantomData<(P, B, LC)>,
}

impl<P, B, LC> Spatial for SpatialDefinition<P, B, LC>
where
    P: Point,
    B: Bounds<Point = P>,
    LC: LeafCompression<B>,
{
    type Bounds = B;
    type Point = P;
    type Compression = LC;
    const COORDINATES: usize = Self::Bounds::COORDINATES;
    const DIMENSIONS: usize = Self::COORDINATES / 2;
}

/// HUSH
pub trait Point: Copy {
    /// HUSH
    fn from_geo(geo: GeoPoint) -> Self;
    /// HUSH
    fn x(&self) -> f64;
    /// HUSH
    fn y(&self) -> f64;
}

/// HUSH
#[derive(Copy, Clone)]
pub struct XY(pub f64, pub f64);

impl Point for XY {
    fn from_geo(geo: GeoPoint) -> Self {
        XY(geo.lon, geo.lat)
    }
    fn x(&self) -> f64 {
        self.0
    }
    fn y(&self) -> f64 {
        self.1
    }
}

/// HUSH
pub trait Bounds: Copy {
    /// HUSH
    type Point: Point;
    /// HUSH
    const COORDINATES: usize;
    /// HUSH
    fn empty() -> Self;
    /// HUSH
    fn get(&self, index: usize) -> f64;
    /// HUSH
    fn set(&mut self, index: usize, value: f64);
    /// HUSH
    fn extend_by_point(&mut self, point: Self::Point);
    /// HUSH
    fn extend_by_line(&mut self, from: Self::Point, to: Self::Point);
}

/// HUSH
#[derive(Copy, Clone)]
pub struct XYBounds<P: Point> {
    bounds: [f64; 4],
    _marker: PhantomData<P>,
}

impl<P: Point> Bounds for XYBounds<P> {
    type Point = P;
    const COORDINATES: usize = 4;
    fn empty() -> Self {
        XYBounds {
            bounds: [f64::MAX, f64::MAX, f64::MIN, f64::MIN],
            _marker: PhantomData,
        }
    }
    fn get(&self, index: usize) -> f64 {
        self.bounds[index]
    }
    fn set(&mut self, index: usize, value: f64) {
        self.bounds[index] = value
    }
    fn extend_by_point(&mut self, point: Self::Point) {
        self.bounds[0] = self.bounds[0].min(point.y()); // min_y
        self.bounds[1] = self.bounds[1].min(point.x()); // min_x
        self.bounds[2] = self.bounds[2].max(point.y()); // max_y
        self.bounds[3] = self.bounds[3].max(point.x()); // max_x
    }
    fn extend_by_line(&mut self, from: Self::Point, to: Self::Point) {
        // No arc bulge in 2D - just extend for both endpoints
        self.extend_by_point(from);
        self.extend_by_point(to);
    }
}

/// XOR compression of stored boudning areas/volumes.
pub struct XORCompression;

/// Bitshuffle+ztd compression of stored bounding areas/volumes.
pub struct BitshuffleCompression;

/// Raw compression indicates that the bounding areas/volumes should be stored uncompressed.
pub struct RawCompression;

/// Compression strategy for leaf pages.
pub trait LeafCompression<B: Bounds> {
    /// HUSH
    const PAGE_SIZE: usize;
    /// HUSH
    fn compress<W: io::Write>(envelopes: &mut [Envelope<B>], write: &mut W) -> io::Result<()>;
    /// HUSH
    fn decompress(data: &[u8]) -> io::Result<Vec<Envelope<B>>>;
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
impl<B: Bounds> LeafCompression<B> for XORCompression {
    const PAGE_SIZE: usize = 512;
    fn compress<W: io::Write>(envelopes: &mut [Envelope<B>], write: &mut W) -> io::Result<()> {
        write.write_all(&(envelopes.len() as u16).to_le_bytes())?;

        if envelopes.is_empty() {
            return Ok(());
        }

        // Find max-spread dimension for sorting
        let mut spreads: Vec<SpreadSurvey> = (0..B::COORDINATES)
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
            .max_by(|(_, a), (_, b)| a.spread().total_cmp(&b.spread()))
            .unwrap();

        // Sort for spatial locality
        envelopes.sort_by(|a, b| a.bounds.get(dimension).total_cmp(&b.bounds.get(dimension)));

        // Doc IDs
        compress(&CompressibleDocId { envelopes }, write)?;

        // Each dimension: length-prefixed compressed coords
        for dim in 0..B::COORDINATES {
            let coords: Vec<f64> = envelopes.iter().map(|e| e.bounds.get(dim)).collect();
            let compressed = compress_f64(&coords);
            write.write_all(&(compressed.len() as u32).to_le_bytes())?;
            write.write_all(&compressed)?;
        }

        Ok(())
    }

    fn decompress(data: &[u8]) -> io::Result<Vec<Envelope<B>>> {
        let count = u16::from_le_bytes(data[0..2].try_into().unwrap()) as usize;
        if count == 0 {
            return Ok(Vec::new());
        }

        let mut envelopes: Vec<Envelope<B>> = Vec::with_capacity(count);

        let mut offset = 2;
        offset += decompress::<u32, _>(&data[offset..], count, |_, doc_id| {
            envelopes.push(Envelope::skeleton(doc_id))
        })?;
        for dim in 0..B::COORDINATES {
            let len = u32::from_le_bytes(data[offset..offset + 4].try_into().unwrap()) as usize;
            offset += 4;
            let coords = decompress_f64(&data[offset..offset + len], count);
            for (i, coord) in coords.into_iter().enumerate() {
                envelopes[i].bounds.set(dim, coord);
            }
            offset += len;
        }

        Ok(envelopes)
    }
}
#[cfg(feature = "zstd-compression")]
impl<B: Bounds> LeafCompression<B> for BitshuffleCompression {
    const PAGE_SIZE: usize = 2048;

    fn compress<W: io::Write>(envelopes: &mut [Envelope<B>], write: &mut W) -> io::Result<()> {
        write.write_all(&(envelopes.len() as u16).to_le_bytes())?;

        if envelopes.is_empty() {
            return Ok(());
        }

        // Find max-spread dimension for sorting
        let mut spreads: Vec<SpreadSurvey> = (0..B::COORDINATES)
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
            .max_by(|(_, a), (_, b)| a.spread().total_cmp(&b.spread()))
            .unwrap();

        // Sort for spatial locality
        envelopes.sort_by(|a, b| a.bounds.get(dimension).total_cmp(&b.bounds.get(dimension)));

        // Element size: doc_id (4) + bounds (COORDINATES * 8)
        let element_size = 4 + B::COORDINATES * 8;
        let padded_count = envelopes.len().div_ceil(8);

        // Collect as bytes, padding is zeros
        let mut data = vec![0u8; padded_count * element_size];
        for (i, envelope) in envelopes.iter().enumerate() {
            let offset = i * element_size;
            data[offset..offset + 4].copy_from_slice(&envelope.doc_id.to_le_bytes());
            for j in 0..B::COORDINATES {
                let bound_offset = offset + 4 + j * 8;
                data[bound_offset..bound_offset + 8]
                    .copy_from_slice(&envelope.bounds.get(j).to_le_bytes());
            }
        }

        let shuffled = bitshuffle(&data, padded_count, element_size);
        let compressed = zstd::bulk::compress(&shuffled, 3)?; // level 3
        write.write_all(&(compressed.len() as u32).to_le_bytes())?;
        write.write_all(&compressed)?;
        write.write_all(&compressed)?;

        Ok(())
    }

    fn decompress(data: &[u8]) -> io::Result<Vec<Envelope<B>>> {
        let count = u16::from_le_bytes(data[0..2].try_into().unwrap()) as usize;
        if count == 0 {
            return Ok(Vec::new());
        }

        let element_size = 4 + B::COORDINATES * 8;
        let padded_count = count.div_ceil(8);

        let compressed_len = u32::from_le_bytes(data[2..6].try_into().unwrap()) as usize;
        let decompressed =
            zstd::bulk::decompress(&data[6..6 + compressed_len], padded_count * element_size)?;
        let unshuffled = bitunshuffle(&decompressed, padded_count, element_size);

        let mut envelopes: Vec<Envelope<B>> = Vec::with_capacity(count);
        for i in 0..count {
            let offset = i * element_size;
            let doc_id = u32::from_le_bytes(unshuffled[offset..offset + 4].try_into().unwrap());
            let mut envelope: Envelope<B> = Envelope::skeleton(doc_id);
            for j in 0..B::COORDINATES {
                let bound_offset = offset + 4 + j * 8;
                let value = f64::from_le_bytes(
                    unshuffled[bound_offset..bound_offset + 8]
                        .try_into()
                        .unwrap(),
                );
                envelope.bounds.set(j, value);
            }
            envelopes.push(envelope);
        }

        Ok(envelopes)
    }
}

impl<B: Bounds> LeafCompression<B> for RawCompression {
    const PAGE_SIZE: usize = 512;
    fn compress<W: io::Write>(envelopes: &mut [Envelope<B>], write: &mut W) -> io::Result<()> {
        write.write_all(&(envelopes.len() as u16).to_le_bytes())?;
        for envelope in envelopes.iter() {
            write.write_all(&envelope.doc_id.to_le_bytes())?;
        }
        for i in 0..B::COORDINATES {
            for envelope in envelopes.iter() {
                write.write_all(&envelope.bounds.get(i).to_le_bytes())?;
            }
        }
        Ok(())
    }
    fn decompress(data: &[u8]) -> io::Result<Vec<Envelope<B>>> {
        let count = u16::from_le_bytes(data[0..2].try_into().unwrap()) as usize;
        let mut envelopes: Vec<Envelope<B>> = Vec::with_capacity(count);

        let doc_id_start = 2;
        for i in 0..count {
            let offset = doc_id_start + i * 4;
            let doc_id = u32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
            envelopes.push(Envelope::skeleton(doc_id));
        }

        let coords_start = doc_id_start + count * 4;
        for dim in 0..B::COORDINATES {
            for i in 0..count {
                let offset = coords_start + (dim * count + i) * 8;
                let value = f64::from_le_bytes(data[offset..offset + 8].try_into().unwrap());
                envelopes[i].bounds.set(dim, value);
            }
        }
        Ok(envelopes)
    }
}
