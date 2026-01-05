//! Quadtree implementation for 2D spatial indexing.
//!
//! This module provides a multi-shape quadtree for Tantivy's spatial indexing, based on the
//! architecture of Google S2's MutableS2ShapeIndex but adapted for planar geometry instead of
//! unit-sphere geometry.

mod builder;
mod cell_id;
mod geometry;
mod index;
mod interior_tracker;
mod merge;
mod predicates;
mod serialization;

pub use builder::{BuilderOptions, QuadtreeIndex, QuadtreeIndexBuilder};
pub use cell_id::{QuadtreeCellId, MAX_LEVEL};
pub use geometry::{Bounds, Interval, Point2D, Rect};
pub use index::{ClippedShape, PaddedCell, QuadtreeCell};
pub use interior_tracker::{contains_tracker_origin, InteriorTracker, ShapeIdSet};
pub use merge::{
    merge, CacheStats, DeleteBitSet, DocIdMap, GeometryCache, GeometryReader,
    InMemoryGeometryReader, InputIterator, LruGeometryCache, MergeOptions, MergeStats,
};
pub use predicates::{
    brute_force_contains_2d, compute_origin_inside_2d, contains_with_edges_2d,
    contains_with_origin_2d, crossing_sign_2d, edge_or_vertex_crossing_2d, orient_2d,
    EdgeCrosser2D,
};
