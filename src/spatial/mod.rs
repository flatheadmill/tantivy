//! Spatial indexing support for geometric data.
//!
//! This module provides BKD (Block KD-Tree) spatial indexing capabilities
//! for Tantivy, enabling efficient spatial queries on geometric data types.

pub mod triangle;

pub use triangle::Triangle;
