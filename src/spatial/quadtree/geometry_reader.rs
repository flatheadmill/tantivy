//! HUSH

// =============================================================================
// GeometryCache Trait
// =============================================================================

use std::collections::HashMap;
use std::io;

use crate::spatial::quadtree::Point2D;

/// Cache statistics for monitoring geometry access patterns.
#[derive(Debug, Default, Clone)]
pub struct CacheStats {
    /// Number of cache hits.
    pub hits: u64,

    /// Number of cache misses.
    pub misses: u64,

    /// Number of evictions from cache.
    pub evictions: u64,

    /// Total bytes read from underlying storage.
    pub bytes_read: u64,
}

impl CacheStats {
    /// Returns the hit rate as a percentage (0.0 - 100.0).
    pub fn hit_rate(&self) -> f64 {
        let total = self.hits + self.misses;
        if total == 0 {
            0.0
        } else {
            100.0 * self.hits as f64 / total as f64
        }
    }
}

/// Cache for decoded geometries during merge.
///
/// Implementations may read from disk, memory, or a hybrid storage.
/// The cache should be sized to the working set of cells being split,
/// not total geometries in the segment.
///
/// # Thread Safety
///
/// GeometryCache implementations are not required to be thread-safe.
/// Each merge operation uses its own cache instance.
pub trait GeometryCache {
    /// Returns the vertices for the given document ID.
    ///
    /// The returned slice contains the polygon vertices in order.
    /// For closed polygons, the first and last vertices may be the same,
    /// or the caller should treat edge N-1 as connecting vertex N-1 to vertex 0.
    ///
    /// # Errors
    ///
    /// Returns an error if the geometry cannot be read (I/O error,
    /// missing document, corrupt data).
    fn get(&mut self, doc_id: u32) -> io::Result<&[Point2D]>;

    /// Prefetch hint for upcoming document IDs.
    ///
    /// The cache may use this to prefetch geometries in the background.
    /// The default implementation does nothing.
    fn prefetch(&mut self, _doc_ids: &[u32]) {}

    /// Returns cache statistics for monitoring.
    fn stats(&self) -> CacheStats {
        CacheStats::default()
    }
}

// =============================================================================
// GeometryReader Trait
// =============================================================================

/// Trait for reading raw geometry data from storage.
///
/// This is the low-level interface that backs a GeometryCache.
/// Implementations might read from:
/// - In-memory vectors
/// - Memory-mapped files
/// - Columnar storage
pub trait GeometryReader {
    /// Reads the geometry for the given document ID.
    ///
    /// Returns the vertices as a new Vec.
    fn read_geometry(&mut self, doc_id: u32) -> io::Result<Vec<Point2D>>;
}

// =============================================================================
// LruGeometryCache
// =============================================================================

/// An LRU cache for geometries backed by a GeometryReader.
///
/// Uses a simple Vec + HashMap implementation for LRU behavior:
/// - HashMap maps doc_id to index in the entries vector
/// - Entries are ordered by access time (most recent at end)
/// - When at capacity, the least recently used entry is evicted
pub struct LruGeometryCache<R: GeometryReader> {
    /// Underlying reader for cache misses.
    reader: R,

    /// Cached geometries: (doc_id, vertices).
    entries: Vec<(u32, Vec<Point2D>)>,

    /// Maps doc_id to index in entries for O(1) lookup.
    index: HashMap<u32, usize>,

    /// Maximum number of cached geometries.
    capacity: usize,

    /// Access statistics.
    stats: CacheStats,
}

impl<R: GeometryReader> LruGeometryCache<R> {
    /// Creates a new LRU cache with the given capacity.
    ///
    /// # Arguments
    ///
    /// * `reader` - The underlying geometry reader
    /// * `capacity` - Maximum number of geometries to cache
    pub fn new(reader: R, capacity: usize) -> Self {
        Self {
            reader,
            entries: Vec::with_capacity(capacity),
            index: HashMap::with_capacity(capacity),
            capacity: capacity.max(1), // At least 1
            stats: CacheStats::default(),
        }
    }

    /// Returns a reference to the underlying reader.
    pub fn reader(&self) -> &R {
        &self.reader
    }

    /// Returns the current number of cached entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns true if the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Clears all cached entries.
    pub fn clear(&mut self) {
        self.entries.clear();
        self.index.clear();
    }

    /// Moves an entry to the end (most recently used position).
    fn touch(&mut self, idx: usize) {
        if idx == self.entries.len() - 1 {
            return; // Already at end
        }

        let entry = self.entries.remove(idx);
        let doc_id = entry.0;

        // Update indices for shifted entries
        for i in idx..self.entries.len() {
            if let Some(old_idx) = self.index.get_mut(&self.entries[i].0) {
                *old_idx = i;
            }
        }

        // Add to end
        self.index.insert(doc_id, self.entries.len());
        self.entries.push(entry);
    }

    /// Evicts the least recently used entry.
    fn evict_lru(&mut self) {
        if self.entries.is_empty() {
            return;
        }

        let (doc_id, _) = self.entries.remove(0);
        self.index.remove(&doc_id);
        self.stats.evictions += 1;

        // Update indices for shifted entries
        for (i, (id, _)) in self.entries.iter().enumerate() {
            self.index.insert(*id, i);
        }
    }
}

impl<R: GeometryReader> GeometryCache for LruGeometryCache<R> {
    fn get(&mut self, doc_id: u32) -> io::Result<&[Point2D]> {
        // Check if already cached
        if let Some(&idx) = self.index.get(&doc_id) {
            self.stats.hits += 1;
            self.touch(idx);
            // After touch, it's at the end
            return Ok(&self.entries.last().unwrap().1);
        }

        self.stats.misses += 1;

        // Evict if at capacity
        if self.entries.len() >= self.capacity {
            self.evict_lru();
        }

        // Load from reader
        let vertices = self.reader.read_geometry(doc_id)?;
        self.stats.bytes_read += (vertices.len() * std::mem::size_of::<Point2D>()) as u64;

        // Insert at end (most recently used)
        let new_idx = self.entries.len();
        self.entries.push((doc_id, vertices));
        self.index.insert(doc_id, new_idx);

        Ok(&self.entries.last().unwrap().1)
    }

    fn prefetch(&mut self, doc_ids: &[u32]) {
        // Simple synchronous prefetch - just warm the cache
        for &doc_id in doc_ids {
            if !self.index.contains_key(&doc_id) {
                let _ = self.get(doc_id); // Ignore errors in prefetch
            }
        }
    }

    fn stats(&self) -> CacheStats {
        self.stats.clone()
    }
}

// =============================================================================
// InMemoryGeometryReader
// =============================================================================

/// A simple in-memory geometry reader for testing.
///
/// Stores all geometries in a HashMap.
pub struct InMemoryGeometryReader {
    geometries: HashMap<u32, Vec<Point2D>>,
}

impl InMemoryGeometryReader {
    /// Creates a new empty reader.
    pub fn new() -> Self {
        Self {
            geometries: HashMap::new(),
        }
    }

    /// Adds a geometry for the given document ID.
    pub fn add(&mut self, doc_id: u32, vertices: Vec<Point2D>) {
        self.geometries.insert(doc_id, vertices);
    }

    /// Returns the number of stored geometries.
    pub fn len(&self) -> usize {
        self.geometries.len()
    }

    /// Returns true if empty.
    pub fn is_empty(&self) -> bool {
        self.geometries.is_empty()
    }
}

impl Default for InMemoryGeometryReader {
    fn default() -> Self {
        Self::new()
    }
}

impl GeometryReader for InMemoryGeometryReader {
    fn read_geometry(&mut self, doc_id: u32) -> io::Result<Vec<Point2D>> {
        self.geometries.get(&doc_id).cloned().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("geometry not found for doc_id {}", doc_id),
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // fn test_bounds() -> Bounds {
    // Bounds::new(0.0, 0.0, 100.0, 100.0)
    // }

    // -------------------------------------------------------------------------
    // CacheStats Tests
    // -------------------------------------------------------------------------

    mod cache_stats_tests {
        use super::*;

        #[test]
        fn test_hit_rate() {
            let mut stats = CacheStats::default();
            assert_eq!(stats.hit_rate(), 0.0);

            stats.hits = 75;
            stats.misses = 25;
            assert!((stats.hit_rate() - 75.0).abs() < 0.01);

            stats.hits = 100;
            stats.misses = 0;
            assert!((stats.hit_rate() - 100.0).abs() < 0.01);
        }
    }

    // -------------------------------------------------------------------------
    // InMemoryGeometryReader Tests
    // -------------------------------------------------------------------------

    mod geometry_reader_tests {
        use super::*;

        #[test]
        fn test_add_and_read() {
            let mut reader = InMemoryGeometryReader::new();

            let vertices = vec![
                Point2D::new(0.0, 0.0),
                Point2D::new(10.0, 0.0),
                Point2D::new(10.0, 10.0),
                Point2D::new(0.0, 10.0),
            ];

            reader.add(42, vertices.clone());

            let result = reader.read_geometry(42).unwrap();
            assert_eq!(result.len(), 4);
            assert_eq!(result[0], vertices[0]);
        }

        #[test]
        fn test_not_found() {
            let mut reader = InMemoryGeometryReader::new();
            let result = reader.read_geometry(999);
            assert!(result.is_err());
        }
    }

    // -------------------------------------------------------------------------
    // LruGeometryCache Tests
    // -------------------------------------------------------------------------

    mod lru_cache_tests {
        use super::*;

        fn make_cache() -> LruGeometryCache<InMemoryGeometryReader> {
            let mut reader = InMemoryGeometryReader::new();
            reader.add(1, vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 1.0)]);
            reader.add(2, vec![Point2D::new(2.0, 2.0), Point2D::new(3.0, 3.0)]);
            reader.add(3, vec![Point2D::new(4.0, 4.0), Point2D::new(5.0, 5.0)]);

            LruGeometryCache::new(reader, 2)
        }

        #[test]
        fn test_cache_hit() {
            let mut cache = make_cache();

            // First access - miss
            let _ = cache.get(1).unwrap();
            assert_eq!(cache.stats().misses, 1);
            assert_eq!(cache.stats().hits, 0);

            // Second access - hit
            let _ = cache.get(1).unwrap();
            assert_eq!(cache.stats().misses, 1);
            assert_eq!(cache.stats().hits, 1);
        }

        #[test]
        fn test_cache_eviction() {
            let mut cache = make_cache();

            // Fill cache (capacity = 2)
            let _ = cache.get(1).unwrap();
            let _ = cache.get(2).unwrap();
            assert_eq!(cache.len(), 2);
            assert_eq!(cache.stats().evictions, 0);

            // Access 3 - should evict 1 (LRU)
            let _ = cache.get(3).unwrap();
            assert_eq!(cache.len(), 2);
            assert_eq!(cache.stats().evictions, 1);

            // Access 1 again - should be a miss (was evicted)
            let _ = cache.get(1).unwrap();
            assert_eq!(cache.stats().misses, 4); // 1, 2, 3, 1 again
        }

        #[test]
        fn test_lru_touch() {
            let mut cache = make_cache();

            // Access 1, 2
            let _ = cache.get(1).unwrap();
            let _ = cache.get(2).unwrap();

            // Touch 1 (makes it most recent)
            let _ = cache.get(1).unwrap();

            // Access 3 - should evict 2 (now LRU)
            let _ = cache.get(3).unwrap();

            // 1 should still be cached
            let old_misses = cache.stats().misses;
            let _ = cache.get(1).unwrap();
            assert_eq!(cache.stats().misses, old_misses); // No new miss
        }
    }
}
