// =============================================================================
// QuadtreeCellId - Z-Order Cell Encoding
// =============================================================================

use crate::spatial::quadtree::{Interval, Point2D, Rect};

/// Maximum level for quadtree cells (0-30).
/// At level 30 with typical geographic bounds, resolution is sub-centimeter.
pub const MAX_LEVEL: u8 = 30;

/// Interleaves bits of two 32-bit values into a 64-bit Z-order code.
///
/// The result has bits of x in even positions and bits of y in odd positions.
fn interleave_bits(x: u32, y: u32) -> u64 {
    let mut x = x as u64;
    let mut y = y as u64;

    // Spread x bits: 0000abcd -> 0a0b0c0d
    x = (x | (x << 16)) & 0x0000_FFFF_0000_FFFF;
    x = (x | (x << 8)) & 0x00FF_00FF_00FF_00FF;
    x = (x | (x << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    x = (x | (x << 2)) & 0x3333_3333_3333_3333;
    x = (x | (x << 1)) & 0x5555_5555_5555_5555;

    // Spread y bits
    y = (y | (y << 16)) & 0x0000_FFFF_0000_FFFF;
    y = (y | (y << 8)) & 0x00FF_00FF_00FF_00FF;
    y = (y | (y << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    y = (y | (y << 2)) & 0x3333_3333_3333_3333;
    y = (y | (y << 1)) & 0x5555_5555_5555_5555;

    // Combine: x in even bits, y in odd bits
    x | (y << 1)
}

/// Deinterleaves a Z-order code back into (x, y) coordinates.
fn deinterleave_bits(z: u64, level: u8) -> (u32, u32) {
    let mask = (1u64 << (2 * level as u32)) - 1;
    let z = z & mask;

    // Extract even bits (x) and odd bits (y)
    let mut x = z & 0x5555_5555_5555_5555;
    let mut y = (z >> 1) & 0x5555_5555_5555_5555;

    // Compact x bits
    x = (x | (x >> 1)) & 0x3333_3333_3333_3333;
    x = (x | (x >> 2)) & 0x0F0F_0F0F_0F0F_0F0F;
    x = (x | (x >> 4)) & 0x00FF_00FF_00FF_00FF;
    x = (x | (x >> 8)) & 0x0000_FFFF_0000_FFFF;
    x = (x | (x >> 16)) & 0x0000_0000_FFFF_FFFF;

    // Compact y bits
    y = (y | (y >> 1)) & 0x3333_3333_3333_3333;
    y = (y | (y >> 2)) & 0x0F0F_0F0F_0F0F_0F0F;
    y = (y | (y >> 4)) & 0x00FF_00FF_00FF_00FF;
    y = (y | (y >> 8)) & 0x0000_FFFF_0000_FFFF;
    y = (y | (y >> 16)) & 0x0000_0000_FFFF_FFFF;

    (x as u32, y as u32)
}

/// A cell ID using Z-order (Morton code) encoding.
///
/// The encoding uses 64 bits:
/// - Bits 63-1: Position (2 bits per level, up to 30 levels = 60 bits)
/// - Bit 0: Always 1 when shifted to mark the level (sentinel)
///
/// The cell ID format follows the pattern used by S2CellId:
/// - The position bits are stored at the high end
/// - A trailing 1-bit marks where the meaningful bits end
/// - The level can be computed by finding the position of the lowest set bit
///
/// Example for level 2 cell at position (1, 2) = binary 01, 10 -> Z-order 0110:
/// ```text
/// bit:  63      60 59                                    4  3  2  1  0
///       |-------|--|------------------------------------|----|----|----|
///       | 0110  |  |           ... zeros ...            |    | 1  |    |
///       position                                        sentinel
/// ```
/// The sentinel 1-bit is at position (62 - 2*level) = 62 - 4 = 58... actually
/// let me recalculate.
///
/// Actually, the encoding is:
/// - Position bits at the top (MSB)
/// - Sentinel 1-bit immediately after position bits
/// - All trailing bits are 0
///
/// For level L: the sentinel is at bit position 62 - 2*L
/// For level 0: position is in bits 63-62, sentinel at bit 61
/// For level 30: position uses all 60 bits (63-4), sentinel at bit 3
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QuadtreeCellId(pub u64);

impl QuadtreeCellId {
    /// The invalid/sentinel cell ID (value 0).
    pub const NONE: Self = Self(0);

    /// The root cell ID at level 0 covering the entire space.
    /// Position = 0, level = 0, so sentinel at bit 62.
    pub const ROOT: Self = Self(1 << 62);

    /// Maximum level for quadtree cells (0-30).
    /// At level 30 with typical geographic bounds, resolution is sub-centimeter.
    pub const MAX_LEVEL: u8 = MAX_LEVEL;

    /// Creates a cell ID from a raw u64 value.
    #[inline]
    pub const fn from_raw(value: u64) -> Self {
        Self(value)
    }

    /// Returns the raw u64 value.
    #[inline]
    pub const fn to_raw(&self) -> u64 {
        self.0
    }

    /// Checks if this is a valid cell ID (non-zero with proper sentinel).
    #[inline]
    pub fn is_valid(&self) -> bool {
        self.0 != 0 && self.0.trailing_zeros() % 2 == 0
    }

    /// Returns the level of this cell (0 = root, MAX_LEVEL = finest).
    ///
    /// The level is determined by the position of the sentinel bit.
    #[inline]
    pub fn level(&self) -> u8 {
        debug_assert!(self.is_valid());
        // Sentinel is at bit (62 - 2*level), so level = (62 - sentinel_pos) / 2
        // trailing_zeros gives us the sentinel position
        let sentinel_pos = self.0.trailing_zeros();
        ((62 - sentinel_pos) / 2) as u8
    }

    /// Returns the first cell at the given level.
    ///
    /// This is the cell at integer coordinates (0, 0) at the specified level,
    /// which is the first cell in Z-order traversal at that level.
    ///
    /// # Arguments
    ///
    /// * `level` - The cell level (0 = root, MAX_LEVEL = finest)
    ///
    /// # Example
    ///
    /// ```ignore
    /// let first_leaf = QuadtreeCellId::begin(QuadtreeCellId::MAX_LEVEL);
    /// assert_eq!(first_leaf.to_ij(), (0, 0));
    /// ```
    #[inline]
    pub fn begin(level: u8) -> Self {
        debug_assert!(level <= MAX_LEVEL);
        Self::from_ij(0, 0, level)
    }

    /// Returns the cell just past the last cell at the given level.
    ///
    /// This is useful for iteration bounds. The returned value is NOT a valid
    /// cell ID but represents the exclusive upper bound of cells at this level.
    ///
    /// # Arguments
    ///
    /// * `level` - The cell level (0 = root, MAX_LEVEL = finest)
    #[inline]
    pub fn end(level: u8) -> Self {
        debug_assert!(level <= MAX_LEVEL);
        // The "end" sentinel is one past the last valid cell
        // At level L, there are 4^L cells, so end is at position 4^L
        // which equals 2^(2*L)
        let sentinel_pos = 62 - 2 * level as u32;
        let max_position = 1u64 << (2 * level as u32);
        Self((max_position << (sentinel_pos + 1)) | (1u64 << sentinel_pos))
    }

    /// Returns the next cell at the same level.
    ///
    /// Returns the cell immediately following this one in Z-order at the same
    /// level. If this is the last cell at the level, returns `end(level)`.
    pub fn next(&self) -> Self {
        let level = self.level();
        let sentinel_pos = 62 - 2 * level as u32;

        // Add 1 to the position bits
        let lsb = 1u64 << (sentinel_pos + 1);
        Self(self.0 + lsb)
    }

    /// Returns the previous cell at the same level.
    ///
    /// Returns the cell immediately preceding this one in Z-order at the same
    /// level. Panics in debug mode if called on `begin(level)`.
    pub fn prev(&self) -> Self {
        let level = self.level();
        let sentinel_pos = 62 - 2 * level as u32;

        debug_assert!(self.0 > Self::begin(level).0, "prev() called on begin()");

        // Subtract 1 from the position bits
        let lsb = 1u64 << (sentinel_pos + 1);
        Self(self.0 - lsb)
    }

    /// Creates a cell ID from normalized coordinates at the given level.
    ///
    /// Coordinates should be in [0, 1) range. Values exactly at 1.0 are
    /// clamped to just below 1.0.
    pub fn from_normalized(x: f64, y: f64, level: u8) -> Self {
        debug_assert!(level <= MAX_LEVEL);

        // Clamp coordinates to valid range
        let x = x.clamp(0.0, 1.0 - f64::EPSILON);
        let y = y.clamp(0.0, 1.0 - f64::EPSILON);

        // Convert to integer coordinates at this level
        let scale = (1u64 << level) as f64;
        let ix = (x * scale) as u32;
        let iy = (y * scale) as u32;

        Self::from_ij(ix, iy, level)
    }

    /// Creates a cell ID from world (x, y) coordinates at the given level.
    ///
    /// This is a convenience method that normalizes coordinates using the provided
    /// bounds before creating the cell ID. The coordinates are first converted to
    /// the [0, 1) normalized space, then mapped to the appropriate cell.
    ///
    /// # Arguments
    ///
    /// * `x` - The x-coordinate in world space
    /// * `y` - The y-coordinate in world space
    /// * `level` - The cell level (0 = root, MAX_LEVEL = finest)
    /// * `bounds` - The global bounds used for coordinate normalization
    ///
    /// # Example
    ///
    /// ```ignore
    /// let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);
    /// let cell = QuadtreeCellId::from_xy(50.0, 50.0, 5, &bounds);
    /// // Creates cell containing the center point at level 5
    /// ```
      pub fn from_xy(x: f64, y: f64, level: u8, bounds: &Rect) -> Self {
          let nx = (x - bounds.x.lo) / (bounds.x.hi - bounds.x.lo);
          let ny = (y - bounds.y.lo) / (bounds.y.hi - bounds.y.lo);
          Self::from_normalized(nx, ny, level)
      }

    /// Creates a cell ID from integer (i, j) coordinates at the given level.
    ///
    /// i is the x-coordinate (column), j is the y-coordinate (row).
    /// Both must be in range [0, 2^level).
    pub fn from_ij(i: u32, j: u32, level: u8) -> Self {
        debug_assert!(level <= MAX_LEVEL);
        let max_coord = 1u32 << level;
        debug_assert!(i < max_coord && j < max_coord);

        // Interleave i and j bits to create Z-order (Morton) code
        let position = interleave_bits(i, j);

        // Simpler approach: position in high bits, sentinel marks end
        let sentinel_pos = 62 - 2 * level as u32;
        let id = (position << (sentinel_pos + 1)) | (1u64 << sentinel_pos);

        Self(id)
    }

    /// Extracts the (i, j) integer coordinates at this cell's level.
    pub fn to_ij(&self) -> (u32, u32) {
        let level = self.level();
        let sentinel_pos = 62 - 2 * level as u32;
        let position = self.0 >> (sentinel_pos + 1);
        deinterleave_bits(position, level)
    }

    /// Returns the center of this cell in normalized [0, 1] coordinates.
    pub fn to_center_normalized(&self) -> Point2D {
        let level = self.level();
        let (i, j) = self.to_ij();
        let scale = (1u64 << level) as f64;

        Point2D {
            x: (i as f64 + 0.5) / scale,
            y: (j as f64 + 0.5) / scale,
        }
    }

    /// Returns the center of this cell in world coordinates.
  pub fn to_center(&self, bounds: &Rect) -> Point2D {
      let norm = self.to_center_normalized();
      Point2D {
          x: bounds.x.lo + norm.x * (bounds.x.hi - bounds.x.lo),
          y: bounds.y.lo + norm.y * (bounds.y.hi - bounds.y.lo),
      }
  }

    /// Returns the bounds of this cell in normalized [0, 1] coordinates.
    pub fn to_bounds_normalized(&self) -> Rect {
        let level = self.level();
        let (i, j) = self.to_ij();
        let scale = (1u64 << level) as f64;

        Rect {
            x: Interval::new(i as f64 / scale, (i + 1) as f64 / scale),
            y: Interval::new(j as f64 / scale, (j + 1) as f64 / scale),
        }
    }

    /// Returns the bounds of this cell in world coordinates.
  pub fn to_bounds(&self, bounds: &Rect) -> Rect {
      let norm = self.to_bounds_normalized();
      let width = bounds.x.hi - bounds.x.lo;
      let height = bounds.y.hi - bounds.y.lo;
      Rect {
          x: Interval::new(
              bounds.x.lo + norm.x.lo * width,
              bounds.x.lo + norm.x.hi * width,
          ),
          y: Interval::new(
              bounds.y.lo + norm.y.lo * height,
              bounds.y.lo + norm.y.hi * height,
          ),
      }
  }

    /// Returns the parent cell (one level up).
    ///
    /// Returns None if this is already the root cell.
    pub fn parent(&self) -> Option<Self> {
        if self.level() == 0 {
            return None;
        }

        // Move sentinel two bits to the right and clear the old position bits
        let sentinel_pos = 62 - 2 * self.level() as u32;
        let new_sentinel_pos = sentinel_pos + 2;
        let mask = !((1u64 << (new_sentinel_pos + 1)) - 1);
        let id = (self.0 & mask) | (1u64 << new_sentinel_pos);

        Some(Self(id))
    }

    /// Returns the four child cells (one level down).
    ///
    /// Returns None if this is already at MAX_LEVEL.
    pub fn children(&self) -> Option<[Self; 4]> {
        let level = self.level();
        if level >= MAX_LEVEL {
            return None;
        }

        let (i, j) = self.to_ij();
        let new_level = level + 1;

        Some([
            Self::from_ij(2 * i, 2 * j, new_level),
            Self::from_ij(2 * i + 1, 2 * j, new_level),
            Self::from_ij(2 * i, 2 * j + 1, new_level),
            Self::from_ij(2 * i + 1, 2 * j + 1, new_level),
        ])
    }

    /// Returns the child at the given index (0-3).
    ///
    /// Index layout (Z-order):
    /// - 0: bottom-left
    /// - 1: bottom-right
    /// - 2: top-left
    /// - 3: top-right
    pub fn child(&self, index: usize) -> Option<Self> {
        debug_assert!(index < 4);
        let level = self.level();
        if level >= MAX_LEVEL {
            return None;
        }

        let (i, j) = self.to_ij();
        let di = (index & 1) as u32;
        let dj = ((index >> 1) & 1) as u32;

        Some(Self::from_ij(2 * i + di, 2 * j + dj, level + 1))
    }

    /// Returns which child (0-3) this cell is of its parent.
    ///
    /// Returns the child index in Z-order:
    /// - 0: bottom-left  (i even, j even)
    /// - 1: bottom-right (i odd, j even)
    /// - 2: top-left     (i even, j odd)
    /// - 3: top-right    (i odd, j odd)
    ///
    /// Returns 0 for the root cell (which has no parent).
    pub fn child_position(&self) -> usize {
        let (i, j) = self.to_ij();
        ((i & 1) | ((j & 1) << 1)) as usize
    }

    /// Returns the minimum cell ID in the range covered by this cell.
    ///
    /// This is the cell ID of the first (smallest) leaf cell contained
    /// within this cell. For leaf cells, returns self.
    pub fn range_min(&self) -> Self {
        // The minimum is this cell's position with the sentinel moved to MAX_LEVEL
        let level = self.level();
        let sentinel_pos = 62 - 2 * level as u32;
        let position = self.0 & !((1u64 << (sentinel_pos + 1)) - 1);

        // Set sentinel at MAX_LEVEL position
        let new_sentinel_pos = 62 - 2 * MAX_LEVEL as u32;
        Self(position | (1u64 << new_sentinel_pos))
    }

    /// Returns the maximum cell ID in the range covered by this cell.
    ///
    /// This is the cell ID of the last (largest) leaf cell contained
    /// within this cell. For leaf cells, returns self.
    pub fn range_max(&self) -> Self {
        // The maximum fills all descendant bits with 1s
        let level = self.level();
        let sentinel_pos = 62 - 2 * level as u32;

        // Keep position bits, fill lower bits with 1s up to sentinel
        let position = self.0 & !((1u64 << (sentinel_pos + 1)) - 1);
        let fill_mask = (1u64 << (sentinel_pos + 1)) - (1u64 << (62 - 2 * MAX_LEVEL as u32));

        // Set sentinel at MAX_LEVEL position
        let new_sentinel_pos = 62 - 2 * MAX_LEVEL as u32;
        Self(position | fill_mask | (1u64 << new_sentinel_pos))
    }

    /// Checks if this cell contains another cell (other is a descendant).
    pub fn contains(&self, other: &Self) -> bool {
        other.0 >= self.range_min().0 && other.0 <= self.range_max().0
    }

    /// Checks if this cell intersects another cell (one contains the other or they're equal).
    pub fn intersects(&self, other: &Self) -> bool {
        self.contains(other) || other.contains(self)
    }

    /// Returns the common ancestor of this cell and another.
    pub fn common_ancestor(&self, other: &Self) -> Self {
        // XOR to find differing bits, then mask to common prefix
        let diff = self.0 ^ other.0;
        if diff == 0 {
            return *self;
        }

        // Find the highest differing bit
        let highest_diff = 63 - diff.leading_zeros();

        // The common ancestor level is where the first difference occurs
        // Each level uses 2 bits starting from bit 62
        let diff_level = (62 - highest_diff) / 2;

        // Get parent at that level
        let mut result = *self;
        while result.level() > diff_level as u8 {
            result = result.parent().unwrap();
        }
        result
    }
}
