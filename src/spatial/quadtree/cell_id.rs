// =============================================================================
// QuadtreeCellId - Z-Order Cell Encoding
// =============================================================================

/// Maximum level for quadtree cells (0-30).
/// At level 30 with typical geographic bounds, resolution is sub-centimeter.
pub const MAX_LEVEL: u8 = 30;

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

        // Shift position to top and add sentinel
        // Position uses 2*level bits, starting at bit 62 (just below bit 63)
        let shift = 62 - 2 * level as u32;
        let id = (position << shift) | (1u64 << shift.saturating_sub(1).max(0));

        // Actually, let's reconsider the encoding:
        // For level L, we have 2*L position bits at the top, then a sentinel 1.
        // The sentinel is at bit (62 - 2*L).
        // But we need position bits in 63..(64-2*L), sentinel at (62-2*L).
        // Hmm, let me think again...

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
    pub fn to_center(&self, bounds: &Bounds) -> Point2D {
        bounds.denormalize(&self.to_center_normalized())
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
    pub fn to_bounds(&self, bounds: &Bounds) -> Rect {
        let norm = self.to_bounds_normalized();
        Rect {
            x: Interval::new(
                bounds.min_x + norm.x.lo * bounds.width(),
                bounds.min_x + norm.x.hi * bounds.width(),
            ),
            y: Interval::new(
                bounds.min_y + norm.y.lo * bounds.height(),
                bounds.min_y + norm.y.hi * bounds.height(),
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
