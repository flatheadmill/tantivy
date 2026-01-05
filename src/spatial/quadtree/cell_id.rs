// =============================================================================
// QuadtreeCellId - Z-Order Cell Encoding
// =============================================================================

use crate::spatial::quadtree::{Bounds, Interval, Point2D, Rect};

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

#[cfg(test)]
mod tests {
    use super::*;
    mod cell_id_tests {
        use super::*;
        use crate::spatial::quadtree::Bounds;

        #[test]
        fn test_root_cell() {
            let root = QuadtreeCellId::ROOT;
            assert!(root.is_valid());
            assert_eq!(root.level(), 0);

            // Root should cover entire space
            let bounds = root.to_bounds_normalized();
            assert!((bounds.x.lo - 0.0).abs() < 1e-10);
            assert!((bounds.x.hi - 1.0).abs() < 1e-10);
            assert!((bounds.y.lo - 0.0).abs() < 1e-10);
            assert!((bounds.y.hi - 1.0).abs() < 1e-10);
        }

        #[test]
        fn test_none_cell() {
            let none = QuadtreeCellId::NONE;
            assert!(!none.is_valid());
        }

        #[test]
        fn test_from_ij_level() {
            // Level 1: 2x2 grid
            let cell_00 = QuadtreeCellId::from_ij(0, 0, 1);
            let cell_10 = QuadtreeCellId::from_ij(1, 0, 1);
            let cell_01 = QuadtreeCellId::from_ij(0, 1, 1);
            let cell_11 = QuadtreeCellId::from_ij(1, 1, 1);

            assert!(cell_00.is_valid());
            assert_eq!(cell_00.level(), 1);
            assert_eq!(cell_10.level(), 1);
            assert_eq!(cell_01.level(), 1);
            assert_eq!(cell_11.level(), 1);

            // Should be ordered by Z-order
            assert!(cell_00 < cell_10);
            assert!(cell_10 < cell_01);
            assert!(cell_01 < cell_11);
        }

        #[test]
        fn test_to_ij_roundtrip() {
            for level in 0..=10 {
                let max = 1u32 << level;
                for i in [0, max / 4, max / 2, max - 1].iter().filter(|&&x| x < max) {
                    for j in [0, max / 4, max / 2, max - 1].iter().filter(|&&x| x < max) {
                        let cell = QuadtreeCellId::from_ij(*i, *j, level);
                        let (ri, rj) = cell.to_ij();
                        assert_eq!(
                            (*i, *j),
                            (ri, rj),
                            "Roundtrip failed for level={}, i={}, j={}",
                            level,
                            i,
                            j
                        );
                    }
                }
            }
        }

        #[test]
        fn test_from_normalized() {
            // Corner points
            let bl = QuadtreeCellId::from_normalized(0.0, 0.0, 2);
            assert_eq!(bl.to_ij(), (0, 0));

            // Near top-right (but not exactly 1.0)
            let tr = QuadtreeCellId::from_normalized(0.999, 0.999, 2);
            assert_eq!(tr.to_ij(), (3, 3));

            // Center should be (2, 2) at level 2
            let center = QuadtreeCellId::from_normalized(0.5, 0.5, 2);
            assert_eq!(center.to_ij(), (2, 2));
        }

        #[test]
        fn test_to_center() {
            let bounds = Bounds::new(0.0, 0.0, 100.0, 100.0);

            // Root cell center should be (50, 50)
            let root_center = QuadtreeCellId::ROOT.to_center(&bounds);
            assert!((root_center.x - 50.0).abs() < 1e-10);
            assert!((root_center.y - 50.0).abs() < 1e-10);

            // Level 1 cell (0,0) should have center at (25, 25)
            let cell = QuadtreeCellId::from_ij(0, 0, 1);
            let center = cell.to_center(&bounds);
            assert!((center.x - 25.0).abs() < 1e-10);
            assert!((center.y - 25.0).abs() < 1e-10);
        }

        #[test]
        fn test_to_bounds() {
            let global = Bounds::new(0.0, 0.0, 100.0, 100.0);

            // Root cell should cover entire bounds
            let root_bounds = QuadtreeCellId::ROOT.to_bounds(&global);
            assert!((root_bounds.x.lo - 0.0).abs() < 1e-10);
            assert!((root_bounds.x.hi - 100.0).abs() < 1e-10);

            // Level 1 cell (0,0) should cover bottom-left quadrant
            let cell = QuadtreeCellId::from_ij(0, 0, 1);
            let cell_bounds = cell.to_bounds(&global);
            assert!((cell_bounds.x.lo - 0.0).abs() < 1e-10);
            assert!((cell_bounds.x.hi - 50.0).abs() < 1e-10);
            assert!((cell_bounds.y.lo - 0.0).abs() < 1e-10);
            assert!((cell_bounds.y.hi - 50.0).abs() < 1e-10);
        }

        #[test]
        fn test_parent() {
            let root = QuadtreeCellId::ROOT;
            assert!(root.parent().is_none());

            let level1 = QuadtreeCellId::from_ij(0, 0, 1);
            let parent = level1.parent().unwrap();
            assert_eq!(parent.level(), 0);

            let level2 = QuadtreeCellId::from_ij(2, 3, 2);
            let parent2 = level2.parent().unwrap();
            assert_eq!(parent2.level(), 1);
            assert_eq!(parent2.to_ij(), (1, 1)); // (2,3) -> (1,1) at parent level
        }

        #[test]
        fn test_children() {
            let root = QuadtreeCellId::ROOT;
            let children = root.children().unwrap();

            assert_eq!(children.len(), 4);
            for child in &children {
                assert_eq!(child.level(), 1);
            }

            // Children should be in Z-order
            assert_eq!(children[0].to_ij(), (0, 0));
            assert_eq!(children[1].to_ij(), (1, 0));
            assert_eq!(children[2].to_ij(), (0, 1));
            assert_eq!(children[3].to_ij(), (1, 1));

            // All children's parent should be root
            for child in &children {
                assert_eq!(child.parent().unwrap(), root);
            }
        }

        #[test]
        fn test_child() {
            let root = QuadtreeCellId::ROOT;

            let c0 = root.child(0).unwrap();
            let c1 = root.child(1).unwrap();
            let c2 = root.child(2).unwrap();
            let c3 = root.child(3).unwrap();

            assert_eq!(c0.to_ij(), (0, 0));
            assert_eq!(c1.to_ij(), (1, 0));
            assert_eq!(c2.to_ij(), (0, 1));
            assert_eq!(c3.to_ij(), (1, 1));
        }

        #[test]
        fn test_max_level_no_children() {
            let max_cell = QuadtreeCellId::from_ij(0, 0, MAX_LEVEL);
            assert!(max_cell.children().is_none());
            assert!(max_cell.child(0).is_none());
        }

        #[test]
        fn test_contains() {
            let root = QuadtreeCellId::ROOT;
            let level1 = QuadtreeCellId::from_ij(0, 0, 1);
            let level2 = QuadtreeCellId::from_ij(0, 0, 2);

            // Root contains everything
            assert!(root.contains(&root));
            assert!(root.contains(&level1));
            assert!(root.contains(&level2));

            // Parent contains child
            assert!(level1.contains(&level2));

            // Child doesn't contain parent
            assert!(!level2.contains(&level1));
            assert!(!level1.contains(&root));

            // Different branches don't contain each other
            let other = QuadtreeCellId::from_ij(1, 1, 1);
            assert!(!level1.contains(&other));
            assert!(!other.contains(&level1));
        }

        #[test]
        fn test_intersects() {
            let cell_a = QuadtreeCellId::from_ij(0, 0, 2);
            let cell_b = QuadtreeCellId::from_ij(0, 0, 3); // child of cell_a
            let cell_c = QuadtreeCellId::from_ij(3, 3, 2); // different cell

            assert!(cell_a.intersects(&cell_b)); // parent-child
            assert!(cell_b.intersects(&cell_a)); // child-parent
            assert!(!cell_a.intersects(&cell_c)); // disjoint
        }

        #[test]
        fn test_range_min_max() {
            let cell = QuadtreeCellId::from_ij(0, 0, 2);

            let min = cell.range_min();
            let max = cell.range_max();

            assert!(min.is_valid());
            assert!(max.is_valid());
            assert_eq!(min.level(), MAX_LEVEL);
            assert_eq!(max.level(), MAX_LEVEL);
            assert!(min < max || min == max);

            // All children should be in range
            if let Some(children) = cell.children() {
                for child in &children {
                    assert!(child.0 >= min.0);
                    assert!(child.0 <= max.0);
                }
            }
        }

        #[test]
        fn test_common_ancestor() {
            let a = QuadtreeCellId::from_ij(0, 0, 4);
            let b = QuadtreeCellId::from_ij(1, 0, 4);

            let ancestor = a.common_ancestor(&b);
            assert!(ancestor.contains(&a));
            assert!(ancestor.contains(&b));

            // Same cell should return itself
            assert_eq!(a.common_ancestor(&a), a);

            // Root should be common ancestor of any two cells in different quadrants
            let c = QuadtreeCellId::from_ij(0, 0, 1);
            let d = QuadtreeCellId::from_ij(1, 1, 1);
            let root_ancestor = c.common_ancestor(&d);
            assert_eq!(root_ancestor, QuadtreeCellId::ROOT);
        }
    }
    mod interleave_tests {
        use super::*;

        #[test]
        fn test_interleave_basic() {
            // (0, 0) -> 0
            assert_eq!(interleave_bits(0, 0), 0);

            // (1, 0) -> 1 (x in bit 0)
            assert_eq!(interleave_bits(1, 0), 0b01);

            // (0, 1) -> 2 (y in bit 1)
            assert_eq!(interleave_bits(0, 1), 0b10);

            // (1, 1) -> 3
            assert_eq!(interleave_bits(1, 1), 0b11);

            // (2, 0) -> 4 (x bit 1 goes to bit 2)
            assert_eq!(interleave_bits(2, 0), 0b0100);

            // (0, 2) -> 8 (y bit 1 goes to bit 3)
            assert_eq!(interleave_bits(0, 2), 0b1000);
        }

        #[test]
        fn test_deinterleave_roundtrip() {
            for level in 1..=10 {
                let max = 1u32 << level;
                for x in [0, 1, max / 2, max - 1] {
                    for y in [0, 1, max / 2, max - 1] {
                        let z = interleave_bits(x, y);
                        let (dx, dy) = deinterleave_bits(z, level);
                        assert_eq!(
                            (x, y),
                            (dx, dy),
                            "Roundtrip failed for level={}, x={}, y={}",
                            level,
                            x,
                            y
                        );
                    }
                }
            }
        }
    }
}
