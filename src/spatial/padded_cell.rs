  pub trait PaddedCell: Clone + Sized {
      type Point: Copy;
      type Rect: Copy;
      type CellId: Copy + Eq + Ord;

      // === Identity ===
      fn id(&self) -> Self::CellId;
      fn level(&self) -> u8;

      // === Hierarchy ===
      /// Child by spatial index (0-3), not curve order.
      fn child(&self, index: usize) -> Option<Self>;

      fn children(&self) -> Option<[Self; 4]> {
          Some([self.child(0)?, self.child(1)?, self.child(2)?, self.child(3)?])
      }

      // === Geometry ===
      fn center(&self) -> Self::Point;
      fn bounds(&self) -> Self::Rect;
      fn padded_bounds(&self) -> Self::Rect;
      fn middle(&mut self) -> Option<Self::Rect>;

      // === Curve Traversal ===
      /// Entry vertex where space-filling curve enters cell.
      fn entry_vertex(&self) -> Self::Point;

      /// Exit vertex where space-filling curve exits cell.
      fn exit_vertex(&self) -> Self::Point;

      /// Curve position (0-3) to spatial child index.
      /// Z-order: identity. Hilbert: orientation-dependent.
      fn curve_pos_to_child(&self, pos: usize) -> usize;

      /// Spatial child index to curve position.
      fn child_to_curve_pos(&self, child: usize) -> usize;

      // === Operations ===
      fn shrink_to_fit(&self, bounds: &Self::Rect) -> Self::CellId;
      fn contains_point(&self, point: &Self::Point) -> bool;
  }
