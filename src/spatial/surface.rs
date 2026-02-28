  pub trait Surface {
      type Point: Copy;
      type Rect: Copy;
      type CellId: Copy + Eq + Ord;
      type Crosser;
      type PaddedCell: PaddedCell<
          Point = Self::Point,
          Rect = Self::Rect,
          CellId = Self::CellId,
      >;

      const MAX_LEVEL: u8;

      // === Cell Construction ===
      fn root_cell(&self, padding: f64) -> Self::PaddedCell;
      fn padded_cell(&self, id: Self::CellId, padding: f64) -> Self::PaddedCell;
      fn cell_to_rect(&self, id: Self::CellId) -> Self::Rect;

      // === Existing methods ===
      fn origin(&self) -> Self::Point;
      fn edge_bound(&self, v0: &Self::Point, v1: &Self::Point) -> Self::Rect;
      fn compute_max_level(&self, v0: &Self::Point, v1: &Self::Point) -> u8;

      fn clip_edge_to_children(
          &self,
          edge_bounds: &Self::Rect,
          v0: &Self::Point,
          v1: &Self::Point,
          cell: &Self::PaddedCell,
      ) -> [[Option<Self::Rect>; 2]; 2];

      fn new_crosser(a: &Self::Point, b: &Self::Point) -> Self::Crosser;
      fn test_crossing(crosser: &mut Self::Crosser, v0: &Self::Point, v1: &Self::Point) -> bool;

      fn cell_begin(level: u8) -> Self::CellId;
      fn cell_range_min(cell: Self::CellId) -> Self::CellId;
      fn cell_next(cell: Self::CellId) -> Self::CellId;

      fn brute_force_contains(&self, point: &Self::Point, vertices: &[Self::Point]) -> bool;
      fn bounds_intersect(a: &Self::Rect, b: &Self::Rect) -> bool;
  fn clip_bound_to_children(
      &self,
      bound: &Self::Rect,
      middle: &Self::Rect,
      v0: &Self::Point,
      v1: &Self::Point,
  ) -> [[Option<Self::Rect>; 2]; 2];
  }
