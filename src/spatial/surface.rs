pub trait Surface {
      type Point: Copy;
      type CellId: Copy + Eq;
      type Crosser;
      const MAX_LEVEL: u8;

    fn origin(&self) -> Self::Point;
      fn origin_outside_bounds(/* bounds info */) -> Self::Point;
      fn new_crosser(a: &Self::Point, b: &Self::Point) -> Self::Crosser;
      fn test_crossing(crosser: &mut Self::Crosser, v0: &Self::Point, v1: &Self::Point) -> bool;
      fn cell_begin(level: u8) -> Self::CellId;
      fn cell_range_min(cell: Self::CellId) -> Self::CellId;
      fn brute_force_contains(&self, point: &Self::Point, vertices: &[Self::Point]) ->
  bool; }
