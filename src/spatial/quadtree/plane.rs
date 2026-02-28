use crate::spatial::{quadtree::{brute_force_contains_2d, Point2D, Rect}, surface::Surface};

 /// Ratio of cell size to edge length that determines when an edge is "long".
  const CELL_SIZE_TO_LONG_EDGE_RATIO: f64 = 1.0;

struct Plane {
      bounds: Rect,
}

impl Surface for Plane {
    type Point = Point2D;
    type CellId = QuadtreeCellId;
    type Rect = Rect;
      const MAX_LEVEL: u8 = 30;
      fn clip_bound_to_children(
          &self,
          bound: &Rect,
          middle: &Rect,
          v0: &Point2D,
          v1: &Point2D,
      ) -> [[Option<Rect>; 2]; 2] {
          let mut result: [[Option<Rect>; 2]; 2] = [[None; 2]; 2];

          if bound.x.hi <= middle.x.lo {
              // Entirely in left children
              self.clip_v_to_children(bound, middle, v0, v1, 0, &mut result);
          } else if bound.x.lo >= middle.x.hi {
              // Entirely in right children
              self.clip_v_to_children(bound, middle, v0, v1, 1, &mut result);
          } else if bound.y.hi <= middle.y.lo {
              // Entirely in bottom children
              result[0][0] = Some(self.clip_u_bound(bound, v0, v1, true, middle.x.hi))
;
              result[1][0] = Some(self.clip_u_bound(bound, v0, v1, false, middle.x.lo)
);
          } else if bound.y.lo >= middle.y.hi {
              // Entirely in top children
              result[0][1] = Some(self.clip_u_bound(bound, v0, v1, true, middle.x.hi))
;
              result[1][1] = Some(self.clip_u_bound(bound, v0, v1, false, middle.x.lo)
);
          } else {
              // Spans multiple quadrants
              let left = self.clip_u_bound(bound, v0, v1, true, middle.x.hi);
              self.clip_v_to_children(&left, middle, v0, v1, 0, &mut result);

              let right = self.clip_u_bound(bound, v0, v1, false, middle.x.lo);
              self.clip_v_to_children(&right, middle, v0, v1, 1, &mut result);
          }

          result
      }

      fn compute_max_level(&self, v0: &Point2D, v1: &Point2D) -> u8 {
          let dx = v1.x - v0.x;
          let dy = v1.y - v0.y;
          let edge_length = (dx * dx + dy * dy).sqrt();

          if edge_length == 0.0 {
              return Self::MAX_LEVEL;
          }

          // The cell edge length at which this edge becomes "long"
          let max_cell_edge = edge_length * CELL_SIZE_TO_LONG_EDGE_RATIO;

          // Find the level where cell edge <= max_cell_edge
          let span = (self.bounds.x.hi - self.bounds.x.lo)
              .max(self.bounds.y.hi - self.bounds.y.lo);

          let mut level = 0u8;
          let mut cell_edge = span;

          while cell_edge > max_cell_edge && level < Self::MAX_LEVEL {
              cell_edge /= 2.0;
              level += 1;
          }

          level
      }
}

pub fn contains_tracker_origin(origin: &Point2D, vertices: &[Point2D]) -> bool {
    if vertices.len() < 3 {
        return false;
    }

    // Use brute force contains check
    brute_force_contains_2d(origin, vertices)
}

 impl Plane {
      /// Clips a bound along the U (X) axis.
      fn clip_u_bound(&self, bound: &Rect, v0: &Point2D, v1: &Point2D, clip_hi: bool, u: f64) -> Rect {
          let mut new_bound = *bound;

          if clip_hi {
              if bound.x.hi <= u { return *bound; }
              new_bound.x.hi = u.min(bound.x.hi);
          } else {
              if bound.x.lo >= u { return *bound; }
              new_bound.x.lo = u.max(bound.x.lo);
          }

          // Interpolate V bound based on where edge crosses U = u
          if (v1.x - v0.x).abs() > f64::EPSILON {
              let t = (u - v0.x) / (v1.x - v0.x);
              if t > 0.0 && t < 1.0 {
                  let v_at_u = v0.y + t * (v1.y - v0.y);
                  let v_clamped = v_at_u.clamp(bound.y.lo, bound.y.hi);
                  let slope_positive = (v1.x > v0.x) == (v1.y > v0.y);

                  if clip_hi == slope_positive {
                      new_bound.y.hi = new_bound.y.hi.min(v_clamped + f64::EPSILON);
                  } else {
                      new_bound.y.lo = new_bound.y.lo.max(v_clamped - f64::EPSILON);
                  }
              }
          }
          new_bound
      }

      /// Clips a bound along the V (Y) axis.
      fn clip_v_bound(&self, bound: &Rect, v0: &Point2D, v1: &Point2D, clip_hi: bool, v: f64) -> Rect {
          let mut new_bound = *bound;

          if clip_hi {
              if bound.y.hi <= v { return *bound; }
              new_bound.y.hi = v.min(bound.y.hi);
          } else {
              if bound.y.lo >= v { return *bound; }
              new_bound.y.lo = v.max(bound.y.lo);
          }

          // Interpolate U bound based on where edge crosses V = v
          if (v1.y - v0.y).abs() > f64::EPSILON {
              let t = (v - v0.y) / (v1.y - v0.y);
              if t > 0.0 && t < 1.0 {
                  let u_at_v = v0.x + t * (v1.x - v0.x);
                  let u_clamped = u_at_v.clamp(bound.x.lo, bound.x.hi);
                  let slope_positive = (v1.x > v0.x) == (v1.y > v0.y);

                  if clip_hi == slope_positive {
                      new_bound.x.hi = new_bound.x.hi.min(u_clamped + f64::EPSILON);
                  } else {
                      new_bound.x.lo = new_bound.x.lo.max(u_clamped - f64::EPSILON);
                  }
              }
          }
          new_bound
      }

      /// Helper: clip along V for a given X index
      fn clip_v_to_children(
          &self, bound: &Rect, middle: &Rect, v0: &Point2D, v1: &Point2D,
          i: usize, result: &mut [[Option<Rect>; 2]; 2],
      ) {
          if bound.y.hi <= middle.y.lo {
              result[i][0] = Some(*bound);
          } else if bound.y.lo >= middle.y.hi {
              result[i][1] = Some(*bound);
          } else {
              result[i][0] = Some(self.clip_v_bound(bound, v0, v1, true, middle.y.hi));
              result[i][1] = Some(self.clip_v_bound(bound, v0, v1, false, middle.y.lo));
          }
      }
  }
