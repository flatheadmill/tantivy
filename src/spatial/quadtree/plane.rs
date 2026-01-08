use crate::spatial::quadtree::{brute_force_contains_2d, Point2D};

pub fn contains_tracker_origin(origin: &Point2D, vertices: &[Point2D]) -> bool {
    if vertices.len() < 3 {
        return false;
    }

    // Use brute force contains check
    brute_force_contains_2d(origin, vertices)
}
