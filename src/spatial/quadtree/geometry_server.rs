use crate::spatial::quadtree::Point2D;
use std::io;

pub struct GeometryServer;

impl GeometryServer {
      pub fn new() -> Self {
          Self
      }

      pub fn get_geometry(&self, _geometry_id: u32) -> io::Result<Vec<Point2D>> {
            Ok(Vec::new())
      }

    pub fn get_doc_id(&self, _geometry_id: u32) -> u32 {
        0
    }
  }
