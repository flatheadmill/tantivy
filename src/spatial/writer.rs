//! HUSH
use std::collections::HashMap;
use std::io;

use crate::schema::Field;
use crate::spatial::envelope::{Bounds, Envelope, Point, Spatial};
use crate::spatial::geometry::Geometry;
use crate::spatial::point::GeoPoint;
use crate::spatial::serializer::SpatialSerializer;
use crate::DocId;

/// HUSH
pub struct SpatialWriter<S: Spatial> {
    /// Map from field to its boudning envelope.
    envelopes_by_field: HashMap<Field, HashMap<DocId, S::Bounds>>,
}

impl<S: Spatial> SpatialWriter<S> {
    /// HUST
    pub fn add_geometry(&mut self, doc_id: DocId, field: Field, geometry: Geometry) {
        let envelopes = self.envelopes_by_field.entry(field).or_default();
        let bounds = envelopes.entry(doc_id).or_insert_with(S::Bounds::empty);
        Self::extend_by_geometry(bounds, geometry);
    }

    fn extend_by_geometry(bounds: &mut S::Bounds, geometry: Geometry) {
        match geometry {
            Geometry::Point(point) => {
                bounds.extend_by_point(S::Point::from_geo(point));
            }
            Geometry::MultiPoint(points) => {
                for point in points {
                    bounds.extend_by_point(S::Point::from_geo(point));
                }
            }
            Geometry::LineString(line) => {
                Self::extend_by_line_string(bounds, &line);
            }
            Geometry::MultiLineString(lines) => {
                for line in lines {
                    Self::extend_by_line_string(bounds, &line);
                }
            }
            Geometry::Polygon(polygon) => {
                for ring in polygon {
                    Self::extend_by_line_string(bounds, &ring);
                }
            }
            Geometry::MultiPolygon(polygons) => {
                for polygon in polygons {
                    for ring in polygon {
                        Self::extend_by_line_string(bounds, &ring);
                    }
                }
            }
            Geometry::GeometryCollection(geometries) => {
                for geometry in geometries {
                    Self::extend_by_geometry(bounds, geometry);
                }
            }
        }
    }

    fn extend_by_line_string(bounds: &mut S::Bounds, line: &[GeoPoint]) {
        let mut iter = line.iter();
        let Some(first) = iter.next() else { return };
        let mut prev = S::Point::from_geo(*first);
        bounds.extend_by_point(prev);
        for geo in iter {
            let point = S::Point::from_geo(*geo);
            bounds.extend_by_line(prev, point);
            prev = point;
        }
    }

    /// Memory usage estimate
    pub fn mem_usage(&self) -> usize {
        self.envelopes_by_field
            .values()
            .map(|envelopes| envelopes.len() * std::mem::size_of::<Envelope<S::Bounds>>())
            .sum()
    }

    /// Serializing our field.
    pub fn serialize(&mut self, mut serializer: SpatialSerializer<S>) -> io::Result<()> {
        for (field, envelope_map) in &mut self.envelopes_by_field {
            let mut envelopes: Vec<Envelope<S::Bounds>> = envelope_map
                .iter()
                .map(|(&doc_id, &bounds)| Envelope::from_bounds(doc_id, bounds))
                .collect();
            serializer.serialize_field(*field, &mut envelopes)?;
        }
        serializer.close()?;
        Ok(())
    }
}

impl<S: Spatial> Default for SpatialWriter<S> {
    /// HUSH
    fn default() -> Self {
        SpatialWriter {
            envelopes_by_field: HashMap::new(),
        }
    }
}
