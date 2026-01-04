//! A point in the geographical coordinate system.

/// A point in the geographical coordinate system.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeoPoint {
    /// Longitude
    pub lon: f64,
    /// Latitude
    pub lat: f64,
}

impl From<GeoPoint> for (f64, f64) {
    fn from(p: GeoPoint) -> (f64, f64) {
        (p.lon, p.lat)
    }
}
