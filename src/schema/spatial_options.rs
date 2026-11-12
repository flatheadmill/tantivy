use std::borrow::Cow;
use std::ops::BitOr;

use serde::{Deserialize, Serialize};

use crate::schema::flags::{SchemaFlagList, StoredFlag};

/// Name of a registered spatial index implementation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SpatialIndexName(Cow<'static, str>);

impl Default for SpatialIndexName {
    fn default() -> Self {
        SpatialIndexName::from_static("")
    }
}

impl SpatialIndexName {
    /// Create from a static string.
    pub const fn from_static(name: &'static str) -> Self {
        SpatialIndexName(Cow::Borrowed(name))
    }

    /// Create from a runtime string.
    pub fn from_name(name: &str) -> Self {
        SpatialIndexName(Cow::Owned(name.to_string()))
    }

    /// The name as a string slice.
    pub fn name(&self) -> &str {
        &self.0
    }
}

/// Define how a spatial field should be handled by tantivy.
///
/// The default has no implementation name. Select one before creating an index.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Default)]
pub struct SpatialOptions {
    #[serde(default)]
    stored: bool,
    #[serde(default)]
    spatial_index: SpatialIndexName,
}

/// The field will be indexed on the unit sphere.
pub const SPHERE: SpatialOptions = SpatialOptions {
    stored: false,
    spatial_index: SpatialIndexName(Cow::Borrowed("sphere")),
};

/// The field will be indexed on the Euclidean plane.
pub const PLANE: SpatialOptions = SpatialOptions {
    stored: false,
    spatial_index: SpatialIndexName(Cow::Borrowed("plane")),
};

impl SpatialOptions {
    /// Returns true if the geometry is to be stored.
    #[inline]
    pub fn is_stored(&self) -> bool {
        self.stored
    }

    /// Returns the name of the spatial index implementation.
    pub fn spatial_index_name(&self) -> &str {
        self.spatial_index.name()
    }

    /// Set the spatial index name.
    pub fn set_spatial_index(mut self, name: &str) -> Self {
        self.spatial_index = SpatialIndexName::from_name(name);
        self
    }
}

impl<T: Into<SpatialOptions>> BitOr<T> for SpatialOptions {
    type Output = SpatialOptions;

    fn bitor(self, other: T) -> SpatialOptions {
        let other = other.into();
        SpatialOptions {
            stored: self.stored | other.stored,
            spatial_index: if !other.spatial_index.name().is_empty() {
                other.spatial_index
            } else {
                self.spatial_index
            },
        }
    }
}

impl From<()> for SpatialOptions {
    fn from(_: ()) -> SpatialOptions {
        SpatialOptions::default()
    }
}

impl<Head, Tail> From<SchemaFlagList<Head, Tail>> for SpatialOptions
where
    Head: Clone,
    Tail: Clone,
    Self: BitOr<Output = Self> + From<Head> + From<Tail>,
{
    fn from(head_tail: SchemaFlagList<Head, Tail>) -> Self {
        Self::from(head_tail.head) | Self::from(head_tail.tail)
    }
}

impl From<StoredFlag> for SpatialOptions {
    fn from(_: StoredFlag) -> SpatialOptions {
        SpatialOptions {
            stored: true,
            spatial_index: SpatialIndexName::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schema::STORED;

    #[test]
    fn test_explicit_name_composition() {
        for (options, name, stored) in [
            (SPHERE | PLANE, "plane", false),
            (PLANE | SPHERE, "sphere", false),
            (STORED | SPHERE, "sphere", true),
            (SPHERE | STORED, "sphere", true),
            (STORED | PLANE | SPHERE, "sphere", true),
            (SpatialOptions::from(STORED), "", true),
            (SpatialOptions::default() | PLANE, "plane", false),
            (PLANE | SpatialOptions::default(), "plane", false),
            (STORED | SPHERE.set_spatial_index("named"), "named", true),
            (SPHERE.set_spatial_index("named") | STORED, "named", true),
        ] {
            assert_eq!(options.spatial_index_name(), name);
            assert_eq!(options.is_stored(), stored);
        }
    }

    #[test]
    fn test_serialized_names_are_explicit() -> serde_json::Result<()> {
        for name in ["sphere", "plane", "named"] {
            let value = serde_json::json!({"stored": true, "spatial_index": name});
            let options: SpatialOptions = serde_json::from_value(value.clone())?;
            assert_eq!(options.spatial_index_name(), name);
            assert_eq!(serde_json::to_value(options)?, value);
        }
        let options: SpatialOptions = serde_json::from_str(r#"{"stored":true}"#)?;
        assert!(options.spatial_index_name().is_empty());
        Ok(())
    }
}
