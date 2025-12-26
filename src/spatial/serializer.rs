//! HUSH
use std::io;
use std::io::Write;
use std::marker::PhantomData;

use crate::directory::{CompositeWrite, WritePtr};
use crate::schema::Field;
use crate::spatial::bvh::write_tree;
use crate::spatial::envelope::{Envelope, Spatial};

/// The fieldnorms serializer is in charge of
/// the serialization of field norms for all fields.
pub struct SpatialSerializer<S: Spatial> {
    composite_write: CompositeWrite,
    _marker: PhantomData<S>,
}

impl<S: Spatial> SpatialSerializer<S> {
    /// Create a composite file from the write pointer.
    pub fn from_write(write: WritePtr) -> io::Result<SpatialSerializer<S>> {
        // just making room for the pointer to header.
        let composite_write = CompositeWrite::wrap(write);
        Ok(SpatialSerializer::<S> {
            composite_write,
            _marker: PhantomData,
        })
    }

    /// Serialize the given field
    pub fn serialize_field(
        &mut self,
        field: Field,
        envelopes: &mut [Envelope<S::Bounds>],
    ) -> io::Result<()> {
        if envelopes.is_empty() {
            return Ok(());
        }
        let write = self.composite_write.for_field(field);
        write_tree::<S>(envelopes, write)?;
        write.flush()?;
        Ok(())
    }

    /// Clean up, flush, and close.
    pub fn close(self) -> io::Result<()> {
        self.composite_write.close()?;
        Ok(())
    }
}
