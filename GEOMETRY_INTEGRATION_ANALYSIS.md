# Geometry Type Integration Analysis

## Current Status
We've added `ValueType::Geometry` (variant 13) to store serialized coordinate data as bytes in the document buffer. The Geometry type is intended to store serialized spatial data, NOT raw `MultiPolygon<f64>`.

## 1. Incorrect MultiPolygon<f64> Usage Locations

### ReferenceValueLeaf::Spatial Usage (INCORRECT - should be Geometry with &[u8])
- `/src/core/json_utils.rs:229` - `ReferenceValueLeaf::Spatial(_) => todo!()`
- `/src/fastfield/writer.rs:192` - `ReferenceValueLeaf::Spatial(_) => todo!()`
- `/src/fastfield/writer.rs:324` - `ReferenceValueLeaf::Spatial(_) => todo!()`
- `/src/schema/document/default_document.rs:263` - `ReferenceValueLeaf::Spatial(multi_polygon)` - writes MultiPolygon
- `/src/schema/document/default_document.rs:636` - `ReferenceValueLeaf::Spatial(_) => todo!()`
- `/src/schema/document/owned_value.rs:83` - `OwnedValue::Spatial(val) => ReferenceValueLeaf::Spatial(val.clone())`
- `/src/schema/document/owned_value.rs:293` - `ReferenceValueLeaf::Spatial(_) => todo!()`
- `/src/schema/document/value.rs:98` - `ReferenceValueLeaf::Spatial(ref val)` - returns MultiPolygon
- `/src/schema/document/value.rs:236` - `ReferenceValueLeaf::Spatial(_) => todo!()`
- `/src/schema/document/se.rs:136` - `ReferenceValueLeaf::Spatial(_) => todo!()`

### OwnedValue::Spatial Usage (INCORRECT - should be Geometry with Vec<u8>)
- `/src/schema/document/default_document.rs:129` - `OwnedValue::Spatial(value)`
- `/src/schema/document/owned_value.rs:83` - converts to ReferenceValueLeaf::Spatial
- `/src/schema/document/owned_value.rs:205` - `OwnedValue::Spatial(_) => todo!()`

### Current Enum Definitions (INCORRECT)
```rust
// value.rs:151
ReferenceValueLeaf::Spatial(MultiPolygon<f64>),  // Should be Geometry(&'a [u8])

// owned_value.rs:54
OwnedValue::Spatial(MultiPolygon<f64>),  // Should be Geometry(Vec<u8>)
```

## 2. Object/Array Pattern Analysis

### How Object and Array Work as Complex Leaf Types

**Storage:**
- Object: Stores `Vec<(String, OwnedValue)>` in OwnedValue
- Array: Stores `Vec<OwnedValue>` in OwnedValue
- In the buffer: Both serialize as `&[ValueAddr]` - a list of addresses to their elements

**Access Pattern:**
```rust
// Object iterator returns key-value pairs
ReferenceValue::Object(ObjectIter) -> Iterator<Item = (&'a str, Self)>

// Array iterator returns values
ReferenceValue::Array(ArrayIter) -> Iterator<Item = Self>
```

**Buffer Serialization (default_document.rs):**
```rust
// Array serialization
ReferenceValue::Array(elements) => {
    let mut addresses = Vec::new();
    for elem in elements {
        let value_addr = self.add_value(elem);
        write_into(&mut addresses, value_addr);
    }
    ValueAddr {
        type_id,
        val_addr: write_bytes_into(&mut self.node_data, &addresses),
    }
}

// Object serialization
ReferenceValue::Object(entries) => {
    let mut addresses = Vec::new();
    for (key, value) in entries {
        let key_addr = self.add_value_leaf(ReferenceValueLeaf::Str(key));
        let value_addr = self.add_value(value);
        write_into(&mut addresses, key_addr);
        write_into(&mut addresses, value_addr);
    }
    ValueAddr {
        type_id,
        val_addr: write_bytes_into(&mut self.node_data, &addresses),
    }
}
```

**Key Insight:** Object and Array are NOT leaf types - they're container types in `ReferenceValue`, not `ReferenceValueLeaf`.

## 3. Correct Geometry Integration Design

### What Should Change:

**ReferenceValueLeaf:**
```rust
pub enum ReferenceValueLeaf<'a> {
    // ... other variants ...
    Geometry(&'a [u8]),  // Serialized geometry data, NOT MultiPolygon<f64>
}
```

**OwnedValue:**
```rust
pub enum OwnedValue {
    // ... other variants ...
    Geometry(Vec<u8>),  // Owned serialized geometry data
}
```

**Value Trait Methods:**
```rust
trait Value<'a> {
    // Add new method
    fn as_geometry(&self) -> Option<&'a [u8]> {
        self.as_leaf().and_then(|leaf| leaf.as_geometry())
    }
}
```

**ReferenceValueLeaf Methods:**
```rust
impl<'a> ReferenceValueLeaf<'a> {
    pub fn as_geometry(&self) -> Option<&'a [u8]> {
        if let Self::Geometry(data) = self {
            Some(data)
        } else {
            None
        }
    }
}
```

## 4. Path to spatial_writer

### Current Flow in segment_writer.rs:
```rust
FieldType::Spatial(_) => {
    for value in values {
        if let Some(multi_polygon) = value.as_spatial() {  // Returns &MultiPolygon<f64>
            self.spatial_writer
                .add_multi_polygon(field, doc_id, multi_polygon);
        }
    }
}
```

### New Flow Should Be:
```rust
FieldType::Spatial(_) => {
    for value in values {
        if let Some(geometry_bytes) = value.as_geometry() {  // Returns &[u8]
            // Deserialize bytes to MultiPolygon<f64> here
            let multi_polygon = deserialize_geometry(geometry_bytes)?;
            self.spatial_writer
                .add_multi_polygon(field, doc_id, &multi_polygon);
        }
    }
}
```

The spatial_writer expects `&MultiPolygon<f64>`, so deserialization happens at the point of use, not storage.

## 5. Is "leaf" JSON-specific?

**Answer: No, "leaf" is NOT JSON-specific.**

The term "leaf" refers to terminal/atomic values in the value tree:
- Scalars: u64, i64, f64, bool, str, bytes, etc.
- Complex atomic types: Date, IpAddr, Facet, PreTokStr
- **Geometry should be a leaf** - it's an atomic unit of spatial data

Object and Array are NOT leaves - they're containers that hold other values (which can be leaves or other containers).

## 6. Required Changes Summary

### A. Update Enums:
1. Change `ReferenceValueLeaf::Spatial(MultiPolygon<f64>)` to `ReferenceValueLeaf::Geometry(&'a [u8])`
2. Change `OwnedValue::Spatial(MultiPolygon<f64>)` to `OwnedValue::Geometry(Vec<u8>)`

### B. Update Methods:
1. Rename `as_spatial()` to `as_geometry()` returning `Option<&'a [u8]>`
2. Update `add_spatial()` to serialize MultiPolygon to bytes before storing

### C. Add Missing Implementation in CompactDocValue:
```rust
// In get_ref_value() method, add:
ValueType::Geometry => {
    let data = self.container.extract_bytes(addr);
    Ok(ReferenceValueLeaf::Geometry(data).into())
}
```

### D. Update default_document.rs:
1. Implement `write_multi_polygon()` to serialize MultiPolygon to bytes
2. Update `add_value_leaf()` to handle `ReferenceValueLeaf::Geometry(&[u8])`

### E. Update segment_writer.rs:
Add deserialization step before calling spatial_writer

### F. Update From implementations:
Fix all the `From` trait implementations to handle Geometry properly

### G. Update serialization (se.rs):
Handle `ReferenceValueLeaf::Geometry` serialization

## Key Design Principle

Geometry is stored as **serialized bytes** in the document buffer, just like how strings and bytes are stored. The MultiPolygon<f64> is only the **input format** that gets serialized on write and deserialized on read. This matches the pattern of other complex types like PreTokenizedString which also serialize their internal structure.