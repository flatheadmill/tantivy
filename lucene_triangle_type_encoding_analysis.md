# Lucene Triangle TYPE Encoding Analysis

## Summary

**The triangle TYPE (POINT/LINE/TRIANGLE) is NOT stored in the BKD tree encoding. It is computed from vertex equality when triangles are decoded.**

## Detailed Analysis

### 1. Triangle Encoding in BKD Tree

When triangles are written to the BKD tree, they are encoded as 7 integers (28 bytes total):

```java
// From ShapeField.encodeTriangle():
// Bytes 0-3:   minY (4 bytes)
// Bytes 4-7:   minX (4 bytes)
// Bytes 8-11:  maxY (4 bytes)
// Bytes 12-15: maxX (4 bytes)
// Bytes 16-19: y (5th coordinate, 4 bytes)
// Bytes 20-23: x (6th coordinate, 4 bytes)
// Bytes 24-27: bits (packed integer, 4 bytes)
```

### 2. Bit Packing Format

The 7th integer (bits) contains:

```java
// Bits 0-2: Reconstruction code (8 possible values, 0-7)
//   - Indicates how to reconstruct the triangle from the 6 stored coordinates
//   - Values: MINY_MINX_MAXY_MAXX_Y_X (0), MINY_MINX_Y_X_MAXY_MAXX (1), etc.
// Bit 3: ab edge flag (1 if edge ab is from original shape)
// Bit 4: bc edge flag (1 if edge bc is from original shape)
// Bit 5: ca edge flag (1 if edge ca is from original shape)
// Bits 6-31: Unused in BKD encoding

bits |= (ab) ? (1 << 3) : 0;
bits |= (bc) ? (1 << 4) : 0;
bits |= (ca) ? (1 << 5) : 0;
```

**NOTE: The TYPE enum is NOT encoded in these bits.**

### 3. TYPE Computation During Decoding

When triangles are read from the BKD tree:

```java
// From ShapeField.decodeTriangle():
public static void decodeTriangle(byte[] t, DecodedTriangle triangle) {
    // 1. Extract coordinates and edge flags from the encoded bytes
    int bits = NumericUtils.sortableBytesToInt(t, 6 * BYTES);
    int tCode = (((1 << 3) - 1) & (bits >> 0)); // Extract bits 0-2

    // 2. Reconstruct the three vertices (aX, aY, bX, bY, cX, cY)
    // based on the reconstruction code

    // 3. Extract edge flags
    ab = (bits & 1 << 3) == 1 << 3;
    bc = (bits & 1 << 4) == 1 << 4;
    ca = (bits & 1 << 5) == 1 << 5;

    // 4. Compute TYPE from vertex equality
    resolveTriangleType(triangle);
}

static void resolveTriangleType(DecodedTriangle triangle) {
    if (triangle.aX == triangle.bX && triangle.aY == triangle.bY) {
        if (triangle.aX == triangle.cX && triangle.aY == triangle.cY) {
            triangle.type = DecodedTriangle.TYPE.POINT;  // All vertices equal
        } else {
            triangle.type = DecodedTriangle.TYPE.LINE;   // a == b, c different
        }
    } else if (triangle.aX == triangle.cX && triangle.aY == triangle.cY) {
        triangle.type = DecodedTriangle.TYPE.LINE;       // a == c, b different
    } else if (triangle.bX == triangle.cX && triangle.bY == triangle.cY) {
        triangle.type = DecodedTriangle.TYPE.LINE;       // b == c, a different
    } else {
        triangle.type = DecodedTriangle.TYPE.TRIANGLE;   // All vertices different
    }
}
```

### 4. TYPE Storage in ShapeDocValues (Different from BKD)

ShapeDocValues (used for doc values, not the BKD index) does explicitly store TYPE:

```java
// ShapeDocValues stores TYPE in a header byte:
if (node.triangle.type == TYPE.POINT) {
    header |= 0x04;  // Bit 2 set for POINT
} else if (node.triangle.type == TYPE.LINE) {
    header |= 0x08;  // Bit 3 set for LINE
}
// No bits set (0x00) means TRIANGLE
```

But this is a separate storage mechanism from the BKD tree encoding.

### 5. Implications for Tantivy Implementation

Since Lucene infers TYPE from vertex equality:

1. **No bits needed for TYPE in the BKD encoding** - TYPE is computed, not stored
2. **Degenerate shapes are detected by vertex comparison**:
   - POINT: All three vertices have identical coordinates
   - LINE: Two vertices are identical (three possible cases)
   - TRIANGLE: All vertices are different
3. **The 7th integer has only 6 bits used** in the BKD encoding:
   - Bits 0-2: Reconstruction code
   - Bits 3-5: Edge flags (ab, bc, ca)
   - Bits 6-31: Available for future use

### 6. Example Encodings

```java
// POINT (all vertices equal)
aX = 100, aY = 200
bX = 100, bY = 200  // Same as a
cX = 100, cY = 200  // Same as a
// Encoded with 6 coordinates, TYPE inferred as POINT during decode

// LINE (two vertices equal)
aX = 100, aY = 200
bX = 150, bY = 250  // Different from a
cX = 100, cY = 200  // Same as a
// Encoded with 6 coordinates, TYPE inferred as LINE during decode

// TRIANGLE (all vertices different)
aX = 100, aY = 200
bX = 150, bY = 250
cX = 120, cY = 180
// Encoded with 6 coordinates, TYPE inferred as TRIANGLE during decode
```

## Conclusion

The triangle TYPE is **inferred from vertex equality** rather than explicitly stored in the BKD tree encoding. This saves bit space and simplifies the encoding, as the TYPE can always be reliably determined from the vertex coordinates themselves.

For Tantivy's implementation, this means:
- No need to reserve bits for TYPE in the packed integer
- TYPE should be computed during decoding by comparing vertices
- The approach is efficient and unambiguous since degenerate cases are well-defined