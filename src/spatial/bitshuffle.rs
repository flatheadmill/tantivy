//! Rust implementation of Kiyoshi Masui's `bitshuffle`. `bitsuffle` is a shuffle algorithm that
//! rearranges typed data for more efficent compression. It takes an array of typed structures and
//! transposes the bits so that the bits for each structure run consecutively by the bit position
//! in the structure. This organization lends itself to LZ4 and ZSTD compression.

/// Transpose an 8x8 bit matrix packed in a u64.
#[cfg(target_endian = "little")]
fn transpose_8x8(mut x: u64) -> u64 {
    let mut t: u64;

    // Swap 1x1 blocks (adjacent bits in alternating positions)
    t = (x ^ (x >> 7)) & 0x00AA00AA00AA00AA;
    x = x ^ t ^ (t << 7);

    // Swap 2x2 blocks
    t = (x ^ (x >> 14)) & 0x0000CCCC0000CCCC;
    x = x ^ t ^ (t << 14);

    // Swap 4x4 blocks
    t = (x ^ (x >> 28)) & 0x00000000F0F0F0F0;
    x = x ^ t ^ (t << 28);

    x
}

#[cfg(target_endian = "big")]
fn transpose_8x8(mut x: u64) -> u64 {
    let mut t: u64;

    // Swap 1x1 blocks
    t = (x ^ (x >> 9)) & 0x0055005500550055;
    x = x ^ t ^ (t << 9);

    // Swap 2x2 blocks
    t = (x ^ (x >> 18)) & 0x0000333300003333;
    x = x ^ t ^ (t << 18);

    // Swap 4x4 blocks
    t = (x ^ (x >> 36)) & 0x000000000F0F0F0F;
    x = x ^ t ^ (t << 36);

    x
}

/// Transpose the bytes of the elements in the array such that the first byte of each element is
/// consecutive, followed by the second bytes of each element and so on.
fn transpose_bytes_within_elements(
    input: &[u8],
    output: &mut [u8],
    element_count: usize,
    element_size: usize,
) {
    for i in 0..element_count {
        for j in 0..element_size {
            output[j * element_count + i] = input[i * element_size + j];
        }
    }
}

/// Transpose the bits in an array of bytes a u64 integer at a time. Take the 8 bytes from the
/// resulting transposed u64 integers and group them by significance. We end up with 8 columns of
/// sizeof(element) / 8.
fn transpose_bits_within_bytes(input: &[u8], output: &mut [u8]) {
    let groups = input.len() / 8;
    for i in 0..groups {
        let x = u64::from_ne_bytes(input[i * 8..i * 8 + 8].try_into().unwrap());
        let bytes = transpose_8x8(x).to_ne_bytes();
        for j in 0..8 {
            output[j * groups + i] = bytes[j];
        }
    }
}

/// We now have 8 columns of sizeof(element) / 8. Because we first ordered the array by columns of
/// bytes where each column contains bytes from the same offset into each element, we know that our
/// 8 columns of grouped bits are going to be subsequently ordered by element byte offset, so that
/// if we have a 16 byte element, the first column of leading bits is going to contain values for
/// the bits 0 and 64 alternating.
fn transpose_bitrow_eight(
    input: &[u8],
    output: &mut [u8],
    element_count: usize,
    element_size: usize,
) {
    let bitrow_size = element_count / 8;
    for i in 0..8 {
        for j in 0..element_size {
            let from = (i * element_size + j) * bitrow_size;
            let to = (j * 8 + i) * bitrow_size;
            output[to..to + bitrow_size].copy_from_slice(&input[from..from + bitrow_size]);
        }
    }
}

/// HUSH
pub fn bitshuffle(input: &[u8], element_count: usize, element_size: usize) -> Vec<u8> {
    let size = element_count * element_size;

    assert!(
        element_count % 8 == 0,
        "element_count must be divisible by 8"
    );

    let mut transposed_bytes = vec![0u8; size];
    transpose_bytes_within_elements(input, &mut transposed_bytes, element_count, element_size);

    let mut transposed_bits = vec![0u8; size];
    transpose_bits_within_bytes(&transposed_bytes, &mut transposed_bits);

    let mut output = vec![0u8; size];
    transpose_bitrow_eight(&transposed_bits, &mut output, element_count, element_size);

    output
}

fn untranspose_bitrow_eight(
    input: &[u8],
    output: &mut [u8],
    element_count: usize,
    element_size: usize,
) {
    let bitrow_size = element_count / 8;
    for i in 0..element_size {
        for j in 0..8 {
            let from = (i * 8 + j) * bitrow_size;
            let to = (j * element_size + i) * bitrow_size;
            output[to..to + bitrow_size].copy_from_slice(&input[from..from + bitrow_size]);
        }
    }
}

fn untranspose_bits_within_bytes(input: &[u8], output: &mut [u8]) {
    let groups = input.len() / 8;

    for i in 0..groups {
        let mut bytes = [0u8; 8];
        for j in 0..8 {
            bytes[j] = input[j * groups + i];
        }
        let transposed = transpose_8x8(u64::from_ne_bytes(bytes));
        output[i * 8..i * 8 + 8].copy_from_slice(&transposed.to_ne_bytes());
    }
}

fn untranspose_bytes_within_elements(
    input: &[u8],
    output: &mut [u8],
    element_count: usize,
    element_size: usize,
) {
    for i in 0..element_count {
        for j in 0..element_size {
            output[i * element_size + j] = input[j * element_count + i];
        }
    }
}

/// HUSH
pub fn bitunshuffle(input: &[u8], element_count: usize, element_size: usize) -> Vec<u8> {
    let size = element_count * element_size;

    assert!(
        element_count % 8 == 0,
        "element_count must be divisible by 8"
    );

    let mut untransposed_bitrows = vec![0u8; size];
    untranspose_bitrow_eight(
        input,
        &mut untransposed_bitrows,
        element_count,
        element_size,
    );

    let mut untransposed_bits = vec![0u8; size];
    untranspose_bits_within_bytes(&untransposed_bitrows, &mut untransposed_bits);

    let mut output = vec![0u8; size];
    untranspose_bytes_within_elements(&untransposed_bits, &mut output, element_count, element_size);

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_bitshuffle_roundtrip() {
        let original: Vec<u8> = (0..64).collect(); // 16 elements of 4 bytes each
        let shuffled = bitshuffle(&original, 16, 4);
        let unshuffled = bitunshuffle(&shuffled, 16, 4);
        assert_eq!(original, unshuffled);
    }
}
