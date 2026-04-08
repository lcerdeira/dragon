/// Variable-length integer encoding (LEB128 / varint).
/// Used for compressing genome path indices and delta-encoded sequences.

/// Encode a u64 as a varint into a byte buffer. Returns number of bytes written.
pub fn encode_varint(mut value: u64, buf: &mut Vec<u8>) -> usize {
    let start = buf.len();
    loop {
        let mut byte = (value & 0x7F) as u8;
        value >>= 7;
        if value > 0 {
            byte |= 0x80;
        }
        buf.push(byte);
        if value == 0 {
            break;
        }
    }
    buf.len() - start
}

/// Decode a varint from a byte slice. Returns (value, bytes_consumed).
pub fn decode_varint(buf: &[u8]) -> (u64, usize) {
    let mut value: u64 = 0;
    let mut shift = 0;
    for (i, &byte) in buf.iter().enumerate() {
        value |= ((byte & 0x7F) as u64) << shift;
        if byte & 0x80 == 0 {
            return (value, i + 1);
        }
        shift += 7;
        if shift >= 64 {
            break;
        }
    }
    (value, buf.len())
}

/// Encode a signed i64 using zigzag encoding + varint.
pub fn encode_zigzag(value: i64, buf: &mut Vec<u8>) -> usize {
    let encoded = ((value << 1) ^ (value >> 63)) as u64;
    encode_varint(encoded, buf)
}

/// Decode a zigzag-encoded varint to i64.
pub fn decode_zigzag(buf: &[u8]) -> (i64, usize) {
    let (encoded, consumed) = decode_varint(buf);
    let value = ((encoded >> 1) as i64) ^ -((encoded & 1) as i64);
    (value, consumed)
}

/// Encode a slice of u64 values as varints.
pub fn encode_varint_slice(values: &[u64]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(values.len() * 2);
    for &v in values {
        encode_varint(v, &mut buf);
    }
    buf
}

/// Decode a varint-encoded buffer into a Vec<u64>.
pub fn decode_varint_slice(buf: &[u8], count: usize) -> Vec<u64> {
    let mut values = Vec::with_capacity(count);
    let mut offset = 0;
    for _ in 0..count {
        if offset >= buf.len() {
            break;
        }
        let (val, consumed) = decode_varint(&buf[offset..]);
        values.push(val);
        offset += consumed;
    }
    values
}

/// Delta-encode a sorted slice of u64, then varint-encode the deltas.
pub fn delta_encode(sorted_values: &[u64]) -> Vec<u8> {
    if sorted_values.is_empty() {
        return Vec::new();
    }
    let mut buf = Vec::with_capacity(sorted_values.len() * 2);
    encode_varint(sorted_values[0], &mut buf);
    for i in 1..sorted_values.len() {
        let delta = sorted_values[i] - sorted_values[i - 1];
        encode_varint(delta, &mut buf);
    }
    buf
}

/// Decode delta + varint encoded values.
pub fn delta_decode(buf: &[u8], count: usize) -> Vec<u64> {
    let deltas = decode_varint_slice(buf, count);
    if deltas.is_empty() {
        return Vec::new();
    }
    let mut values = Vec::with_capacity(deltas.len());
    values.push(deltas[0]);
    for i in 1..deltas.len() {
        values.push(values[i - 1] + deltas[i]);
    }
    values
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_varint_small() {
        let mut buf = Vec::new();
        encode_varint(42, &mut buf);
        let (val, consumed) = decode_varint(&buf);
        assert_eq!(val, 42);
        assert_eq!(consumed, 1);
    }

    #[test]
    fn test_varint_large() {
        let mut buf = Vec::new();
        let large = 1_000_000_000u64;
        encode_varint(large, &mut buf);
        let (val, consumed) = decode_varint(&buf);
        assert_eq!(val, large);
        assert!(consumed > 1);
    }

    #[test]
    fn test_zigzag() {
        for v in [-1i64, 0, 1, -100, 100, i64::MIN, i64::MAX] {
            let mut buf = Vec::new();
            encode_zigzag(v, &mut buf);
            let (decoded, _) = decode_zigzag(&buf);
            assert_eq!(decoded, v);
        }
    }

    #[test]
    fn test_delta_encode_decode() {
        let values = vec![10, 20, 25, 100, 200];
        let encoded = delta_encode(&values);
        let decoded = delta_decode(&encoded, values.len());
        assert_eq!(decoded, values);
    }

    #[test]
    fn test_varint_slice() {
        let values = vec![1, 128, 16384, 0, u64::MAX];
        let encoded = encode_varint_slice(&values);
        let decoded = decode_varint_slice(&encoded, values.len());
        assert_eq!(decoded, values);
    }
}
