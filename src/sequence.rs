use crate::ZhuntError;
use std::borrow::Cow;
use std::fs;
use std::io::{BufReader, Read};
use std::path::Path;

const READ_BUFFER_BYTES: usize = 1024 * 1024;
const INVALID_BASE: u8 = u8::MAX;
const NORMALIZED_BASES: [u8; 256] = normalized_base_table();

const fn normalized_base_table() -> [u8; 256] {
    let mut table = [INVALID_BASE; 256];
    table[b'a' as usize] = b'a';
    table[b'A' as usize] = b'a';
    table[b't' as usize] = b't';
    table[b'T' as usize] = b't';
    table[b'g' as usize] = b'g';
    table[b'G' as usize] = b'g';
    table[b'c' as usize] = b'c';
    table[b'C' as usize] = b'c';
    table[b'n' as usize] = b'n';
    table[b'N' as usize] = b'n';
    table
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CircularSequence {
    pub(crate) bases: Vec<u8>,
    wrap_prefix: Vec<u8>,
}

impl CircularSequence {
    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    pub(crate) fn window(&self, start: usize, len: usize) -> Cow<'_, [u8]> {
        if start + len <= self.bases.len() {
            Cow::Borrowed(&self.bases[start..start + len])
        } else {
            let base_tail = &self.bases[start..];
            let wrap_len = len - base_tail.len();
            debug_assert!(wrap_len <= self.wrap_prefix.len());

            let mut window = Vec::with_capacity(len);
            window.extend_from_slice(base_tail);
            window.extend_from_slice(&self.wrap_prefix[..wrap_len]);
            Cow::Owned(window)
        }
    }

    fn from_bases(bases: Vec<u8>, circular_nucleotides: usize) -> Result<Self, ZhuntError> {
        if bases.is_empty() {
            return Err(ZhuntError::EmptySequence);
        }

        let mut wrap_prefix = Vec::with_capacity(circular_nucleotides);
        for index in 0..circular_nucleotides {
            wrap_prefix.push(bases[index % bases.len()]);
        }

        Ok(Self { bases, wrap_prefix })
    }
}

pub fn read_sequence_file<P: AsRef<Path>>(
    path: P,
    circular_nucleotides: usize,
) -> Result<CircularSequence, ZhuntError> {
    let file = fs::File::open(path)?;
    parse_sequence_reader(
        BufReader::with_capacity(READ_BUFFER_BYTES, file),
        circular_nucleotides,
    )
}

pub fn parse_sequence_reader<R: Read>(
    mut reader: R,
    circular_nucleotides: usize,
) -> Result<CircularSequence, ZhuntError> {
    let mut bases = Vec::new();
    let mut buffer = vec![0_u8; READ_BUFFER_BYTES];

    loop {
        let bytes_read = reader.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        for &byte in &buffer[..bytes_read] {
            let base = NORMALIZED_BASES[byte as usize];
            if base != INVALID_BASE {
                bases.push(base);
            }
        }
    }

    CircularSequence::from_bases(bases, circular_nucleotides)
}

pub fn parse_sequence_bytes(
    bytes: &[u8],
    circular_nucleotides: usize,
) -> Result<CircularSequence, ZhuntError> {
    let mut bases = Vec::new();
    for &byte in bytes {
        let base = NORMALIZED_BASES[byte as usize];
        if base != INVALID_BASE {
            bases.push(base);
        }
    }
    CircularSequence::from_bases(bases, circular_nucleotides)
}
