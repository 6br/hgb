/* Making it possible to store any types of header including bam */

use bam::header;
use std::io::{Result, Read, Write, Error, ErrorKind};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

/// Local headers are categorized as types. 
#[derive(Clone)]
pub enum HeaderType {
    None,
    BAM(header::Header)
}

impl HeaderType {
    /// Save to stream.
    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        match self {
            HeaderType::None => stream.write_i32::<LittleEndian>(0 as i32),
            HeaderType::BAM(header) => {
                stream.write_i32::<LittleEndian>(1 as i32)?;
                header.write_bam(stream)
            }
        }
    }
    /// Load from stream.
    pub fn from_stream<R: Read>(stream: &mut R) -> Result<HeaderType> {
        let header_type = stream.read_i32::<LittleEndian>()?;
        match header_type {
            0 => {Ok(HeaderType::None)}
            1 => {
                let header = header::Header::from_bam(stream).map_err(|e| Error::new(e.kind(), format!("Failed to read local BAI header: {}", e)))?;
                Ok(HeaderType::BAM(header))
            }
            _ => Err(Error::new(ErrorKind::InvalidData, "Invalid header type id"))
        }
    }
}

/// GHB/GHI Header.
/// 
/// You can modify it by pushing new entry using [push_entry](#method.push_entry),
/// You cannot remove lines.
#[derive(Clone)]
pub struct Header {
    global_header : header::Header, // Need to be replaced.
    headers: Vec<HeaderType>
}

impl Header {
    /// Create header from stream
    pub fn new_from_stream<R: Read>(stream: &mut R) -> Result<Self> {
        let mut header = Header::new();
        header.from_stream(stream)?;
        Ok(header)
    }
    /// Returns the length of the reference with `ref_id` (0-based). as u64.
    /// TODO() Support u64; now raise errors when you pass u32.
    /// Returns None if there is no such reference
    pub fn reference_len(&self, id: u64) -> Option<u64> {
        self.global_header.reference_len(id as u32).map(|e| e as u64)
    }
    /// Returns reference id from its name, if possible.
    pub fn reference_id(&self, ref_name: &str) -> Option<u64> {
        self.global_header.reference_id(ref_name).map(|e| e as u64)
    }
    /// Returns the name of the reference with `ref_id` (0-based).
    /// TODO() Support u64; now raise errors when you pass u32.
    /// Returns None if there is no such reference
    pub fn reference_name(&self, ref_id: u64) -> Option<&str> {
        self.global_header.reference_name(ref_id as u32)
    }
    /// Returns the number of reference sequences in the BAM file.
    pub fn n_references(&self) -> usize {
        self.global_header.n_references()
    }
    pub fn push_entry(&mut self, header_entry: header::HeaderEntry) -> std::result::Result<(), String> {
        self.global_header.push_entry(header_entry)
    }
    pub fn new() -> Self {
        Header{global_header: header::Header::new(), headers:vec![]}
    }

    pub fn set_entry(&mut self, entries: Vec<header::HeaderEntry>) {
        for i in entries {
            self.global_header.push_entry(i).unwrap();
        }
    }

    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        self.global_header.write_bam(stream)?;
        let n_samples = self.headers.len() as i32;
        stream.write_i32::<LittleEndian>(n_samples)?;
        for i in &self.headers {
            i.to_stream(stream)?;
        }
        Ok(())
    }

    pub fn from_stream<R: Read>(&mut self, stream: &mut R) -> Result<bool> {
        let global_header = header::Header::from_bam(stream).map_err(|e| Error::new(e.kind(), format!("Failed to read BAI header: {}", e)))?;
        let n_samples = stream.read_i32::<LittleEndian>()? as usize;
        let mut headers = Vec::with_capacity(n_samples);
        for _i in 0..n_samples {
            headers.push(HeaderType::from_stream(stream)?);
        }
        self.global_header = global_header;
        self.headers = headers;
        Ok(true)
    }
}