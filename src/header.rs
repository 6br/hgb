/* Making it possible to store any types of header including bam */

use bam::header;
use std::io::{Result, Read, Write, Seek, Error, ErrorKind};
use std::io::ErrorKind::InvalidData;
use std::collections::BTreeMap;
// use std::path::Path;
// use std::fs::File;
// use std::marker::PhantomData;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

pub enum HeaderType {
    None,
    BAM(header::Header)
}

impl HeaderType {
    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        match self {
            HeaderType::None => stream.write_i32::<LittleEndian>(0 as i32),
            HeaderType::BAM(header) => {
                stream.write_i32::<LittleEndian>(1 as i32)?;
                header.write_bam(stream)
            }
        }
    }
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

pub struct Header {
    global_header : header::Header, // Need to be replaced.
    headers: Vec<HeaderType>
}

impl Header {
    pub fn new_from_stream<R: Read>(stream: &mut R) -> Result<Self> {
        let mut header = Header::new();
        header.from_stream(stream)?;
        Ok(header)
    }
    pub fn reference_len(&self, id: u32) -> Option<u32> {
        self.global_header.reference_len(id)
    }
    pub fn reference_id(&self, ref_name: &str) -> Option<u32> {
        self.global_header.reference_id(ref_name)
    }
    pub fn reference_name(&self, ref_id: u32) -> Option<&str> {
        self.global_header.reference_name(ref_id)
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
        for i in 0..n_samples {
            headers.push(HeaderType::from_stream(stream)?);
        }
        self.global_header = global_header;
        self.headers = headers;
        Ok(true)
    }
}