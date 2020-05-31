/* Making it possible to store any types of header including bam */

use bam::header;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use csv::Reader;
use header::HeaderLine;
use std::io::{Error, ErrorKind, Read, Result, Write};
/// Local headers are categorized as types.
#[derive(Clone)]
pub enum HeaderType {
    None,
    BAM(header::Header),
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
            0 => Ok(HeaderType::None),
            1 => {
                let header = header::Header::from_bam(stream).map_err(|e| {
                    Error::new(e.kind(), format!("Failed to read local BAI header: {}", e))
                })?;
                Ok(HeaderType::BAM(header))
            }
            _ => Err(Error::new(ErrorKind::InvalidData, "Invalid header type id")),
        }
    }
}

/// GHB/GHI Header.
///
/// You can modify it by pushing new entry using [push_entry](#method.push_entry),
/// You cannot remove lines.
#[derive(Clone)]
pub struct Header {
    global_header: header::Header, // Need to be replaced.
    headers: Vec<HeaderType>,
    names: Vec<String>,
//     bam_readers: BTreeMap<u64, IndexedRe>
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
        self.global_header
            .reference_len(id as u32)
            .map(|e| e as u64)
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
    /// Returns reference names.
    pub fn reference_names(&self) -> &[String] {
        &self.global_header.reference_names()
    }
    /// Pushes a new header entry.
    ///
    /// Returns an error if the same reference appears twice or @SQ line has an incorrect format.
    pub fn push_entry(
        &mut self,
        header_entry: header::HeaderEntry,
    ) -> std::result::Result<(), String> {
        self.global_header.push_entry(header_entry)
    }
    /// Keep a bam Header  on the local header
    pub fn set_local_header(&mut self, header: &bam::Header, name: &str, index: usize) -> () {
        if let None = self.headers.get(index) {
            self.headers.resize(index + 1, HeaderType::None);
        }
        if let None = self.names.get(index) {
            self.names.resize(index + 1, "".to_string());
        }
        self.headers[index] = HeaderType::BAM(header.clone());
        self.names[index] = name.to_string();
    }
    /// Load chromosome from chrom.sizes.
    pub fn set_header_from_sizes<R: Read>(
        &mut self,
        reader: &mut Reader<R>,
    ) -> std::result::Result<(), Box<dyn std::error::Error>> {
        for result in reader.records() {
            let record = result?;
            self.global_header
                .push_entry(header::HeaderEntry::ref_sequence(
                    record[0].to_string(),
                    record[1].parse::<u32>()?,
                ))?;
        }
        Ok(())
    }
    pub fn get_local_header(&self, index: usize) -> Option<&HeaderType> {
        self.headers.get(index)
    }
    pub fn get_name(&self, index: usize) -> Option<&String> {
        self.names.get(index)
    }
    pub fn get_local_bam_header(&self, index: usize) -> Option<&bam::Header> {
        self.headers.get(index).and_then(|f| match &f {
            HeaderType::None => None,
            HeaderType::BAM(a) => Some(a),
        })
    }
    /// Transfer the reference into global header.
    /// If it intends to override existing line, the new line would be ignored.
    pub fn transfer(&mut self, header: &bam::Header) -> () {
        for i in header.lines() {
            match i {
                HeaderLine::Entry(a) => {
                    let _ = self.global_header.push_entry(a.clone());
                }
                HeaderLine::Comment(_) => {}
            }
        }
    }
    /// Creates an empty header.
    pub fn new() -> Self {
        Header {
            global_header: header::Header::new(),
            headers: vec![],
            names: vec![], //BTreeMap::new(),
        }
    }

    /// Pushes a new header entry.
    ///
    /// Do not insert if the same reference appears twice or @SQ line has an incorrect format.
    pub fn set_entry(&mut self, entries: Vec<header::HeaderEntry>) {
        for i in entries {
            self.global_header.push_entry(i).unwrap();
        }
    }

    /// Writes header in an uncompressed binary format.
    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        self.global_header.write_bam(stream)?;
        let n_samples = self.headers.len() as i32;
        stream.write_i32::<LittleEndian>(n_samples)?;
        for i in &self.headers {
            i.to_stream(stream)?;
        }
        for item in &self.names {
            stream.write_u64::<LittleEndian>(item.len() as u64)?;
            stream.write_all(&item.as_bytes())?;
        }
        Ok(())
    }

    /// Parse uncompressed header.
    pub fn from_stream<R: Read>(&mut self, stream: &mut R) -> Result<bool> {
        let global_header = header::Header::from_bam(stream)
            .map_err(|e| Error::new(e.kind(), format!("Failed to read global header: {}", e)))?;
        let n_samples = stream.read_i32::<LittleEndian>()? as usize;
        let mut headers = Vec::with_capacity(n_samples);
        for _i in 0..n_samples {
            headers.push(HeaderType::from_stream(stream)?);
        }
        let mut names = Vec::with_capacity(n_samples);
        for _i in 0..n_samples {
            let len = stream.read_u64::<LittleEndian>()? as usize;
            let mut buf = vec![0u8; len];
            stream.read_exact(&mut buf)?;
            let converted: String = String::from_utf8(buf.to_vec()).unwrap();
            names.push(converted);
        }
        self.global_header = global_header;
        self.headers = headers;
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::Header;
    use std::fs::File;

    #[test]
    fn chrom_sizes() {
        let mut header = Header::new();
        let file = File::open(r#"test/hg38.chrom.sizes"#).unwrap();
        let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);
        let _result = header.set_header_from_sizes(&mut rdr).unwrap();
        assert_eq!(header.reference_names().len(), 454);
        assert_eq!(header.reference_len(0).unwrap(), 242_193_529);
    }
}
