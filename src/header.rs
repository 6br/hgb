/* Making it possible to store any types of header including bam */

use bam::header

#[derive(Clone, PartialEq, Eq, Debug)]
enum HeaderType {
    BAM(header::Header)
}

pub struct Header {
    global_header : header::Header // Need to be replaced.
    headers: Vec<HeaderType>
}

impl Header {
    pub fn new(global_header: header::Header) -> Self {
        Header{global_header: Header::new(), headers:vec![]}
    }

    pub fn set_entry(&mut self, entries: Vec<HeaderEntry>) {
        for i in entries {
            self.global_header.push_entry(i).unwrap();
        }
    }

    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        global_header.write_bam(&mut writer)?
        let n_samples = self.headess.len() as i32;
        stream.write_i32::<LittleEndian>(n_samples)?;
        for i in headers {
            i.to_steam(&mut_writer)?
        }
    }

    pub fn from_stream<R: Read>(&mut self, stream: &mut R) -> Result<bool> {
        let global_header = Header::from_bam(&mut index_reader).map_err(|e| Error::new(e.kind(), format!("Failed to read BAI header: {}", e)))?;
        let n_samples = stream.read_i32::<LittleEndian>()? as usize;
        let mut chunks = Vec::with_capacity(n_chunks);
        for i in 0..n_samples {
            chunks.push(HeaderType::from_stream(stream))?
        }
        self.
    }
}