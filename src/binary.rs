//! GHB reader and writer.
//!
//! Contains a [GHB reader](struct.GhbReader.html) and [GHB writer](struct.GhbWriter.html).
//! You can construct them as
//!
//! ```rust
//! let reader = GhbReader::from_path("in.ghb").unwrap();
//! let writer = GhbWriter::from_path("out.ghb", reader.header().clone()).unwrap();
//! ```
//!
//! The reader implements [RecordReader](../trait.RecordReader.html) trait,
//! and the writer implements [RecordWriter](../trait.RecordWriter.html) trait. See them for
//! more information.

use std::io::{Write, BufWriter, Result, Seek, SeekFrom, BufReader, BufRead};
use std::fs::File;
use std::path::Path;
use crate::ColumnarSet;
use bam::header::Header;
use crate::index::{VirtualOffset, Chunk};
use crate::range::Record;
//use crate::InvertedRecordWriter;
use super::{InvertedRecordWriter, ChunkWriter, InvertedRecord};

/// Builder of the [SamWriter](struct.SamWriter.html).
pub struct GhbWriterBuilder {
    write_header: bool,
}

impl GhbWriterBuilder {
    pub fn new() -> Self {
        Self {
            write_header: true,
        }
    }

    /// The option to write or skip header when creating the SAM writer (writing by default).
    pub fn write_header(&mut self, write: bool) -> &mut Self {
        self.write_header = write;
        self
    }

    /// Creates a SAM writer from a file and a header.
    pub fn from_path<P: AsRef<Path>>(&mut self, path: P, header: Header)
            -> Result<GhbWriter<BufWriter<File>>> {
        let stream = BufWriter::new(File::create(path)?);
        self.from_stream(stream, header)
    }

    /// Creates a SAM writer from a stream and a header. Preferably the stream should be wrapped
    /// in a buffer writer, such as `BufWriter`.
    pub fn from_stream<W: Write>(&mut self, mut stream: W, header: Header) -> Result<GhbWriter<W>> {
        if self.write_header {
            header.write_bam(&mut stream)?;
        }
        Ok(GhbWriter { stream, header })
    }
}

/// Writes records in G format.
///
/// Can be created as
/// ```rust
/// let writer = SamWriter::from_path("out.sam", header).unwrap();
/// ```
/// or using a [builder](struct.SamWriterBuilder.html)
/// ```rust
/// let writer = SamWriter::build()
///     .write_header(false)
///     .from_path("out.sam", header).unwrap();
/// ```
///
/// You can clone a [header](../header/struct.Header.html) from SAM/BAM reader or
/// create one yourself.
///
/// You need to import [RecordWriter](../trait.RecordWriter.html)
/// to write [records](../record/struct.Record.html):
/// ```rust
/// use bam::RecordWriter;
/// let mut writer = bam::SamWriter::from_path("out.sam", header).unwrap();
/// let mut record = bam::Record::new();
/// // Filling the record.
/// writer.write(&record).unwrap();
/// ```
pub struct GhbWriter<W: Write> {
    stream: W,
    header: Header,
}

impl GhbWriter<BufWriter<File>> {
    /// Create a [builder](struct.SamWriterBuilder.html).
    pub fn build() -> GhbWriterBuilder {
        GhbWriterBuilder::new()
    }

    /// Creates a SAM writer from a path and a header.
    pub fn from_path<P: AsRef<Path>>(path: P, header: Header) -> Result<Self> {
        GhbWriterBuilder::new().from_path(path, header)
    }
}

impl<W: Write> GhbWriter<W> {
    /// Creates a SAM writer from a stream and a header. Preferably the stream should be wrapped
    /// in a buffer writer, such as `BufWriter`.
    pub fn from_stream(stream: W, header: Header) -> Result<Self> {
        GhbWriterBuilder::new().from_stream(stream, header)
    }

    /// Returns [header](../header/struct.Header.html).
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Flushes contents to output.
    pub fn flush(&mut self) -> Result<()> {
        self.stream.flush()
    }

    /// Consumes the writer and returns inner stream.
    pub fn take_stream(mut self) -> W {
        let _ignore = self.stream.flush();
        self.stream
    }
}

impl<W: Write+Seek, T: ColumnarSet> ChunkWriter<T> for GhbWriter<W> {
    /// Writes a single record in SAM format.
    fn write(&mut self, record: &Record<T>) -> Result<Chunk> {
        let start = self.stream.seek(SeekFrom::Current(0)).map(|a| VirtualOffset::from_raw(a))?;
        record.to_stream(&mut self.stream)?;
        let stop = self.stream.seek(SeekFrom::Current(0)).map(|a| VirtualOffset::from_raw(a))?;
        // Ok(Chunk{start: start, end: stop, sample_id: record.sample_id(), file_id:0})
        Ok(Chunk::new(record.sample_id(), 0, start, stop))
    }

    fn finish(&mut self) -> Result<()> {
        self.flush()
    }

    fn flush(&mut self) -> Result<()> {
        self.flush()
    }
}

impl<W: Write+Seek> InvertedRecordWriter for GhbWriter<W> {
    /// Writes a single record in SAM format.
    fn write(&mut self, record: &InvertedRecord) -> Result<VirtualOffset> {
        record.to_stream(&mut self.stream)?; // , &self.header)
        self.stream.seek(SeekFrom::Current(0)).map(|a| VirtualOffset::from_raw(a))
    }

    fn finish(&mut self) -> Result<()> {
        self.flush()
    }

    fn flush(&mut self) -> Result<()> {
        self.flush()
    }
}

/// Reads records from SAM format.
///
/// Can be opened as
/// ```rust
/// let reader = SamReader::from_path("in.sam").unwrap();
/// ```
///
/// You can iterate over records:
/// ```rust
/// let mut reader = bam::SamReader::from_path("in.sam").unwrap();
/// for record in reader {
///     let record = record.unwrap();
///     // Do something.
/// }
/// ```
/// You can use [RecordReader](../trait.RecordReader.html) trait to read records without excess
/// allocation.
pub struct GhbReader<R: BufRead> {
    stream: R,
    header: Header,
    buffer: String,
}
/*
impl GhbReader<BufReader<File>> {
    /// Opens SAM reader from `path`.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let stream = BufReader::new(File::open(path)?);
        SamReader::from_stream(stream)
    }
}

impl<R: BufRead> GhbReader<R> {
    /// Opens SAM reader from a buffered stream.
    pub fn from_stream(mut stream: R) -> Result<Self> {
        let mut header = Header::new();
        let mut buffer = String::new();
        loop {
            buffer.clear();
            if stream.read_line(&mut buffer)? == 0 {
                break;
            };
            if buffer.starts_with('@') {
                header.push_line(buffer.trim_end())?;
            } else {
                break;
            }
        }
        Ok(SamReader { stream, header, buffer })
    }

    /// Returns [header](../header/struct.Header.html).
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Consumes the reader and returns inner stream.
    pub fn take_stream(self) -> R {
        self.stream
    }
}

impl<R: BufRead> RecordReader for GhbReader<R> {
    fn read_into(&mut self, record: &mut InvertedRecord) -> Result<bool> {
        if self.buffer.is_empty() {
            return Ok(false);
        }
        let res = match record.fill_from_sam(self.buffer.trim(), &self.header) {
            Ok(()) => Ok(true),
            Err(e) => {
                record.clear();
                Err(e)
            },
        };
        self.buffer.clear();
        match self.stream.read_line(&mut self.buffer) {
            Ok(_) => res,
            Err(e) => res.or(Err(e)),
        }
    }

    /// Does nothing, as SAM readers are single-thread.
    fn pause(&mut self) {}
}

/// Iterator over records.
impl<R: BufRead> Iterator for SamReader<R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Record::new();
        match self.read_into(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}
*/