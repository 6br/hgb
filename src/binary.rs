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

use super::{ChunkWriter, InvertedRecord, InvertedRecordWriter};
use crate::header::Header;
use crate::index::{Chunk, VirtualOffset};
use crate::range::Record;
use crate::ColumnarSet;
use std::cmp::max;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Result, Seek, SeekFrom, Write};
use std::path::Path;

/// Builder of the [GhbWriter](struct.GhbWriter.html).
pub struct GhbWriterBuilder {
    write_header: bool,
}

impl GhbWriterBuilder {
    pub fn new() -> Self {
        Self { write_header: true }
    }

    /// The option to write or skip header when creating the GHB writer (writing by default).
    pub fn write_header(&mut self, write: bool) -> &mut Self {
        self.write_header = write;
        self
    }

    /// Creates a GHB writer from a file and a header.
    pub fn from_path<P: AsRef<Path>>(
        &mut self,
        path: P,
        header: Header,
    ) -> Result<GhbWriter<BufWriter<File>>> {
        let stream = BufWriter::new(File::create(path)?);
        self.from_stream(stream, header)
    }

    /// Creates a GHB writer from a stream and a header. Preferably the stream should be wrapped
    /// in a buffer writer, such as `BufWriter`.
    pub fn from_stream<W: Write>(&mut self, mut stream: W, header: Header) -> Result<GhbWriter<W>> {
        if self.write_header {
            header.to_stream(&mut stream)?;
        }
        stream.flush()?;
        Ok(GhbWriter { stream, header })
    }
}

/// Writes records in GHB format.
///
/// Can be created as
/// ```rust
/// let writer = GhbWriter::from_path("out.ghb", header).unwrap();
/// ```
/// or using a [builder](struct.GhbWriterBuilder.html)
/// ```rust
/// let writer = GhbWriter::build()
///     .write_header(false)
///     .from_path("out.ghb", header).unwrap();
/// ```
///
/// You can clone a [header](../header/struct.Header.html) from GHB/GHI reader or
/// create one yourself.
pub struct GhbWriter<W: Write> {
    stream: W,
    header: Header,
}

impl GhbWriter<BufWriter<File>> {
    /// Create a [builder](struct.GhbWriterBuilder.html).
    pub fn build() -> GhbWriterBuilder {
        GhbWriterBuilder::new()
    }

    /// Creates a GHB writer from a path and a header.
    pub fn from_path<P: AsRef<Path>>(path: P, header: Header) -> Result<Self> {
        GhbWriterBuilder::new().from_path(path, header)
    }
}

impl<W: Write> GhbWriter<W> {
    /// Creates a GHB writer from a stream and a header. Preferably the stream should be wrapped
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

impl<W: Write + Seek> ChunkWriter for GhbWriter<W> {
    /// Writes a single record in GHB format.
    fn write(&mut self, record: &Record) -> Result<Chunk> {
        let start = self
            .stream
            .seek(SeekFrom::Current(0))
            .map(|a| VirtualOffset::from_raw(a))?;
        record.to_stream(&mut self.stream)?;
        let stop = self
            .stream
            .seek(SeekFrom::Current(0))
            .map(|a| VirtualOffset::from_raw(a))?;
        Ok(Chunk::new(record.sample_id(), 0, start, stop))
    }

    fn finish(&mut self) -> Result<()> {
        self.flush()
    }

    fn flush(&mut self) -> Result<()> {
        self.flush()
    }
}

impl<W: Write + Seek> InvertedRecordWriter for GhbWriter<W> {
    /// Writes a single record in GHB format.
    fn write(&mut self, record: &InvertedRecord) -> Result<VirtualOffset> {
        record.to_stream(&mut self.stream)?;
        self.stream
            .seek(SeekFrom::Current(0))
            .map(|a| VirtualOffset::from_raw(a))
    }

    fn finish(&mut self) -> Result<()> {
        self.flush()
    }

    fn flush(&mut self) -> Result<()> {
        self.flush()
    }
}

/// Reads records from GHB format.
///
/// Can be opened as
/// ```rust
/// let reader = GhbReader::from_path("in.ghb").unwrap();
/// ```
///
/// You can iterate over records:
/// ```rust
/// let mut reader = binary::GhbReader::from_path("in.ghb").unwrap();
/// for record in reader {
///     let record = record.unwrap();
///     // Do something.
/// }
/// ```
/// You can use [RecordReader](../trait.RecordReader.html) trait to read records without excess
/// allocation.
/// See: https://stackoverflow.com/questions/54791718/whats-the-difference-between-a-traits-generic-type-and-a-generic-associated-ty/54792178
pub struct GhbReader<R: Read + Seek> {
    // _marker: std::marker::PhantomData<fn() -> T>, // https://in-neuro.hatenablog.com/entry/2019/01/22/220639
    stream: R,
    header: Header,
    chunks: Vec<Chunk>,
    index: usize,
    started: bool,
    offset: u64,
}

impl GhbReader<BufReader<File>> {
    /// Opens GHB reader from `path`.
    pub fn from_path<P: AsRef<Path>>(path: P, header: Header) -> Result<Self> {
        let stream = BufReader::new(File::open(path)?);
        GhbReader::from_stream(stream, header)
    }
}

impl<R: Read + Seek> GhbReader<R> {
    /// Opens GHB reader from a buffered stream.
    pub fn from_stream(mut stream: R, header: Header) -> Result<Self> {
        // let header = Header::from_bam(&mut stream)?;
        // let header = Header::new();
        // let _marker = std::marker::PhantomData;
        let chunks = Vec::new();
        let offset = stream.seek(SeekFrom::Current(0))?;
        Ok(GhbReader {
            stream,
            header,
            chunks,
            index: 0 as usize,
            started: false,
            offset,
        })
    }

    /// Returns [header](../header/struct.Header.html).
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Consumes the reader and returns inner stream.
    pub fn take_stream(self) -> R {
        self.stream
    }

    pub fn set_header(mut self, header: Header) {
        self.header = header;
    }

    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }

    pub fn fill_from_binary(&mut self, record: &mut Record) -> Result<bool> {
        if let Some(new_offset) = self.next_offset() {
            if new_offset != self.offset {
                self.stream.seek(SeekFrom::Start(new_offset))?;
                self.offset = new_offset;
            }
            record.fill_from_bam(&mut self.stream)
        // self.offset = TODO(offset should be updated, but how?)
        } else {
            Ok(false)
            //Err(Error::new(ErrorKind::Other, "The end of stream"))
        }
    }

    fn next_offset(&mut self) -> Option<u64> {
        if self.index >= self.chunks.len() {
            return None;
        }
        if !self.started {
            self.started = true;
            return Some(self.chunks[0].start().raw());
        }

        let curr_offset = VirtualOffset::new(self.offset, 0);
        // When not full query:
        if curr_offset > self.chunks[self.index].end()
            || VirtualOffset::MAX != self.chunks[self.index].end()
        {
            self.index += 1;
        }
        //}
        if self.index >= self.chunks.len() {
            None
        } else {
            Some(max(self.offset, self.chunks[self.index].start().raw()))
        }
    }

    pub fn set_chunks<I: IntoIterator<Item = Chunk>>(&mut self, chunks: I) {
        self.chunks.clear();
        self.chunks.extend(chunks);
        for i in 1..self.chunks.len() {
            if self.chunks[i - 1].intersect(&self.chunks[i]) {
                panic!(
                    "Cannot set chunks: chunk {:?} intersects chunk {:?}",
                    self.chunks[i - 1],
                    self.chunks[i]
                );
            } else if self.chunks[i - 1] >= self.chunks[i] {
                panic!(
                    "Cannot set chunks: chunks are unordered: {:?} >= {:?}",
                    self.chunks[i - 1],
                    self.chunks[i]
                );
            }
        }
        self.index = 0;
        self.started = false;
    }
}
