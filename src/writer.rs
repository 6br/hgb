//! GHB writer.

use std::io::{Write, Result};
use std::fs::File;
use std::path::Path;

use bam::bgzip;
use bam::header::Header;
use super::IndexWriter;
//use crate::range::InvertedRecord;
//use crate::ColumnarSet;
use crate::index::Index;

/// Builder of the [GhbWriter](struct.GhbWriter.html).
pub struct GhbWriterBuilder {
    write_header: bool,
    level: u8,
    additional_threads: u16,
}


impl GhbWriterBuilder {
    pub fn new() -> Self {
        Self {
            write_header: true,
            level: 6,
            additional_threads: 0,
        }
    }

    /// The option to write or skip header when creating the BAM writer (writing by default).
    pub fn write_header(&mut self, write: bool) -> &mut Self {
        self.write_header = write;
        self
    }

    /// Specify compression level from 0 to 9 (6 by default).
    pub fn compression_level(&mut self, level: u8) -> &mut Self {
        assert!(level <= 9, "Maximal compression level is 9");
        self.level = level;
        self
    }

    /// Specify the number of additional threads.
    /// Additional threads are used to compress blocks, while the
    /// main thread reads the writes to a file/stream.
    /// If `additional_threads` is 0 (default), the main thread
    /// will compress blocks itself.
    pub fn additional_threads(&mut self, additional_threads: u16) -> &mut Self {
        self.additional_threads = additional_threads;
        self
    }

    /// Creates a BAM writer from a file and a header.
    pub fn from_path<P: AsRef<Path>>(&mut self, path: P, header: Header)
            -> Result<GhbWriter<File>> {
        let stream = File::create(path)?;
        self.from_stream(stream, header)
    }

    /// Creates a BAM writer from a stream and a header.
    pub fn from_stream<W: Write>(&mut self, stream: W, header: Header) -> Result<GhbWriter<W>> {
        let mut writer = bgzip::Writer::build()
            .additional_threads(self.additional_threads)
            .compression_level(self.level)
            .from_stream(stream);
        if self.write_header {
            header.write_bam(&mut writer)?;
        }
        writer.flush_contents()?;
        Ok(GhbWriter { writer, header })
    }
}

/// Bam writer. Can be created using [from_path](#method.from_path) or using
/// [GhbWriterBuilder](struct.GhbWriterBuilder.html).
///
/// Use [RecordWriter](../trait.RecordWriter.html) trait to write records.
pub struct GhbWriter<W: Write> {
    writer: bgzip::Writer<W>,
    header: Header,
}

impl GhbWriter<File> {
    /// Creates a [GhbWriterBuilder](struct.GhbWriterBuilder.html).
    pub fn build() -> GhbWriterBuilder {
        GhbWriterBuilder::new()
    }

    /// Creates a new `GhbWriter` from a path and header.
    pub fn from_path<P: AsRef<Path>>(path: P, header: Header) -> Result<Self> {
        GhbWriterBuilder::new().from_path(path, header)
    }
}

impl<W: Write> GhbWriter<W> {
    /// Creates a new `GhbWriter` from a stream and header.
    pub fn from_stream(stream: W, header: Header) -> Result<Self> {
        GhbWriterBuilder::new().from_stream(stream, header)
    }

    /// Returns BAM header.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Consumes the writer and returns inner stream.
    pub fn take_stream(self) -> W {
        self.writer.take_stream()
    }


    /// Pauses multi-thread writer until the next write operation. Does nothing to a single-thread writer.
    ///
    /// Use with caution: pausing and unpausing takes some time. Additionally, blocks that are compressed
    /// at the moment will finish compressing, but will not be written.
    /// All other blocks in the queue will not be compressed nor written.
    ///
    /// To compress and write all remaining blocks you can call [flush](#method.flush) before calling `pause`.
    pub fn pause(&mut self) {
        self.writer.pause();
    }
}

impl<W: Write> IndexWriter for GhbWriter<W> {
    fn write(&mut self, record: &Index) -> Result<()> {
        record.to_stream(&mut self.writer)?;
        self.writer.end_context();
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.writer.finish()
    }

    fn flush(&mut self) -> Result<()> {
        self.writer.flush()
    }
}