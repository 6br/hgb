
use std::path::{Path, PathBuf};
use crate::index::Index;
use bam::bgzip;
use std::io::{Result, Error, Read, Seek};
use std::fs::File;
use std::io::ErrorKind::{InvalidInput};
use bam::header::Header;

/// Defines how to react to a BAI index being younger than BAM file.
///
/// # Variants
/// * `Error` - [IndexedReader](struct.IndexedReader.html) will not be constructed if the BAI
/// index is was modified earlier than the BAM file. `io::Error` will be raised.
/// * `Ignore` - does nothing if the index is younger than the BAM file.
/// * `Warn` - calls a function `Fn(&str)` and continues constructing
/// [IndexedReader](struct.IndexedReader.html);
pub enum ModificationTime {
    Error,
    Ignore,
    Warn(Box<dyn Fn(&str)>),
}

impl ModificationTime {
    fn check<T: AsRef<Path>, U: AsRef<Path>>(&self, bam_path: T, bai_path: U) -> Result<()> {
        let bam_modified = bam_path.as_ref().metadata().and_then(|metadata| metadata.modified());
        let bai_modified = bai_path.as_ref().metadata().and_then(|metadata| metadata.modified());
        let bam_younger = match (bam_modified, bai_modified) {
            (Ok(bam_time), Ok(bai_time)) => bai_time < bam_time,
            _ => false, // Modification time not available.
        };
        if !bam_younger {
            return Ok(());
        }

        match &self {
            ModificationTime::Ignore => {},
            ModificationTime::Error => return Err(Error::new(InvalidInput,
                "the binary file is younger than the index")),
            ModificationTime::Warn(box_fun) =>
                box_fun("the binary file is younger than the index index"),
        }
        Ok(())
    }

    /// Create a warning strategy `ModificationTime::Warn`.
    pub fn warn<F: Fn(&str) + 'static>(warning: F) -> Self {
        ModificationTime::Warn(Box::new(warning))
    }
}

/// [IndexedReader](struct.IndexedReader.html) builder. Allows to specify paths to BAM and BAI
/// files, as well as the number of threads
/// and an option to ignore or warn BAI modification time check.
pub struct IndexedReaderBuilder {
    bai_path: Option<PathBuf>,
    modification_time: ModificationTime,
    additional_threads: u16,
}

impl IndexedReaderBuilder {
    /// Creates a new [IndexedReader](struct.IndexedReader.html) builder.
    pub fn new() -> Self {
        Self {
            bai_path: None,
            modification_time: ModificationTime::Error,
            additional_threads: 0,
        }
    }

    /// Sets a path to a BAI index. By default, it is `{bam_path}.bai`.
    /// Overwrites the last value, if any.
    pub fn bai_path<P: AsRef<Path>>(&mut self, path: P) -> &mut Self {
        self.bai_path = Some(path.as_ref().to_path_buf());
        self
    }

    /// By default, [IndexedReader::from_path](struct.IndexedReader.html#method.from_path) and
    /// [IndexedReaderBuilder::from_path](struct.IndexedReaderBuilder.html#method.from_path)
    /// returns an `io::Error` if the last modification of the BAI index was earlier
    /// than the last modification of the BAM file.
    ///
    /// Enum [ModificationTime](enum.ModificationTime.html) contains options to skip
    /// this check or raise a warning instead of returning an error.
    pub fn modification_time(&mut self, modification_time: ModificationTime) -> &mut Self {
        self.modification_time = modification_time;
        self
    }

    /// Sets the number of additional threads.
    ///
    /// Additional threads are used to decompress bgzip blocks, while the
    /// main thread reads the blocks from a file/stream.
    /// If `additional_threads` is 0 (default), the main thread will decompress blocks itself.
    pub fn additional_threads(&mut self, additional_threads: u16) -> &mut Self {
        self.additional_threads = additional_threads;
        self
    }

    /// Creates a new [IndexedReader](struct.IndexedReader.html) from `bam_path`.
    /// If BAI path was not specified, the functions tries to open `{bam_path}.bai`.
    pub fn from_path<P: AsRef<Path>>(&self, bam_path: P) -> Result<IndexedReader<File>> {
        let bam_path = bam_path.as_ref();
        let bai_path = self.bai_path.as_ref().map(PathBuf::clone)
            .unwrap_or_else(|| PathBuf::from(format!("{}.ghi", bam_path.display())));
        self.modification_time.check(&bam_path, &bai_path)?;

        let reader = bgzip::SeekReader::from_path(bam_path, self.additional_threads)
            .map_err(|e| Error::new(e.kind(), format!("Failed to open BAM file: {}", e)))?;

        let index = Index::from_path(bai_path)
            .map_err(|e| Error::new(e.kind(), format!("Failed to open BAI index: {}", e)))?;
        IndexedReader::new(reader, index)
    }

    /// Creates a new [IndexedReader](struct.IndexedReader.html) from two streams.
    /// BAM stream should support random access, while BAI stream does not need to.
    /// `check_time` and `bai_path` values are ignored.
    pub fn from_streams<R: Read + Seek, T: Read>(&self, bam_stream: R, bai_stream: T)
            -> Result<IndexedReader<R>> {
        let reader = bgzip::SeekReader::from_stream(bam_stream, self.additional_threads)
            .map_err(|e| Error::new(e.kind(), format!("Failed to read BAM stream: {}", e)))?;

        let index = Index::from_stream(bai_stream)
            .map_err(|e| Error::new(e.kind(), format!("Failed to read BAI index: {}", e)))?;
        IndexedReader::new(reader, index)
    }
}

pub struct IndexedReader<R: Read + Seek> {
    reader: bgzip::SeekReader<R>,
    header: Header,
    index: Index,
}

impl IndexedReader<File> {
    /// Creates [IndexedReaderBuilder](struct.IndexedReaderBuilder.html).
    pub fn build() -> IndexedReaderBuilder {
        IndexedReaderBuilder::new()
    }

    /// Opens bam file from `path`. Bai index will be loaded from `{path}.bai`.
    ///
    /// Same as `Self::build().from_path(path)`.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        IndexedReaderBuilder::new().from_path(path)
    }
}

impl<R: Read + Seek> IndexedReader<R> {
    fn new(mut reader: bgzip::SeekReader<R>, index: Index) -> Result<Self> {
        reader.make_consecutive();
        let header = Header::from_bam(&mut reader)?;
        Ok(Self { reader, header, index })
    }

}