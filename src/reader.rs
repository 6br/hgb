use crate::binary::GhbReader;
use crate::checker_index::Index;
use crate::header::Header;
use crate::index::{Chunk, Region, VirtualOffset};
use crate::range::Record;
use crate::ChunkReader;
use bam::bgzip;
use std::fs::File;
use std::io::ErrorKind::{self, InvalidInput};
use std::io::{BufReader, Error, Read, Result, Seek};
use std::path::{Path, PathBuf};

/// Defines how to react to a GHI index being younger than GHB file.
///
/// # Variants
/// * `Error` - [IndexedReader](struct.IndexedReader.html) will not be constructed if the GHI
/// index is was modified earlier than the GHB file. `io::Error` will be raised.
/// * `Ignore` - does nothing if the index is younger than the GHB file.
/// * `Warn` - calls a function `Fn(&str)` and continues constructing
/// [IndexedReader](struct.IndexedReader.html);
pub enum ModificationTime {
    Error,
    Ignore,
    Warn(Box<dyn Fn(&str)>),
}

impl ModificationTime {
    fn check<T: AsRef<Path>, U: AsRef<Path>>(&self, ghb_path: T, ghi_path: U) -> Result<()> {
        let ghb_modified = ghb_path
            .as_ref()
            .metadata()
            .and_then(|metadata| metadata.modified());
        let ghi_modified = ghi_path
            .as_ref()
            .metadata()
            .and_then(|metadata| metadata.modified());
        let ghb_younger = match (ghb_modified, ghi_modified) {
            (Ok(ghb_time), Ok(ghi_time)) => ghi_time < ghb_time,
            _ => false, // Modification time not available.
        };
        if !ghb_younger {
            return Ok(());
        }

        match &self {
            ModificationTime::Ignore => {}
            ModificationTime::Error => {
                return Err(Error::new(
                    InvalidInput,
                    "the binary file is younger than the index",
                ))
            }
            ModificationTime::Warn(box_fun) => {
                box_fun("the binary file is younger than the index index")
            }
        }
        Ok(())
    }

    /// Create a warning strategy `ModificationTime::Warn`.
    pub fn warn<F: Fn(&str) + 'static>(warning: F) -> Self {
        ModificationTime::Warn(Box::new(warning))
    }
}

/// [IndexedReader](struct.IndexedReader.html) builder. Allows to specify paths to GHB and GHI
/// files, as well as the number of threads
/// and an option to ignore or warn GHI modification time check.
pub struct IndexedReaderBuilder {
    ghi_path: Option<PathBuf>,
    modification_time: ModificationTime,
    additional_threads: u16,
}

impl IndexedReaderBuilder {
    /// Creates a new [IndexedReader](struct.IndexedReader.html) builder.
    pub fn new() -> Self {
        Self {
            ghi_path: None,
            modification_time: ModificationTime::Error,
            additional_threads: 0,
        }
    }

    /// Sets a path to GHI index. By default, it is `{ghb_path}.ghi`.
    /// Overwrites the last value, if any.
    pub fn ghi_path<P: AsRef<Path>>(&mut self, path: P) -> &mut Self {
        self.ghi_path = Some(path.as_ref().to_path_buf());
        self
    }

    /// By default, [IndexedReader::from_path](struct.IndexedReader.html#method.from_path) and
    /// [IndexedReaderBuilder::from_path](struct.IndexedReaderBuilder.html#method.from_path)
    /// returns an `io::Error` if the last modification of the GHI index was earlier
    /// than the last modification of the GHB file.
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

    /// Creates a new [IndexedReader](struct.IndexedReader.html) from `ghb_path`.
    /// If GHI path was not specified, the functions tries to open `{ghb_path}.ghi`.
    pub fn from_path<P: AsRef<Path>, R: Read + Seek>(
        &self,
        ghb_path: P,
    ) -> Result<IndexedReader<BufReader<File>>> {
        let ghb_path = ghb_path.as_ref();
        let ghi_path = self
            .ghi_path
            .as_ref()
            .map(PathBuf::clone)
            .unwrap_or_else(|| PathBuf::from(format!("{}.ghi", ghb_path.display())));
        self.modification_time.check(&ghb_path, &ghi_path)?;

        let mut index_reader = bgzip::SeekReader::from_path(ghi_path, self.additional_threads)
            .map_err(|e| Error::new(e.kind(), format!("Failed to open GHI file: {}", e)))?;
        index_reader.make_consecutive();

        let header = Header::new_from_stream(&mut index_reader)
            .map_err(|e| Error::new(e.kind(), format!("Failed to read GHI header: {}", e)))?;

        let index = Index::from_stream(&mut index_reader)
            .map_err(|e| Error::new(e.kind(), format!("Failed to read GHI file: {}", e)))?;

        let reader = GhbReader::from_path(ghb_path, header, self.additional_threads)
            .map_err(|e| Error::new(e.kind(), format!("Failed to open GHB file: {}", e)))?;

        IndexedReader::new(reader, index)
    }

    /// Creates a new [IndexedReader](struct.IndexedReader.html) from two streams.
    /// GHB stream should support random access, while GHI stream does not need to.
    /// `check_time` and `ghi_path` values are ignored.
    pub fn from_streams<R: Read + Seek>(
        &self,
        bam_stream: R,
        bai_stream: R,
    ) -> Result<IndexedReader<R>> {
        let mut index_reader = bgzip::SeekReader::from_stream(bai_stream, self.additional_threads)
            .map_err(|e| Error::new(e.kind(), format!("Failed to read GHI index: {}", e)))?;
        index_reader.make_consecutive();

        let header = Header::new_from_stream(&mut index_reader)
            .map_err(|e| Error::new(e.kind(), format!("Failed to read GHI header: {}", e)))?;

        let index = Index::from_stream(&mut index_reader)?;

        let reader = GhbReader::from_stream(bam_stream, header, self.additional_threads)
            .map_err(|e| Error::new(e.kind(), format!("Failed to read GHB stream: {}", e)))?;

        IndexedReader::new(reader, index)
    }
}

pub struct IndexedReader<R: Read + Seek> {
    //    _marker: std::marker::PhantomData<T>,
    reader: GhbReader<R>,
    index: Index,
}

impl IndexedReader<BufReader<File>> {
    /// Creates [IndexedReaderBuilder](struct.IndexedReaderBuilder.html).
    pub fn build() -> IndexedReaderBuilder {
        IndexedReaderBuilder::new()
    }

    pub fn reference_id(&self, chrom: &str) -> Option<u64> {
        self.reader.header().reference_id(chrom)
    }

    pub fn reference_name(&self, chrom: u64) -> Option<&str> {
        self.reader.header().reference_name(chrom)
    }

    /// Opens GHB file from `path`. Ghi index will be loaded from `{path}.ghi`.
    ///
    /// Same as `Self::build().from_path(path)`.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        IndexedReaderBuilder::new().from_path::<P, BufReader<File>>(path)
    }

    /// Opens GHB file from `path`. Ghi index will be loaded from `{path}.ghi`.
    ///
    /// Same as `Self::build().from_path(path)`.
    pub fn from_path_with_additional_threads<P: AsRef<Path>>(
        path: P,
        threads: u16,
    ) -> Result<Self> {
        IndexedReaderBuilder::new()
            .additional_threads(threads)
            .from_path::<P, BufReader<File>>(path)
    }
}

impl<R: Read + Seek> IndexedReader<R> {
    fn new(reader: GhbReader<R>, index: Index) -> Result<Self> {
        // reader.make_consecutive();
        // let _marker = std::marker::PhantomData;
        // let header = reader.header()// Header::from_bam(&mut index_reader)?;
        Ok(Self { reader, index })
    }

    /// Returns [index](../index/struct.Index.html).
    pub fn index(&self) -> &Index {
        &self.index
    }

    /// Returns [header](../header/struct.Header.html).
    pub fn header(&self) -> &Header {
        &self.reader.header()
    }

    pub fn take_stream(self) -> R {
        self.reader.take_stream()
    }

    /// Returns an iterator over records aligned to the [reference region](struct.Region.html).
    pub fn fetch<'a>(&'a mut self, region: &Region) -> Result<RegionViewer<'a, R>> {
        self.fetch_by(region, |_| true)
    }

    /// Returns an iterator over records aligned to the [reference region](struct.Region.html).
    ///
    /// Records will be filtered by `predicate`. It helps to slightly reduce fetching time,
    /// as some records will be removed without allocating new memory and without calculating
    /// alignment length.
    pub fn fetch_by<'a, F>(
        &'a mut self,
        region: &Region,
        predicate: F,
    ) -> Result<RegionViewer<'a, R>>
    where
        F: 'static + Fn(&Record) -> bool,
    {
        match self.header().reference_len(region.ref_id()) {
            None => {
                return Err(Error::new(
                    InvalidInput,
                    format!(
                        "Failed to fetch records: out of bounds reference {}",
                        region.ref_id()
                    ),
                ))
            }
            Some(len) if len < region.end() => {
                return Err(Error::new(
                    InvalidInput,
                    format!(
                        "Failed to fetch records: end > reference length ({} > {})",
                        region.end(),
                        len
                    ),
                ))
            }
            _ => {}
        }

        let chunks = self
            .index
            .fetch_chunks(region.ref_id(), region.start(), region.end());
        //eprintln!("{:?}", chunks);

        self.reader.set_chunks(chunks);

        Ok(RegionViewer {
            parent: self,
            //            start: region.start() as i32,
            //            end: region.end() as i32,
            predicate: Box::new(predicate),
        })
    }

    /// Returns an iterator over records aligned to the [reference region](struct.Region.html).
    pub fn chunk<'a>(&'a mut self, region: Vec<Chunk>) -> Result<RegionViewer<'a, R>> {
        self.chunk_by(region, |_| true)
    }

    /// Returns an iterator over records aligned to the [reference region](struct.Region.html).
    ///
    /// Records will be filtered by `predicate`. It helps to slightly reduce fetching time,
    /// as some records will be removed without allocating new memory and without calculating
    /// alignment length.
    pub fn chunk_by<'a, F>(
        &'a mut self,
        chunks: Vec<Chunk>,
        predicate: F,
    ) -> Result<RegionViewer<'a, R>>
    where
        F: 'static + Fn(&Record) -> bool,
    {
        self.reader.set_chunks(chunks);

        Ok(RegionViewer {
            parent: self,
            //            start: region.start() as i32,
            //            end: region.end() as i32,
            predicate: Box::new(predicate),
        })
    }

    /// Returns an iterator over all records from the start of the GHB file.
    pub fn full<'a>(&'a mut self) -> RegionViewer<'a, R> {
        self.full_by(|_| true)
    }

    /// Returns an iterator over all records from the start of the GHB file.
    ///
    /// Records will be filtered by `predicate`, which allows to skip some records without allocating new memory.
    pub fn full_by<'a, F>(&'a mut self, predicate: F) -> RegionViewer<'a, R>
    where
        F: 'static + Fn(&Record) -> bool,
    {
        //if let Some(offset) = self.index.start_offset() {
        self.reader.set_chunks(vec![Chunk::new(
            0,
            0,
            VirtualOffset::new(0, 0),
            VirtualOffset::MAX,
        )]);
        //self.reader.set_chunks(self.index.chunks());
        //}
        RegionViewer {
            parent: self,
            //                start: std::i32::MIN,
            //                end: std::i32::MAX,
            predicate: Box::new(predicate),
        }
    }
}

/// Iterator over records in a specific region.
/// Implements [RecordReader](../trait.RecordReader.html) trait.
///
/// If possible, create a single record using [Record::new](../record/struct.Record.html#method.new)
/// and then use [read_into](../trait.RecordReader.html#method.read_into) instead of iterating,
/// as it saves time on allocation.
pub struct RegionViewer<'a, R: Read + Seek> {
    parent: &'a mut IndexedReader<R>,
    // start: u64,
    // end: u64,
    predicate: Box<dyn Fn(&Record) -> bool>,
}

impl<'a, R: Read + Seek> RegionViewer<'a, R> {
    /// Returns [header](../header/struct.Header.html).
    pub fn header(&self) -> &Header {
        self.parent.header()
    }

    /// Returns [index](../index/struct.Index.html).
    pub fn index(&self) -> &Index {
        self.parent.index()
    }
}

impl<'a, R: Read + Seek> ChunkReader for RegionViewer<'a, R> {
    fn read_into(&mut self, record: &mut Record) -> Result<bool> {
        loop {
            let res = self.parent.reader.fill_from_binary(record);
            if !res.as_ref().unwrap_or(&false) {
                record.clear();
                return res;
            }
            if !(self.predicate)(&record) {
                continue;
            }
            return Ok(true);
        }
    }

    fn pause(&mut self) {
        // self.parent.pause();
    }
}

/// Iterator over records.
impl<'a, R: Read + Seek> Iterator for RegionViewer<'a, R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Record::new();
        match self.read_into(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => None,
            Err(e) => match e.kind() {
                ErrorKind::UnexpectedEof => None, // The end of file; normal.
                _ => Some(Err(e)),
            }, // Some(Err(e)),
        }
    }
}
