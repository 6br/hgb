
use std::path::{Path, PathBuf};
use crate::index::Index;
use bam::bgzip;
use std::io::{Result, Error, Read, Seek, BufReader};
use std::fs::File;
use std::io::ErrorKind::{InvalidInput, InvalidData};
use bam::header::Header;
use crate::range::Record;
// use crate::ColumnarSet;
use crate::index;
use crate::ChunkReader;
use crate::binary::GhbReader;

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
    pub fn from_path<P: AsRef<Path>, R: Read + Seek>(&self, bam_path: P) -> Result<IndexedReader<BufReader<File>>> {
        let bam_path = bam_path.as_ref();
        let bai_path = self.bai_path.as_ref().map(PathBuf::clone)
            .unwrap_or_else(|| PathBuf::from(format!("{}.ghi", bam_path.display())));
        self.modification_time.check(&bam_path, &bai_path)?;

        let mut index_reader = bgzip::SeekReader::from_path(bai_path, self.additional_threads)
            .map_err(|e| Error::new(e.kind(), format!("Failed to open BAI file: {}", e)))?;
        index_reader.make_consecutive();

        let header = Header::from_bam(&mut index_reader).map_err(|e| Error::new(e.kind(), format!("Failed to read BAI header: {}", e)))?;

        let index = Index::from_stream(&mut index_reader)
        .map_err(|e| Error::new(e.kind(), format!("Failed to read BAI file: {}", e)))?;
    
        let reader = GhbReader::from_path(bam_path, header)
        .map_err(|e| Error::new(e.kind(), format!("Failed to open BAM file: {}", e)))?;


        IndexedReader::new(reader, index)
    }

    /// Creates a new [IndexedReader](struct.IndexedReader.html) from two streams.
    /// BAM stream should support random access, while BAI stream does not need to.
    /// `check_time` and `bai_path` values are ignored.
    pub fn from_streams<R: Read + Seek>(&self, bam_stream: R, bai_stream: R)
            -> Result<IndexedReader<R>> {

        let mut index_reader = bgzip::SeekReader::from_stream(bai_stream, self.additional_threads)
            .map_err(|e| Error::new(e.kind(), format!("Failed to read BAI index: {}", e)))?;
        index_reader.make_consecutive();

        let header = Header::from_bam(&mut index_reader).map_err(|e| Error::new(e.kind(), format!("Failed to read BAI header: {}", e)))?;

        let index = Index::from_stream(&mut index_reader)?;

        let reader = GhbReader::from_stream(bam_stream, header)
        .map_err(|e| Error::new(e.kind(), format!("Failed to read BAM stream: {}", e)))?;


        IndexedReader::new(reader, index)
    }
}

pub struct IndexedReader<R: Read + Seek> {
//    _marker: std::marker::PhantomData<T>,
    reader: GhbReader<R>,
    index: Index
}

impl IndexedReader<BufReader<File>> {
    /// Creates [IndexedReaderBuilder](struct.IndexedReaderBuilder.html).
    pub fn build() -> IndexedReaderBuilder {
        IndexedReaderBuilder::new()
    }

    /// Opens bam file from `path`. Bai index will be loaded from `{path}.bai`.
    ///
    /// Same as `Self::build().from_path(path)`.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        IndexedReaderBuilder::new().from_path::<P,BufReader<File>>(path)
    }
    
}

impl<R: Read + Seek> IndexedReader<R> {
    fn new(mut reader: GhbReader<R>, index: Index) -> Result<Self> {
        // reader.make_consecutive();
        // let _marker = std::marker::PhantomData;
        // let header = reader.header()// Header::from_bam(&mut index_reader)?;
        Ok(Self { reader, index })
    }

    pub fn index(&self) -> &Index {
        &self.index
    }

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

    pub fn fetch_by<'a, F>(&'a mut self, region: &Region, predicate: F) -> Result<RegionViewer<'a, R>>
    where F: 'static + Fn(&Record) -> bool
    {

        match self.header().reference_len(region.ref_id()) {
            None => return Err(Error::new(InvalidInput,
                format!("Failed to fetch records: out of bounds reference {}", region.ref_id()))),
            Some(len) if len < region.end() => return Err(Error::new(InvalidInput,
                format!("Failed to fetch records: end > reference length ({} > {})", region.end(), len))),
            _ => {},
        }

        let chunks = self.index.fetch_chunks(region.ref_id(), region.start() as i32, region.end() as i32);
        self.reader.set_chunks(chunks);

        Ok(RegionViewer {
            parent: self,
            start: region.start() as i32,
            end: region.end() as i32,
            predicate: Box::new(predicate),
        })
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
    // parent_stream: &'a mut GhbReader<R>,
    start: i32,
    end: i32,
    predicate: Box<dyn Fn(&Record) -> bool>,
}

impl<'a, R: Read + Seek> RegionViewer<'a, R> {
    /// Returns [header](../header/struct.Header.html).
    pub fn header(&self) -> &Header {
        self.parent.header()
    }

    /// Returns [BAI index](../index/struct.Index.html).
    pub fn index(&self) -> &Index {
        self.parent.index()
    }
/*
    pub fn take_stream(self) -> R {
        self.parent.take_stream()
    }*/
}



/// Genomic coordinates, used in [struct.IndexedReader.html#method.fetch] and [struct.IndexedReader.html#method.pileup].
/// `ref_id` is 0-based, `start-end` is 0-based half-open interval.
#[derive(Clone, Debug)]
pub struct Region {
    ref_id: u32,
    start: u32,
    end: u32,
}

impl Region {
    /// Creates new region. `ref_id` is 0-based, `start-end` is 0-based half-open interval.
    pub fn new(ref_id: u32, start: u32, end: u32) -> Region {
        assert!(start <= end, "Region: start should not be greater than end ({} > {})", start, end);
        Region { ref_id, start, end }
    }

    pub fn ref_id(&self) -> u32 {
        self.ref_id
    }

    pub fn start(&self) -> u32 {
        self.start
    }

    pub fn end(&self) -> u32 {
        self.end
    }

    pub fn len(&self) -> u32 {
        self.end - self.start
    }

    pub fn set_ref_id(&mut self, ref_id: u32) {
        self.ref_id = ref_id;
    }

    pub fn set_start(&mut self, start: u32) {
        assert!(start <= self.end, "Region: start should not be greater than end ({} > {})", start, self.end);
        self.start = start;
    }

    pub fn set_end(&mut self, end: u32) {
        assert!(self.start <= end, "Region: start should not be greater than end ({} > {})", self.start, end);
        self.end = end;
    }

    pub fn contains(&self, ref_id: u32, pos: u32) -> bool {
        self.ref_id == ref_id && self.start <= pos && pos < self.end
    }
}


impl<'a, R: Read + Seek> ChunkReader for RegionViewer<'a, R> {
    fn read_into(&mut self, record: &mut Record) -> Result<bool> {
        loop {
            // let res = record.fill_from_bam(&mut self.parent.take_stream());
            let res = self.parent.reader.fill_from_binary(record);
            if !res.as_ref().unwrap_or(&false) {
                record.clear();
                return res;
            }
            // Reads are sorted, so no more reads would be in the region.
            /*
            if record.start() >= self.end {
                record.clear();
                return Ok(false);
            }
            */
            if !(self.predicate)(&record) {
                continue;
            }
//            let record_bin = record.calculate_bin();
            /*
            if record_bin > index::MAX_BIN {
                record.clear();
                return Err(Error::new(InvalidData, "Read has BAI bin bigger than max possible value"));
            }

            let (min_start, max_end) = index::bin_to_region(record_bin);
            if min_start >= self.start && max_end <= self.end {
                return Ok(true);
            }
            let record_end = record.calculate_end();
            
            if record.flag().is_mapped() && record_end < record.start() {
                record.clear();
                return Err(Error::new(InvalidData, "Corrupted record: aln_end < aln_start"));
            }
            if record.flag().is_mapped() {
                if record_end > self.start {
                    return Ok(true);
                }
            } else if record.start() >= self.start {
                return Ok(true);
            }*/
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
            Err(e) => Some(Err(e)),
        }
    }
}