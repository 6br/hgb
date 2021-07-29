//#![feature(associated_type_bounds)]
//#![feature(rustc_private)]
//extern crate libc;
extern crate csv;
extern crate log;

extern crate bam;
extern crate bitpacking;
extern crate byteorder;
extern crate enum_map;
extern crate regex;
extern crate serde_json;
extern crate twobit;

//pub mod alignment;
pub mod bed;
pub mod binary;
//pub mod buffer;
pub mod builder;
pub mod checker_index;
pub mod color;
pub mod compression;
pub mod dump;
pub mod gff;
pub mod header;
pub mod index;
//pub mod server;
pub mod range;
pub mod reader;
pub mod simple_buffer;

pub mod twopass_alignment;
pub mod vis;
//pub mod vis_orig;
pub mod writer;

use bam::IndexedReader;
use checker_index::Index;
use genomic_range::StringRegion;
use io::Seek;
use range::InvertedRecord;
use range::{Format, Record};
use std::io;
use std::{
    collections::BTreeMap,
    io::{Read, Result, Write},
    str::FromStr,
};

pub type Frequency = BTreeMap<u64, Vec<(u64, u32, char)>>;

/// A trait for writing records.
pub trait ChunkWriter<R: Read + Seek> {
    /// Writes a single record.
    fn write(
        &mut self,
        record: &Record,
        //         threads: u16,
        bam_reader: Option<&mut IndexedReader<R>>,
    ) -> io::Result<index::Chunk>;

    /// Finishes the stream, same as `std::mem::drop(writer)`, but can return an error.
    fn finish(&mut self) -> io::Result<()>;

    /// Flushes contents.
    fn flush(&mut self) -> io::Result<()>;
}

/// A trait for reading records.
pub trait ChunkReader: Iterator<Item = io::Result<Record>> {
    /// Writes the next record into `record`. It allows to skip excessive memory allocation.
    /// If there are no more records to iterate over, the function returns `false`.
    ///
    /// If the function returns an error, the record is cleared.
    fn read_into(&mut self, record: &mut Record) -> io::Result<bool>;

    /// Pauses multi-thread reader until the next read operation. Does nothing to a single-thread reader.
    ///
    /// Use with caution: pausing and unpausing takes some time.
    fn pause(&mut self);
}

/// A trait for writing records. (Deprecated.)
pub trait InvertedRecordWriter {
    /// Writes a single record.
    fn write(&mut self, record: &InvertedRecord) -> io::Result<index::VirtualOffset>;

    /// Finishes the stream, same as `std::mem::drop(writer)`, but can return an error.
    fn finish(&mut self) -> io::Result<()>;

    /// Flushes contents.
    fn flush(&mut self) -> io::Result<()>;
}

/// A trait for writing records.
pub trait IndexWriter {
    /// Writes a single record.
    fn write(&mut self, record: &Index) -> io::Result<()>;

    /// Finishes the stream, same as `std::mem::drop(writer)`, but can return an error.
    fn finish(&mut self) -> io::Result<()>;

    /// Flushes contents.
    fn flush(&mut self) -> io::Result<()>;
}

/// A trait for Column-convertable dataset.
pub trait ColumnarSet {
    fn new() -> Self;

    fn to_stream<W: Write, R: Read + Seek>(
        &self,
        stream: &mut W,
        threads: u16,
        bam_reader: Option<&mut IndexedReader<R>>,
    ) -> Result<()>;

    // Returns a consumed offset.
    fn from_stream<U: Read>(&mut self, stream: &mut U, threads: u16) -> Result<u64>;
}

/// A trait for Formattable data structure.
pub trait Builder {
    fn to_format(self) -> Format;
}

pub trait Annotation {
    fn start(&self) -> u64;
    fn end(&self) -> u64;
    fn name(&self) -> &str;
}

/// Visualization Presets
#[derive(Debug)]
pub enum VisPreset {
    Auto,       // Switch presets Depends on range length.
    Base,       // all cigars and insertions, compatible with IGV
    Gene,       // Only length and split alignment
    OnlySV,     // Split alignment enabled mode
    Qual,       // Quality mode
    Chromosome, // Only coverages.
}

impl VisPreset {
    pub fn augment(range: StringRegion) -> VisPreset {
        let diff = range.end() - range.start();
        if diff <= 10000 {
            VisPreset::Base
        } else if diff <= 100000 {
            VisPreset::Gene
        } else {
            VisPreset::Chromosome
        }
    }
}

// Implement the trait
impl FromStr for VisPreset {
    type Err = &'static str;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "gene" => Ok(VisPreset::Gene),
            "auto" => Ok(VisPreset::Auto),
            "base" => Ok(VisPreset::Base),
            "chrom" => Ok(VisPreset::Chromosome),
            "sv" => Ok(VisPreset::OnlySV),
            "qual" => Ok(VisPreset::Qual),
            _ => Ok(VisPreset::Auto),
        }
    }
}
#[derive(Debug, Clone)]
pub struct Vis {
    pub range: StringRegion,
    //list: Vec<(u64, bam::Record)>,
    pub annotation: Vec<(u64, bed::Record)>,
    pub freq: BTreeMap<u64, Vec<(u64, u32, char)>>,
    pub compressed_list: Vec<(u64, usize)>,
    pub index_list: Vec<usize>,
    pub prev_index: usize,
    pub supplementary_list: Vec<(Vec<u8>, usize, usize, i32, i32)>,
    pub prefetch_max: u64,
}

impl Vis {
    pub fn new(
        range: StringRegion,
        //list: Vec<(u64, bam::Record)>
        annotation: Vec<(u64, bed::Record)>,
        freq: BTreeMap<u64, Vec<(u64, u32, char)>>,
        compressed_list: Vec<(u64, usize)>,
        index_list: Vec<usize>,
        prev_index: usize,
        supplementary_list: Vec<(Vec<u8>, usize, usize, i32, i32)>,
        prefetch_max: u64,
    ) -> Vis {
        Vis {
            range,
            //list,
            annotation,
            freq,
            compressed_list,
            index_list,
            prev_index,
            supplementary_list,
            prefetch_max,
        }
    }
}

#[derive(Clone)]
pub struct VisOrig<'a> {
    range: StringRegion,
    list: Vec<(u64, bam::Record)>,
    annotation: Vec<(u64, bed::Record)>,
    frequency: BTreeMap<u64, Vec<(u64, u32, char)>>,
    compressed_list: &'a Vec<(u64, usize)>,
    index_list: Vec<usize>,
    prev_index: usize,
    supplementary_list: &'a Vec<(Vec<u8>, usize, usize, i32, i32)>,
}

impl<'a> VisOrig<'a> {
    pub fn new(
        range: StringRegion,
        list: Vec<(u64, bam::Record)>,
        annotation: Vec<(u64, bed::Record)>,
        frequency: BTreeMap<u64, Vec<(u64, u32, char)>>,
        compressed_list: &'a Vec<(u64, usize)>,
        index_list: Vec<usize>,
        prev_index: usize,
        supplementary_list: &'a Vec<(Vec<u8>, usize, usize, i32, i32)>,
    ) -> Self {
        VisOrig {
            range,
            list,
            annotation,
            frequency,
            compressed_list,
            index_list,
            prev_index,
            supplementary_list,
        }
    }

    pub fn convert(&self) -> VisRef {
        VisRef {
            range: self.range.clone(),
            list: &self.list,
            annotation: &self.annotation,
            frequency: &self.frequency,
            compressed_list: &self.compressed_list,
            index_list: &self.index_list,
            prev_index: self.prev_index,
            supplementary_list: &self.supplementary_list,
        }
    }
}

#[derive(Clone)]
pub struct VisRef<'a> {
    range: StringRegion,
    list: &'a Vec<(u64, bam::Record)>,
    annotation: &'a Vec<(u64, bed::Record)>,
    frequency: &'a BTreeMap<u64, Vec<(u64, u32, char)>>,
    compressed_list: &'a Vec<(u64, usize)>,
    index_list: &'a Vec<usize>,
    prev_index: usize,
    supplementary_list: &'a Vec<(Vec<u8>, usize, usize, i32, i32)>,
}

impl<'a> VisRef<'a> {
    pub fn new(
        range: StringRegion,
        list: &'a Vec<(u64, bam::Record)>,
        annotation: &'a Vec<(u64, bed::Record)>,
        frequency: &'a BTreeMap<u64, Vec<(u64, u32, char)>>,
        compressed_list: &'a Vec<(u64, usize)>,
        index_list: &'a std::vec::Vec<usize>,
        prev_index: usize,
        supplementary_list: &'a std::vec::Vec<(std::vec::Vec<u8>, usize, usize, i32, i32)>,
    ) -> Self {
        VisRef {
            range,
            list,
            annotation,
            frequency,
            compressed_list,
            index_list,
            prev_index,
            supplementary_list,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::bed;
    use crate::binary;
    use crate::builder::InvertedRecordBuilder;
    use crate::header::Header;
    use crate::range::{Format, InvertedRecordEntire};
    use crate::reader::IndexedReader;
    use crate::writer::GhiWriter;
    use crate::IndexWriter;
    use crate::{index::Region, range::Set};
    use std::fs::File;

    #[test]
    fn full_works() {
        let path = "./test/test.bed";
        let reader = bed::Reader::from_file(path).unwrap();
        //let set = InvertedRecordBuilderSet::new(reader, 0 as u64);
        let mut header_2 = Header::new();
        let set: Set<InvertedRecordBuilder, File> =
            Set::<InvertedRecordBuilder, File>::new(reader, 1_u64, &mut header_2).unwrap();

        //        let set_vec = vec![set];
        // println!("{:?}", set);

        let set_vec = vec![set];
        //        let entire = InvertedRecordEntire::new(set_vec);
        let mut entire: InvertedRecordEntire<File> = InvertedRecordEntire::new_from_set(set_vec);

        let header = Header::new();
        let mut writer = binary::GhbWriter::build()
            .write_header(false)
            .from_path("./test/test2.ghb", header)
            .unwrap();
        let index = entire.write_binary(&mut writer).unwrap();
        writer.flush().unwrap();
        // println!("{}", index);

        let mut header2 = Header::new();
        entire.write_header(&mut header2);
        let mut index_writer = GhiWriter::build()
            .write_header(true)
            .from_path("./test/test2.ghb.ghi", header2)
            .unwrap();
        let _result = index_writer.write(&index);
        let _result = index_writer.flush();

        let mut reader2 = IndexedReader::from_path("./test/test2.ghb").unwrap();
        println!("{}", reader2.index());

        let viewer = reader2.full();

        let bed_records = viewer
            .into_iter()
            .flat_map(|t| {
                t.map(|f| {
                    if let Format::Range(rec) = f.data() {
                        // println!("debug {:#?}", rec.to_record(chrom));
                        rec.to_record("null")
                    } else {
                        return vec![];
                    }
                })
                .unwrap()
            })
            .collect::<Vec<bed::Record>>();
        assert_eq!(bed_records.len(), 10);
    }

    #[test]
    fn fetch_works() {
        // let _example = b"1\t5\t5000\tname1\t0.5\n1\t5\t5000\tname1\t0.5";
        let path = "./test/test.bed";
        let reader = bed::Reader::from_file(path).unwrap();
        let mut header2 = Header::new();
        let set: Set<InvertedRecordBuilder, File> =
            Set::<InvertedRecordBuilder, File>::new(reader, 1_u64, &mut header2).unwrap();

        let set_vec = vec![set];
        let mut entire: InvertedRecordEntire<File> = InvertedRecordEntire::new_from_set(set_vec);

        let header = Header::new();
        let mut writer = binary::GhbWriter::build()
            .write_header(false)
            .from_path("./test/test.ghb", header)
            .unwrap();
        let index = entire.write_binary(&mut writer).unwrap();
        writer.flush().unwrap();

        // let mut header2 = Header::new();
        // entire.write_header(&mut header2);
        let mut index_writer = GhiWriter::build()
            .write_header(true)
            .from_path("./test/test.ghb.ghi", header2)
            .unwrap();
        let _result = index_writer.write(&index);
        assert_eq!(_result.ok(), Some(()));
        let _result = index_writer.flush();

        let mut reader2 = IndexedReader::from_path("./test/test.ghb").unwrap();
        // println!("{}", reader2.index());

        // assert_eq!(format!("{}", index), format!("{}", reader2.index()));
        // assert_eq!(&index, reader2.index());
        let chrom = "2";
        let chrom_id = reader2.reference_id(&chrom).unwrap();
        assert_eq!(chrom_id, 1);

        let viewer = reader2
            .fetch(&Region::new(chrom_id, 17_000, 17_500))
            .unwrap();

        let records = viewer
            .into_iter()
            .flat_map(|t| {
                t.map(|f| {
                    if let Format::Range(rec) = f.data() {
                        // println!("debug {:#?}", rec);
                        rec.to_record(chrom)
                    } else {
                        return vec![];
                    }
                })
                .unwrap()
            })
            .collect::<Vec<bed::Record>>();
        println!("Records: {:?}", records);

        let example = "2\t16382\t16385\tbin4682\t20\t-\n2\t16388\t31768\tbin4683\t20\t-\n";

        let mut buf = Vec::new();
        {
            let mut writer = bed::Writer::new(&mut buf);
            for i in records {
                writer.write(&i).ok().unwrap();
            }
        }
        assert_eq!(example, String::from_utf8(buf).unwrap().as_str());
        // let records_from = reader.records().into_iter().flat_map(|t| t).collect::<Vec<bed::Record>>();
        //assert_eq!(records[0], records_from[0]);
        //assert_eq!(records[1], records_from[1]);
        //assert_eq!(records.len(), records_from.len());

        let viewer = reader2.fetch(&Region::new(0, 17_000, 17_500)).unwrap();

        for record in viewer {
            let _record = record.unwrap();
            // println!("Record: {:?}", record);
        }

        let viewer = reader2.fetch(&Region::new(1, 1, 3)).unwrap();
        // println!("{}", viewer.index());
        for record in viewer {
            // println!("A");
            let _record = record.unwrap();
            // println!("Record: {:?}", record);
        }
    }
}
