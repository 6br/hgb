#![feature(associated_type_bounds)]
#![feature(rustc_private)]
extern crate log;
extern crate libc;

extern crate serde_json;
extern crate bitpacking;
extern crate regex;
extern crate byteorder;
extern crate bam;

#[no_mangle]
pub mod alignment;
pub mod binary;
pub mod builder;
pub mod compression;
pub mod header;
pub mod index;
pub mod range;
pub mod reader;
pub mod writer;


use std::io;
use range::InvertedRecord;
use range::Record;
use index::Index;
use std::io::{Result, Write, Read};

/// A trait for writing records.
pub trait ChunkWriter {
    /// Writes a single record.
    fn write(&mut self, record: &Record) -> io::Result<index::Chunk>;

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

    fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()>;

    fn from_stream<R: Read>(&mut self, stream: &mut R) -> Result<bool>;
}


#[cfg(test)]
mod tests {
    use crate::header::{Header};
    use bio::io::bed;
    use crate::binary;
    use crate::writer::GhiWriter;
    use crate::reader::IndexedReader;
    use crate::range::{InvertedRecordEntire, Format, Record};
    use crate::IndexWriter;
    use crate::builder::InvertedRecordBuilderSet;
    use crate::index::Region;

    #[test]
    fn full_works() {
        let path = "./test/test.bed";
        let reader = bed::Reader::from_file(path).unwrap();
        let set = InvertedRecordBuilderSet::new(reader, 0 as u64);

        let set_vec = vec![set];
        let entire = InvertedRecordEntire::new(set_vec);
        let header = Header::new();
        let mut writer = binary::GhbWriter::build()
            .write_header(false)
            .from_path("./test/test2.ghb", header).unwrap();
        let index = entire.write_binary(&mut writer).unwrap();
        writer.flush().unwrap();
        println!("{}", index);


        let mut header2 = Header::new();
        entire.write_header(&mut header2);
        let mut index_writer = GhiWriter::build().write_header(true).from_path("./test/test2.ghb.ghi", header2).unwrap();
        let _result = index_writer.write(&index);
        let _result = index_writer.flush();

        let mut reader2 = IndexedReader::from_path("./test/test2.ghb").unwrap();
        println!("{}", reader2.index());

        let viewer = reader2.full();
        let records = viewer.into_iter().scan((), |_,x| x.ok()).collect::<Vec<Record>>();
        println!("Records: {:?}", records);

        assert_eq!(records.len(), 10);
    }

    #[test]
    fn it_works() {
        // let _example = b"1\t5\t5000\tname1\t0.5\n1\t5\t5000\tname1\t0.5";
        let path = "./test/test.bed";
        let reader = bed::Reader::from_file(path).unwrap();
        let set = InvertedRecordBuilderSet::new(reader, 0 as u64);

        let set_vec = vec![set];
        let entire = InvertedRecordEntire::new(set_vec);

        let header = Header::new();
        let mut writer = binary::GhbWriter::build()
            .write_header(false)
            .from_path("./test/test.ghb", header).unwrap();
        let index = entire.write_binary(&mut writer).unwrap();
        writer.flush().unwrap();


        let mut header2 = Header::new();
        entire.write_header(&mut header2);
        let mut index_writer = GhiWriter::build().write_header(true).from_path("./test/test.ghb.ghi", header2).unwrap();
        let _result = index_writer.write(&index);
        assert_eq!(_result.ok(), Some(()));
        let _result = index_writer.flush();

        let mut reader2 = IndexedReader::from_path("./test/test.ghb").unwrap();
        println!("{}", reader2.index());
        

        //assert_eq!(format!("{}", index), format!("{}", reader2.index()));
        // assert_eq!(&index, reader2.index());
        let chrom = "2";
        let chrom_id = reader2.reference_id(&chrom).unwrap();
        println!("{}", chrom_id);
        let viewer = reader2.fetch(&Region::new(chrom_id, 17_000, 17_500)).unwrap();
        // println!("{}", viewer.index());
        /*
        for record in viewer {
            // println!("A");
            let record = record.unwrap();
            println!("Record: {:?}", record);
        }*/
        let records = viewer.into_iter().flat_map(|t| t.map(|f| 
            if let Format::Range(rec) = f.data() {
                // println!("debug {:#?}", rec.to_record(chrom));
                return rec.to_record(chrom)
            } else {
                return vec![]
            }
        ).unwrap()).collect::<Vec<bed::Record>>();
        println!("Records: {:?}", records);

        let example = "2\t16382\t16385\tbin4682\t20\t-\n2\t16388\t31768\tbin4683\t20\t-\n";
        // let mut test_reader = bed::Reader::new(&example[..]);
        let mut buf = Vec::new();
        {
            let mut writer = bed::Writer::new(&mut buf);
            for  i in records {
                writer.write(&i).ok().unwrap();
            } 
        }
        assert_eq!(
            example,
            String::from_utf8(buf).unwrap().as_str()
        );
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