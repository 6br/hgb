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
pub mod index;
pub mod range;
pub mod alignment;
pub mod writer;
pub mod reader;
pub mod binary;
pub mod builder;


use std::io;
use range::InvertedRecord;
use range::Record;
use index::Index;
use std::io::{Result, Write, Read};

/// A trait for writing BAM/SAM records.
pub trait ChunkWriter {
    /// Writes a single record.
    fn write(&mut self, record: &Record) -> io::Result<index::Chunk>;

    /// Finishes the stream, same as `std::mem::drop(writer)`, but can return an error.
    fn finish(&mut self) -> io::Result<()>;

    /// Flushes contents.
    fn flush(&mut self) -> io::Result<()>;
}

/// A trait for reading BAM/SAM records.
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


/// A trait for writing BAM/SAM records. (Deprecated.)
pub trait InvertedRecordWriter {
    /// Writes a single record.
    fn write(&mut self, record: &InvertedRecord) -> io::Result<index::VirtualOffset>;

    /// Finishes the stream, same as `std::mem::drop(writer)`, but can return an error.
    fn finish(&mut self) -> io::Result<()>;

    /// Flushes contents.
    fn flush(&mut self) -> io::Result<()>;
}

pub trait InvertedRecordReader: Iterator<Item = io::Result<InvertedRecord>> {
    /// Writes the next record into `record`. It allows to skip excessive memory allocation.
    /// If there are no more records to iterate over, the function returns `false`.
    ///
    /// If the function returns an error, the record is cleared.
    fn read_into(&mut self, record: &mut InvertedRecord) -> io::Result<bool>;

    /// Pauses multi-thread reader until the next read operation. Does nothing to a single-thread reader.
    ///
    /// Use with caution: pausing and unpausing takes some time.
    fn pause(&mut self);
}

/// A trait for writing BAM/SAM records.
pub trait IndexWriter {
    /// Writes a single record.
    fn write(&mut self, record: &Index) -> io::Result<()>;

    /// Finishes the stream, same as `std::mem::drop(writer)`, but can return an error.
    fn finish(&mut self) -> io::Result<()>;

    /// Flushes contents.
    fn flush(&mut self) -> io::Result<()>;
}

pub trait ColumnarSet {
    fn new() -> Self;

    //fn clear(&mut self);
    fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()>;

    fn from_stream<R: Read>(&mut self, stream: &mut R) -> Result<bool>;
}

#[no_mangle]
pub extern "C" fn hello_rust() -> *const u8 {
    "Hello, world!\0".as_ptr()
}


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        
        assert_eq!(2 + 2, 4);
    }
}
