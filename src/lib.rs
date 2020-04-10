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
pub mod binary;
pub mod builder;

use std::io;
use range::InvertedRecord;
use std::io::{Result, Write};

/// A trait for writing BAM/SAM records.
pub trait InvertedRecordWriter {
    /// Writes a single record.
    fn write(&mut self, record: &InvertedRecord) -> io::Result<()>;

    /// Finishes the stream, same as `std::mem::drop(writer)`, but can return an error.
    fn finish(&mut self) -> io::Result<()>;

    /// Flushes contents.
    fn flush(&mut self) -> io::Result<()>;
}

pub trait ColumnarSet {
    fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()>;

    /*fn from_stream<R: Read>(mut stream: R) -> Result<InvertedRecord>;*/
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
