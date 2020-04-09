use bio::io::bed;
use std::io::{Result, Read, Write};
use std::collections::HashMap;
use std::cell::RefCell;
use crate::index::region_to_bin_2;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Chunk<T> {
    length: u32,
    sample_id: u32,
    sample_file_id: u32, // When we split files of inverted data structure
    // format_type: u32, // Enum
    data: T,
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordReference {
    bins: HashMap<u32, InvertedRecordBuilder> // Mutex?
}

impl InvertedRecordReference {
    pub fn new() -> Self {
        InvertedRecordReference{
            bins:  HashMap::new()
        }
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordSet {
    chrom: HashMap<String, InvertedRecordReference> // Mutex?
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct InvertedRecord {
    start: Vec<u64>,
    end: Vec<u64>,
    name: Vec<String>, 
    // aux: Vec<String> // We should update better data format.
}

impl InvertedRecord {
    pub fn to_stream<W: Write>(&self, f: &mut W) -> Result<()> {
        stream.write_u64::<LittleEndian>(self.start.len)?;
        stream.write_u64::<LittleEndian>(3)?;
        stream.write_u8::<LittleEndian>(1)?; // u64
        stream.write_u8::<LittleEndian>(1)?; // u64 
        stream.write_u8::<LittleEndian>(0)?; // String

        for i in self.start {
            stream.write_u64::<LittleEndian>(i)?;
        }
        for i in self.end {
            stream.write_u64::<LittleEndian>(i)?;
        }
        for i in self.name {
            stream.write_u64::<LittleEndian>(i.len as u64)?;
            stream.write_all(i)?
        }

        Ok(())
    }

    pub fn from_stream<R: Read>(mut stream: R) -> Result<InvertedRecord> {
        let n_header = stream.read_u64::<LittleEndian>()?;
        let n_items = stream.read_u64::<LittleEndian>()?;
        for i in 0..n_header {
            stream.read_u8::<LittleEndian>()?;
        }
        let mut start = Vec::with_capacity(n_items);
        let mut end = Vec::with_capacity(n_items);
        let mut name = Vec::with_capacity(n_items);
        for i in 0..n_items {
            start.push(stream.read_u64::<LittleEndian>()?);
        }
        for i in 0..n_items {
            end.push(stream.read_u64::<LittleEndian>()?);
        }
        for i in 0..n_items {
            unsafe {

            }
            let mut 
            name.push(stream.read_exact()?);
        }
        return Ok(InvertedRecord{start: start, end: end, name: name})
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct InvertedRecordBuilder {
    start: RefCell<Vec<u64>>,
    end: RefCell<Vec<u64>>,
    name: RefCell<Vec<String>> //Just a tab-separated string here.
}

impl InvertedRecordBuilder {
    pub fn new() -> Self {
        InvertedRecordBuilder{
            start: RefCell::new(Vec::new()),
            end: RefCell::new(Vec::new()),
            name: RefCell::new(Vec::new())
        }
    }
    pub fn add(&self, start: u64, end: u64, name: String) {
        self.start.borrow_mut().push(start);
        self.end.borrow_mut().push(end);
        self.name.borrow_mut().push(name);

    }
}

impl InvertedRecordSet {
    pub fn new<R:Read>(mut reader: bed::Reader<R>) -> Self {
        let mut inverted_record_set = HashMap::new();
        /*let mut start: Vec<u64> = vec![];
        let mut end: Vec<u64> = vec![];
        let mut aux: Vec<String> = vec![];*/
        for record in reader.records() {
            let rec = record.ok().expect("Error reading record.");
            let chrom = inverted_record_set.entry(rec.chrom().to_string()).or_insert(InvertedRecordReference::new());
            let stat = chrom.bins.entry(region_to_bin_2(rec.start(), rec.end())).or_insert(InvertedRecordBuilder::new());
            // rec.chrom;
            stat.add(rec.start(), rec.end(), rec.name().unwrap_or("").to_string())
            //start.push(rec.start());
            //end.push(rec.end());
            //aux.push(rec.name().unwrap_or("").to_string());
            // aux.push(rec.name().join("\t"));
        }
        return InvertedRecordSet{chrom: inverted_record_set} //{start: start, end: end, name: aux}
    }
}

#[cfg(test)]
mod tests {
    use bio::io::bed;
    use super::InvertedRecordSet;
    #[test]
    fn it_works() {
        let example = b"1\t5\t5000\tname1\t0.5\n1\t5\t5000\tname1\t0.5";
        let mut reader = bed::Reader::new(&example[..]);
        let set = InvertedRecordSet::new(reader);
        println!("{:?}", set);

    }
}
