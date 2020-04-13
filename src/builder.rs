use bio::io::bed;
use crate::index::region_to_bin_3;
use std::cell::RefCell;
//use crate::range::InvertedRecord;
use std::collections::HashMap;
use std::io::{Read};

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordReference {
    pub bins: HashMap<u32, InvertedRecordBuilder> // Mutex?
}

impl InvertedRecordReference {
    pub fn new() -> Self {
        InvertedRecordReference{
            bins:  HashMap::new()
        }
    }
}


#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct InvertedRecordBuilder {
    pub start: RefCell<Vec<u64>>,
    pub end: RefCell<Vec<u64>>,
    pub name: RefCell<Vec<String>> //Just a tab-separated string here.
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

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordSet {
    pub sample_id: u32,
    pub chrom: HashMap<String, InvertedRecordReference> // Mutex? // Is it better to use u32 here?
}

impl InvertedRecordSet {
    pub fn new<R:Read>(mut reader: bed::Reader<R>, sample_id: u32) -> Self {
        let mut inverted_record_set = HashMap::new();
        /*let mut start: Vec<u64> = vec![];
        let mut end: Vec<u64> = vec![];
        let mut aux: Vec<String> = vec![];*/
        for record in reader.records() {
            println!("{:?}", record);
            let rec = record.ok().expect("Error reading record.");
            let chrom = inverted_record_set.entry(rec.chrom().to_string()).or_insert(InvertedRecordReference::new());
            let stat = chrom.bins.entry(region_to_bin_3(rec.start(), rec.end())).or_insert(InvertedRecordBuilder::new());
            // rec.chrom;
            stat.add(rec.start(), rec.end(), rec.name().unwrap_or("").to_string())
            //start.push(rec.start());
            //end.push(rec.end());
            //aux.push(rec.name().unwrap_or("").to_string());
            // aux.push(rec.name().join("\t"));
        }
        return InvertedRecordSet{sample_id: sample_id, chrom: inverted_record_set} //{start: start, end: end, name: aux}
    }
}
#[cfg(test)]
mod tests {
    use bam::header::{Header, HeaderEntry};
    use bio::io::bed;
    use std::io;
    use crate::binary;
    use crate::writer;
    use crate::reader::IndexedReader;
    use crate::range::InvertedRecordEntire;
    use crate::IndexWriter;
    use super::InvertedRecordSet;

    #[test]
    fn it_works() {
        let _example = b"1\t5\t5000\tname1\t0.5\n1\t5\t5000\tname1\t0.5";
        let path = "./test/test.bed";
        let mut reader = bed::Reader::from_file(path).unwrap();
        let set = InvertedRecordSet::new(reader, 0 as u32);
        println!("{:?}", set);
        let set_vec = vec![set];
        let entire = InvertedRecordEntire::new(set_vec);
        // let output = io::BufWriter::new(io::stdout());
        let header = Header::new();
        let mut writer = binary::GhbWriter::build()
            .write_header(false)
            .from_path("./test/test.ghb", header).unwrap();
        let index = entire.write_binary(&mut writer).unwrap();
        println!("{}", index);
        let header2 = Header::new();
        let mut writer = writer::GhbWriter::build().write_header(false).from_path("./test/test.ghb.ghi", header2).unwrap();
        writer.write(&index);

        let reader2 = IndexedReader::from_path("./test/test.ghb").unwrap();
        println!("{}", reader2.index());

        for record in reader2.fetch(&bam::Region::new(1, 1_000, 1_500)).unwrap() {
            let record = record.unwrap();
            println!("{}", record);
        }
        assert_eq!(write!(index), write!(reader2.index()));
    }
}