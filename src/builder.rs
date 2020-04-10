use bio::io::bed;
use crate::index::region_to_bin_3;
use std::cell::RefCell;
//use crate::range::InvertedRecord;
use std::collections::HashMap;
use std::io::{Read};

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

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordSet {
    sample_id: u32,
    chrom: HashMap<String, InvertedRecordReference> // Mutex? // Is it better to use u32 here?
}

impl InvertedRecordSet {
    pub fn new<R:Read>(mut reader: bed::Reader<R>, sample_id: u32) -> Self {
        let mut inverted_record_set = HashMap::new();
        /*let mut start: Vec<u64> = vec![];
        let mut end: Vec<u64> = vec![];
        let mut aux: Vec<String> = vec![];*/
        for record in reader.records() {
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
    use bio::io::bed;
    use super::InvertedRecordSet;
    #[test]
    fn it_works() {
        let example = b"1\t5\t5000\tname1\t0.5\n1\t5\t5000\tname1\t0.5";
        let mut reader = bed::Reader::new(&example[..]);
        let set = InvertedRecordSet::new(reader, 0 as u32);
        println!("{:?}", set);

    }
}