use bio::io::bed;
use crate::index::region_to_bin_3;
use std::cell::RefCell;
//use crate::range::InvertedRecord;
use std::collections::{HashMap, BTreeMap};
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
    pub name: RefCell<Vec<String>>, //Just a tab-separated string here.
    pub aux: RefCell<Vec<Vec<String>>>,
}

impl InvertedRecordBuilder {
    pub fn new() -> Self {
        InvertedRecordBuilder{
            start: RefCell::new(Vec::new()),
            end: RefCell::new(Vec::new()),
            name: RefCell::new(Vec::new()),
            aux: RefCell::new(Vec::new()),
        }
    }
    pub fn add(&self, start: u64, end: u64, name: String, aux: Vec<String>) {
        self.start.borrow_mut().push(start);
        self.end.borrow_mut().push(end);
        self.name.borrow_mut().push(name);
        self.aux.borrow_mut().push(aux);
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordSet {
    pub sample_id: u64,
    pub chrom: BTreeMap<String, InvertedRecordReference> // Mutex?
}

impl InvertedRecordSet {
    pub fn new<R:Read>(mut reader: bed::Reader<R>, sample_id: u64) -> Self {
        let mut inverted_record_set = BTreeMap::new();
        /*let mut start: Vec<u64> = vec![];
        let mut end: Vec<u64> = vec![];
        let mut aux: Vec<String> = vec![];*/
        for record in reader.records() {
            println!("{:?}", record);
            let rec = record.ok().expect("Error reading record.");
            let chrom = inverted_record_set.entry(rec.chrom().to_string()).or_insert(InvertedRecordReference::new());
            let stat = chrom.bins.entry(region_to_bin_3(rec.start(), rec.end())).or_insert(InvertedRecordBuilder::new());
            // rec.chrom;
            let mut aux = vec![];
            let mut n = 4; // Ignore name field
            while let Some(item) = rec.aux(n){
                aux.push(item.to_string());
                n += 1;
            }
            stat.add(rec.start(), rec.end(), rec.name().unwrap_or("").to_string(), aux)
            //start.push(rec.start());
            //end.push(rec.end());
            //aux.push(rec.name().unwrap_or("").to_string());
            // aux.push(rec.name().join("\t"));
        }
        return InvertedRecordSet{sample_id: sample_id, chrom: inverted_record_set} //{start: start, end: end, name: aux}
    }
}
