use bio::io::bed;
use std::cell::RefCell;
use std::collections::{HashMap, BTreeMap};
use std::io::Read;
use crate::{alignment::Set, index::region_to_bin_3, header::Header};
use bam::header::HeaderEntry;
use crate::range::InvertedRecord;
use crate::alignment::Bins;
use crate::range::Format;
use crate::Builder;

#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct InvertedRecordBuilder {
    pub start: RefCell<Vec<u64>>,
    pub end: RefCell<Vec<u64>>,
    pub name: RefCell<Vec<String>>, 
    pub aux: RefCell<Vec<Vec<String>>>, //Just a tab-separated string here.
}

impl Bins<InvertedRecordBuilder> {
    pub fn new() -> Self {
        Bins {
            bins: HashMap::new()
        }
    }
}


/// Builder of the [InvertedRecord](struct.InvertedRecord.html).
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

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordBuilderSet {
    pub sample_id: u64,
    pub chrom: BTreeMap<String, InvertedRecordReference> // Mutex?
}

impl InvertedRecordBuilderSet {
    pub fn new<R:Read>(mut reader: bed::Reader<R>, sample_id: u64) -> Self {
        let mut inverted_record_set = BTreeMap::new();

        for record in reader.records() {
            let rec = record.ok().expect("Error reading record.");
            let chrom = inverted_record_set.entry(rec.chrom().to_string()).or_insert(InvertedRecordReference::new());
            let stat = chrom.bins.entry(region_to_bin_3(rec.start(), rec.end())).or_insert(InvertedRecordBuilder::new());

            let mut aux = vec![];
            let mut n = 4; // Ignore name field
            while let Some(item) = rec.aux(n){
                aux.push(item.to_string());
                n += 1;
            }
            stat.add(rec.start(), rec.end(), rec.name().unwrap_or("").to_string(), aux)
        }
        return InvertedRecordBuilderSet{sample_id: sample_id, chrom: inverted_record_set}
    }
}

impl Set<InvertedRecordBuilder> {
    pub fn new<R:Read>(mut reader: bed::Reader<R>, sample_id: u64, header: &mut Header) -> Result<Self, String> {
        let mut inverted_record_set = BTreeMap::new();

        for record in reader.records() {
            let rec = record.ok().expect("Error reading record.");
            let mut chrom_id = header.reference_id(rec.chrom());
            if let None = chrom_id {
                let chrom_item = HeaderEntry::ref_sequence(rec.chrom().to_string(), u32::max_value());
                header.push_entry(chrom_item)?;
                chrom_id = header.reference_id(rec.chrom());
            }
            
            let chrom = inverted_record_set.entry(chrom_id.unwrap()).or_insert(Bins::<InvertedRecordBuilder>::new());
            let stat = chrom.bins.entry(region_to_bin_3(rec.start(), rec.end())).or_insert(InvertedRecordBuilder::new());


            let mut aux = vec![];
            let mut n = 4; // Ignore name field
            while let Some(item) = rec.aux(n){
                aux.push(item.to_string());
                n += 1;
            }
            stat.add(rec.start(), rec.end(), rec.name().unwrap_or("").to_string(), aux)
        }
        Ok(Set::<InvertedRecordBuilder>{sample_id: sample_id, chrom: inverted_record_set, unmapped: InvertedRecordBuilder::new()})
    }
}

impl Builder for InvertedRecordBuilder {
    fn to_format(self) -> Format {
        Format::Range(self.to_record())
    }
}

impl InvertedRecordBuilder {
    pub fn to_record(self) -> InvertedRecord {
        InvertedRecord::from_builder(&self)
        //InvertedRecord { data: self.take() }
    }
}