use std::cell::RefCell;
use crate::compression::{IntegerEncode, StringEncode, DeltaVByte, VByte, Deflate, integer_encode_wrapper, string_encode, integer_decode, string_decode};
use bam::record::Record;
// use sugars::refcell;

/// BAM-Compatible Alignment Inverted-Record
#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct Alignment {
    start: DeltaVByte,
    end: VByte,
    // ref_id: u64,
    mate_ref_id: VByte,
    mate_start: VByte, // If the data is almost absent, how to encode them?
    cigar: Vec<u8>,
    seq: Vec<u8>, // Should sequence store as difference encoding?
    quals: Vec<u8>,
    tags: Vec<u8>,
}

/// BAM-Compatible Alignment Record
#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct AlignmentBuilder {
    pub ref_id: Option<u64>,
    pub mate_ref_id: RefCell<Vec<i64>>, // We need to know whether cell is enough.
    pub start: RefCell<Vec<u64>>,
    pub end: RefCell<Vec<u64>>,
    pub mate_start: RefCell<Vec<i32>>,
//    pub bin: RefCell<Vec<u32>> // We have different binning system.
    pub cigar: RefCell<Vec<u8>>, //We regard them as a string at first.
    pub seq: RefCell<Vec<u8>>, //We regard them as a string at first.
    pub quals: RefCell<Vec<u8>>, //We regard them as a string at first.
    pub tags: RefCell<Vec<Vec<String>>>, //Just a tab-separated string here.
}

impl AlignmentBuilder {
    pub fn new(ref_id: Option<u64>) -> Self {
        AlignmentBuilder{ ref_id, mate_ref_id: RefCell::new(vec![]), start: RefCell::new(vec![]), end: RefCell::new(vec![]), mate_start: RefCell::new(vec![]), cigar: RefCell::new(vec![]), seq: RefCell::new(vec![]), quals: RefCell::new(vec![]), tags: RefCell::new(vec![])}
    }
    pub fn add(&self, alignment: &Record) {
        self.start.borrow_mut().push(alignment.start() as u64);
        self.end.borrow_mut().push(alignment.end() as u64);
        self.mate_start.borrow_mut().push(alignment.mate_start());
        self.cigar.borrow_mut().push(alignment.cigar());

    }
}

pub impl Alignment {
    pub fn from_builder(builder: &AlignmentBuilder) -> Alignment {
        Alignment {
            start: integer_encode_wrapper(&builder.start.borrow(), true),
            end: integer_encode_wrapper(&builder.end.borrow(), false),
            mate_ref_id: (),
            mate_start: (),
            cigar: (),
            seq: (),
            quals: (),
            tags: (), 
        }
    }

    pub fn to_record(self) -> Vec<Record> {

    }
}