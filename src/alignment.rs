use std::cell::RefCell;
use bam::record::Record;
use std::io::{Result, Read, Write, Seek, Error};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};


/// BAM-Compatible Alignment Inverted-Record
#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct Alignment {
    data: Vec<u8> //xz compressed records
}


/// BAM-Compatible Alignment Record
#[derive(Clone, Debug)]
pub struct AlignmentBuilder {
    alignments: RefCell<Vec<Record>>,
}

impl AlignmentBuilder {
    pub fn new(ref_id: Option<u64>) -> Self {
        AlignmentBuilder{ 
            alignments: RefCell::new(vec![])
        }
    }
    pub fn add(&self, alignment: Record) {
        self.alignments.borrow_mut().push(alignment);
    }
}

impl Alignment {
    pub fn from_builder(builder: &AlignmentBuilder) -> Result<Alignment> {
        let binary = vec![];
        binary.write_u64::<LittleEndian>(self.alignments.borrow_mut().len());
        for i in builder.alignments.borrow_mut().into_iter() {
            i.write_bam(&binary)?;
        }
        Ok(Alignment { alignments: binary})
    }

    pub fn to_record(self) -> Result<Vec<Record>> {
        let len = self.data.read_u64::<LittleEndian>()?;
        let mut records = Vec::with_capacity(len);
        let mut record = record::Record::new();
        record.fill_from_bam(&mut self.alignments)?
        
    }
}
