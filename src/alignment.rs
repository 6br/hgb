use std::cell::RefCell;
use bam::RecordWriter;
use bam::RecordReader;
use bam::record::Record;
use std::io::Result;
use byteorder::{LittleEndian, WriteBytesExt, ReadBytesExt};
use byteorder::ByteOrder;
use crate::ColumnarSet;

/// BAM-Compatible Alignment Inverted-Record
#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct Alignment {
    data: Vec<u8> //bgzip compressed records
}

/// BAM-Compatible Alignment Record
#[derive(Clone, Debug)]
pub struct AlignmentBuilder {
    alignments: RefCell<Vec<Record>>,
}

impl AlignmentBuilder {
    pub fn new() -> Self {
        AlignmentBuilder{ 
            alignments: RefCell::new(vec![])
        }
    }
    pub fn add(&self, alignment: Record) {
        self.alignments.borrow_mut().push(alignment);
    }
}

impl Alignment {
    /*
    pub fn new_from_builder(builder: &AlignmentBuilder) -> Result<Alignment> {
        let mut binary = vec![];
        binary.write_u64::<LittleEndian>(builder.alignments.borrow_mut().len() as u64);
        let header = bam::Header::new();
        let mut writer = bam::bam_writer::BamWriterBuilder::new().write_header(false).from_stream(&mut binary, header)?;
        for i in builder.alignments.borrow().iter() {
            // i.write_bam(&mut binary)?;
            writer.write(&i)?;
        }
        Ok(Alignment { data: binary })
    }
    */

    pub fn from_builder(mut self, builder: &AlignmentBuilder) -> Result<()> {
        self.data.write_u64::<LittleEndian>(builder.alignments.borrow_mut().len() as u64)?;
        let header = bam::Header::new();
        let mut writer = bam::bam_writer::BamWriterBuilder::new().write_header(false).from_stream(&mut self.data, header)?;
        for i in builder.alignments.borrow().iter() {
            // i.write_bam(&mut binary)?;
            writer.write(&i)?;
        }
        Ok(())
    }


    pub fn to_record(self) -> Result<Vec<Record>> {
        // let len = (&mut self.data).read_u64::<LittleEndian>()?;
        let len = LittleEndian::read_u64(self.data.as_ref());
        let records = Vec::with_capacity(len as usize);
        let mut reader = bam::BamReader::from_stream(self.data.as_ref() as &[u8], 4).unwrap();
        for _i in 0..len as usize {
            let mut record = Record::new();
            // record.fill_from_bam(&mut self.data)?;
            reader.read_into(&mut record)?;
            // records.push(record);
        }
        Ok(records)
    }
}

impl ColumnarSet for Alignment {
    fn new() -> Alignment {
        Alignment { data: vec![] }
    }
    fn to_stream<W: std::io::Write>(&self, stream: &mut W) -> Result<()> {
        stream.write_u64::<LittleEndian>(self.data.len() as u64)?;
        stream.write_all(&self.data)?;
        Ok(())
    }
    fn from_stream<R: std::io::Read>(&mut self, stream: &mut R) -> Result<bool> {
        let n_vec = stream.read_u64::<LittleEndian>()?;
        self.data = vec![];
        for _i in 0..n_vec {
            self.data.push(stream.read_u8()?);
        }
        Ok(true)
    }
    
}