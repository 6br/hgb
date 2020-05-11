use std::cell::RefCell;
use bam::RecordWriter;
use bam::RecordReader;
use bam::record::Record;
use std::{collections::{BTreeMap, HashMap}, io::{Read, Result, BufWriter}};
use byteorder::{LittleEndian, WriteBytesExt, ReadBytesExt};
use byteorder::ByteOrder;
use crate::{index::region_to_bin_3, ColumnarSet, range::Format, Builder};

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

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Bins<T> {
    pub bins: HashMap<u32, T> // Bin id is regarded as u32 now.
}

impl Bins<AlignmentBuilder> {
    pub fn new() -> Self {
        Bins {
            bins: HashMap::new()
        }
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Set<T> {
    pub sample_id: u64,
    pub chrom: BTreeMap<u64, Bins<T>>, // Mutex?
    pub unmapped: T,
}

impl Set<AlignmentBuilder> {
    pub fn new<R: Read>(reader: bam::BamReader<R>, sample_id: u64) -> Self {
        let mut chrom = BTreeMap::new();
        let unmapped = AlignmentBuilder::new();

        for record in reader {
            let rec = record.ok().expect("Error reading record.");
            if rec.ref_id() >= 0 {
                let bin = chrom.entry(rec.ref_id() as u64).or_insert(Bins::<AlignmentBuilder>::new());
                if rec.start() > 0 && rec.calculate_end() > 0 {                
                    let stat = bin.bins.entry(region_to_bin_3(rec.start() as u64, rec.calculate_end() as u64)).or_insert(AlignmentBuilder::new());
                    stat.add(rec);
                }
            } else if rec.ref_id() == -1 {
                unmapped.add(rec);
            } else {
                // return Err(io::Error::new(InvalidData, "Reference id < -1"));
                panic!("Reference id < -1");
            }
        }
        Set::<AlignmentBuilder> { sample_id, chrom, unmapped}
    }
}

impl AlignmentBuilder {
    pub fn to_record(&self) -> Alignment {
        //let mut aln = Alignment::new();
        //aln.from_builder(self).unwrap();
        //aln
        Alignment::new_from_builder(self).unwrap()
    }
}
impl Builder for AlignmentBuilder {
    fn to_format(&self) -> Format {
        Format::Alignment(self.to_record())
    }
}

impl Alignment {
    
    pub fn new_from_builder(builder: &AlignmentBuilder) -> Result<Alignment> {
        let mut binary = vec![];
        binary.write_u64::<LittleEndian>(builder.alignments.borrow_mut().len() as u64)?;
        let header = bam::Header::new();
        let output = BufWriter::new(binary);
        let mut writer = bam::bam_writer::BamWriterBuilder::new().write_header(false).from_stream(output, header)?;
        for i in builder.alignments.borrow().iter() {
            // i.write_bam(&mut binary)?;
            writer.write(&i)?;
        }
        writer.finish()?;
        let inner = writer.take_stream().into_inner()?;

        Ok(Alignment {data: inner})
        
    }
    

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

#[cfg(test)]
mod tests {
    use crate::header::{Header};
    
    use crate::binary;
    use crate::writer::GhiWriter;
    use crate::reader::IndexedReader;
    use crate::range::{InvertedRecordEntire, Format, Record};
    use crate::alignment::Set;
    #[test]
    fn bam_works() {
        let bam_path = "./test/index_test.bam";
        let mut reader = bam::BamReader::from_path(bam_path, 4).unwrap();
        let set = Set::new(reader, 0 as u64);

        let set_vec = vec![set];
        let entire = InvertedRecordEntire::new(set_vec);
        let header = Header::new();
        let mut writer = binary::GhbWriter::build()
            .write_header(false)
            .from_path("./test/test_bam.ghb", header).unwrap();
        let index = entire.write_binary(&mut writer).unwrap();
        writer.flush().unwrap();
        println!("{}", index);
    }
}