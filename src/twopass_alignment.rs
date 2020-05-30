use crate::{checker_index::Reference, range::Format, Builder, ColumnarSet};
use crate::{
    header::Header,
    index::Region,
    range::{Bins, Set},
};
use bam::bgzip::ReadBgzip;
use bam::record::Record;
use bam::RecordReader;
use bam::{
    index::{Chunk, VirtualOffset},
    RecordWriter,
};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use std::{
    collections::{BTreeMap, HashMap},
    io::{Read, Result},
};
/*
#[derive(Eq, PartialEq, Debug)]
pub enum Alignment {
    Offset(String, Vec<Chunk>),
    Object(Vec<Record>)
}
*/
#[derive(Eq, PartialEq, Debug)]
pub struct Alignment {
    pub data: Vec<Chunk>, // bgzip compressed records
    reader: String,       // TODO()Path should be a reference to bam::Reader.
}

/// BAM-Compatible Alignment Record
pub struct AlignmentBuilder {
    alignments: Vec<Chunk>,
    reader: String,
}

impl AlignmentBuilder {
    pub fn new() -> Self {
        AlignmentBuilder {
            alignments: vec![],
            reader: "".to_string(),
        }
    }
    pub fn add(&mut self, alignment: Chunk, reader: &str) {
        self.alignments.push(alignment);
        self.reader = reader.to_string();
    }

    pub fn take(self) -> Vec<Chunk> {
        self.alignments
    }
}

impl Bins<AlignmentBuilder> {
    pub fn new() -> Self {
        Bins {
            bins: HashMap::new(),
            reference: Reference::new_with_bai_half_overlapping(),
        }
    }
    pub fn new_from_reference(reference: Reference) -> Self {
        Bins {
            bins: HashMap::new(),
            reference,
        }
    }
}

impl Set<AlignmentBuilder> {
    pub fn new<R: Read>(
        mut reader: bam::BamReader<R>,
        sample_id: u64,
        header: &mut Header,
    ) -> Self {
        let mut chrom = BTreeMap::new();
        let mut unmapped = AlignmentBuilder::new();
        let mut prev = reader.reader.current().and_then(|t| t.offset()).unwrap();
        let mut rec = Record::new();

        //for record in reader {
        while let Ok(true) = reader.read_into(&mut rec) {
            //.ok().expect("Error reading record.");
            if rec.ref_id() >= 0 {
                let chrom_len = header.reference_len(rec.ref_id() as u64).unwrap();
                let reference = Reference::new_from_len(chrom_len);
                let bin = chrom
                    .entry(rec.ref_id() as u64)
                    .or_insert(Bins::<AlignmentBuilder>::new_from_reference(reference));
                if rec.start() > 0 && rec.calculate_end() > 0 {
                    let bin_id = bin.reference.region_to_bin(Region::new(
                        rec.ref_id() as u64,
                        rec.start() as u64,
                        rec.calculate_end() as u64,
                    ));
                    let stat = bin
                        .bins
                        .entry(bin_id as u32)
                        .or_insert(AlignmentBuilder::new());
                    let end: u64 = reader.reader.current().and_then(|t| t.offset()).unwrap();
                    stat.add(
                        Chunk::new(VirtualOffset::from_raw(prev), VirtualOffset::from_raw(end)),
                        header.get_name(sample_id as usize).unwrap(),
                    );
                    prev = end;
                }
            } else if rec.ref_id() == -1 {
                let end: u64 = reader.reader.current().and_then(|t| t.offset()).unwrap();
                unmapped.alignments.push(Chunk::new(
                    VirtualOffset::from_raw(prev),
                    VirtualOffset::from_raw(end),
                ));
            } else {
                // return Err(io::Error::new(InvalidData, "Reference id < -1"));
                panic!("Reference id < -1");
            }
        }
        Set::<AlignmentBuilder> {
            sample_id,
            chrom,
            unmapped,
        }
    }
}

impl AlignmentBuilder {
    pub fn to_record(self) -> Alignment {
        Alignment {
            reader: self.reader.clone(),
            data: self.take(),
        }
    }
}
/*
impl Builder for AlignmentBuilder {
    fn to_format(self) -> Format {
        Format::Alignment(self.to_record())
    }
}*/
impl ColumnarSet for Alignment {
    fn new() -> Alignment {
        Alignment {
            data: vec![],
            reader: "".to_string(),
        }
    }
    fn to_stream<W: std::io::Write>(&self, stream: &mut W) -> Result<()> {
        stream.write_u64::<LittleEndian>(self.data.len() as u64)?;
        let header = bam::Header::new();
        let mut writer = bam::bam_writer::BamWriterBuilder::new()
            .additional_threads(4)
            .write_header(false)
            .from_stream(stream, header)?;
        // TODO() Inject the number of threads.
        let mut reader = bam::IndexedReader::from_path(&self.reader)?;
        let viewer = reader.chunk(self.data.clone());
        for i in viewer {
            writer.write(&i?)?;
        }
        writer.finish()?;
        Ok(())
    }
    fn from_stream<R: Read>(&mut self, stream: &mut R) -> Result<bool> {
        todo!("Unimplemented")
        /*
        let len = stream.read_u64::<LittleEndian>()?;

        let mut reader =
            bam::BamReader::from_stream_no_header(stream, bam::Header::new(), 0).unwrap();

        for _i in 0..len as usize {
            let mut record = Record::new();
            reader.read_into(&mut record)?;
            self.data.push(record);
        }
        let mut _record = Record::new();
        // For consuming all bgzip blocks.
        let _ = reader.next();
        let _ = reader.next();

        Ok(true)
        */
    }
}
