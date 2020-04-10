use std::io::{Result, Read, Write};
use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use crate::index::Index;
use crate::ColumnarSet;
use crate::builder::{InvertedRecordBuilder, InvertedRecordSet};

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Chunk<T: ColumnarSet> {
    /*length: u32,*/
    sample_id: u32,
    sample_file_id: u32, // When we split files of inverted data structure
    // format_type: u32, // Enum
    data: T,
}

impl<T: ColumnarSet> Chunk<T> {
    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        stream.write_u32::<LittleEndian>(self.sample_id)?;
        stream.write_u32::<LittleEndian>(self.sample_file_id)?;
        self.data.to_stream(stream)?;
        Ok(())
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordChromosome {
    bins: HashMap<u32, Vec<Chunk<InvertedRecord>>> // Mutex?
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordEntire {
    chrom: HashMap<String, InvertedRecordChromosome>, // Mutex?
}

impl InvertedRecordEntire {
    pub fn new(sample_file_list: Vec<InvertedRecordSet>) -> InvertedRecordEntire {
        let mut inverted_record = HashMap::new();
        for (sample_file_id, set)  in sample_file_list.iter().enumerate() {
            for (name, chromosome) in &set.chrom {
                let chrom = inverted_record.entry(name.clone()).or_insert(InvertedRecordChromosome{bins: HashMap::new() }); // TODO() DO NOT USE CLONE
                for (bin_id, bin) in &chromosome.bins {
                    let chunks = chrom.bins.entry(*bin_id).or_insert(Vec::new());
                    chunks.push(Chunk::<InvertedRecord>{sample_id: set.sample_id, sample_file_id: sample_file_id as u32, data: InvertedRecord::new(bin)})
                }
            }    
        }
        return InvertedRecordEntire{chrom: inverted_record}
    }

    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<Index> {
        /* Unimplemented */
        let index = Index{references: vec![]};
        for (_name, chromosome) in &self.chrom {
            for (_bin_id, bin) in &chromosome.bins {
                for chunk in bin {
                    chunk.to_stream(stream)?;
                }
            }
        }

        Ok(index)
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct InvertedRecord {
    start: Vec<u64>,
    end: Vec<u64>,
    name: Vec<String>, 
    // aux: Vec<String> // We should update better data format.
}

impl ColumnarSet for InvertedRecord {
    fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        stream.write_u64::<LittleEndian>(self.start.len() as u64)?;
        stream.write_u64::<LittleEndian>(3 as u64)?;
        stream.write_u16::<LittleEndian>(1 as u16)?; // u64
        stream.write_u16::<LittleEndian>(1 as u16)?; // u64 
        stream.write_u16::<LittleEndian>(0 as u16)?; // String

        for i in &self.start {
            stream.write_u64::<LittleEndian>(*i)?;
        }
        for i in &self.end {
            stream.write_u64::<LittleEndian>(*i)?;
        }
        for i in &self.name {
            stream.write_u64::<LittleEndian>(i.len() as u64)?;
            stream.write_all(i.as_bytes())?
        }

        Ok(())
    }
}

impl InvertedRecord {
    pub fn new(builder: &InvertedRecordBuilder) -> InvertedRecord {
        /* TODO() Use Bit packing for converting RefCell to Just vector, however now we just clone */
        let start = builder.start.clone().into_inner();
        let end = builder.end.clone().into_inner();
        let name = builder.name.clone().into_inner();
        InvertedRecord{start: start, end: end, name: name}
    }

    pub fn from_stream<R: Read>(mut stream: R) -> Result<InvertedRecord> {
        let n_header = stream.read_u64::<LittleEndian>()?;
        let n_items = stream.read_u64::<LittleEndian>()?;
        for _i in 0..n_header {
            let _ = stream.read_u16::<LittleEndian>()?;
        }
        let mut start = Vec::with_capacity(n_items as usize);
        let mut end = Vec::with_capacity(n_items as usize);
        let mut name = Vec::with_capacity(n_items as usize);
        for _i in 0..n_items {
            start.push(stream.read_u64::<LittleEndian>()?);
        }
        for _i in 0..n_items {
            end.push(stream.read_u64::<LittleEndian>()?);
        }
        for _i in 0..n_items {
            let size = stream.read_u64::<LittleEndian>()?;
            let mut raw = Vec::with_capacity(size as usize);
            stream.read_exact(&mut raw)?;
            name.push(String::from_utf8(raw).unwrap());
        }
        return Ok(InvertedRecord{start: start, end: end, name: name})
    }

    /// Loads index from path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<InvertedRecord> {
        let f = File::open(&path)?;
        InvertedRecord::from_stream(f)
    }
}


