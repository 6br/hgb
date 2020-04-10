use bio::io::bed;
use std::io::{Result, Read, Write};
use std::collections::HashMap;
use std::cell::RefCell;
use std::path::Path;
use std::fs::File;
use crate::index::region_to_bin_2;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Chunk<T> {
    length: u32,
    sample_id: u32,
    sample_file_id: u32, // When we split files of inverted data structure
    // format_type: u32, // Enum
    data: T,
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
                    chunks.push(Chunk::<InvertedRecord>{length: 0, sample_id: set.sample_id, sample_file_id: sample_file_id as u32, data: InvertedRecord::new(bin)})
                }
            }    
        }
        return InvertedRecordEntire{chrom: inverted_record}
    }

    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<Index> {
        /* Unimplemented */
        let index = Index::new();
        for (name, chromosome) in &self.chrom　{
            for (bin_id, bin) in &chromosome.bins {
                
            }
        }

        Ok(())
    }
}

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
pub struct InvertedRecord {
    start: Vec<u64>,
    end: Vec<u64>,
    name: Vec<String>, 
    // aux: Vec<String> // We should update better data format.
}

impl InvertedRecord {
    pub fn new(builder: &InvertedRecordBuilder) -> InvertedRecord {
        /* TODO() Use Bit packing for converting RefCell to Just vector, however now we just clone */
        let start = builder.start.clone().into_inner();
        let end = builder.end.clone().into_inner();
        let name = builder.name.clone().into_inner();
        InvertedRecord{start: start, end: end, name: name}
    }

    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
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
            let stat = chrom.bins.entry(region_to_bin_2(rec.start(), rec.end())).or_insert(InvertedRecordBuilder::new());
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
        let set = InvertedRecordSet::new(reader);
        println!("{:?}", set);

    }
}
