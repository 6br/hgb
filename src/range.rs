use std::io::{Result, Read, Write, Seek, Error};
use std::io::ErrorKind::InvalidData;
use std::collections::BTreeMap;
use bio::io::bed;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use crate::index::{Index, Reference, Bin};
use crate::{ColumnarSet, ChunkWriter};
use crate::binary::GhbWriter;
use crate::compression::{IntegerEncode, StringEncode, DeltaVByte, VByte, Deflate, integer_encode_wrapper, string_encode, integer_decode, string_decode};
use bam::header::HeaderEntry;
use crate::header::Header;
use crate::alignment::{Set, Alignment};
use crate::{Builder, builder::{InvertedRecordBuilder, InvertedRecordBuilderSet}};

type Chromosome = Vec<HeaderEntry>;

#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Format {
    Default(Default),
    Range(InvertedRecord),
    Alignment(Alignment),
}

impl Format {
    pub fn id(format: &Format) -> u32 {
        match format {
            Format::Default(_) => 0,
            Format::Range(_) => 1,
            Format::Alignment(_) => 2
        }
    } 
}

/// GHB Record.
///
///https://qiita.com/mhgp/items/41a75915413aec781fe0
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Record {
    sample_id: u64,
    sample_file_id: u32, // When we split files of inverted data structure
    format: u32, // Enum
    data: Format,
}

impl Record {
    /// Creates an empty record.
    pub fn new() -> Self {
        Record {
            sample_id: 0,
            sample_file_id: 0,
            format: 0,
            data: Format::Default(Default{})
        }
    }
    pub fn data(self) -> Format {
        self.data
    }
    /// Clears the record.
    pub fn clear(&mut self) {
        self.sample_id = 0;
        self.sample_file_id = 0;
        self.format = 0;
        self.data = Format::Default(Default{});
    }
    /// Returns 0-based sample id.
    pub fn sample_id(&self) -> u64 {
        return self.sample_id;
    }
    /// Returns 0-based sample file id.
    pub fn sample_file_id(&self) -> u32 {
        return self.sample_file_id;
    }

    /// Writes a record in binary format.
    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        stream.write_u64::<LittleEndian>(self.sample_id)?;
        stream.write_u32::<LittleEndian>(self.sample_file_id)?;
        stream.write_u32::<LittleEndian>(self.format)?;

        match &self.data {
            Format::Default(data) => data.to_stream(stream)?,
            Format::Range(data) => data.to_stream(stream)?,
            Format::Alignment(data) => data.to_stream(stream)?
        }
        Ok(())
    }

    /// Fills the record from a `stream` of uncompressed binary contents.
    pub fn from_stream<R: Read, T: ColumnarSet>(&self, stream: &mut R) -> Result<Self> {
        let sample_id = stream.read_u64::<LittleEndian>()?;
        let sample_file_id = stream.read_u32::<LittleEndian>()?;
        let format = stream.read_u32::<LittleEndian>()?;
        
        let data = match format {
            0 => {
                let mut data = Default::new();
                data.from_stream(stream)?;
                Format::Default(data)
            },
            1 => {
                let mut record = InvertedRecord::new();
                record.from_stream(stream)?;
                Format::Range(record)
            },
            2 => {
                let mut record = Alignment::new();
                record.from_stream(stream)?;
                Format::Alignment(record)
            }
            _ => panic!("Panic!")
        };
        return Ok(Record{sample_id, sample_file_id,format, data})
    }

    /// Fills the record from a `stream` of uncompressed BAM contents.
    pub(crate) fn fill_from_bam<R: Read>(&mut self, stream: &mut R) -> Result<bool> {
        self.sample_id = stream.read_u64::<LittleEndian>()?;
        self.sample_file_id = stream.read_u32::<LittleEndian>()?;
        self.format = stream.read_u32::<LittleEndian>()?;
        let data = match self.format {
            0 => {
                let mut data = Default::new();
                data.from_stream(stream)?;
                Format::Default(data)
            },
            1 => {
                let mut record = InvertedRecord::new();
                record.from_stream(stream)?;
                Format::Range(record)
            },
            2 => {
                let mut record = Alignment::new();
                record.from_stream(stream)?;
                Format::Alignment(record)
            }
            _ => panic!("Panic!")
        };
        self.data = data;
        return Ok(true);
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordChromosome {
    bins: BTreeMap<u32, Vec<Record>> // Mutex?
}

#[derive(Clone, Debug)]
pub struct InvertedRecordEntire {
    chrom: Vec<InvertedRecordChromosome>, // Mutex?
    unmapped: Vec<Record>,
    chrom_table: Chromosome,
    sample_file_id_max: usize,
}


impl InvertedRecordEntire {
    pub fn add<T: Builder>(&mut self, sample_file: Set<T>) {
        let sample_file_id = self.sample_file_id_max;
        self.sample_file_id_max += 1;
        for (id, chromosome) in sample_file.chrom {
            for (bin_id, bin) in chromosome.bins {
                if let None = self.chrom.get(id as usize) {
                    self.chrom.resize((id + 1) as usize,  InvertedRecordChromosome{bins: BTreeMap::new() });
                }
                let chunks = self.chrom[id as usize].bins.entry(bin_id).or_insert(Vec::new());
                let data = bin.to_format();
                chunks.push(
                    Record{
                        sample_id: sample_file.sample_id, sample_file_id: sample_file_id as u32, 
                        format: Format::id(&data), data: data
                    }
                )
            }
        }
        let data = sample_file.unmapped.to_format();
        let unmapped =
            Record{
                sample_id: sample_file.sample_id, sample_file_id: sample_file_id as u32, 
                format: Format::id(&data), data: data
            };
        self.unmapped.push(unmapped);
    }
    pub fn new_from_set<T: Builder>(sample_file_list: Vec<Set<T>>) -> Self {
        let mut inverted_record = vec![];
        let chrom_table = vec![];
        let sample_len = sample_file_list.len();
        let mut unmapped_list = Vec::with_capacity(sample_len);
        for (sample_file_id, set) in sample_file_list.into_iter().enumerate() {
            for (_name, chromosome) in set.chrom {
                let mut chrom = InvertedRecordChromosome{bins: BTreeMap::new() };

                for (bin_id, bin) in chromosome.bins {
                    let chunks = chrom.bins.entry(bin_id).or_insert(Vec::new());
                    let data = bin.to_format();
                    chunks.push(
                        Record{
                            sample_id: set.sample_id, sample_file_id: sample_file_id as u32, 
                            format: Format::id(&data), data: data
                        })
                }
                inverted_record.push(chrom);
            }    
            let data = set.unmapped.to_format();
            let unmapped =
                Record{
                    sample_id: set.sample_id, sample_file_id: sample_file_id as u32, 
                    format: Format::id(&data), data: data
                };
                unmapped_list.push(unmapped);
        }

        InvertedRecordEntire{chrom: inverted_record, unmapped: unmapped_list, chrom_table, sample_file_id_max: sample_len}
    }
    pub fn chrom_table(self) -> Chromosome {
        self.chrom_table
    }
    pub fn write_header(self, header: &mut Header) {
        for i in self.chrom_table {
            header.push_entry(i).unwrap();
        }
    }
    pub fn new(sample_file_list: Vec<InvertedRecordBuilderSet>) -> InvertedRecordEntire {
        let mut inverted_record = vec![];
        let mut chrom_table = vec![];
        for (sample_file_id, set) in sample_file_list.iter().enumerate() {
            for (name, chromosome) in &set.chrom {
                let mut chrom = InvertedRecordChromosome{bins: BTreeMap::new() };

                let chrom_item = HeaderEntry::ref_sequence(name.clone(), i32::max_value() as u32);
                chrom_table.push(chrom_item);
                for (bin_id, bin) in &chromosome.bins {
                    let chunks = chrom.bins.entry(*bin_id).or_insert(Vec::new());
                    chunks.push(Record{sample_id: set.sample_id, sample_file_id: sample_file_id as u32, format: 1, data: Format::Range(InvertedRecord::from_builder(bin))})
                }
                inverted_record.push(chrom);
            }    
        }
        return InvertedRecordEntire{chrom: inverted_record, unmapped: vec![], chrom_table,sample_file_id_max: sample_file_list.len()}
    }

    pub fn write_binary<W: Write+Seek>(&self, writer: &mut GhbWriter<W>) -> Result<Index> {
        let mut reference = vec![];
        for chromosome in &self.chrom {
            let mut bins = BTreeMap::new();
            for (bin_id, bin) in &chromosome.bins {
                let mut records = vec![];
                for chunk in bin {
                    let record = writer.write(&chunk)?;
                    records.push(record);
                }
                bins.insert(*bin_id, Bin::new(*bin_id, records));
            }
            reference.push(Reference::new(bins));
        }
        Ok(Index::new(reference))
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct Default {

}

impl ColumnarSet for Default {
    fn new() -> Self {
        Default{}
    }
    fn to_stream<W: Write>(&self, _stream: &mut W) -> Result<()> {
        Err(Error::new(InvalidData, format!("No data.")))
    }

    fn from_stream<R: Read>(&mut self, _stream: &mut R) -> Result<bool> {
        Err(Error::new(InvalidData, format!("No data.")))
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct InvertedRecord {
    start: DeltaVByte, //Vec<u8>, //IntegerEncode::DeltaVByte // Vec<u64>,
    end: VByte, //Vec<u8>, //IntegerEncode::VByte,　// Vec<u64>
    name: Deflate, //Vec<u8>, //StringEncode::Deflate,  // Vec<String>
    aux: Deflate //Vec<u8>, //StringEncode::Deflate // We should update better data format.
}

impl ColumnarSet for InvertedRecord {
    fn new() -> InvertedRecord {
        let start = vec![];
        let end = vec![];
        let name = vec![];
        let aux = vec![];
        InvertedRecord{start: start, end: end, name: name, aux: aux}
    }
    
    fn from_stream<R: Read>(&mut self, stream: &mut R) -> Result<bool> {
        let _n_items = stream.read_u64::<LittleEndian>()?;
        let n_header = stream.read_u64::<LittleEndian>()?;
        for _i in 0..n_header {
            let _a = stream.read_u16::<LittleEndian>()?;
        }
        self.start = vec![];
        self.end = vec![];
        self.name = vec![];

/*        let mut start = Vec::with_capacity(n_items as usize);
        let mut end = Vec::with_capacity(n_items as usize);
        let mut name = Vec::with_capacity(n_items as usize);*/
        let mut n_items = stream.read_u64::<LittleEndian>()?;
        for _i in 0..n_items {
//            self.start.push(stream.read_u64::<LittleEndian>()?);
            self.start.push(stream.read_u8()?);
        }
        n_items = stream.read_u64::<LittleEndian>()?;
        for _i in 0..n_items {
//            self.end.push(stream.read_u64::<LittleEndian>()?);
            self.end.push(stream.read_u8()?);
        }
        n_items = stream.read_u64::<LittleEndian>()?;
        for _i in 0..n_items {
//            let size = stream.read_u64::<LittleEndian>()?;
//            let mut raw = Vec::with_capacity(size as usize);
//            stream.take(size).read_to_end(&mut raw);
//            self.name.push(String::from_utf8(raw).unwrap());
            self.name.push(stream.read_u8()?);
        }
        n_items = stream.read_u64::<LittleEndian>()?;
        for _i in 0..n_items {
//            let size = stream.read_u64::<LittleEndian>()?;
//            let mut raw = Vec::with_capacity(size as usize);
//            stream.take(size).read_to_end(&mut raw);
//            self.aux.push(String::from_utf8(raw).unwrap());
            self.aux.push(stream.read_u8()?);
        }

        Ok(true) //Ok(InvertedRecord{start: start, end: end, name: name})
    }
    
    fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        stream.write_u64::<LittleEndian>(self.start.len() as u64)?; //n_item
        stream.write_u64::<LittleEndian>(3 as u64)?; // n_header
        stream.write_u16::<LittleEndian>(1 as u16)?; // u64
        stream.write_u16::<LittleEndian>(1 as u16)?; // u64 
        stream.write_u16::<LittleEndian>(0 as u16)?; // String

        stream.write_u64::<LittleEndian>(self.start.len() as u64)?;
        for i in &self.start {
            stream.write_u8(*i)?;
        }
        stream.write_u64::<LittleEndian>(self.end.len() as u64)?;
        for i in &self.end {
            stream.write_u8(*i)?;
        }
        stream.write_u64::<LittleEndian>(self.name.len() as u64)?;
        for i in &self.name {
//            stream.write_u64::<LittleEndian>(i.len() as u64)?;
//            stream.write_all(i.as_bytes())?
            stream.write_u8(*i)?;
        }
        stream.write_u64::<LittleEndian>(self.aux.len() as u64)?;
        for i in &self.aux {
//            stream.write_u64::<LittleEndian>(i.len() as u64)?;
//            stream.write_all(i.as_bytes())?
            stream.write_u8(*i)?;
        }
        Ok(())
    }
}


impl InvertedRecord {
    /// Construct from InvertedRecordBuilder.
    pub fn from_builder(builder: &InvertedRecordBuilder) -> InvertedRecord {
        let start = integer_encode_wrapper(&builder.start.borrow(), true);
        let end = integer_encode_wrapper(&builder.end.borrow(), false);
        let name = string_encode(&builder.name.borrow());
        let aux_raw = builder.aux.clone().into_inner().into_iter().map(|t| t.join("\t").to_owned()).collect();
        let aux = string_encode(&aux_raw);
        InvertedRecord{start: start, end: end, name: name, aux: aux}
    }
/*
    pub fn from_builder(builder: &InvertedRecordBuilder) -> InvertedRecord {
        /* TODO() Use Bit packing for converting RefCell to Just vector, however now we just clone */
        let start = builder.start.clone().into_inner();
        let end = builder.end.clone().into_inner();
        let name = builder.name.clone().into_inner();
        let aux = builder.aux.clone().into_inner().into_iter().map(|t| t.join("\t").to_owned()).collect();
        InvertedRecord{start: start, end: end, name: name, aux: aux}
    }

    pub fn from_builder_packing(builder: &InvertedRecordBuilder) -> InvertedRecord {
        /* TODO() Use Bit packing for converting RefCell to Just vector, however now we just clone */
        // let start = 
        let bitpacker = BitPacker4x::new();
        let num_bits: u8 = bitpacker.num_bits(&builder.start.into_inner());
        let mut start = vec![0u8; 4 * BitPacker4x::BLOCK_LEN];
        let _compressed_len = bitpacker.compress(&builder.start.into_inner(), &mut start[..], num_bits);
        // assert_eq!((num_bits as usize) *  BitPacker4x::BLOCK_LEN / 8, compressed_len);

        let start = builder.start.clone().into_inner();
        let end = builder.end.clone().into_inner();
        let name = builder.name.clone().into_inner();
        let aux = builder.aux.clone().into_inner().into_iter().map(|t| t.join("\t").to_owned()).collect();
        InvertedRecord{start: start, end: end, name: name, aux: aux}
    }
*/

    /// Output as Record.
    pub fn to_record(self, chromosome: &str) -> Vec<bed::Record> {
        let mut records = vec![];
        let start = integer_decode(IntegerEncode::DeltaVByte(self.start));
        let end = integer_decode(IntegerEncode::VByte(self.end));
        let name = string_decode(&StringEncode::Deflate(self.name));
        let aux = string_decode(&StringEncode::Deflate(self.aux));
        let len = start.len();
        for i in 0..len {
            let mut rec = bed::Record::new();
            rec.set_chrom(chromosome);
            rec.set_start(start[i]);
            rec.set_end(end[i]);
            rec.set_name(&name[i]);
            for aux in aux[i].split("\t") {
                rec.push_aux(aux);
            }
            records.push(rec)
        }
        records
    }
}

