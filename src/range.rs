use crate::binary::GhbWriter;
use crate::checker_index::{Index, Reference};
use crate::compression::{
    integer_decode, integer_encode_wrapper, string_decode, string_encode, Deflate, DeltaVByte,
    IntegerEncode, StringEncode, VByte,
};
use crate::header::Header;
use crate::index::Bin;
use crate::twopass_alignment::Alignment;
use crate::{builder::InvertedRecordBuilder, Builder};
use crate::{ChunkWriter, ColumnarSet};
use bam::{header::HeaderEntry, IndexedReader};
use bio::io::bed;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use log::debug;
use std::collections::{BTreeMap, HashMap};
use std::io::ErrorKind::InvalidData;
use std::io::{Error, Read, Result, Seek, Write};
use std::str::FromStr;

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
            Format::Alignment(_) => 2,
        }
    }
}

// Implement the trait
impl FromStr for Format {
    type Err = &'static str;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "default" => Ok(Format::Default(Default {})),
            "alignment" => Ok(Format::Alignment(Alignment::new())),
            "range" => Ok(Format::Range(InvertedRecord {
                aux: vec![],
                end: vec![],
                start: vec![],
                name: vec![],
            })),
            _ => Err("no match"),
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
    format: u32,         // Enum
    data: Format,
}

impl Record {
    /// Creates an empty record.
    pub fn new() -> Self {
        Record {
            sample_id: 0,
            sample_file_id: 0,
            format: 0,
            data: Format::Default(Default {}),
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
        self.data = Format::Default(Default {});
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
    pub fn to_stream<W: Write, R: Read + Seek>(
        &self,
        stream: &mut W,
        threads: u16,
        bam_reader: Option<&mut IndexedReader<R>>,
    ) -> Result<()> {
        stream.write_u64::<LittleEndian>(self.sample_id)?;
        stream.write_u32::<LittleEndian>(self.sample_file_id)?;
        stream.write_u32::<LittleEndian>(self.format)?;

        match &self.data {
            Format::Default(data) => data.to_stream(stream, threads, bam_reader)?,
            Format::Range(data) => data.to_stream(stream, threads, bam_reader)?,
            Format::Alignment(data) => data.to_stream(stream, threads, bam_reader)?,
        }
        Ok(())
    }

    /// Fills the record from a `stream` of uncompressed binary contents.
    pub fn from_stream<U: Read, T: ColumnarSet>(
        &self,
        stream: &mut U,
        threads: u16,
    ) -> Result<Self> {
        let sample_id = stream.read_u64::<LittleEndian>()?;
        let sample_file_id = stream.read_u32::<LittleEndian>()?;
        let format = stream.read_u32::<LittleEndian>()?;

        let data = match format {
            0 => {
                let mut data = Default::new();
                data.from_stream(stream, threads)?;
                Format::Default(data)
            }
            1 => {
                let mut record = InvertedRecord::new();
                record.from_stream(stream, threads)?;
                Format::Range(record)
            }
            2 => {
                let mut record = Alignment::new();
                record.from_stream(stream, threads)?;
                Format::Alignment(record)
            }
            _ => panic!("Panic!"),
        };
        return Ok(Record {
            sample_id,
            sample_file_id,
            format,
            data,
        });
    }

    /// Fills the record from a `stream` of uncompressed BAM contents.
    pub(crate) fn fill_from_bam<R: Read>(&mut self, stream: &mut R, threads: u16) -> Result<bool> {
        self.sample_id = stream.read_u64::<LittleEndian>()?;
        self.sample_file_id = stream.read_u32::<LittleEndian>()?;
        self.format = stream.read_u32::<LittleEndian>()?;
        let data = match self.format {
            0 => {
                let mut data = Default::new();
                data.from_stream(stream, threads)?;
                Format::Default(data)
            }
            1 => {
                let mut record = InvertedRecord::new();
                record.from_stream(stream, threads)?;
                Format::Range(record)
            }
            2 => {
                let mut record = Alignment::new();
                record.from_stream(stream, threads)?;
                Format::Alignment(record)
            }
            _ => return Err(Error::new(InvalidData, "Invalid format record")),
        };
        self.data = data;
        return Ok(true);
    }
}

#[derive(Clone, Debug)]
pub struct InvertedRecordChromosome {
    bins: BTreeMap<u32, Vec<Record>>, // Mutex?,
    reference: Reference,
}

#[derive(Clone, Debug)]
pub struct Bins<T> {
    pub bins: HashMap<u32, T>, // Bin id is regarded as u32 now.
    pub reference: Reference,
}

/// Set is a data structure for storing entire inverted data structure typed T.
/// #[derive(Clone)]
pub struct Set<T, R: Read + Seek> {
    pub sample_id: u64,
    pub chrom: BTreeMap<u64, Bins<T>>, // Mutex?
    pub unmapped: T,
    pub bam_reader: Option<IndexedReader<R>>,
}

pub struct InvertedRecordEntire<R: Read + Seek> {
    chrom: Vec<InvertedRecordChromosome>, // Mutex?
    unmapped: Vec<Record>,
    chrom_table: Chromosome,
    sample_file_id_max: usize,
    bam_reader: BTreeMap<u64, IndexedReader<R>>,
}

impl<R: Read + Seek> InvertedRecordEntire<R> {
    pub fn new() -> InvertedRecordEntire<R> {
        InvertedRecordEntire {
            chrom: vec![],
            unmapped: vec![],
            chrom_table: vec![],
            sample_file_id_max: 0,
            bam_reader: BTreeMap::new(),
        }
    }
    /// Add new inverted record from Set.
    pub fn add<T: Builder>(&mut self, sample_file: Set<T, R>) {
        let sample_file_id = self.sample_file_id_max;
        self.sample_file_id_max += 1;
        for (id, chromosome) in sample_file.chrom {
            for (bin_id, bin) in chromosome.bins {
                if let None = self.chrom.get(id as usize) {
                    self.chrom.resize(
                        (id + 1) as usize,
                        InvertedRecordChromosome {
                            bins: BTreeMap::new(),
                            reference: Reference::new_with_bai_half_overlapping(),
                        },
                    );
                }
                let chunks = self.chrom[id as usize]
                    .bins
                    .entry(bin_id)
                    .or_insert(Vec::new());
                let data = bin.to_format();
                chunks.push(Record {
                    sample_id: sample_file.sample_id,
                    sample_file_id: sample_file_id as u32,
                    format: Format::id(&data),
                    data: data,
                })
            }
        }
        let data = sample_file.unmapped.to_format();
        let unmapped = Record {
            sample_id: sample_file.sample_id,
            sample_file_id: sample_file_id as u32,
            format: Format::id(&data),
            data: data,
        };
        self.unmapped.push(unmapped);
        if let Some(bam_reader) = sample_file.bam_reader {
            self.bam_reader
                .insert(sample_file.sample_id as u64, bam_reader);
        }
    }
    pub fn add_reader(&mut self, index: u64, reader: IndexedReader<R>) {
        self.bam_reader.insert(index, reader);
    }
    /// Create an initial inverted record from Set.
    pub fn new_from_set<T: Builder>(sample_file_list: Vec<Set<T, R>>) -> Self {
        let mut inverted_record = vec![];
        let chrom_table = vec![];
        let sample_len = sample_file_list.len();
        let mut unmapped_list = Vec::with_capacity(sample_len);
        let mut bam_readers = BTreeMap::new();

        for (sample_file_id, set) in sample_file_list.into_iter().enumerate() {
            debug!("{:?}", sample_file_id);
            if let Some(bam_reader) = set.bam_reader {
                bam_readers.insert(sample_file_id as u64, bam_reader);
            }
            for (_name, chromosome) in set.chrom {
                let mut chrom = InvertedRecordChromosome {
                    bins: BTreeMap::new(),
                    reference: chromosome.reference,
                };

                for (bin_id, bin) in chromosome.bins {
                    let chunks = chrom.bins.entry(bin_id).or_insert(Vec::new());
                    let data = bin.to_format();
                    debug!("{:?} {:?} {:?}", bin_id, set.sample_id, sample_file_id);
                    chunks.push(Record {
                        sample_id: set.sample_id,
                        sample_file_id: sample_file_id as u32,
                        format: Format::id(&data),
                        data: data,
                    })
                }
                if let None = inverted_record.get(_name as usize) {
                    inverted_record.resize(
                        (_name + 1) as usize,
                        InvertedRecordChromosome {
                            bins: BTreeMap::new(),
                            reference: Reference::new_with_bai_half_overlapping(),
                        },
                    );
                }
                inverted_record[_name as usize] = chrom;
            }
            let data = set.unmapped.to_format();
            let unmapped = Record {
                sample_id: set.sample_id,
                sample_file_id: sample_file_id as u32,
                format: Format::id(&data),
                data: data,
            };
            unmapped_list.push(unmapped);
        }

        InvertedRecordEntire {
            chrom: inverted_record,
            unmapped: unmapped_list,
            chrom_table,
            sample_file_id_max: sample_len,
            bam_reader: bam_readers,
        }
    }
    pub fn chrom_table(self) -> Chromosome {
        self.chrom_table
    }
    pub fn write_header(self, header: &mut Header) {
        for i in self.chrom_table {
            header.push_entry(i).unwrap();
        }
    }

    pub fn write_binary<W: Write + Seek>(&mut self, writer: &mut GhbWriter<W>) -> Result<Index> {
        let mut references = vec![];
        for chromosome in &self.chrom {
            let mut reference = chromosome.reference.clone();
            // let mut bins = BTreeMap::new();
            for (bin_id, bin) in &chromosome.bins {
                let mut records = vec![];
                for chunk in bin {
                    //println!("{:?} {:?} {:?}", bin_id, chunk.sample_id, self.bam_reader.keys());
                    let record =
                        writer.write(&chunk, (self.bam_reader).get_mut(&chunk.sample_id))?;
                    //debug!("{:?} {:?}", bin_id, record);
                    records.push(record);
                }
                reference.update(*bin_id as usize, Bin::new(*bin_id, records));
            }
            debug!("{:?}", reference);
            references.push(reference);
        }
        Ok(Index::new(references))
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct Default {}

impl ColumnarSet for Default {
    fn new() -> Self {
        Default {}
    }
    fn to_stream<W: Write, R: Read + Seek>(
        &self,
        _stream: &mut W,
        _threads: u16,
        _bam_reader: Option<&mut IndexedReader<R>>,
    ) -> Result<()> {
        Err(Error::new(InvalidData, format!("No data.")))
    }

    fn from_stream<U: Read>(&mut self, _stream: &mut U, _threads: u16) -> Result<bool> {
        Err(Error::new(InvalidData, format!("No data.")))
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct InvertedRecord {
    start: DeltaVByte, //Vec<u8>, //IntegerEncode::DeltaVByte // Vec<u64>,
    end: VByte,        //Vec<u8>, //IntegerEncode::VByte,　// Vec<u64>
    name: Deflate,     //Vec<u8>, //StringEncode::Deflate,  // Vec<String>
    aux: Deflate,      //Vec<u8>, //StringEncode::Deflate // We should update better data format.
}

impl ColumnarSet for InvertedRecord {
    fn new() -> InvertedRecord {
        let start = vec![];
        let end = vec![];
        let name = vec![];
        let aux = vec![];
        InvertedRecord {
            start: start,
            end: end,
            name: name,
            aux: aux,
        }
    }

    fn from_stream<U: Read>(&mut self, stream: &mut U, _threads: u16) -> Result<bool> {
        let _n_items = stream.read_u64::<LittleEndian>()?;
        let n_header = stream.read_u64::<LittleEndian>()?;
        for _i in 0..n_header {
            let _a = stream.read_u16::<LittleEndian>()?;
        }
        self.start = vec![];
        self.end = vec![];
        self.name = vec![];

        let mut n_items = stream.read_u64::<LittleEndian>()?;
        for _i in 0..n_items {
            self.start.push(stream.read_u8()?);
        }
        n_items = stream.read_u64::<LittleEndian>()?;
        for _i in 0..n_items {
            self.end.push(stream.read_u8()?);
        }
        n_items = stream.read_u64::<LittleEndian>()?;
        for _i in 0..n_items {
            self.name.push(stream.read_u8()?);
        }
        n_items = stream.read_u64::<LittleEndian>()?;
        for _i in 0..n_items {
            self.aux.push(stream.read_u8()?);
        }

        Ok(true)
    }

    fn to_stream<W: Write, R: Read + Seek>(
        &self,
        stream: &mut W,
        _threads: u16,
        _bam_reader: Option<&mut IndexedReader<R>>,
    ) -> Result<()> {
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
            stream.write_u8(*i)?;
        }
        stream.write_u64::<LittleEndian>(self.aux.len() as u64)?;
        for i in &self.aux {
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
        let aux_raw = builder
            .aux
            .clone()
            .into_inner()
            .into_iter()
            .map(|t| t.join("\t").to_owned())
            .collect();
        let aux = string_encode(&aux_raw);
        InvertedRecord {
            start: start,
            end: end,
            name: name,
            aux: aux,
        }
    }
    /*
        // Output as Record with F(u64) -> &str
        pub fn to_record_f<F>(self, to_str: F) -> Vec<bed::Record> where
        F: Fn(u64) -> &str
        {
            let chromosome =
            to_record(self, chromosome);
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
