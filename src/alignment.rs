use bam::RecordWriter;
use bam::RecordReader;
use bam::record::Record;
use std::{collections::{BTreeMap, HashMap}, io::{Read, Result}};
use byteorder::{LittleEndian, WriteBytesExt, ReadBytesExt};
use crate::{index::region_to_bin_3, ColumnarSet, range::Format, Builder};

/// BAM-Compatible Alignment Inverted-Record
#[derive(Clone, Debug)]
pub struct Alignment {
    data: Vec<Record> //bgzip compressed records
}

impl PartialEq for Alignment {
    fn eq(&self, other: &Self) -> bool {
        if self.data.len() != other.data.len() {
            return false
        }
        for i in 0..self.data.len() {
            if self.data[i].name() != other.data[i].name() {
                return false
            }
        }
        true
    }
}

impl Eq for Alignment {

}

/// BAM-Compatible Alignment Record
pub struct AlignmentBuilder {
    alignments: Vec<Record>,
}

impl AlignmentBuilder {
    pub fn new() -> Self {
        AlignmentBuilder{ 
            alignments: vec![]
        }
    }
    pub fn add(&mut self, alignment: Record) {
        self.alignments.push(alignment);
    }

    pub fn take(self) -> Vec<Record> {
        self.alignments
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
        let mut unmapped = AlignmentBuilder::new();

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
    pub fn to_record(self) -> Alignment {
        //let mut aln = Alignment::new();
        //aln.from_builder(self).unwrap();
        //aln
        Alignment { data: self.take() }
    }
}
impl Builder for AlignmentBuilder {
    fn to_format(self) -> Format {
        Format::Alignment(self.to_record())
    }
}
/*
impl AlignmentOld { 
    pub fn new_from_builder(builder: &AlignmentBuilder) -> Result<AlignmentOld> {
        let mut binary = vec![];
        binary.write_u64::<LittleEndian>(builder.alignments.borrow_mut().len() as u64)?;
        let header = bam::Header::new();
        let output = BufWriter::new(binary);
        let mut writer = bam::bam_writer::BamWriterBuilder::new().write_header(false).from_stream(output, header)?;
        for i in builder.alignments.borrow().iter() {
            // i.write_bam(&mut binary)?;
            writer.write(&i)?;
        }
        writer.flush()?;
        writer.finish()?;
        let inner = writer.take_stream().into_inner()?;

        Ok(AlignmentOld{data: inner})
        
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
*/
impl ColumnarSet for Alignment {
    fn new() -> Alignment {
        Alignment {data: vec![]}
    }
    fn to_stream<W: std::io::Write>(&self, stream: &mut W) -> Result<()> {
        stream.write_u64::<LittleEndian>(self.data.len() as u64)?;
        let header = bam::Header::new();
        let mut writer = bam::bam_writer::BamWriterBuilder::new().write_header(false).from_stream(stream, header)?;
        for i in self.data.iter() {
            // i.write_bam(&mut binary)?;
            writer.write(&i)?;
        }
        writer.flush()?;
        writer.finish()?;
        Ok(())
    }
    fn from_stream<R: Read>(&mut self, stream: &mut R) -> Result<bool> {
        let len = stream.read_u64::<LittleEndian>()?;
        // let records = Vec::with_capacity(len as usize);
        let mut reader = bam::BamReader::from_stream(stream, 4).unwrap();
        for _i in 0..len as usize {
            let mut record = Record::new();
            // record.fill_from_bam(&mut self.data)?;
            reader.read_into(&mut record)?;
            self.data.push(record);
        }
        Ok(true)
    }
    
}
/*
impl ColumnarSet for AlignmentOld {
    fn new() -> AlignmentOld {
        AlignmentOld { data: vec![] }
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
*/
#[cfg(test)]
mod tests {
    use crate::header::{Header};
    use crate::index::Region;
    use crate::binary;
    //use crate::writer::GhiWriter;
    use crate::reader::IndexedReader;
    use crate::range::{InvertedRecordEntire, Format};
    use crate::{builder::{InvertedRecordBuilder}, alignment::Set};
    use bio::io::bed;
    use super::AlignmentBuilder;
    #[test]
    fn bam_works() {
        let mut header = Header::new();

        let bam_path = "./test/index_test.bam";
        let reader = bam::BamReader::from_path(bam_path, 4).unwrap();
        let set = Set::<AlignmentBuilder>::new(reader, 0 as u64);

        let path = "./test/test.bed";
        let reader = bed::Reader::from_file(path).unwrap();
        let set2: Set<InvertedRecordBuilder> = Set::<InvertedRecordBuilder>::new(reader, 1 as u64, &mut header).unwrap();
        println!("{:?}", header.reference_id("1"));
        println!("{:?}", set2);

        //let set_vec = vec![set, set2];
        let set_vec = vec![set2];
        let mut entire = InvertedRecordEntire::new_from_set(set_vec);
        entire.add(set);
        let mut writer = binary::GhbWriter::build()
            .write_header(false)
            .from_path("./test/test_bam.ghb", header).unwrap();
        let index = entire.write_binary(&mut writer).unwrap();
        writer.flush().unwrap();
        println!("{}", index);

        let mut reader2 = IndexedReader::from_path("./test/test.ghb").unwrap();
        println!("{}", reader2.index());

        let chrom = "2";
        let chrom_id = reader2.reference_id(&chrom).unwrap();
        assert_eq!(chrom_id, 1);
        let viewer = reader2.fetch(&Region::new(chrom_id, 17_000, 17_500)).unwrap();
        let example = "2\t16382\t16385\tbin4682\t20\t-\n2\t16388\t31768\tbin4683\t20\t-\n";
        let records = viewer.into_iter().flat_map(|t| t.map(|f| 
            if let Format::Range(rec) = f.data() {
                // println!("debug {:#?}", rec.to_record(chrom));
                return rec.to_record(chrom)
            } else {
                return vec![]
            }
        ).unwrap()).collect::<Vec<bed::Record>>();
        println!("Records: {:?}", records);
        let mut buf = Vec::new();
        {
            let mut writer = bed::Writer::new(&mut buf);
            for i in records {
                writer.write(&i).ok().unwrap();
            } 
        }
        assert_eq!(
            example,
            String::from_utf8(buf).unwrap().as_str()
        );
    }
}