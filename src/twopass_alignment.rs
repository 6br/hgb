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
    io::{Read, Result, Seek, Write},
};

#[derive(Clone, Debug)]
pub enum Alignment {
    Offset(Vec<Chunk>),
    Object(Vec<Record>),
}

impl PartialEq for Alignment {
    fn eq(&self, other: &Self) -> bool {
        use Alignment::*;
        match (self, other) {
            (Offset(a), Offset(b)) => a == b,
            (Object(data), Object(other)) => {
                if data.len() != other.len() {
                    return false;
                }
                for i in 0..data.len() {
                    if data[i].name() != other[i].name() {
                        return false;
                    }
                }
                true
            }
            _ => false,
        }
    }
}

impl Eq for Alignment {}

/// BAM-Compatible Alignment Record
pub struct AlignmentBuilder {
    alignments: Vec<Chunk>,
    //    reader: String,
}

impl AlignmentBuilder {
    pub fn new() -> Self {
        AlignmentBuilder {
            alignments: vec![],
            //            reader: "".to_string(),
        }
    }
    pub fn add(&mut self, alignment: Chunk) {
        self.alignments.push(alignment);
        //        self.reader = reader.to_string();
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

impl<R: Read + Seek> Set<AlignmentBuilder, R> {
    pub fn new(
        //mut reader: bam::BamReader<R>,
        mut reader: bam::IndexedReader<R>,
        sample_id: u64,
        header: &mut Header,
    ) -> Self {
        let mut chrom = BTreeMap::new();
        let mut unmapped = AlignmentBuilder::new();
        // let mut prev = reader.reader.current().and_then(|t| t.offset()).unwrap();
        let mut prev = reader.index().start_offset().unwrap();
        let mut rec = Record::new();
        let mut viewer = reader.full();
        let mut prev_next_offset = 0; //viewer.parent.reader.reader.next_offset().unwrap();
                                      //for record in reader {
        while let Ok(true) = viewer.read_into(&mut rec) {
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
                    let end_offset: u64 = viewer
                        .parent
                        .reader
                        .current()
                        .and_then(|t| t.offset())
                        .unwrap();
                    let contents_offset = viewer.parent.reader.contents_offset();
                    let end = VirtualOffset::from_raw(end_offset << 16 | contents_offset as u64);
                    let next_offset = viewer.parent.reader.reader.next_offset().unwrap();

                    assert!(end > prev);
                    stat.add(
                        Chunk::new(prev, end),
                        //  header.get_name(sample_id as usize).unwrap(),
                    );
                    //TODO() Chunks should be merged if the two chunks are neighbor.
                    if prev.block_offset() != end.block_offset()
                        && end.block_offset() != prev_next_offset
                        && end.contents_offset() == 0
                    {
                        // println!("a: {} {} {} {}", prev, end, contents_offset, next_offset);
                        prev = VirtualOffset::new(next_offset, 0);
                    } else {
                        prev = end;
                    }
                    prev_next_offset = next_offset;
                }
            } else if rec.ref_id() == -1 {
                let end_offset: u64 = viewer
                    .parent
                    .reader
                    .current()
                    .and_then(|t| t.offset())
                    .unwrap();
                let content_offset = viewer.parent.reader.contents_offset();
                let end = VirtualOffset::from_raw(end_offset << 16 | content_offset as u64);
                unmapped.alignments.push(Chunk::new(prev, end));
                prev = end;
            } else {
                // return Err(io::Error::new(InvalidData, "Reference id < -1"));
                panic!("Reference id < -1");
            }
        }
        Set::<AlignmentBuilder, R> {
            sample_id,
            chrom,
            unmapped,
            bam_reader: Some(reader),
        }
    }
}

impl AlignmentBuilder {
    pub fn to_record(self) -> Alignment {
        //        Alignment {
        //            reader: self.reader.clone(),
        //            data: self.take(),
        //        }
        Alignment::Offset(self.take())
    }
}

impl Builder for AlignmentBuilder {
    fn to_format(self) -> Format {
        Format::Alignment(self.to_record())
    }
}
impl ColumnarSet for Alignment {
    fn new() -> Alignment {
        Alignment::Object(vec![])
    }
    fn to_stream<W: Write, R: Read + Seek>(
        &self,
        stream: &mut W,
        threads: u16,
        bam_reader: Option<&mut bam::IndexedReader<R>>,
    ) -> Result<()> {
        let len = match self {
            Alignment::Offset(vec) => vec.len(),
            Alignment::Object(vec) => vec.len(),
        };
        stream.write_u64::<LittleEndian>(len as u64)?;
        let header = bam::Header::new();
        let mut writer = bam::bam_writer::BamWriterBuilder::new()
            .additional_threads(threads)
            .write_header(false)
            .from_stream(stream, header)?;
        // TODO() Inject the number of threads.
        // let mut reader = bam::IndexedReader::from_path(&self.reader)?;
        let _ = match self {
            Alignment::Offset(data) => {
                let viewer = bam_reader.unwrap().chunk(data.clone());
                for i in viewer {
                    writer.write(&i?)?;
                    //TODO() Don't have to parse in this implementation.
                }
            }
            Alignment::Object(vec) => {
                //let viewer = bam_reader.unwrap().chunk(self.data.clone());
                for i in vec {
                    writer.write(&i)?;
                    //TODO() Don't have to parse in this implementation.
                }
            }
        };

        writer.finish()?;
        Ok(())
    }
    fn from_stream<R: Read>(&mut self, stream: &mut R, threads: u16) -> Result<bool> {
        let len = stream.read_u64::<LittleEndian>()?;

        let mut reader =
            bam::BamReader::from_stream_no_header(stream, bam::Header::new(), threads).unwrap();
        let mut data = Vec::with_capacity(len as usize);
        for _i in 0..len as usize {
            let mut record = Record::new();
            reader.read_into(&mut record)?;
            data.push(record);
        }
        let alignment = Alignment::Object(data);
        let _old = std::mem::replace(self, alignment);
        let mut _record = Record::new();
        // For consuming all bgzip blocks.
        let _ = reader.next();
        let _ = reader.next();

        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::{Alignment, AlignmentBuilder};
    use crate::bam::RecordWriter;
    use crate::binary;
    use crate::header::Header;
    use crate::index::Region;
    use crate::range::{Format, InvertedRecordEntire};
    use crate::reader::IndexedReader;
    use crate::writer::GhiWriter;
    use crate::IndexWriter;
    use crate::{builder::InvertedRecordBuilder, twopass_alignment::Set};
    use bam::record::Record;
    use bio::io::bed;
    use std::{
        fs::File,
        io::{BufRead, BufReader, Write},
        path::Path,
        process::Command,
        time::Instant,
    };

    #[test]
    fn two_bam_works() {
        let bam_path = "./test/index_test.bam";
        // let reader = bam::BamReader::from_path(bam_path, 4).unwrap();
        let reader2 = bam::IndexedReader::from_path(bam_path).unwrap();
        // println!("{}", reader2.index());

        let bam_header = reader2.header();
        let mut header = Header::new();
        header.transfer(bam_header);
        header.set_local_header(bam_header, bam_path, 0);
        {
            let set = Set::<AlignmentBuilder, File>::new(reader2, 0 as u64, &mut header);
            let bam_path2 = "./test/test-in.bam";
            let reader = bam::IndexedReader::from_path(bam_path2).unwrap();
            let set2: Set<AlignmentBuilder, File> =
                Set::<AlignmentBuilder, File>::new(reader, 1 as u64, &mut header);

            assert_eq!(None, header.reference_id("1"));
            assert_eq!(Some(1), header.reference_id("chr1"));
            assert_eq!(Some(2), header.reference_id("chr2"));

            let dummy_header = Header::new();
            let set_vec = vec![set];
            let mut entire: InvertedRecordEntire<File> =
                InvertedRecordEntire::new_from_set(set_vec);
            // println!("{:?}", entire);
            entire.add(set2);
            // entire.add_reader(0, reader2);
            let mut writer = binary::GhbWriter::build()
                .write_header(false)
                .from_path("./test/test_bam.ghb", dummy_header)
                .unwrap();
            let index = entire.write_binary(&mut writer).unwrap();
            writer.flush().unwrap();

            entire.write_header(&mut header);
            let mut index_writer = GhiWriter::build()
                .write_header(true)
                .from_path("./test/test_bam.ghb.ghi", header)
                .unwrap();
            let _result = index_writer.write(&index);
            assert_eq!(_result.ok(), Some(()));
            let _result = index_writer.flush();
        }
    }
    #[test]
    fn bam_works() {
        let bam_path = "./test/index_test.bam";
        // let reader = bam::BamReader::from_path(bam_path, 4).unwrap();
        let reader2 = bam::IndexedReader::from_path(bam_path).unwrap();
        // println!("{}", reader2.index());

        let bam_header = reader2.header();
        let mut header = Header::new();
        header.transfer(bam_header);
        header.set_local_header(bam_header, bam_path, 0);
        {
            let set = Set::<AlignmentBuilder, File>::new(reader2, 0 as u64, &mut header);

            let example =
                b"chr2\t16382\t16385\tbin4682\t20\t-\nchr2\t16388\t31768\tbin4683\t20\t-\n";
            let reader = bed::Reader::new(&example[..]);
            let set2: Set<InvertedRecordBuilder, File> =
                Set::<InvertedRecordBuilder, File>::new(reader, 1 as u64, &mut header).unwrap();

            assert_eq!(None, header.reference_id("1"));
            assert_eq!(Some(1), header.reference_id("chr1"));
            assert_eq!(Some(2), header.reference_id("chr2"));

            let dummy_header = Header::new();
            let set_vec = vec![set2];
            let mut entire: InvertedRecordEntire<File> =
                InvertedRecordEntire::new_from_set(set_vec);
            // println!("{:?}", entire);
            entire.add(set);
            // entire.add_reader(0, reader2);
            let mut writer = binary::GhbWriter::build()
                .write_header(false)
                .from_path("./test/test_bam.ghb", dummy_header)
                .unwrap();
            let index = entire.write_binary(&mut writer).unwrap();
            writer.flush().unwrap();

            entire.write_header(&mut header);
            let mut index_writer = GhiWriter::build()
                .write_header(true)
                .from_path("./test/test_bam.ghb.ghi", header)
                .unwrap();
            let _result = index_writer.write(&index);
            assert_eq!(_result.ok(), Some(()));
            let _result = index_writer.flush();
        }

        let mut reader2 = IndexedReader::from_path("./test/test_bam.ghb").unwrap();
        // println!("{}", reader2.index().references()[2]);

        let chrom = "chr2";
        let chrom_id = reader2.reference_id(&chrom).unwrap();

        assert_eq!(chrom_id, 2);

        {
            let chrom_name = reader2.reference_name(0).unwrap();
            assert_eq!(chrom_name, "chrM");
        }
        let viewer = reader2
            .fetch(&Region::new(chrom_id, 17_000, 17_500))
            .unwrap();
        let example = "chr2\t16382\t16385\tbin4682\t20\t-\nchr2\t16388\t31768\tbin4683\t20\t-\n";
        let records = viewer
            .into_iter()
            .flat_map(|t| {
                t.map(|f| {
                    if let Format::Range(rec) = f.data() {
                        return rec.to_record(chrom);
                    } else {
                        return vec![];
                    }
                })
                .unwrap()
            })
            .collect::<Vec<bed::Record>>();
        let mut buf = Vec::new();
        {
            let mut writer = bed::Writer::new(&mut buf);
            for i in records {
                writer.write(&i).ok().unwrap();
            }
        }
        assert_eq!(example, String::from_utf8(buf).unwrap().as_str());

        /* check if chr1 paired end read is rescued */
        let chrom_1 = reader2.reference_id("chr1").unwrap();
        let viewer = reader2
            .fetch(&Region::new(chrom_1, 470_000, 471_500))
            .unwrap();
        let records = viewer
            .into_iter()
            .flat_map(|t| {
                t.map(|f|
            // println!("debug {:#?}", t.to_record(chrom));
            if let Format::Alignment(Alignment::Object(rec)) = f.data() {
                return rec
            } else {
                return vec![]
            }
        ).unwrap()
            })
            .collect::<Vec<Record>>();
        assert_eq!(records.len(), 2);

        /* Check if bam can reconstruct, except for unmapped reads */
        let sam_output = format!("./test/index_output.sam");
        let mut count = 0;
        let output1 = format!("./test/index_test.sam");
        let header = reader2.header().get_local_bam_header(0).unwrap().clone();
        // let mut writer = bam::BamWriter::from_path(bam_output, header).unwrap();
        let mut writer = bam::SamWriter::from_path(&sam_output, header).unwrap();
        let viewer2 = reader2.full();
        viewer2.into_iter().for_each(|t| {
            t.map(|f| match f.data() {
                Format::Alignment(Alignment::Object(record)) => {
                    for i in record {
                        writer.write(&i).unwrap();
                        count += 1;
                    }
                    // return rec.to_record(c)
                    ()
                }
                _ => (),
            })
            .unwrap()
        });
        writer.finish().unwrap();

        let mut log = File::create("./test/compare_alignment2_pass.log").unwrap();

        let timer = Instant::now();
        let mut child = Command::new("samtools")
            .args(&["view", "-h", "--no-PG", "-F", "4"])
            .arg(&bam_path)
            .args(&["-o", &output1])
            .spawn()
            .expect("Failed to run samtools view");
        let ecode = child.wait().expect("Failed to wait on samtools view");
        assert!(ecode.success());
        writeln!(log, "        samtools view:  {:?}", timer.elapsed()).unwrap();
        writeln!(log, "        total {} records", count).unwrap();

        compare_sam_files(&output1, &sam_output, &mut log);
    }

    fn compare_sam_files<P: AsRef<Path>, T: AsRef<Path>, W: Write>(
        filename1: P,
        filename2: T,
        log: &mut W,
    ) {
        let filename1 = filename1.as_ref();
        let filename2 = filename2.as_ref();
        let mut file1 = BufReader::new(File::open(&filename1).unwrap());
        let mut file2 = BufReader::new(File::open(&filename2).unwrap());
        let mut line1 = String::new();
        let mut line2 = String::new();

        for i in 1_usize.. {
            line1.clear();
            line2.clear();
            match (file1.read_line(&mut line1), file2.read_line(&mut line2)) {
                (Ok(x), Ok(y)) => {
                    if x == 0 && y != 0 {
                        writeln!(
                            log,
                            "Comparing files {} and {}",
                            filename1.display(),
                            filename2.display()
                        )
                        .unwrap();
                        writeln!(log, "Samtools output: {}", line2.trim()).unwrap();
                        writeln!(log, "Samtools output is longer").unwrap();
                        panic!("Samtools output is longer");
                    } else if x != 0 && y == 0 {
                        writeln!(
                            log,
                            "Comparing files {} and {}",
                            filename1.display(),
                            filename2.display()
                        )
                        .unwrap();
                        writeln!(log, "Crate output:    {}", line1.trim()).unwrap();
                        writeln!(log, "Crate output is longer").unwrap();
                        panic!("Crate output is longer");
                    } else if x == 0 && y == 0 {
                        break;
                    }
                }
                (Ok(_), Err(e)) => panic!("Could not read samtools output: {:?}", e),
                (Err(e), Ok(_)) => panic!("Could not read crate output: {:?}", e),
                (Err(e1), Err(e2)) => panic!("Could not read both outputs: {:?}, {:?}", e1, e2),
            }
            if line1 != line2 {
                writeln!(
                    log,
                    "Comparing files {} and {}",
                    filename1.display(),
                    filename2.display()
                )
                .unwrap();
                writeln!(log, "Crate output:    {}", line1.trim()).unwrap();
                writeln!(log, "Samtools output: {}", line2.trim()).unwrap();
                writeln!(log, "Outputs do not match on line {}", i).unwrap();
                panic!("Outputs do not match on line {}", i);
            }
        }
    }
}
