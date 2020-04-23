use bio::io::bed;
use crate::index::region_to_bin_3;
use std::cell::RefCell;
//use crate::range::InvertedRecord;
use std::collections::{HashMap, BTreeMap};
use std::io::{Read};

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordReference {
    pub bins: HashMap<u32, InvertedRecordBuilder> // Mutex?
}

impl InvertedRecordReference {
    pub fn new() -> Self {
        InvertedRecordReference{
            bins:  HashMap::new()
        }
    }
}


#[derive(Clone, PartialEq, Eq, PartialOrd, Debug)]
pub struct InvertedRecordBuilder {
    pub start: RefCell<Vec<u64>>,
    pub end: RefCell<Vec<u64>>,
    pub name: RefCell<Vec<String>>, //Just a tab-separated string here.
    pub aux: RefCell<Vec<Vec<String>>>,
}

impl InvertedRecordBuilder {
    pub fn new() -> Self {
        InvertedRecordBuilder{
            start: RefCell::new(Vec::new()),
            end: RefCell::new(Vec::new()),
            name: RefCell::new(Vec::new()),
            aux: RefCell::new(Vec::new()),
        }
    }
    pub fn add(&self, start: u64, end: u64, name: String, aux: Vec<String>) {
        self.start.borrow_mut().push(start);
        self.end.borrow_mut().push(end);
        self.name.borrow_mut().push(name);
        self.aux.borrow_mut().push(aux);
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct InvertedRecordSet {
    pub sample_id: u64,
    pub chrom: BTreeMap<String, InvertedRecordReference> // Mutex?
}

impl InvertedRecordSet {
    pub fn new<R:Read>(mut reader: bed::Reader<R>, sample_id: u64) -> Self {
        let mut inverted_record_set = BTreeMap::new();
        /*let mut start: Vec<u64> = vec![];
        let mut end: Vec<u64> = vec![];
        let mut aux: Vec<String> = vec![];*/
        for record in reader.records() {
            println!("{:?}", record);
            let rec = record.ok().expect("Error reading record.");
            let chrom = inverted_record_set.entry(rec.chrom().to_string()).or_insert(InvertedRecordReference::new());
            let stat = chrom.bins.entry(region_to_bin_3(rec.start(), rec.end())).or_insert(InvertedRecordBuilder::new());
            // rec.chrom;
            let mut aux = vec![];
            let mut n = 4; // Ignore name field
            while let Some(item) = rec.aux(n){
                aux.push(item.to_string());
                n += 1;
            }
            stat.add(rec.start(), rec.end(), rec.name().unwrap_or("").to_string(), aux)
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
    use crate::header::{Header};
    use bio::io::bed;
    use std::io;
    use crate::binary;
    use crate::writer;
    use crate::reader::IndexedReader;
    use crate::range::{InvertedRecordEntire, InvertedRecord, Format, Record};
    use crate::IndexWriter;
    use crate::reader;
    use crate::index::Region;
    use super::InvertedRecordSet;

    #[test]
    fn full_works() {
        let path = "./test/test.bed";
        let mut reader = bed::Reader::from_file(path).unwrap();
        let set = InvertedRecordSet::new(reader, 0 as u64);

        let set_vec = vec![set];
        let entire = InvertedRecordEntire::new(set_vec);
        let header = Header::new();
        let mut writer = binary::GhbWriter::build()
            .write_header(false)
            .from_path("./test/test2.ghb", header).unwrap();
        let index = entire.write_binary(&mut writer).unwrap();
        writer.flush().unwrap();
        println!("{}", index);


        let mut header2 = Header::new();
        entire.write_header(&mut header2);
        let mut index_writer = writer::GhbWriter::build().write_header(true).from_path("./test/test2.ghb.ghi", header2).unwrap();
        index_writer.write(&index);
        index_writer.flush();

        let mut reader2 = IndexedReader::from_path("./test/test2.ghb").unwrap();
        println!("{}", reader2.index());

        let chrom = "1";
        let chrom_id = reader2.reference_id(&chrom).unwrap();
        let viewer = reader2.full();
        let records = viewer.into_iter().scan((), |_,x| x.ok()).collect::<Vec<Record>>();
        println!("Records: {:?}", records);

        assert_eq!(records.len(), 10);
    }

    #[test]
    fn it_works() {
        // let _example = b"1\t5\t5000\tname1\t0.5\n1\t5\t5000\tname1\t0.5";
        let path = "./test/test.bed";
        let mut reader = bed::Reader::from_file(path).unwrap();
        let set = InvertedRecordSet::new(reader, 0 as u64);
        // println!("{:?}", set);
        let set_vec = vec![set];
        let entire = InvertedRecordEntire::new(set_vec);
        // let output = io::BufWriter::new(io::stdout());
        let header = Header::new();
        let mut writer = binary::GhbWriter::build()
            .write_header(false)
            .from_path("./test/test.ghb", header).unwrap();
        let index = entire.write_binary(&mut writer).unwrap();
        writer.flush().unwrap();
        println!("{}", index);

        let mut header2 = Header::new();
        entire.write_header(&mut header2);
        let mut index_writer = writer::GhbWriter::build().write_header(true).from_path("./test/test.ghb.ghi", header2).unwrap();
        index_writer.write(&index);
        index_writer.flush();

        let mut reader2 = IndexedReader::from_path("./test/test.ghb").unwrap();
        println!("{}", reader2.index());
        

        //assert_eq!(format!("{}", index), format!("{}", reader2.index()));
        // assert_eq!(&index, reader2.index());
        let chrom = "2";
        let chrom_id = reader2.reference_id(&chrom).unwrap();
        println!("{}", chrom_id);
        let viewer = reader2.fetch(&Region::new(chrom_id, 17_000, 17_500)).unwrap();
        // println!("{}", viewer.index());
        /*
        for record in viewer {
            // println!("A");
            let record = record.unwrap();
            println!("Record: {:?}", record);
        }*/
        let records = viewer.into_iter().flat_map(|t| t.map(|f| 
            if let Format::Range(rec) = f.data() {
                // println!("debug {:#?}", rec.to_record(chrom));
                return rec.to_record(chrom)
            } else {
                return vec![]
            }
        ).unwrap()).collect::<Vec<bed::Record>>();
        println!("Records: {:?}", records);

        let example = "2\t16382\t16385\tbin4682\t20\t-\n2\t16388\t31768\tbin4683\t20\t-\n";
        // let mut test_reader = bed::Reader::new(&example[..]);
        let mut buf = Vec::new();
        {
            let mut writer = bed::Writer::new(&mut buf);
            for  i in records {
                writer.write(&i).ok().unwrap();
            } 
        }
        assert_eq!(
            example,
            String::from_utf8(buf).unwrap().as_str()
        );
        // let records_from = reader.records().into_iter().flat_map(|t| t).collect::<Vec<bed::Record>>();
        //assert_eq!(records[0], records_from[0]);
        //assert_eq!(records[1], records_from[1]);
        //assert_eq!(records.len(), records_from.len());

        let viewer = reader2.fetch(&Region::new(0, 17_000, 17_500)).unwrap();
        // println!("{}", viewer.index());
        for record in viewer {
            // println!("A");
            let record = record.unwrap();
            // println!("Record: {:?}", record);
        }

        let viewer = reader2.fetch(&Region::new(1, 1, 3)).unwrap();
        // println!("{}", viewer.index());
        for record in viewer {
            // println!("A");
            let record = record.unwrap();
            // println!("Record: {:?}", record);
        }
    }
}