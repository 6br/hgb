use bam;
use bam::{Record, RecordWriter};
use clap::ArgMatches;
use genomic_range::StringRegion;

use log::{debug, info};
// use rayon::prelude::*;
use crate::server::server;
use ghi::bed;
use ghi::binary::GhbWriter;
use ghi::builder::InvertedRecordBuilder;
use ghi::header::Header;
use ghi::index::{Chunk, Region, VirtualOffset};
use ghi::range::Default;
use ghi::range::{Format, InvertedRecordEntire, Set};
use ghi::twopass_alignment::{Alignment, AlignmentBuilder};
use ghi::vis::{bam_record_vis, RecordIter};
use ghi::writer::GhiWriter;
use ghi::{gff, reader::IndexedReader, IndexWriter};
use io::{BufReader, Error, ErrorKind, Write};
use itertools::EitherOrBoth::{Both, Left};
use itertools::Itertools;
use std::{
    collections::{BTreeMap, HashMap},
    env::Args,
    fs::File,
    io,
    path::Path,
};

pub fn bam_vis(
    matches: &ArgMatches,
    args: Vec<String>,
    threads: u16,
) -> Result<(), Box<dyn std::error::Error>> {
    // let output_path = matches.value_of("OUTPUT").unwrap();
    let min_read_len = matches
        .value_of("min-read-length")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(0u32);

    if let Some(ranges) = matches.values_of("range") {
        let ranges: Vec<&str> = ranges.collect();
        let prefetch_ranges: Vec<&str> = matches
            .values_of("prefetch-range")
            .and_then(|t| Some(t.collect()))
            .unwrap_or(vec![]);
        /*let ranged_zip = if let Some(prefetch_ranges) = prefetch_ranges {
            ranges.into_iter().zip(prefetch_ranges)
        } else {
            ranges.into_iter().zip(ranges.clone())
        };*/
        for eob in ranges.into_iter().zip_longest(prefetch_ranges) {
            let (range, prefetch_str) = match eob {
                Both(a, b) => (a, b),
                Left(a) => (a, a),
                _ => panic!("Range is not specified."),
            };
            let string_range = StringRegion::new(range).unwrap();
            let prefetch_range = StringRegion::new(prefetch_str).unwrap();
            /*let prefetch_range = if let Some(prefetch_ranges) = prefetch_ranges {
                StringRegion::new(prefetch_ranges[index])?
            } else {
                string_range
            };*/
            if let Some(bam_files) = matches.values_of("bam") {
                let bam_files: Vec<&str> = bam_files.collect();
                let mut list: Vec<(u64, Record)> = vec![];
                println!("Input file: {:?}", bam_files);
                for (index, bam_path) in bam_files.iter().enumerate() {
                    println!("Loading {}", bam_path);
                    // let reader = bam::BamReader::from_path(bam_path, threads).unwrap();
                    let mut reader2 = bam::IndexedReader::build()
                        .additional_threads(threads - 1)
                        .from_path(bam_path)?;

                    // Here all threads can be used, but I suspect that runs double
                    //reader2.fetch()
                    let ref_id = reader2
                        .header()
                        .reference_id(prefetch_range.path.as_ref())
                        .ok_or(Error::new(ErrorKind::Other, "Invalid reference id."))?;
                    let viewer = reader2.fetch(&bam::bam_reader::Region::new(
                        ref_id,
                        prefetch_range.start as u32,
                        prefetch_range.end as u32,
                    ))?;
                    for record in viewer {
                        let record = record?;
                        if !record.flag().is_secondary() && record.query_len() > min_read_len {
                            list.push((index as u64, record));
                        }
                    }
                }
                let mut ann = vec![];
                let mut idx = bam_files.len();
                let mut freq = BTreeMap::new();
                if let Some(freq_files) = matches.values_of("frequency") {
                    // let bed_files: Vec<_> = matches.values_of("bed").unwrap().collect();
                    // frequency bed file needs to be (start, score).
                    let freq_files: Vec<&str> = freq_files.collect();
                    for (_idx, bed_path) in freq_files.iter().enumerate() {
                        info!("Loading {}", bed_path);
                        let mut reader = bed::Reader::from_file(bed_path)?;
                        let mut values = vec![];
                        for record in reader.records() {
                            let record = record?;
                            if record.end() > prefetch_range.start()
                                && record.start() < prefetch_range.end()
                                && record.chrom() == prefetch_range.path
                            {
                                values.push((
                                    record.start(),
                                    record
                                        .score()
                                        .and_then(|t| t.parse::<u32>().ok())
                                        .unwrap_or(0),
                                ));
                            }
                        }
                        freq.insert(idx as u64, values);
                        idx += 1;
                    }
                }

                if let Some(bed_files) = matches.values_of("bed") {
                    let bed_files: Vec<&str> = bed_files.collect();
                    for (_idx, bed_path) in bed_files.iter().enumerate() {
                        info!("Loading {}", bed_path);
                        let mut reader = bed::Reader::from_file(bed_path)?;
                        for record in reader.records() {
                            let record = record?;
                            if record.end() > prefetch_range.start()
                                && record.start() < prefetch_range.end()
                                && record.chrom() == prefetch_range.path
                            {
                                ann.push((idx as u64, record));
                            }
                        }
                        idx += 1;
                    }
                }
                if let Some(gff_files) = matches.values_of("gff3") {
                    let gff_files: Vec<&str> = gff_files.collect();
                    for (_idx, gff_path) in gff_files.iter().enumerate() {
                        info!("Loading {}", gff_path);
                        let mut reader =
                            gff::Reader::from_file(gff_path, gff::GffType::GFF3).unwrap();
                        for gff_record in reader.records() {
                            let gff = gff_record?;
                            if *gff.end() > prefetch_range.start()
                                && *gff.start() < prefetch_range.end()
                                && gff.seqname() == prefetch_range.path
                            {
                                let mut record = bed::Record::new();
                                record.set_chrom(gff.seqname());
                                record.set_start(*gff.start());
                                record.set_end(*gff.end());
                                record.set_name(&gff.attributes()["gene_id"]);
                                record.set_score(&gff.score().unwrap_or(0).to_string());
                                if let Some(strand) = gff.strand() {
                                    record.push_aux(strand.strand_symbol()); // Strand
                                }
                                ann.push((idx as u64, record));
                            }
                        }
                        idx += 1;
                    }
                }
                bam_record_vis_pre_calculate(
                    matches,
                    string_range,
                    prefetch_range,
                    &args,
                    list,
                    ann,
                    freq,
                    threads,
                    |idx| Some(bam_files[idx]),
                )?;
                /*
                if matches.is_present("web") {
                    server(
                        matches.clone(),
                        string_range,
                        prefetch_range,
                        args.clone(),
                        list,
                        ann,
                        freq,
                        threads,
                    )?;
                } else {
                    bam_record_vis(matches, string_range, list, ann, freq, |idx| {
                        Some(bam_files[idx])
                    })?;
                }*/
            }
        }
    }
    Ok(())
}

pub fn build(matches: &ArgMatches, threads: u16) -> () {
    let mut header = Header::new();
    let mut alignment_transfer = false;
    let output_path = matches.value_of("OUTPUT").unwrap();

    if let Some(o) = matches.value_of("chrom") {
        info!("Loading chromosome sizes.");
        let file = File::open(o).unwrap();
        let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);
        let _result = header.set_header_from_sizes(&mut rdr).unwrap();
    } else {
        alignment_transfer = true;
    }

    // let mut set_vec = vec![];
    let mut i = 0;

    let mut records: InvertedRecordEntire<BufReader<File>> =
        InvertedRecordEntire::<BufReader<File>>::new();

    if let Some(bam_files) = matches.values_of("bam") {
        let bam_files: Vec<&str> = bam_files.collect();
        println!("Input file: {:?}", bam_files);
        for bam_path in bam_files.iter() {
            println!("Loading {}", bam_path);
            // let reader = bam::BamReader::from_path(bam_path, threads).unwrap();
            let reader2 = bam::IndexedReader::build()
                .additional_threads(threads - 1)
                .from_path(bam_path)
                .unwrap();
            // Here all threads can be used, but I suspect that runs double
            let bam_header = reader2.header();
            if alignment_transfer {
                header.transfer(bam_header);
            }
            header.set_local_header(bam_header, bam_path, i);
            let set = Set::<AlignmentBuilder, BufReader<File>>::new(reader2, i as u64, &mut header);
            i += 1;
            records.add(set);
        }
    }

    if let Some(bed_files) = matches.values_of("bed") {
        // let bed_files: Vec<_> = matches.values_of("bed").unwrap().collect();
        let bed_files: Vec<&str> = bed_files.collect();
        for bed_path in bed_files {
            info!("Loading {}", bed_path);
            let reader = bed::Reader::from_file(bed_path).unwrap();
            let set: Set<InvertedRecordBuilder, BufReader<File>> =
                Set::<InvertedRecordBuilder, BufReader<File>>::new(reader, i as u64, &mut header)
                    .unwrap();
            header.set_local_header(&bam::Header::new(), bed_path, i);
            i += 1;
            records.add(set);
        }
    }

    let dummy_header = Header::new();
    let mut writer = GhbWriter::build()
        .write_header(false)
        .additional_threads(threads - 1)
        .from_path(output_path, dummy_header)
        .unwrap();
    let index = records.write_binary(&mut writer).unwrap();
    writer.flush().unwrap();
    let output_index_path = format!("{}.ghi", output_path);
    let mut index_writer = GhiWriter::build()
        .write_header(true)
        .additional_threads(threads - 1)
        .from_path(output_index_path, header)
        .unwrap();
    let _result = index_writer.write(&index);
    assert_eq!(_result.ok(), Some(()));
    let _result = index_writer.flush();
}

pub fn query(matches: &ArgMatches, threads: u16) -> () {
    if let Some(o) = matches.value_of("INPUT") {
        let mut reader: IndexedReader<BufReader<File>> =
            IndexedReader::from_path_with_additional_threads(o, threads - 1).unwrap();
        if let Some(ranges) = matches.values_of("range") {
            let ranges: Vec<&str> = ranges.collect();
            for range in ranges {
                eprintln!("{}", range);
                let closure = |x: &str| reader.reference_id(x);
                let string_range = StringRegion::new(range).unwrap();
                let reference_name = &string_range.path;

                let range = Region::convert(&string_range, closure).unwrap();
                let viewer = reader.fetch(&range).unwrap();

                let sample_ids_opt: Option<Vec<u64>> = matches
                    .values_of("id")
                    //.unwrap()
                    .and_then(|a| Some(a.map(|t| t.parse::<u64>().unwrap()).collect()));
                let sample_id_cond = sample_ids_opt.is_some();
                let sample_ids = sample_ids_opt.unwrap_or(vec![]);
                let filter = matches.is_present("filter");

                let format_type_opt = matches.value_of_t::<Format>("type");
                let format_type_cond = format_type_opt.is_ok();
                let format_type = format_type_opt.unwrap_or(Format::Default(Default {}));
                let out = std::io::stdout();
                let out_writer = match matches.value_of("output") {
                    Some(x) => {
                        let path = Path::new(x);
                        Box::new(File::create(&path).unwrap()) as Box<dyn Write>
                    }
                    None => Box::new(out.lock()) as Box<dyn Write>,
                };
                let mut output = io::BufWriter::with_capacity(1048576, out_writer);

                let header = viewer.header().clone();

                let _ = viewer.into_iter().for_each(|t| {
                    //eprintln!("{:?}", t.clone().unwrap());
                    let f = t.unwrap();
                    //if let Ok(f) = t {
                    if !sample_id_cond || sample_ids.iter().any(|&i| i == f.sample_id()) {
                        let sample_id = f.sample_id();
                        let data = f.data();
                        if !format_type_cond
                            || std::mem::discriminant(&format_type) == std::mem::discriminant(&data)
                        {
                            match data {
                                Format::Range(rec) => {
                                    let mut writer = bed::Writer::new(&mut output);
                                    for i in rec.to_record(&reference_name) {
                                        writer.write(&i).unwrap();
                                    }
                                }
                                Format::Alignment(Alignment::Object(rec)) => {
                                    for i in rec {
                                        if !filter
                                            || (i.calculate_end() as u64 > range.start()
                                                && range.end() > i.start() as u64)
                                        {
                                            let _result = i
                                                .write_sam(
                                                    &mut output,
                                                    header
                                                        .get_local_header(sample_id as usize)
                                                        .unwrap()
                                                        .bam_header(),
                                                )
                                                .unwrap();
                                        }
                                    }
                                }
                                _ => {}
                            }
                        }
                    }
                    //}
                });
            }
        }
    }
}

pub fn decompose(matches: &ArgMatches, _threads: u16) -> () {
    if let Some(i) = matches.value_of("INPUT") {
        //if let Some(o) = matches.value_of("OUTPUT") {
        let id = matches
            .value_of("id")
            .and_then(|t| t.parse::<u64>().ok())
            .unwrap();
        let header = matches.is_present("header");
        let mut reader = IndexedReader::from_path_with_additional_threads(i, 0).unwrap();
        // Decomposer doesn't know the record boundary, so it can't parallelize.
        let out = std::io::stdout();
        let out_writer = match matches.value_of("output") {
            Some(x) => {
                let path = Path::new(x);
                Box::new(File::create(&path).unwrap()) as Box<dyn Write>
            }
            None => Box::new(out.lock()) as Box<dyn Write>,
        };
        let mut output = io::BufWriter::with_capacity(1048576, out_writer);

        if let Some(header_type) = reader.header().get_local_header(id as usize) {
            if header {
                header_type.to_text(&mut output).unwrap();
            }
        // header.write_text(&mut writer);
        } else {
            println!("There is no header of id {}", id);
        }

        if header {
            return;
        }

        // todo!("Implement later; Now just returns only header.");

        let viewer = reader.full();
        let header_data = viewer.header().clone();

        let _ = viewer.into_iter().for_each(|t| {
            //eprintln!("{:?}", t);
            t.map(|f| {
                if f.sample_id() == id {
                    match f.data() {
                        Format::Range(rec) => {
                            let mut writer = bed::Writer::new(&mut output);
                            for i in rec.to_record("null") {
                                writer.write(&i).unwrap();
                            }
                        }
                        Format::Alignment(Alignment::Object(rec)) => {
                            for i in rec {
                                //let _result = i.write_bam(&mut output).unwrap();
                                let _result = i
                                    .write_sam(
                                        &mut output,
                                        header_data
                                            .get_local_header(id as usize)
                                            .unwrap()
                                            .bam_header(),
                                    )
                                    .unwrap();
                            }
                        }
                        _ => {}
                    }
                }
            })
            .unwrap()
        });
    }
    //}
}

pub fn split(matches: &ArgMatches, threads: u16) -> () {
    if let Some(bam_path) = matches.value_of("INPUT") {
        //todo!("Implement a web server using actix-web.");
        let mut reader2 = bam::IndexedReader::build()
            .additional_threads(threads - 1)
            .from_path(bam_path)
            .unwrap();
        // Here all threads can be used, but I suspect that runs double
        let bam_header = reader2.header();

        let out = std::io::stdout();
        let out_writer = match matches.value_of("output") {
            Some(x) => {
                let path = Path::new(x);
                Box::new(File::create(&path).unwrap()) as Box<dyn Write>
            }
            None => Box::new(out.lock()) as Box<dyn Write>,
        };
        let output = io::BufWriter::with_capacity(1048576, out_writer);
        let clevel = matches
            .value_of("compression")
            .and_then(|a| a.parse::<u8>().ok())
            .unwrap_or(6u8);
        let header = bam_header;
        let mut _writer = bam::bam_writer::BamWriterBuilder::new()
            .additional_threads(threads - 1)
            .compression_level(clevel)
            .write_header(true)
            .from_stream(output, header.clone())
            .unwrap();

        let viewer = reader2.full();
        for _record in viewer {}
    }
}

pub fn bam_query(matches: &ArgMatches, threads: u16) -> () {
    if let Some(o) = matches.value_of("INPUT") {
        let mut reader: IndexedReader<BufReader<File>> =
            IndexedReader::from_path_with_additional_threads(o, threads - 1).unwrap();
        if let Some(ranges) = matches.values_of("range") {
            let ranges: Vec<&str> = ranges.collect();
            for range in ranges {
                eprintln!("{}", range);
                let closure = |x: &str| reader.reference_id(x);
                let string_range = StringRegion::new(range).unwrap();
                let _reference_name = &string_range.path;

                let range = Region::convert(&string_range, closure).unwrap();
                let viewer = reader.fetch(&range).unwrap();
                let clevel = matches
                    .value_of("compression")
                    .and_then(|a| a.parse::<u8>().ok())
                    .unwrap_or(6u8);
                let sample_ids_opt: Option<Vec<u64>> = matches
                    .values_of("id")
                    //.unwrap()
                    .and_then(|a| Some(a.map(|t| t.parse::<u64>().unwrap()).collect()));
                let sample_id_cond = sample_ids_opt.is_some();
                let sample_ids = sample_ids_opt.unwrap_or(vec![]);
                let filter = matches.is_present("filter");

                let format_type_opt = matches.value_of_t::<Format>("type");
                let format_type_cond = format_type_opt.is_ok();
                let format_type = format_type_opt.unwrap_or(Format::Default(Default {}));
                let out = std::io::stdout();
                let out_writer = match matches.value_of("output") {
                    Some(x) => {
                        let path = Path::new(x);
                        Box::new(File::create(&path).unwrap()) as Box<dyn Write>
                    }
                    None => Box::new(out.lock()) as Box<dyn Write>,
                };
                let output = io::BufWriter::with_capacity(1048576, out_writer);

                let header = viewer
                    .header()
                    .get_local_header(*sample_ids.get(0).unwrap_or(&0) as usize)
                    .unwrap()
                    .bam_header();
                let mut writer = bam::bam_writer::BamWriterBuilder::new()
                    .additional_threads(threads - 1)
                    .compression_level(clevel)
                    .write_header(true)
                    .from_stream(output, header.clone())
                    .unwrap();

                let _ = viewer.into_iter().for_each(|t| {
                    //eprintln!("{:?}", t.clone().unwrap());
                    let f = t.unwrap();
                    if !sample_id_cond || sample_ids.iter().any(|&i| i == f.sample_id()) {
                        let _sample_id = f.sample_id();
                        let data = f.data();
                        if !format_type_cond
                            || std::mem::discriminant(&format_type) == std::mem::discriminant(&data)
                        {
                            match data {
                                Format::Range(_rec) => {
                                    /*let mut writer = bed::Writer::new(&mut output);
                                    for i in rec.to_record(&reference_name) {
                                        writer.write(&i).unwrap();
                                    }*/
                                }
                                Format::Alignment(Alignment::Object(rec)) => {
                                    for i in rec {
                                        if !filter
                                            || (i.calculate_end() as u64 > range.start()
                                                && range.end() > i.start() as u64)
                                        {
                                            writer.write(&i).unwrap();
                                        }
                                    }
                                }
                                _ => {}
                            }
                        }
                    }
                    //}
                });
            }
        }
    }
}

pub fn vis_query(
    matches: &ArgMatches,
    args: Vec<String>,
    threads: u16,
) -> Result<(), Box<dyn std::error::Error>> {
    let min_read_len = matches
        .value_of("min-read-length")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(0u32);
    if let Some(o) = matches.value_of("INPUT") {
        let mut reader: IndexedReader<BufReader<File>> =
            IndexedReader::from_path_with_additional_threads(o, threads - 1).unwrap();

        if let Some(ranges) = matches.values_of("range") {
            let ranges: Vec<&str> = ranges.collect();
            let prefetch_ranges: Vec<&str> = matches
                .values_of("prefetch-range")
                .and_then(|t| Some(t.collect()))
                .unwrap_or(vec![]);
            /*let ranged_zip = if let Some(prefetch_ranges) = prefetch_ranges {
                ranges.into_iter().zip(prefetch_ranges)
            } else {
                ranges.into_iter().zip(ranges.clone())
            };*/
            for eob in ranges.into_iter().zip_longest(prefetch_ranges) {
                let (range, prefetch_str) = match eob {
                    Both(a, b) => (a, b),
                    Left(a) => (a, a),
                    _ => panic!("Range is not specified."),
                };
                let string_range = StringRegion::new(range).unwrap();
                let prefetch_range = StringRegion::new(prefetch_str).unwrap();

                eprintln!("Retrieving: {} (prefetch: {})", range, prefetch_range);
                // let start = Instant::now();
                let closure = |x: &str| reader.reference_id(x);

                let _reference_name = &string_range.path;

                let range = Region::convert(&prefetch_range, closure).unwrap();
                let viewer = reader.fetch(&range).unwrap();

                let sample_ids_opt: Option<Vec<u64>> = matches
                    .values_of("id")
                    //.unwrap()
                    .and_then(|a| Some(a.map(|t| t.parse::<u64>().unwrap()).collect()));
                let sample_id_cond = sample_ids_opt.is_some();
                let sample_ids = sample_ids_opt.unwrap_or(vec![]);
                let filter = matches.is_present("filter");

                let format_type_opt = matches.value_of_t::<Format>("type");
                let format_type_cond = format_type_opt.is_ok();
                let format_type = format_type_opt.unwrap_or(Format::Default(Default {}));
                /*let out = std::io::stdout();
                let out_writer = match matches.value_of("output") {
                    Some(x) => {
                        let path = Path::new(x);
                        Box::new(File::create(&path).unwrap()) as Box<dyn Write>
                    }
                    None => Box::new(out.lock()) as Box<dyn Write>,
                };
                let output = io::BufWriter::with_capacity(1048576, out_writer);*/

                /*let end0 = start.elapsed();
                eprintln!(
                    "{}.{:03} sec.",
                    end0.as_secs(),
                    end0.subsec_nanos() / 1_000_000
                );*/
                let mut list = vec![];
                let mut ann = vec![];
                //let mut list2 = vec![];
                //let mut samples = BTreeMap::new();
                let _ = viewer.into_iter().for_each(|t| {
                    //eprintln!("{:?}", t.clone().unwrap());
                    let f = t.unwrap();
                    if !sample_id_cond || sample_ids.iter().any(|&i| i == f.sample_id()) {
                        let sample_id = f.sample_id();
                        let data = f.data();
                        if !format_type_cond
                            || std::mem::discriminant(&format_type) == std::mem::discriminant(&data)
                        {
                            match data {
                                Format::Range(rec) => {
                                    for i in rec.to_record(&prefetch_range.path) {
                                        if !filter
                                            || (i.end() as u64 > range.start()
                                                && range.end() > i.start() as u64)
                                        {
                                            ann.push((sample_id, i))
                                        }
                                    }
                                }
                                Format::Alignment(Alignment::Object(rec)) => {
                                    for i in rec {
                                        if !filter
                                            || (i.calculate_end() as u64 > range.start()
                                                && range.end() > i.start() as u64)
                                        {
                                            if !i.flag().is_secondary()
                                                && i.query_len() > min_read_len
                                            {
                                                list.push((sample_id, i));
                                            }
                                        }
                                    }
                                }
                                _ => {}
                            }
                        }
                    }
                });
                bam_record_vis_pre_calculate(
                    matches,
                    string_range,
                    prefetch_range,
                    &args,
                    list,
                    ann,
                    BTreeMap::new(),
                    threads,
                    |idx| reader.header().get_name(idx).and_then(|t| Some(t.as_str())),
                )?;
                /*
                if matches.is_present("web") {
                    server(
                        matches.clone(),
                        string_range,
                        prefetch_range,
                        args.clone(),
                        list,
                        ann,
                        BTreeMap::new(),
                        threads,
                    )?;
                } else {
                    bam_record_vis(matches, string_range, list, ann, BTreeMap::new(), |idx| {
                        reader.header().get_name(idx).and_then(|t| Some(t.as_str()))
                    })?;
                }*/
            }
        }
    }
    Ok(())
}

pub fn bin(matches: &ArgMatches, threads: u16) -> () {
    if let Some(o) = matches.value_of("INPUT") {
        let mut reader: IndexedReader<BufReader<File>> =
            IndexedReader::from_path_with_additional_threads(o, threads - 1).unwrap();
        let closure = |x: &str| reader.reference_id(x);
        let mut chunks = vec![];
        if let Some(range) = matches.value_of("range") {
            let string_range = StringRegion::new(range).unwrap();
            let range = Region::convert(&string_range, closure).unwrap();

            let index = matches
                .value_of("layer")
                .and_then(|a| a.parse::<usize>().ok())
                .unwrap_or(0usize);
            {
                for (i, slice) in reader.index().references()[range.ref_id() as usize]
                    .region_to_bins(range)
                    .enumerate()
                {
                    if index == i {
                        if let Some(bin) = slice.slice {
                            for chunk in bin {
                                chunks.extend(chunk.chunks().iter().copied());
                            }
                        }
                    }
                    /*let chunk = reader.index().references()[range.ref_id() as usize].bins()[bin_id]
                        .chunks()
                        .clone();
                    let mut res = Vec::new();
                    for i in chunks {
                        res.push(i);
                    }
                    chunks.extend(res);*/
                }
                // chunks = vec![chunks[index]];
            }
        } else {
            if let Some(bin_id) = matches.value_of("range") {
                let _ref_id = matches
                    .value_of("ref_id")
                    .and_then(|t| t.parse::<usize>().ok())
                    .unwrap();
                let _bin_id = bin_id.parse::<usize>().unwrap();
            /*let chunk = reader.index().references()[ref_id].bins()[bin_id].chunks();
            let mut res = Vec::new();
            for i in chunks {
                res.push(i);
            }
            chunks.extend(res);*/
            } else {
                let start = matches
                    .value_of("bin")
                    .and_then(|t| t.parse::<u64>().ok())
                    .unwrap();
                let end = matches
                    .value_of("end")
                    .and_then(|t| t.parse::<u64>().ok())
                    .unwrap();
                chunks = vec![Chunk::new(
                    0,
                    0,
                    VirtualOffset::from_raw(start),
                    VirtualOffset::from_raw(end),
                )];
            }
        }
        let viewer = reader.chunk(chunks).unwrap();

        let sample_ids_opt: Option<Vec<u64>> = matches
            .values_of("id")
            .and_then(|a| Some(a.map(|t| t.parse::<u64>().unwrap()).collect()));
        let sample_id_cond = sample_ids_opt.is_some();
        let sample_ids = sample_ids_opt.unwrap_or(vec![]);

        let format_type_opt = matches.value_of_t::<Format>("type");
        let format_type_cond = format_type_opt.is_ok();
        let format_type = format_type_opt.unwrap_or(Format::Default(Default {}));

        debug!("{:?} {:?} {:?}", sample_id_cond, sample_ids, format_type);
        let mut output = io::BufWriter::with_capacity(1048576, io::stdout());
        let header = viewer.header().clone();

        let _ = viewer.into_iter().for_each(|t| {
            debug!("{:?}", t);
            if let Ok(f) = t {
                debug!("{:?}", f);
                if !sample_id_cond || sample_ids.iter().any(|&i| i == f.sample_id()) {
                    let sample_id = f.sample_id();
                    let data = f.data();
                    debug!("{:?}", data);
                    if !format_type_cond
                        || std::mem::discriminant(&format_type) == std::mem::discriminant(&data)
                    {
                        match data {
                            Format::Range(rec) => {
                                let mut writer = bed::Writer::new(&mut output);
                                for i in rec.to_record("null") {
                                    writer.write(&i).unwrap();
                                }
                            }
                            Format::Alignment(Alignment::Object(rec)) => {
                                for i in rec {
                                    let _result = i
                                        .write_sam(
                                            &mut output,
                                            header
                                                .get_local_header(sample_id as usize)
                                                .unwrap()
                                                .bam_header(),
                                        )
                                        .unwrap();
                                }
                            }
                            _ => {}
                        }
                    }
                }
            }
        });
    }
}

pub fn bam_record_vis_pre_calculate<'a, F>(
    matches: &ArgMatches,
    range: StringRegion,
    prefetch_range: StringRegion,
    args: &Vec<String>,
    mut list: Vec<(u64, Record)>,
    mut ann: Vec<(u64, bed::Record)>,
    mut freq: BTreeMap<u64, Vec<(u64, u32)>>,
    threads: u16,
    lambda: F,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: Fn(usize) -> Option<&'a str>,
{
    /*let threads = matches
    .value_of("threads")
    .and_then(|t| t.parse::<u16>().ok())
    .unwrap_or(1u16);*/
    let pileup = matches.is_present("pileup");
    let split_only = matches.is_present("only-split-alignment");
    let sort_by_name = matches.is_present("sort-by-name");
    let packing = !matches.is_present("no-packing");
    let split = matches.is_present("split-alignment");
    let max_coverage = matches
        .value_of("max-coverage")
        .and_then(|a| a.parse::<u32>().ok());
    // Calculate coverage; it won't work on sort_by_name
    // let mut frequency = BTreeMap::new(); // Vec::with_capacity();

    if pileup {
        list.iter().group_by(|elt| elt.0).into_iter().for_each(|t| {
            let mut line =
                Vec::with_capacity((prefetch_range.end - prefetch_range.start + 1) as usize);
            for column in bam::Pileup::with_filter(&mut RecordIter::new(t.1), |record| {
                record.flag().no_bits(1796)
            }) {
                let column = column.unwrap();
                /*println!("Column at {}:{}, {} records", column.ref_id(),
                column.ref_pos() + 1, column.entries().len());*/
                // Should we have sparse occurrence table?
                // eprintln!("{:?} {:?}",  range.path, lambda(column.ref_id() as usize).unwrap_or(&column.ref_id().to_string()));
                // lambda(column.ref_id() as usize).unwrap_or(&column.ref_id().to_string())
                // == range.path
                // &&
                if prefetch_range.start <= column.ref_pos() as u64
                    && column.ref_pos() as u64 <= prefetch_range.end
                {
                    line.push((column.ref_pos() as u64, column.entries().len() as u32));
                }
            }
            // eprintln!("{:?}", line);
            freq.insert(t.0, line);
        });
    }
    ann.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
    list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));

    if sort_by_name {
        list.sort_by(|a, b| {
            a.0.cmp(&b.0)
                .then(a.1.name().cmp(&b.1.name()))
                .then(a.1.start().cmp(&b.1.start()))
        });
    }

    // Packing for each genome
    let mut prev_index = 0;
    let mut last_prev_index = 0;
    //let mut compressed_list = BTreeMap::<u64, usize>::new();
    let mut compressed_list = vec![];
    let mut index_list = Vec::with_capacity(list.len());
    let mut supplementary_list = vec![];
    if split {
        let mut end_map = HashMap::new();

        let new_list = {
            list.sort_by(|a, b| {
                a.0.cmp(&b.0)
                    .then(a.1.name().cmp(&b.1.name()))
                    .then(a.1.start().cmp(&b.1.start()))
            });
            list.clone()
        };
        new_list
            .iter()
            .group_by(|elt| elt.0)
            .into_iter()
            .for_each(|t| {
                let sample_id = t.0.clone();
                t.1.group_by(|elt| elt.1.name()).into_iter().for_each(|s| {
                    let items: Vec<&(u64, Record)> = s.1.into_iter().collect();
                    if items.len() > 1 {
                        let last: &(u64, Record) =
                            items.iter().max_by_key(|t| t.1.calculate_end()).unwrap();
                        end_map.insert(
                            (sample_id, s.0),
                            (
                                items[0].1.calculate_end(),
                                last.1.start(),
                                last.1.calculate_end(),
                                items.len(),
                            ),
                        );
                    }
                })
                //group.into_iter().for_each(|t| {})
            });

        if sort_by_name {
            if false {
                // sort_by_cigar {
                list.sort_by(|a, b| {
                    a.0.cmp(&b.0).then(a.1.name().cmp(&b.1.name())).then(
                        (a.1.cigar().soft_clipping(!a.1.flag().is_reverse_strand())
                            + a.1.cigar().hard_clipping(!a.1.flag().is_reverse_strand()))
                        .cmp(
                            &((b.1.cigar().soft_clipping(!b.1.flag().is_reverse_strand()))
                                + (b.1.cigar().hard_clipping(!b.1.flag().is_reverse_strand()))),
                        ),
                    )
                });
            } else {
                list.sort_by(|a, b| {
                    /*
                    a.0.cmp(&b.0).then(a.1.name().cmp(&b.1.name())).then(
                        /*a.1.cigar()
                        .soft_clipping(!a.1.flag().is_reverse_strand())
                        .cmp(&b.1.cigar().soft_clipping(!a.1.flag().is_reverse_strand())),*/
                        a.1.aligned_query_start().cmp(&b.1.aligned_query_start()),
                    )*/
                    a.0.cmp(&b.0).then(a.1.name().cmp(&b.1.name())).then(
                        (a.1.cigar().soft_clipping(true) + a.1.cigar().hard_clipping(true)).cmp(
                            &((b.1.cigar().soft_clipping(true))
                                + (b.1.cigar().hard_clipping(true))),
                        ),
                    )
                });
            }
        } else {
            list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
        }

        if split_only {
            list = list
                .into_iter()
                .filter(|(sample_id, record)| end_map.contains_key(&(*sample_id, record.name())))
                .collect();
        }

        list.iter().group_by(|elt| elt.0).into_iter().for_each(|t| {
            // let mut heap = BinaryHeap::<(i64, usize)>::new();
            let mut packing_vec = vec![0u64];
            let mut name_index = HashMap::new();
            prev_index += 1;
            let sample_id = t.0;
            (t.1).enumerate().for_each(|(e, k)| {
                let end = if !packing {
                    range.end() as i32
                } else if let Some(end) = end_map.get(&(sample_id, k.1.name())) {
                    end.2
                } else {
                    k.1.calculate_end()
                };

                let mut index = if sort_by_name {
                    prev_index += 1;
                    e
                } else if let Some(index) = name_index.get(k.1.name()) {
                    *index
                } else if let Some(index) = packing_vec
                    .iter_mut()
                    .enumerate()
                    .find(|(_, item)| **item < k.1.start() as u64)
                {
                    *index.1 = end as u64;
                    index.0
                } else {
                    packing_vec.push(end as u64);
                    prev_index += 1;
                    packing_vec.len() - 1
                };
                if let Some(end) = end_map.get(&(sample_id, k.1.name())) {
                    if None == name_index.get(k.1.name()) {
                        supplementary_list.push((
                            k.1.name(),
                            index + last_prev_index,
                            index + last_prev_index + end.3,
                            end.0,
                            end.1,
                        ));
                    }
                }
                if let Some(max_cov) = max_coverage {
                    if index > max_cov as usize {
                        index = max_cov as usize;
                    }
                }
                index_list.push(index + last_prev_index);
                name_index.insert(k.1.name(), index);
            });
            if let Some(max_cov) = max_coverage {
                prev_index = max_cov as usize + last_prev_index;
            }
            compressed_list.push((t.0, prev_index));
            last_prev_index = prev_index;
        });
    } else {
        if packing {
            list.iter().group_by(|elt| elt.0).into_iter().for_each(|t| {
                // let mut heap = BinaryHeap::<(i64, usize)>::new();
                let mut packing = vec![0u64];
                prev_index += 1;
                (t.1).for_each(|k| {
                    let mut index = if let Some(index) = packing
                        .iter_mut()
                        .enumerate()
                        .find(|(_, item)| **item < k.1.start() as u64)
                    {
                        //packing[index.0] = k.1.calculate_end() as u64;
                        *index.1 = k.1.calculate_end() as u64;
                        index.0
                    } else {
                        packing.push(k.1.calculate_end() as u64);
                        prev_index += 1;
                        packing.len() - 1
                        //prev_index - 1
                    }; /*
                       let index: usize = if heap.peek() != None
                           && -heap.peek().unwrap().0 < k.1.start() as i64
                       {
                           let hp = heap.pop().unwrap();
                           // let index = hp.1;
                           heap.push((-k.1.calculate_end() as i64, hp.1));
                           hp.1
                       } else {
                           let index = prev_index;
                           prev_index += 1;
                           heap.push((-k.1.calculate_end() as i64, index));
                           index
                       };*/
                    //let index =
                    if let Some(max_cov) = max_coverage {
                        if index > max_cov as usize {
                            index = max_cov as usize;
                            prev_index = max_cov as usize + last_prev_index;
                        }
                    }
                    index_list.push(index + last_prev_index);
                    // eprintln!("{:?}", packing);
                    //(index, (k.0, k.1))
                }); //.collect::<Vec<(usize, (u64, Record))>>()
                    // compressed_list.push(prev_index);
                    //compressed_list.insert(t.0, prev_index);
                    //prev_index += 1;
                if let Some(max_cov) = max_coverage {
                    prev_index = max_cov as usize + last_prev_index;
                }
                compressed_list.push((t.0, prev_index));
                //eprintln!("{:?} {:?} {:?}", compressed_list, packing, index_list);
                last_prev_index = prev_index;
                //(t.0, ((t.1).0, (t.1).1))
                // .collect::<&(u64, Record)>(). // collect::<Vec<(usize, (u64, Record))>>
            });
        } else {
            // Now does not specify the maximal length by max_coverage.
            index_list = (0..list.len()).collect();

            // list.sort_by(|a, b| a.0.cmp(&b.0));
            // eprintln!("{}", list.len());
            list.iter().group_by(|elt| elt.0).into_iter().for_each(
                |(sample_sequential_id, sample)| {
                    let count = sample.count();
                    prev_index += count;
                    // compressed_list.push(prev_index);
                    // compressed_list.insert(sample_sequential_id, prev_index);
                    compressed_list.push((sample_sequential_id, prev_index));
                },
            )
        }
    }
    eprintln!("{:?}", compressed_list);

    Ok(if matches.is_present("web") {
        server(
            matches.clone(),
            range,
            prefetch_range,
            args.clone(),
            list,
            ann,
            freq,
            compressed_list,
            index_list,
            prev_index,
            supplementary_list,
            threads,
        )?;
    } else {
        bam_record_vis(
            matches,
            range,
            list,
            &ann,
            &freq,
            &compressed_list,
            &index_list,
            prev_index,
            &supplementary_list,
            lambda,
        )?;
    })
}
