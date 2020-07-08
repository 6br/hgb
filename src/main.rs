extern crate log;

use bam;
use bam::{
    record::{
        tags::{IntegerType, StringType, TagValue},
        Cigar,
    },
    Record, RecordWriter,
};
use clap::{App, Arg, ArgMatches};
use env_logger;
use genomic_range::StringRegion;
use ghi::bed;
use itertools::Itertools;
use log::{debug, info};
use plotters::coord::ReverseCoordTranslate;
use plotters::prelude::Palette;
use plotters::prelude::*;
use rayon::prelude::*;
use std::{
    collections::{BTreeMap, BinaryHeap, HashMap},
    fs::{self, File},
    io,
    path::Path,
    time::Instant,
};

use ghi::binary::GhbWriter;
use ghi::builder::InvertedRecordBuilder;
use ghi::header::Header;
use ghi::index::{Chunk, Region, VirtualOffset};
use ghi::range::Default;
use ghi::range::{Format, InvertedRecordEntire, Set};
use ghi::twopass_alignment::{Alignment, AlignmentBuilder};
use ghi::writer::GhiWriter;
use ghi::{reader::IndexedReader, IndexWriter};
use io::{BufReader, Write};

fn main() {
    env_logger::init();
    let matches = App::new("GHB/GHI annotation/alignment database")
        .version("0.1")
        .author("6br. ")
        .about("Does awesome things")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple(true)
                .about("Sets the level of verbosity"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .default_value("1")
                .about("Sets the number of threads"),
        )
        .subcommand(
            App::new("build")
                .about("construct GHB and GHI index from bam/")
                .arg(
                    Arg::new("bam")
                        .short('a')
                        .takes_value(true)
                        .multiple(true)
                        .about("sorted bam"),
                )
                .arg(
                    Arg::new("bed")
                        .short('b')
                        .takes_value(true)
                        .multiple(true)
                        .about("sorted bed"),
                )
                .arg(
                    Arg::new("chrom")
                        .short('c')
                        .takes_value(true)
                        .about("chroms.sizes"),
                )
                .arg(Arg::new("force").short('f').about("Outputs only header"))
                .arg(
                    Arg::new("OUTPUT")
                        .about("Sets the output file to use")
                        .takes_value(true)
                        .required(true)
                        .index(1),
                ), /*.arg(Arg::new("gff")
                   .short('g')
                   .multiple(true)
                   .about("sorted gff"))*/
        )
        .subcommand(
            App::new("query")
                .about("construct GHB and GHI index from bam/")
                .arg(
                    Arg::new("range")
                        .short('r')
                        .takes_value(true)
                        .multiple(true)
                        .about("sorted bam"),
                )
                .arg(
                    Arg::new("bed")
                        .short('L')
                        .takes_value(true)
                        .about("sorted bam"),
                )
                .arg(
                    Arg::new("compression")
                        .short('c')
                        .takes_value(true)
                        .about("compression level"),
                )
                .arg(
                    Arg::new("type")
                        .short('t')
                        .takes_value(true)
                        // .default_value("default")
                        .possible_values(&["alignment", "range", "default"])
                        .about("annotation type to fetch"),
                )
                .arg(
                    Arg::new("id")
                        .short('i')
                        .takes_value(true)
                        .multiple(true)
                        .about("annotation sample to fetch"),
                )
                .arg(Arg::new("filter").short('f').about("Pre-filter"))
                .arg(Arg::new("binary").short('b').about("Binary"))
                .arg(Arg::new("vis").short('v').about("Binary"))
                .arg(Arg::new("no-cigar").short('n').about("Binary"))
                .arg(Arg::new("packing").short('p').about("Binary"))
                .arg(Arg::new("legend").short('l').about("Legend"))
                .arg(Arg::new("quality").short('q').about("Quality"))
                .arg(Arg::new("x").short('x').takes_value(true).about("x"))
                .arg(Arg::new("y").short('y').takes_value(true).about("y"))
                .arg(Arg::new("split-alignment").short('s').about("No-md"))
                .arg(Arg::new("no-insertion").short('z').about("No-md"))
                // .arg(Arg::new("insertion").short('').about("No-md"))
                .arg(
                    Arg::new("graph")
                        .short('g')
                        .takes_value(true)
                        .about("graph"),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .takes_value(true)
                        .about("Output format"),
                )
                .arg(
                    Arg::new("INPUT")
                        .about("Sets the input file to use")
                        .required(true)
                        .index(1),
                ),
        )
        .subcommand(
            App::new("decompose")
                .about("extract bam/bed from GHB and GHI index: (*) Now it needs to ")
                .arg(
                    Arg::new("id")
                        .short('i')
                        .takes_value(true)
                        //            .multiple(true)
                        .about("annotation sample to fetch (alignment | annotation)"),
                )
                .arg(Arg::new("header").short('z').about("Outputs only header"))
                .arg(
                    Arg::new("INPUT")
                        .about("Sets the input file to use")
                        .required(true)
                        .index(1),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .takes_value(true)
                        .about("Output format"),
                ),
        )
        .subcommand(
            App::new("bin")
                .about("extract bam/bed from GHB and GHI index: (*) Now it needs to ")
                .arg(
                    Arg::new("id")
                        .short('i')
                        .takes_value(true)
                        //            .multiple(true)
                        .about("annotation sample to fetch (alignment | annotation)"),
                )
                .arg(
                    Arg::new("INPUT")
                        .about("Sets the input file to use")
                        .required(true)
                        .index(1),
                )
                .arg(
                    Arg::new("bin")
                        .short('b')
                        .takes_value(true)
                        //            .multiple(true)
                        .about("annotation sample to fetch (alignment | annotation)"),
                )
                .arg(
                    Arg::new("start")
                        .short('s')
                        .takes_value(true)
                        //            .multiple(true)
                        .about("annotation sample to fetch (alignment | annotation)"),
                )
                .arg(
                    Arg::new("end")
                        .short('e')
                        .takes_value(true)
                        //            .multiple(true)
                        .about("annotation sample to fetch (alignment | annotation)"),
                )
                .arg(
                    Arg::new("ref_id")
                        .short('c')
                        .takes_value(true)
                        //            .multiple(true)
                        .about("annotation sample to fetch (alignment | annotation)"),
                )
                .arg(
                    Arg::new("range")
                        .short('r')
                        .takes_value(true)
                        .multiple(true)
                        .about("sorted bam"),
                )
                .arg(
                    Arg::new("layer")
                        .short('l')
                        .takes_value(true)
                        .about("sorted bam"),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .takes_value(true)
                        .about("Output format"),
                ),
        )
        .subcommand(
            App::new("server")
                .about("start the web server.")
                .arg(Arg::new("host").short('i').about("host"))
                .arg(Arg::new("port").short('p').about("port"))
                .arg(
                    Arg::new("INPUT")
                        .about("Sets the input file to use")
                        .required(true)
                        .index(1),
                ),
        )
        .subcommand(
            App::new("vis")
                .about("construct GHB and GHI index from bam/")
                .arg(
                    Arg::new("range")
                        .short('r')
                        .takes_value(true)
                        .multiple(true)
                        .about("sorted bam"),
                )
                .arg(
                    Arg::new("bed")
                        .short('L')
                        .takes_value(true)
                        .about("sorted bam"),
                )
                .arg(
                    Arg::new("type")
                        .short('t')
                        .takes_value(true)
                        // .default_value("default")
                        .possible_values(&["alignment", "range", "default"])
                        .about("annotation type to fetch"),
                )
                .arg(
                    Arg::new("id")
                        .short('i')
                        .takes_value(true)
                        .multiple(true)
                        .about("annotation sample to fetch"),
                )
                .arg(Arg::new("filter").short('f').about("Pre-filter"))
                .arg(Arg::new("binary").short('b').about("Binary"))
                .arg(Arg::new("no-cigar").short('n').about("Binary"))
                .arg(Arg::new("packing").short('p').about("Binary"))
                .arg(Arg::new("legend").short('l').about("Legend"))
                .arg(Arg::new("quality").short('q').about("Quality"))
                .arg(Arg::new("x").short('x').takes_value(true).about("x"))
                .arg(Arg::new("y").short('y').takes_value(true).about("y"))
                .arg(Arg::new("split-alignment").short('s').about("No-md"))
                .arg(Arg::new("no-insertion").short('z').about("No-md"))
                .arg(
                    Arg::new("graph")
                        .short('g')
                        .takes_value(true)
                        .about("graph"),
                )
                // .arg(Arg::new("insertion").short('').about("No-md"))
                .arg(
                    Arg::new("output")
                        .short('o')
                        .takes_value(true)
                        .about("Output format"),
                )
                .arg(
                    Arg::new("bam")
                        .short('a')
                        .takes_value(true)
                        .multiple(true)
                        .about("sorted bam"),
                ),
        )
        .subcommand(
            App::new("split")
                .about("start the web server.")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .takes_value(true)
                        .about("Output format"),
                )
                .arg(
                    Arg::new("INPUT")
                        .about("Sets the input file to use")
                        .required(true)
                        .index(1),
                ),
        )
        .get_matches();

    let threads = matches
        .value_of("threads")
        .and_then(|t| t.parse::<u16>().ok())
        .unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads as usize - 1)
        .build_global()
        .unwrap();

    if let Some(ref matches) = matches.subcommand_matches("build") {
        build(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("query") {
        match matches.is_present("binary") {
            true => bam_query(matches, threads),
            false => match matches.is_present("vis") {
                true => vis_query(matches, threads).unwrap(),
                false => query(matches, threads),
            },
        }
    } else if let Some(ref matches) = matches.subcommand_matches("decompose") {
        decompose(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("split") {
        split(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("bin") {
        bin(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("vis") {
        bam_vis(matches, threads).unwrap();
    }
}

fn bam_vis(matches: &ArgMatches, threads: u16) -> Result<(), Box<dyn std::error::Error>> {
    // let output_path = matches.value_of("OUTPUT").unwrap();

    if let Some(ranges) = matches.values_of("range") {
        let ranges: Vec<&str> = ranges.collect();
        for range in ranges {
            let string_range = StringRegion::new(range).unwrap();
            if let Some(bam_files) = matches.values_of("bam") {
                let bam_files: Vec<&str> = bam_files.collect();
                let mut list: Vec<(u64, Record)> = vec![];
                println!("Input file: {:?}", bam_files);
                for (index, bam_path) in bam_files.iter().enumerate() {
                    println!("Loading {}", bam_path);
                    // let reader = bam::BamReader::from_path(bam_path, threads).unwrap();
                    let mut reader2 = bam::IndexedReader::build()
                        .additional_threads(threads - 1)
                        .from_path(bam_path)
                        .unwrap();

                    // Here all threads can be used, but I suspect that runs double
                    //reader2.fetch()
                    let ref_id = reader2.header().reference_id(string_range.path.as_ref());
                    //                    let closure = |x: &str| reader.reference_id(x);
                    //                   let string_range = StringRegion::new(range).unwrap();
                    //                    let reference_name = &string_range.path;

                    //                    let range = Region::convert(&string_range, closure).unwrap();
                    let viewer = reader2
                        .fetch(&bam::bam_reader::Region::new(
                            ref_id.unwrap(),
                            string_range.start as u32,
                            string_range.end as u32,
                        ))
                        .unwrap();
                    for record in viewer {
                        list.push((index as u64, record.unwrap()));
                    }
                }
                bam_record_vis(matches, string_range, list, |idx| Some(bam_files[idx]))?;
            }
        }
    }
    Ok(())
}

fn build(matches: &ArgMatches, threads: u16) -> () {
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
                Set::<InvertedRecordBuilder, BufReader<File>>::new(reader, 1 as u64, &mut header)
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

fn query(matches: &ArgMatches, threads: u16) -> () {
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

fn decompose(matches: &ArgMatches, _threads: u16) -> () {
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

fn split(matches: &ArgMatches, threads: u16) -> () {
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
        let mut writer = bam::bam_writer::BamWriterBuilder::new()
            .additional_threads(threads - 1)
            .compression_level(clevel)
            .write_header(true)
            .from_stream(output, header.clone())
            .unwrap();

        let viewer = reader2.full();
        for record in viewer {}
    }
}

fn bam_query(matches: &ArgMatches, threads: u16) -> () {
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

fn vis_query(matches: &ArgMatches, threads: u16) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(o) = matches.value_of("INPUT") {
        let mut reader: IndexedReader<BufReader<File>> =
            IndexedReader::from_path_with_additional_threads(o, threads - 1).unwrap();
        if let Some(ranges) = matches.values_of("range") {
            let ranges: Vec<&str> = ranges.collect();
            for range in ranges {
                eprintln!("{}", range);
                let start = Instant::now();
                let closure = |x: &str| reader.reference_id(x);
                let string_range = StringRegion::new(range).unwrap();
                let _reference_name = &string_range.path;

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
                /*let out = std::io::stdout();
                let out_writer = match matches.value_of("output") {
                    Some(x) => {
                        let path = Path::new(x);
                        Box::new(File::create(&path).unwrap()) as Box<dyn Write>
                    }
                    None => Box::new(out.lock()) as Box<dyn Write>,
                };
                let output = io::BufWriter::with_capacity(1048576, out_writer);*/

                let end0 = start.elapsed();
                eprintln!(
                    "{}.{:03} sec.",
                    end0.as_secs(),
                    end0.subsec_nanos() / 1_000_000
                );
                let mut list = vec![];
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
                                            if !i.flag().is_secondary() {
                                                list.push((sample_id, i));
                                            }
                                            //list2.push((sample_id, true));
                                            //samples.insert(sample_id, true);
                                            //list.insert(sample_id, i);
                                        }
                                    }
                                }
                                _ => {}
                            }
                        }
                    }
                });
                bam_record_vis(matches, string_range, list, |idx| {
                    reader.header().get_name(idx).and_then(|t| Some(t.as_str()))
                })?;
            }
        }
    }
    Ok(())
}

fn bam_record_vis<'a, F>(
    matches: &ArgMatches,
    range: StringRegion,
    mut list: Vec<(u64, Record)>,
    lambda: F,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: Fn(usize) -> Option<&'a str>,
{
    let no_cigar = matches.is_present("no-cigar");
    let output = matches.value_of("output").unwrap();
    let packing = matches.is_present("packing");
    let quality = matches.is_present("quality");
    let legend = matches.is_present("legend");
    let insertion = !matches.is_present("no-insertion");
    let split = matches.is_present("split-alignment");
    let x = matches
        .value_of("x")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(1280u32);
    let y = matches
        .value_of("x")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(15u32);
    let graph = matches
        .value_of("graph")
        //.and_then(|a| fs::read_to_string(a).ok())
        .and_then(|a| {
            Some(
                csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .has_headers(false)
                    .from_reader(File::open(a).ok()?),
            )
        });
    //.and_then(|a| Some(a.split("\n").map(|t| t.split("\t").collect()).collect()));
    list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
    /*
    let iterator = if packing {
        // (0..).zip(list.iter())
        let mut prev_index = 0;
        list.iter().group_by(|elt| elt.0).into_iter().map(|t| {
            let mut heap = BinaryHeap::<(u64, usize)>::new();
            (t.1).map(|k| {
                let index: usize = if heap.peek() != None && heap.peek().unwrap().0 > k.1.start() as u64 {
                    let hp = heap.pop().unwrap();
                    // let index = hp.1;
                    heap.push((k.1.calculate_end() as u64, hp.1));
                    hp.1
                } else {
                    let index = prev_index;
                    prev_index += 1;
                    heap.push((k.1.calculate_end() as u64, index));
                    index
                };
                //let index =
                (index, (k.0, k.1))
            }) //.collect::<Vec<(usize, (u64, Record))>>()
            //(t.0, ((t.1).0, (t.1).1))
            // .collect::<&(u64, Record)>(). // collect::<Vec<(usize, (u64, Record))>>
            }
        ).flatten().collect::<Vec<(usize, (u64, Record))>>().into_iter()
    } else {
        list.into_iter().enumerate().collect::<Vec<(usize, (u64, Record))>>().into_iter()
    };*/

    // Packing for each genome
    let mut prev_index = 0;
    let mut last_prev_index = 0;
    //let mut compressed_list = BTreeMap::<u64, usize>::new();
    let mut compressed_list = vec![];
    let mut index_list = Vec::with_capacity(list.len());
    if packing {
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
                                /*items
                                .get_mut(0)
                                .unwrap()
                                .1
                                .tags_mut()
                                .push_num(b"PE", last.1.calculate_end());*/
                                // end_map.insert((t.0.clone(), s.0.clone()), last.1.calculate_end().clone());
                                end_map.insert((sample_id, s.0), last.1.calculate_end());
                            }
                        })
                        //group.into_iter().for_each(|t| {})
                    });
            
            list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));

            list.iter().group_by(|elt| elt.0).into_iter().for_each(|t| {
                // let mut heap = BinaryHeap::<(i64, usize)>::new();
                let mut packing = vec![0u64];
                let mut name_index = HashMap::new();
                prev_index += 1;
                let sample_id = t.0;
                (t.1).for_each(|k| {
                    let end = if let Some(end) = end_map.get(&(sample_id, k.1.name())) {
                        *end
                    } else {
                        k.1.calculate_end()
                    };

                    let index = if let Some(index) = name_index.get(k.1.name()) {
                        *index
                    } else if let Some(index) = packing
                        .iter_mut()
                        .enumerate()
                        .find(|(_, item)| **item < k.1.start() as u64)
                    {
                        *index.1 = end as u64;
                        index.0
                    } else {
                        packing.push(end as u64);
                        prev_index += 1;
                        packing.len() - 1
                    };
                    index_list.push(index + last_prev_index);
                    name_index.insert(k.1.name(), index);
                });
                compressed_list.push((t.0, prev_index));
                last_prev_index = prev_index;
            });
        } else {
            list.iter().group_by(|elt| elt.0).into_iter().for_each(|t| {
                // let mut heap = BinaryHeap::<(i64, usize)>::new();
                let mut packing = vec![0u64];
                prev_index += 1;
                (t.1).for_each(|k| {
                    let index = if let Some(index) = packing
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
                    index_list.push(index + last_prev_index);
                    // eprintln!("{:?}", packing);
                    //(index, (k.0, k.1))
                }); //.collect::<Vec<(usize, (u64, Record))>>()
                    // compressed_list.push(prev_index);
                    //compressed_list.insert(t.0, prev_index);
                    //prev_index += 1;
                compressed_list.push((t.0, prev_index));
                //eprintln!("{:?} {:?} {:?}", compressed_list, packing, index_list);
                last_prev_index = prev_index;
                //(t.0, ((t.1).0, (t.1).1))
                // .collect::<&(u64, Record)>(). // collect::<Vec<(usize, (u64, Record))>>
            });
        }
    } else {
        index_list = (0..list.len()).collect();
        // list.sort_by(|a, b| a.0.cmp(&b.0));
        // eprintln!("{}", list.len());
        list.iter()
            .group_by(|elt| elt.0)
            .into_iter()
            .for_each(|(sample_sequential_id, sample)| {
                let count = sample.count();
                prev_index += count;
                // compressed_list.push(prev_index);
                // compressed_list.insert(sample_sequential_id, prev_index);
                compressed_list.push((sample_sequential_id, prev_index));
            })
    }
    eprintln!("{:?}", compressed_list);
    /*
    let axis = if let Some(graph) = graph {
        graph.into_iter().collect().unique_by(|s| s[0]).count() // group_by(|elt| elt[0]).count() as usize
    } else {
        0usize
    };*/
    let axis = if let Some(mut graph) = graph {
        // let mut vec = vec![];
        let mut btree = BTreeMap::new();
        for i in graph.records() {
            let k = i?;
            //btree.insert(k[0].to_string(), (k[1].parse::<u64>()?, k[2].parse::<u64>()?));
            //eprintln!("{:?}", k);
            let item = btree.entry(k[0].to_string()).or_insert(vec![]);
            item.push((k[1].parse::<u64>()?, k[2].parse::<u64>()?));
            //vec.push((k[0].to_string(), k[1].parse::<u64>()?, k[2].parse::<u64>()?))
        }
        //btree
        btree
    //        graph.into_iter().collect().unique_by(|s| s[0]).count() // group_by(|elt| elt[0]).count() as usize
    } else {
        BTreeMap::new()
        //vec![]
    };
    let axis_count = axis.len() + 1; //into_iter().unique_by(|s| s.0).count();

    let root = BitMapBackend::new(
        output,
        (x, 40 + (prev_index as u32 + axis_count as u32) * y),
    )
    .into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(10, 10, 10, 10);
    // After this point, we should be able to draw construct a chart context
    // let areas = root.split_by_breakpoints([], compressed_list);
    let mut chart = ChartBuilder::on(&root)
        // Set the caption of the chart
        .caption(format!("{}", range), ("sans-serif", 20).into_font())
        // Set the size of the label region
        .x_label_area_size(20)
        .y_label_area_size(40)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_ranged((range.start() - 1)..(range.end() + 1), 0..(1 + prev_index))?;
    // Then we can draw a mesh
    chart
        .configure_mesh()
        // We can customize the maximum number of labels allowed for each axis
        //.x_labels(5)
        .y_labels(1)
        // We can also change the format of the label text
        .x_label_formatter(&|x| format!("{:.3}", x))
        .draw()?;

    let mut node_id_dict: BTreeMap<u64, (u64, u64)> = BTreeMap::new();
    let mut prev_pos = std::u64::MAX;
    let mut prev_node_id = 0u64;

    // We assume that axis[&range.path] includes all node id to be serialized.
    if let Some(axis) = axis.get(&range.path) {
        axis.iter()
            //.filter(|t| t.0 == range.path)
            .for_each(|(node_id, pos)| {
                if prev_pos > *pos {
                } else {
                    node_id_dict.insert(prev_node_id, (prev_pos, *pos));
                }
                prev_pos = *pos;
                prev_node_id = *node_id;
            });
        //node_id_dict[&prev_node_id] = (prev_pos, range.end());
        node_id_dict.insert(prev_node_id, (prev_pos, range.end()));
    }

    // Draw graph axis if there is graph information.
    axis.into_iter()
        //.group_by(|elt| elt.0.as_str())
        //.into_iter()
        .enumerate()
        .for_each(|(key, value)| {
            let path_name = value.0.clone();
            value.1.into_iter().for_each(|(node_id, _pos)| {
                let points = node_id_dict.get(&node_id); // [&node_id];
                if let Some(points) = points {
                    if points.1 > range.start() && points.0 < range.end() {
                        let start = if points.0 > range.start() {
                            points.0
                        } else {
                            range.start()
                        };
                        let end = if points.1 < range.end() {
                            points.1
                        } else {
                            range.end()
                        };
                        let stroke = Palette99::pick(node_id as usize);
                        let mut bar2 = Rectangle::new(
                            [(start, prev_index + key), (end, prev_index + key + 1)],
                            stroke.stroke_width(2),
                        );
                        bar2.set_margin(1, 1, 0, 0);
                        // prev_index += 1;
                        chart
                            .draw_series(LineSeries::new(
                                vec![(start, prev_index + key + 1), (end, prev_index + key + 1)],
                                stroke.stroke_width(y / 2),
                            ))
                            .unwrap()
                            .label(format!("{}: {}", path_name, node_id))
                            .legend(move |(x, y)| {
                                Rectangle::new([(x - 5, y - 5), (x + 5, y + 5)], stroke.filled())
                            });
                        //Some(bar2)
                    }
                }
            })
        });
    /*
    chart.draw_series(
        axis.into_iter()
            //.group_by(|elt| elt.0.as_str())
            //.into_iter()
            .enumerate()
            .map(|(key, value)| {
                value.1.into_iter().map(|(node_id, _pos)| {
                    let points = node_id_dict[&node_id];
                    if points.1 > range.start() && points.0 < range.end() {
                        let start = if points.0 > range.start() {
                            points.0
                        } else {
                            range.start()
                        };
                        let end = if points.1 < range.end() {
                            points.1
                        } else {
                            range.end()
                        };
                        let stroke = Palette99::pick(node_id as usize);
                        let mut bar2 = Rectangle::new(
                            [(start, prev_index + key), (end, prev_index + key + 1)],
                            stroke.stroke_width(2),
                        );
                        bar2.set_margin(1, 1, 0, 0);
                        // prev_index += 1;
                        Some(bar2)
                    } else {
                        None
                    }
                })
            })
            .flatten()
            .filter_map(|s| s)
    )?;
    */

    if legend {
        //list2.sort_by(|a, b| a.0.cmp(&b.0));
        // eprintln!("{}", list.len());
        // let mut prev_index = 0;

        for (sample_sequential_id, sample) in compressed_list.iter()
        // list2.into_iter().group_by(|elt| elt.0).into_iter()
        {
            // Check that the sum of each group is +/- 4.
            // assert_eq!(4, group.iter().fold(0_i32, |a, b| a + b).abs());
            let count = *sample; //.count();
            let idx = *sample_sequential_id as usize;
            // let idx = sample.next().0;
            chart
                .draw_series(LineSeries::new(
                    vec![(range.start() - 1, count), (range.end() + 1, count)],
                    &Palette99::pick(idx),
                ))?
                .label(format!("{}", lambda(idx).unwrap_or(&idx.to_string())))
                .legend(move |(x, y)| {
                    Rectangle::new(
                        [(x - 5, y - 5), (x + 5, y + 5)],
                        Palette99::pick(idx).filled(),
                    )
                });
            //prev_index += count;
        }
    }
    // For each sample:
    /*
        chart.draw_series(list2.into_iter().group_by(|elt| elt.0).into_iter().map(
            |(sample_sequential_id, sample)| {
                let count = sample.count();
                let stroke = Palette99::pick(sample_sequential_id as usize);
                let mut bar2 = Rectangle::new(
                    [
                        (range.start(), prev_index),
                        (range.end(), prev_index + count),
                    ],
                    stroke.stroke_width(100),
                );
                bar2.set_margin(1, 0, 0, 0);
                prev_index += count;
                bar2
            },
        ))?;
    */

    // For each alignment:

    let series = {
        //list.into_iter().enumerate().map(|(index, data)| {
        let mut bars = vec![];
        index_list.into_iter().zip(list).for_each(|(index, data)| {
            //chart.draw_series(index_list.into_par_iter().zip(list).map(|(index, data)| {
            //for (index, data) in list.iter().enumerate() {
            let bam = data.1;
            let color = if bam.flag().is_reverse_strand() {
                CYAN
            } else {
                RED
            };
            let stroke = Palette99::pick(data.0 as usize); //.unwrap(); //if data.0 % 2 == 0 { CYAN } else { GREEN };
            let start = if bam.start() as u64 > range.start() {
                bam.start() as u64
            } else {
                range.start()
            };
            let end = if bam.calculate_end() as u64 > range.end() {
                range.end()
            } else {
                bam.calculate_end() as u64
            };
            /*chart
            .draw_series(LineSeries::new(vec![(start, index), (end, index)], &color))?;*/
            let mut bar = Rectangle::new([(start, index), (end, index + 1)], color.filled());
            bar.set_margin(1, 1, 0, 0);

            // eprintln!("{:?}", [(start, index), (end, index + 1)]);
            bars.push(bar);
            //let mut bars =  //, bar2];
            if split {
                match bam.tags().get(b"SA") {
                    Some(TagValue::String(array_view, StringType::String)) => {
                        // assert!(array_view.int_type() == IntegerType::U32);
                        let current_left_clip =
                            bam.cigar().soft_clipping(!bam.flag().is_reverse_strand());
                        let sastr = String::from_utf8_lossy(array_view);
                        let sa: Vec<Vec<&str>> =
                            sastr.split(';').map(|t| t.split(',').collect()).collect();
                        let sa_left_clip: Vec<u32> = sa
                            .into_iter()
                            .filter(|t| t.len() > 2)
                            .map(|t| {
                                let strand = t[2];
                                //let cigar = Cigar::from_raw(t[3]).soft_clipping(strand == "+");
                                let mut cigar = Cigar::new();
                                cigar.extend_from_text(t[3].bytes()).ok()?;
                                //eprintln!("{:?}", cigar.soft_clipping(strand == "+"));

                                Some(cigar.soft_clipping(strand == "+"))
                            })
                            .filter_map(|t| t)
                            .collect();
                        let is_smaller = sa_left_clip
                            .iter()
                            .find(|t| t < &&current_left_clip)
                            .is_some();
                        let is_larger = sa_left_clip
                            .iter()
                            .find(|t| t > &&current_left_clip)
                            .is_some();
                        let color = RGBColor(120, 85, 43);
                        if ((is_smaller && !bam.flag().is_reverse_strand())
                            || (is_larger && bam.flag().is_reverse_strand()))
                            && bam.start() as u64 > range.start()
                        {
                            // split alignment on left
                            let mut bar = Rectangle::new(
                                [(start, index), (start + 1, index + 1)],
                                color.stroke_width(3), //.filled(),
                            );
                            bar.set_margin(0, 0, 0, 0);
                            bars.push(bar)
                        }
                        if ((is_larger && !bam.flag().is_reverse_strand())
                            || (is_smaller && bam.flag().is_reverse_strand()))
                            && (bam.calculate_end() as u64) < range.end()
                        {
                            // split alignment on right
                            let mut bar = Rectangle::new(
                                [(end - 1, index), (end, index + 1)],
                                color.stroke_width(2), //.filled(),
                            );
                            bar.set_margin(0, 0, 0, 0);
                            bars.push(bar)
                        }
                        /*eprintln!(
                            "SA = {}, {:?},
                             {}, {}, {}",
                            String::from_utf8_lossy(array_view),
                            sa_left_clip,
                            is_smaller,
                            is_larger,
                            bam.flag().is_reverse_strand()
                        );*/
                    }
                    // Some(TagValue::)
                    Some(_) => {} // panic!("Unexpected type"),
                    _ => {}
                }
            }
            if legend {
            } else {
                let mut bar2 =
                    Rectangle::new([(start, index), (end, index + 1)], stroke.stroke_width(1));
                bar2.set_margin(1, 1, 0, 0);
                //vec![bar,bar2]
                bars.push(bar2);
            };
            if !no_cigar {
                let mut prev_ref = bam.start() as u64;
                let mut prev_pixel_ref = start;
                let left_top = chart.as_coord_spec().translate(&(start, index));
                let right_bottom = chart.as_coord_spec().translate(&(end, index + 1));

                //if let Ok(mut a) = bam.alignment_entries() {
                match (quality, bam.alignment_entries()) {
                    (false, Ok(mut a)) => {
                        for i in left_top.0 + 1..right_bottom.0 {
                            /*let k = {
                                let from = chart.into_coord_trans();
                                from((i, left_top.1))
                            };*/
                            let k = chart.as_coord_spec().reverse_translate((i, left_top.1));
                            let mut color = None;
                            if let Some(k) = k {
                                while k.0 > prev_ref {
                                    let entry = a.next();
                                    if let Some(entry) = entry {
                                        if entry.is_insertion() {
                                            if prev_ref >= range.start() as u64 && insertion {
                                                //let mut bar = Rectangle::new([(prev_ref, index), (prev_ref+1, index + 1)], MAGENTA.stroke_width(1));
                                                //eprintln!("{:?}", [(prev_ref, index), (prev_ref + 1, index + 1)]);
                                                //bar.set_margin(0, 0, 0, 3);
                                                //bars.push(bar);
                                                //prev_ref = 0;
                                                color = Some(MAGENTA);
                                            }
                                        } else if entry.is_deletion() {
                                            prev_ref = entry.ref_pos_nt().unwrap().0 as u64;
                                            if prev_ref > range.end() as u64 {
                                                break;
                                            }
                                            if prev_ref >= range.start() as u64 {
                                                //let mut bar = Rectangle::new([(prev_ref , index), (prev_ref + 1, index + 1)], WHITE.filled());
                                                //bar.set_margin(2, 2, 0, 0);
                                                //eprintln!("{:?}", [(prev_ref, index), (prev_ref + 1, index + 1)]);
                                                //bars.push(bar);
                                                color = Some(WHITE);
                                            }
                                        } else if entry.is_seq_match() {
                                            prev_ref = entry.ref_pos_nt().unwrap().0 as u64;
                                            if prev_ref > range.end() as u64 {
                                                break;
                                            }
                                        } else {
                                            /* Mismatch */
                                            prev_ref = entry.ref_pos_nt().unwrap().0 as u64;

                                            if prev_ref > range.end() as u64 {
                                                break;
                                            }
                                            if prev_ref >= range.start() as u64 {
                                                let record_nt = entry.record_pos_nt().unwrap().1;
                                                color = match record_nt as char {
                                                    'A' => Some(GREEN), //RED,
                                                    'C' => Some(BLUE),  // BLUE,
                                                    'G' => Some(YELLOW),
                                                    'T' => Some(RED), //GREEN,
                                                    _ => Some(BLACK),
                                                };

                                                /*let mut bar =
                                                Rectangle::new([(prev_ref as u64, index), (prev_ref as u64 + 1, index + 1)], color.filled());
                                                bar.set_margin(2, 2, 0, 0);
                                                //eprintln!("{:?}", [(prev_ref, index), (prev_ref + 1, index + 1)]);
                                                bars.push(bar);*/
                                            }
                                        }
                                    }
                                }
                                if prev_ref >= range.start() as u64 {
                                    if let Some(color) = color {
                                        let mut bar = Rectangle::new(
                                            [
                                                (prev_pixel_ref as u64, index),
                                                (prev_ref as u64, index + 1),
                                            ],
                                            color.filled(),
                                        );
                                        bar.set_margin(2, 2, 0, 0);
                                        /*if prev_pixel_ref < start || end < prev_ref {
                                            eprintln!("{:?}", [(prev_pixel_ref, index), (prev_ref, index + 1)]);
                                        }*/
                                        bars.push(bar);
                                    }
                                    prev_pixel_ref = k.0;
                                }
                            }
                            /*if let Some((record_pos, record_nt)) = entry.record_pos_nt() {
                                print!("{} {}", record_pos, record_nt as char);
                            } else {
                                print!("-");
                            }
                            print!(", ");
                            if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
                                println!("{} {}", ref_pos, ref_nt as char);
                            } else {
                                println!("-");
                            }*/
                        }
                    }
                    _ => {
                        for entry in bam.aligned_pairs() {
                            match entry {
                                // (Seq_idx, ref_idx)
                                (Some(_record), Some(reference)) => {
                                    if prev_ref > range.start() as u64 && quality {
                                        if let Some(qual) =
                                            bam.qualities().raw().get(_record as usize)
                                        {
                                            //                                            eprintln!("{:?}", RGBColor(*qual*5, *qual*5, *qual*5));
                                            let mut bar = Rectangle::new(
                                                [
                                                    (reference as u64, index),
                                                    (reference as u64 + 1, index + 1),
                                                ],
                                                RGBColor(*qual * 5, *qual * 5, *qual * 5).filled(),
                                            );
                                            bar.set_margin(2, 2, 0, 0);
                                            bars.push(bar);
                                        }
                                    }
                                    prev_ref = reference as u64;
                                    if reference > range.end() as u32 {
                                        break;
                                    }
                                }
                                (None, Some(reference)) => {
                                    //Deletion
                                    if reference > range.start() as u32 && !quality {
                                        let mut bar = Rectangle::new(
                                            [
                                                (reference as u64, index),
                                                (reference as u64 + 1, index + 1),
                                            ],
                                            WHITE.filled(),
                                        );
                                        bar.set_margin(2, 2, 0, 0);
                                        prev_ref = reference as u64;
                                        bars.push(bar);
                                    }
                                    if reference >= range.end() as u32 {
                                        break;
                                    }
                                }
                                (Some(_record), None) => {
                                    //Insertion
                                    if prev_ref > range.start() as u64 && insertion {
                                        let mut bar = Rectangle::new(
                                            [(prev_ref, index), (prev_ref + 1, index + 1)],
                                            MAGENTA.stroke_width(1),
                                        );
                                        // eprintln!("{:?}", [(prev_ref, index), (prev_ref + 1, index + 1)]);
                                        bar.set_margin(0, 0, 0, 5);
                                        bars.push(bar);
                                        prev_ref = 0;
                                    }
                                    // eprintln!("{}", prev_ref)
                                }
                                _ => {}
                            }
                        }
                    }
                }
            }
            //}).flatten().collect::<Vec<_>>())?;
        });
        bars
    };
    chart.draw_series(series)?;

    if legend {
        chart
            .configure_series_labels()
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .draw()?;
    }
    Ok(())
}

fn bin(matches: &ArgMatches, threads: u16) -> () {
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
