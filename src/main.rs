extern crate log;

use bam::RecordWriter;
use clap::{App, Arg, ArgMatches};
use env_logger;
use genomic_range::StringRegion;
use ghi::bed;
use log::{debug, info};
use std::{fs::File, io, path::Path};

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
                    Arg::new("type")
                        .short('t')
                        .takes_value(true)
                        // .default_value("default")
                        .possible_values(& ["alignment", "range", "default"])
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
        .get_matches();

    let threads = matches
        .value_of("threads")
        .and_then(|t| t.parse::<u16>().ok())
        .unwrap();

    if let Some(ref matches) = matches.subcommand_matches("build") {
        build(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("query") {
        match matches.is_present("binary") {
            true => bam_query(matches, threads),
            false => query(matches, threads),
        }
    } else if let Some(ref matches) = matches.subcommand_matches("decompose") {
        decompose(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("server") {
        server(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("bin") {
        bin(matches, threads);
    }
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
        if let Some(range) = matches.value_of("range") {
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
            let mut output = io::BufWriter::new(out_writer);

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
        let mut output = io::BufWriter::new(out_writer);

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

fn server(_matches: &ArgMatches, _threads: u16) -> () {
    //todo!("Implement a web server using actix-web.");
}

fn bam_query(matches: &ArgMatches, threads: u16) -> () {
    if let Some(o) = matches.value_of("INPUT") {
        let mut reader: IndexedReader<BufReader<File>> =
            IndexedReader::from_path_with_additional_threads(o, threads - 1).unwrap();
            if let Some(ranges) = matches.values_of("range") {
                let ranges: Vec<&str> = ranges.collect();
                for range in ranges {
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
                let out = std::io::stdout();
                let out_writer = match matches.value_of("output") {
                    Some(x) => {
                        let path = Path::new(x);
                        Box::new(File::create(&path).unwrap()) as Box<dyn Write>
                    }
                    None => Box::new(out.lock()) as Box<dyn Write>,
                };
                let output = io::BufWriter::new(out_writer);

                let header = viewer
                    .header()
                    .get_local_header(*sample_ids.get(0).unwrap_or(&0) as usize)
                    .unwrap()
                    .bam_header();
                let mut writer = bam::bam_writer::BamWriterBuilder::new()
                    .additional_threads(threads - 1)
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
                                            writer.write(&i).unwrap();
                                            //i.write_bam(&mut writer).unwrap();
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

fn bin(matches: &ArgMatches, threads: u16) -> () {
    if let Some(o) = matches.value_of("INPUT") {
        let mut reader: IndexedReader<BufReader<File>> =
            IndexedReader::from_path_with_additional_threads(o, threads - 1).unwrap();
        let closure = |x: &str| reader.reference_id(x);
        let mut chunks = vec![];
        if let Some(range) = matches.value_of("range") {
            let string_range = StringRegion::new(range).unwrap();
            let range = Region::convert(&string_range, closure).unwrap();
            {
                let _bin_id =
                    reader.index().references()[range.ref_id() as usize].region_to_bin(range);
                /*let chunk = reader.index().references()[range.ref_id() as usize].bins()[bin_id]
                    .chunks()
                    .clone();
                let mut res = Vec::new();
                for i in chunks {
                    res.push(i);
                }
                chunks.extend(res);*/
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
        let mut output = io::BufWriter::new(io::stdout());
        let header = viewer.header().clone();
        /*let mut writer = bam::BamWriter::build()
        .write_header(true)
        .from_stream(output, reader.header().clone()).unwrap();*/

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
