extern crate log;

use bio::io::bed;
use clap::{App, Arg, ArgMatches};
use log::info;
use std::fs::File;

use ghi::alignment::AlignmentBuilder;
use ghi::binary::GhbWriter;
use ghi::builder::InvertedRecordBuilder;
use ghi::header::Header;
use ghi::index::Region;
use ghi::range::{Format, InvertedRecordEntire, Set};
use ghi::writer::GhiWriter;
use ghi::{reader::IndexedReader, IndexWriter};

fn main() {
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
                        .possible_values(&["alignment", "range"])
                        .about("annotation type to fetch"),
                )
                .arg(
                    Arg::new("id")
                        .short('i')
                        .takes_value(true)
                        .multiple(true)
                        .about("annotation sample to fetch"),
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
                .arg(
                    Arg::new("INPUT")
                        .about("Sets the input file to use")
                        .required(true)
                        .index(1),
                )
                .arg(
                    Arg::new("OUTPUT")
                        .about("Sets the output file to use")
                        .required(true)
                        .index(1),
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
        query(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("decompose") {
        decompose(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("server") {
        server(matches, threads);
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

    let bam_files: Vec<_> = matches.values_of("bam").unwrap().collect();
    let mut set_vec = vec![];

    for (i, bam_path) in bam_files.iter().enumerate() {
        info!("Loading {}", bam_path);
        let reader = bam::BamReader::from_path(bam_path, threads).unwrap();
        let bam_header = reader.header();
        if alignment_transfer {
            header.transfer(bam_header).unwrap();
        }
        header.set_local_header(bam_header, i);
        let set = Set::<AlignmentBuilder>::new(reader, 0 as u64, &mut header);
        set_vec.push(set);
    }
    let mut entire = InvertedRecordEntire::new_from_set(set_vec);

    let bed_files: Vec<_> = matches.values_of("bed").unwrap().collect();

    for bed_path in bed_files {
        info!("Loading {}", bed_path);
        let reader = bed::Reader::from_file(bed_path).unwrap();
        let set: Set<InvertedRecordBuilder> =
            Set::<InvertedRecordBuilder>::new(reader, 1 as u64, &mut header).unwrap();

        entire.add(set);
    }

    let dummy_header = Header::new();
    let mut writer = GhbWriter::build()
        .write_header(false)
        .from_path(output_path, dummy_header)
        .unwrap();
    let index = entire.write_binary(&mut writer).unwrap();
    writer.flush().unwrap();
    let output_index_path = format!("{}.ghi", output_path);
    let mut index_writer = GhiWriter::build()
        .write_header(true)
        .from_path(output_index_path, header)
        .unwrap();
    let _result = index_writer.write(&index);
    assert_eq!(_result.ok(), Some(()));
    let _result = index_writer.flush();
    println!("Annotation saved!");
    // Some(())
}

fn query(matches: &ArgMatches, threads: u16) -> () {
    if let Some(o) = matches.value_of("INPUT") {
        let mut reader = IndexedReader::from_path(o).unwrap();
        if let Some(range) = matches.value_of("range") {
            let closure = |x: &str| reader.reference_id(x);
            let range = Region::parse(range, closure).unwrap();

            let viewer = reader.fetch(&range).unwrap();

            let sample_id = matches.value_of("id").unwrap();
            let format_type = matches.value_of_t("type").unwrap();

            match format_type {
                Format::Alignment(alignment) => {}
                Format::Range(range) => {}
                _ => println!("Format not matched."),
            }

            let records = viewer.into_iter().flat_map(|t| {
                t.map(|f| {
                    if let Format::Range(rec) = f.data() {
                        rec.to_record(range.ref_id())
                    } else if let Format::Alignment(rec) = f.data() {
                        rec.to_record(range.ref_id())
                    }
                })
            });
        }
    }
}

fn decompose(matches: &ArgMatches, threads: u16) -> () {
    if let Some(i) = matches.value_of("INPUT") {
        if let Some(o) = matches.value_of("OUTPUT") {
            let id = matches
                .value_of("id")
                .and_then(|t| t.parse::<usize>().ok())
                .unwrap();
            let mut reader = IndexedReader::from_path(i).unwrap();
            let mut writer = File::open(o).unwrap();
            if let Some(header) = reader.header().get_local_header(id) {
                header.to_stream(&mut writer);
            } else {
                println!("There is no header of id {}", id);
            }
            todo!("Implement later; Now just returns only header.");
            /*
            let viewer = reader.full();
            viewer.into_iter().flat_map(|t|

            )*/
        }
    }
}

fn server(matches: &ArgMatches, threads: u16) -> () {
    //todo!("Implement a web server using actix-web.");
}
