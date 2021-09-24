use bam::{self, record::tags::TagValue};
use bam::{Record, RecordWriter};
use clap::ArgMatches;
use genomic_range::StringRegion;

#[cfg(feature = "web")]
use crate::buffered_server;
#[cfg(feature = "web")]
use crate::rest_server;
#[cfg(feature = "web")]
use crate::server::server;

use ghi::bed;
use ghi::binary::GhbWriter;
use ghi::builder::InvertedRecordBuilder;
use ghi::header::Header;
use ghi::index::{Chunk, Region, VirtualOffset};
use ghi::range::Default;
use ghi::range::{Format, InvertedRecordEntire, Set};
use ghi::twopass_alignment::{Alignment, AlignmentBuilder};
use ghi::vis::{bam_record_vis_orig, RecordIter};
use ghi::writer::GhiWriter;
use ghi::{
    gff, reader::IndexedReader, simple_buffer::ChromosomeBuffer, Frequency, IndexWriter, VisOrig,
};
use io::{BufReader, Error, ErrorKind, Write};
use itertools::EitherOrBoth::{Both, Left};
use itertools::Itertools;
use log::{debug, info, warn};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    convert::TryInto,
    fs::File,
    io,
    path::Path,
    sync::{Arc, Mutex},
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
    let neighbor = matches
        .value_of("neighbor")
        .and_then(|a| a.parse::<u64>().ok())
        .unwrap_or(0u64);
    let no_bits = matches
        .value_of("no-bits")
        .and_then(|t| t.parse::<u16>().ok())
        .unwrap_or(1796u16);
    let full_length = matches.is_present("full-length");
    let separated_by_tag = matches.occurrences_of("separated-by-tag") != 0;
    let separated_by_tag_vec = matches.value_of("separated-by-tag");
    let separated_by_tag_offset = matches
        .value_of("separated-by-tag-offset")
        .and_then(|a| a.parse::<usize>().ok());
    let bam_interval = if separated_by_tag {
        separated_by_tag_offset.unwrap()
    } else {
        1
    };
    let labels: Option<Vec<&str>> = matches.values_of("labels").map(|t| t.collect());

    if let Some(bam_files) = matches.values_of("bam") {
        let mut bam_files: Vec<&str> = bam_files.collect();
        let mut bam_readers = bam_files
            .iter()
            .map(|bam_path| {
                bam::IndexedReader::build()
                    .additional_threads(threads - 1)
                    .from_path(bam_path)
                    .unwrap()
            })
            .collect::<Vec<_>>();

        let mut ranges: Vec<String> = vec![];
        if let Some(bed_range) = matches.value_of("bed-range") {
            let mut reader = bed::Reader::from_file(bed_range)?;
            let mut values = vec![];
            for record in reader.records() {
                let record = record?;
                values.push(record);
            }
            let ranges_tmp: Vec<String> = if let Some(ranges_str) = matches.value_of("range") {
                //let ranges_tmp: Vec<String> =
                //    ranges_str.into_iter().map(|t| t.to_string()).collect();
                let string_range = StringRegion::new(&ranges_str).unwrap();
                values
                    .into_iter()
                    .filter(|record| {
                        record.chrom() == string_range.path
                            && (record.end() > string_range.start()
                                && string_range.end() > record.start())
                    })
                    .map(|record| {
                        format!(
                            "{}:{}-{}",
                            record.chrom(),
                            record.start() - neighbor,
                            record.end() + neighbor
                        )
                    })
                    .collect()
            } else {
                values
                    .into_iter()
                    .map(|record| {
                        format!(
                            "{}:{}-{}",
                            record.chrom(),
                            record.start() - neighbor,
                            record.end() + neighbor
                        )
                    })
                    .collect()
            };
            ranges.extend(ranges_tmp);
        } else if let Some(ranges_str) = matches.values_of("range") {
            let ranges_tmp: Vec<String> = ranges_str.into_iter().map(|t| t.to_string()).collect();
            ranges.extend(ranges_tmp);
        }
        let mut precursor = Vec::with_capacity(ranges.len());
        /*let prefetch_ranges: Vec<String> = matches
        .values_of("prefetch-range")
        .and_then(|t| Some(t.map(|t| t.to_string()).collect()))
        .unwrap_or(vec![]);*/
        let prefetch_ranges: Vec<&str> = matches
            .values_of("prefetch-range")
            .map(|t| t.collect())
            .unwrap_or_default();
        /*let ranged_zip = if let Some(prefetch_ranges) = prefetch_ranges {
            ranges.into_iter().zip(prefetch_ranges)
        } else {
            ranges.into_iter().zip(ranges.clone())
        };*/
        for eob in ranges.into_iter().zip_longest(prefetch_ranges) {
            let (range, prefetch_str) = match eob {
                Both(a, b) => (a, b.to_string()),
                Left(a) => (a.clone(), a),
                _ => panic!("Range is not specified."),
            };

            /*let prefetch_range = if let Some(prefetch_ranges) = prefetch_ranges {
                StringRegion::new(prefetch_ranges[index])?
            } else {
                string_range
            };*/
            let mut list: Vec<(u64, Record)> = vec![];
            println!("Input file: {:?}", bam_files);
            let reader2 = &mut bam_readers[0];
            let lambda = |range: String| {
                eprintln!("{:?}", range);
                let ref_id = reader2
                    .header()
                    .reference_id(&range)
                    .ok_or_else(|| Error::new(ErrorKind::Other, "Invalid reference id."))?;
                let end = reader2
                    .header()
                    .reference_len(ref_id)
                    .ok_or_else(|| Error::new(ErrorKind::Other, "Invalid reference id."))?;
                Ok(StringRegion {
                    path: range,
                    start: 1,
                    end: end as u64,
                })
            };
            let string_range = StringRegion::new(&range)
                .or_else(
                    |_| -> std::result::Result<
                        genomic_range::StringRegion,
                        Box<dyn std::error::Error>,
                    > { lambda(range) },
                )
                .unwrap();
            let prefetch_range = StringRegion::new(&prefetch_str)
                .or_else(
                    |_| -> std::result::Result<
                        genomic_range::StringRegion,
                        Box<dyn std::error::Error>,
                    > { lambda(prefetch_str) },
                )
                .unwrap();
            for (index, reader2) in &mut bam_readers.iter_mut().enumerate() {
                //println!("Loading {}", bam_path);
                // let reader = bam::BamReader::from_path(bam_path, threads).unwrap();
                //let mut reader2 = bam::IndexedReader::build()
                //    .additional_threads(threads - 1)
                //    .from_path(bam_path)?;

                // Here all threads can be used, but I suspect that runs double
                //reader2.fetch()
                let ref_id = reader2
                    .header()
                    .reference_id(prefetch_range.path.as_ref())
                    .ok_or_else(|| Error::new(ErrorKind::Other, "Invalid reference id."))?;
                let ref_len = reader2
                    .header()
                    .reference_len(ref_id)
                    .ok_or_else(|| Error::new(ErrorKind::Other, "Invalid reference length."))?;
                let start = if prefetch_range.start == 0 {
                    warn!("Invalid reference length: start == 0",);
                    1
                } else {
                    prefetch_range.start as u32
                };
                let end = if prefetch_range.end as u32 > ref_len {
                    warn!(
                        "Invalid reference length: end > reference length ({} > {})",
                        prefetch_range.end as u32, ref_len
                    );
                    ref_len
                } else {
                    prefetch_range.end as u32
                };
                let viewer = reader2.fetch(&bam::bam_reader::Region::new(ref_id, start, end))?;
                for record in viewer {
                    let record = record?;
                    if record.flag().no_bits(no_bits)
                        && record.query_len() > min_read_len
                        && (!full_length
                            || record.start() <= prefetch_range.start as i32
                                && record.calculate_end() >= prefetch_range.end as i32)
                    {
                        let idx = if separated_by_tag {
                            if let Some(colored_by_str) = separated_by_tag_vec {
                                if colored_by_str.is_empty() {
                                    let track = if record.flag().is_reverse_strand() {
                                        1
                                    } else {
                                        0
                                    };
                                    index * bam_interval + track
                                } else {
                                    let tag: &[u8; 2] = colored_by_str.as_bytes().try_into().expect("colored by tag with unexpected length: tag name must be two characters.");
                                    if let Some(TagValue::Int(tag_id, _)) = record.tags().get(tag) {
                                        index * bam_interval + tag_id as usize
                                    } else {
                                        index * bam_interval
                                    }
                                }
                            } else {
                                index * bam_interval
                            }
                        } else {
                            index * bam_interval
                        };
                        list.push((idx as u64, record));
                    }
                }
            }
            //    (string_range, prefetch_range)
            //};

            let mut ann = vec![];
            let mut idx = bam_files.len() * bam_interval;
            let mut freq = BTreeMap::new();
            if let Some(freq_files) = matches.values_of("frequency") {
                // let bed_files: Vec<_> = matches.values_of("bed").unwrap().collect();
                // frequency bed file needs to be (start, score).
                let mut freq_files: Vec<&str> = freq_files.collect();
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
                                    .and_then(|t| t.parse::<f32>().ok())
                                    .map(|t| t as u32)
                                    .unwrap_or_else(|| {
                                        record
                                            .name()
                                            .and_then(|t| t.parse::<f32>().ok())
                                            .map(|t| t as u32)
                                            .unwrap_or(0)
                                    }),
                                '*',
                            ));
                        }
                    }
                    freq.insert(idx as u64, values);
                    idx += 1;
                }
                // Need to change immutable
                bam_files.append(&mut freq_files);
            }
            //eprintln!("{:?}", freq);

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
                // bam_files.append(&mut bed_files);
            }
            if let Some(gff_files) = matches.values_of("gff3") {
                let gff_files: Vec<&str> = gff_files.collect();
                for (_idx, gff_path) in gff_files.iter().enumerate() {
                    info!("Loading {}", gff_path);
                    let mut reader = gff::Reader::from_file(gff_path, gff::GffType::GFF3).unwrap();
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
                // bam_files.append(&mut gff_files);
            }
            precursor.push(VisPrecursor::new(
                string_range,
                prefetch_range,
                list,
                ann,
                freq,
            ));
        }

        bam_record_vis_pre_calculate(matches, &args, precursor, threads, |idx| {
            //            if separated_by_tag {
            //                bam_files
            //                    .get(idx / bam_interval)
            //                    .and_then(|t| Some(format!("{}_{}", *t, idx % bam_interval).as_str()))
            //            } else {
            if let Some(labels) = &labels {
                labels.get(idx / bam_interval).copied()
            } else {
                bam_files.get(idx / bam_interval).copied()
            }
            //            }
        })?;
    }
    Ok(())
}

pub fn build(matches: &ArgMatches, threads: u16) {
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

pub fn query(matches: &ArgMatches, threads: u16) {
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
                    .map(|a| a.map(|t| t.parse::<u64>().unwrap()).collect());
                let sample_id_cond = sample_ids_opt.is_some();
                let sample_ids = sample_ids_opt.unwrap_or_default();
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

pub fn decompose(matches: &ArgMatches, _threads: u16) {
    if let Some(i) = matches.value_of("INPUT") {
        //if let Some(o) = matches.value_of("OUTPUT") {
        let header = matches.is_present("header");
        let formatted_header = matches.is_present("formatted-header");
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
        let id = matches.value_of("id").and_then(|t| t.parse::<u64>().ok());
        if let Some(id) = id {
            if let Some(header_type) = reader.header().get_local_header(id as usize) {
                if header {
                    header_type.to_text(&mut output).unwrap();
                } else if formatted_header {
                    header_type.to_tsv(&mut output).unwrap();
                }
            // header.write_text(&mut writer);
            } else {
                println!("There is no header of id {}", id);
            }
        } else {
            // Output global header
            if header {
                reader
                    .header()
                    .get_global_header()
                    .write_text(&mut output)
                    .unwrap();
            } else if formatted_header {
                reader.header().to_tsv(&mut output).unwrap();
            }
        }

        if header || formatted_header {
            return;
        }

        // todo!("Implement later; Now just returns only header.");

        let viewer = reader.full();
        let header_data = viewer.header().clone();

        let _ = viewer.into_iter().for_each(|t| {
            //eprintln!("{:?}", t);
            t.map(|f| {
                if f.sample_id() == id.unwrap() {
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
                                            .get_local_header(id.unwrap() as usize)
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

pub fn split(matches: &ArgMatches, threads: u16) {
    if let Some(bam_path) = matches.value_of("INPUT") {
        //todo!("Implement a web server using actix-web.");
        let mut reader2 = bam::IndexedReader::build()
            .additional_threads(threads - 1)
            .from_path(bam_path)
            .unwrap();
        // Here all threads can be used, but I suspect that runs double
        let bam_header = reader2.header();
        let header = matches.is_present("header");
        let formatted_header = matches.is_present("formatted-header");

        let out = std::io::stdout();
        let out_writer = match matches.value_of("output") {
            Some(x) => {
                let path = Path::new(x);
                Box::new(File::create(&path).unwrap()) as Box<dyn Write>
            }
            None => Box::new(out.lock()) as Box<dyn Write>,
        };
        let mut output = io::BufWriter::with_capacity(1048576, out_writer);

        if header {
            bam_header.write_text(&mut output).unwrap();
            return;
        } else if formatted_header {
            for (name, len) in bam_header
                .reference_names()
                .iter()
                .zip(bam_header.reference_lengths())
            {
                writeln!(&mut output, "{}\t{}", name, len).unwrap();
            }
            return;
        }
        let output_secondary_unmapped = matches.is_present("secondary-unmapped");

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
        let max_coverage = matches
            .value_of("max-coverage")
            .and_then(|a| a.parse::<u32>().ok());

        //let viewer = reader2.full();
        //for _record in viewer {}
        for (id, len) in header.clone().reference_lengths().iter().enumerate() {
            let viewer = reader2
                .fetch(&bam::Region::new(id as u32, 1, *len))
                .unwrap();
            let mut list = vec![];
            for record in viewer {
                let record = record.unwrap();
                if !record.flag().is_secondary() {
                    list.push((0, record));
                } else if output_secondary_unmapped {
                    writer.write(&record).unwrap();
                }
            }

            // Sort
            list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
            let mut prev_index = 0;
            let mut last_prev_index = 0;
            let mut index_list = Vec::with_capacity(list.len());
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
                    };
                    if let Some(max_cov) = max_coverage {
                        if index > max_cov as usize {
                            index = max_cov as usize;
                            prev_index = max_cov as usize + last_prev_index;
                        }
                    }
                    index_list.push(index + last_prev_index);
                    // eprintln!("{:?}", packing);
                });
                if let Some(max_cov) = max_coverage {
                    prev_index = max_cov as usize + last_prev_index;
                }
                last_prev_index = prev_index;
            });

            for ((_, mut record), index) in list.into_iter().zip(index_list) {
                record.tags_mut().push_num(b"YY", index as i32);
                writer.write(&record).unwrap();
            }
        }
    }
}

pub fn bam_query(matches: &ArgMatches, threads: u16) {
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
                    .map(|a| a.map(|t| t.parse::<u64>().unwrap()).collect());
                let sample_id_cond = sample_ids_opt.is_some();
                let sample_ids = sample_ids_opt.unwrap_or_default();
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
pub fn bench_query(matches: &ArgMatches, _args: Vec<String>, threads: u16) {
    if let Some(o) = matches.value_of("INPUT") {
        let mut reader: IndexedReader<BufReader<File>> =
            IndexedReader::from_path_with_additional_threads(o, threads - 1).unwrap();
        if let Some(ranges) = matches.values_of("range") {
            let ranges: Vec<&str> = ranges.collect();
            let prefetch_ranges: Vec<&str> = matches
                .values_of("prefetch-range")
                .map(|t| t.collect())
                .unwrap_or_default();

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
                let _viewer = reader.fetch(&range).unwrap();

                let mut buffer: ChromosomeBuffer = ChromosomeBuffer::new(reader, matches.clone());
                let mut list = vec![];
                let mut list_btree = BTreeSet::new();
                buffer.retrieve(&string_range, &mut list, &mut list_btree);
                let new_vis = buffer.vis(&string_range, &mut list, &mut list_btree);
                println!("{}", new_vis.unwrap().prefetch_max);
                break;
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
    let neighbor = matches
        .value_of("neighbor")
        .and_then(|a| a.parse::<u64>().ok())
        .unwrap_or(0u64);
    let no_bits = matches
        .value_of("no-bits")
        .and_then(|t| t.parse::<u16>().ok())
        .unwrap_or(1796u16);
    if let Some(o) = matches.value_of("INPUT") {
        let mut reader: IndexedReader<BufReader<File>> =
            IndexedReader::from_path_with_additional_threads(o, threads - 1)
                .map_err(|e| Error::new(e.kind(), format!("Failed to read GHB/GHI file: {}", e)))
                .unwrap();
        let mut ranges: Vec<String> = vec![];
        if let Some(bed_range) = matches.value_of("bed-range") {
            let mut reader = bed::Reader::from_file(bed_range)?;
            let mut values = vec![];
            for record in reader.records() {
                let record = record?;
                values.push(record);
                //                let range =
                //                    format!("{}:{}-{}", record.chrom(), record.start(), record.end()).to_string();
                //                ranges.push(&range.as_ref());
            }
            let ranges_tmp: Vec<String> = values
                .into_iter()
                .map(|record| {
                    format!(
                        "{}:{}-{}",
                        record.chrom(),
                        record.start() - neighbor,
                        record.end() + neighbor
                    )
                })
                .collect();
            ranges.extend(ranges_tmp);
        }
        if let Some(ranges_str) = matches.values_of("range") {
            let ranges_tmp: Vec<String> = ranges_str.into_iter().map(|t| t.to_string()).collect();
            ranges.extend(ranges_tmp);
        }
        let mut precursor = Vec::with_capacity(ranges.len());
        let prefetch_ranges: Vec<&str> = matches
            .values_of("prefetch-range")
            .map(|t| t.collect())
            .unwrap_or_default();

        /*let ranged_zip = if let Some(prefetch_ranges) = prefetch_ranges {
            ranges.into_iter().zip(prefetch_ranges)
        } else {
            ranges.into_iter().zip(ranges.clone())
        };*/
        for eob in ranges.into_iter().zip_longest(prefetch_ranges) {
            let (range, prefetch_str) = match eob {
                Both(a, b) => (a, b.to_string()),
                Left(a) => (a.clone(), a),
                _ => panic!("Range is not specified."),
            };
            let string_range = StringRegion::new(&range).unwrap();
            let prefetch_range = StringRegion::new(&prefetch_str).unwrap();

            eprintln!("Retrieving: {} (prefetch: {})", range, prefetch_range);
            // let start = Instant::now();
            let closure = |x: &str| reader.reference_id(x);

            let _reference_name = &string_range.path;

            let range = Region::convert(&prefetch_range, closure).unwrap();
            let viewer = reader.fetch(&range).unwrap();
            if matches.is_present("rest") {
                let buffer: ChromosomeBuffer = ChromosomeBuffer::new(reader, matches.clone());
                rest_server::server(
                    matches.clone(),
                    string_range,
                    prefetch_range,
                    args,
                    buffer,
                    threads,
                )?;
                return Ok(());
            } else if matches.is_present("whole-chromosome") && matches.is_present("web") {
                let buffer: ChromosomeBuffer = ChromosomeBuffer::new(reader, matches.clone());
                buffered_server::server(
                    matches.clone(),
                    string_range,
                    prefetch_range,
                    args,
                    buffer,
                    threads,
                )?;
                return Ok(());
            }

            let sample_ids_opt: Option<Vec<u64>> = matches
                .values_of("id")
                .map(|a| a.map(|t| t.parse::<u64>().unwrap()).collect());
            let sample_id_cond = sample_ids_opt.is_some();
            let sample_ids = sample_ids_opt.unwrap_or_default();
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
                                    if (!filter
                                        || (i.calculate_end() as u64 > range.start()
                                            && range.end() > i.start() as u64))
                                        && i.flag().no_bits(no_bits)
                                        && i.query_len() > min_read_len
                                    {
                                        list.push((sample_id, i));
                                    }
                                }
                            }
                            _ => {}
                        }
                    }
                }
            });
            precursor.push(VisPrecursor::new(
                string_range,
                prefetch_range,
                list,
                ann,
                BTreeMap::new(),
            ));
        }
        bam_record_vis_pre_calculate(matches, &args, precursor, threads, |idx| {
            reader.header().get_name(idx).map(|t| t.as_str())
        })?;
    }
    Ok(())
}

pub fn bin(matches: &ArgMatches, threads: u16) {
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
        } else if let Some(bin_id) = matches.value_of("range") {
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
        let viewer = reader.chunk(chunks).unwrap();

        let sample_ids_opt: Option<Vec<u64>> = matches
            .values_of("id")
            .map(|a| a.map(|t| t.parse::<u64>().unwrap()).collect());
        let sample_id_cond = sample_ids_opt.is_some();
        let sample_ids = sample_ids_opt.unwrap_or_default();

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

pub struct VisPrecursor {
    range: StringRegion,
    prefetch_range: StringRegion,
    list: Arc<Mutex<Vec<(u64, Record)>>>,
    annotation: Arc<Mutex<Vec<(u64, bed::Record)>>>,
    frequency: Arc<Mutex<Frequency>>,
    //list: Cell<Vec<(u64, Record)>>,
    //annotation: Cell<Vec<(u64, bed::Record)>>,
    //frequency: Cell<BTreeMap<u64, Vec<(u64, u32, char)>>>,
    index_list: Arc<Mutex<Vec<usize>>>,
    //suppl_list: Mutex<Vec<(Vec<u8>, usize, usize, usize, usize)>>
}

impl VisPrecursor {
    pub fn new(
        range: StringRegion,
        prefetch_range: StringRegion,
        list: Vec<(u64, Record)>,
        annotation: Vec<(u64, bed::Record)>,
        frequency: Frequency,
    ) -> Self {
        VisPrecursor {
            range,
            prefetch_range,
            list: Arc::new(Mutex::new(list)),
            annotation: Arc::new(Mutex::new(annotation)),
            frequency: Arc::new(Mutex::new(frequency)),
            index_list: Arc::new(Mutex::new(vec![])),
            //suppl_list: Mutex::new(vec![]),
        }
    }
}

pub fn bam_record_vis_pre_calculate<'a, F>(
    matches: &ArgMatches,
    args: &[String],
    vis: Vec<VisPrecursor>,
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
    //let range = &vis[0].range.clone();
    //let prefetch_range = &vis[0].prefetch_range.clone();

    let pileup = matches.is_present("pileup");
    let split_only = matches.is_present("only-split-alignment");
    let split_exclude = matches.is_present("exclude-split-alignment");
    let sort_by_name = matches.is_present("sort-by-name");
    let packing = !matches.is_present("no-packing");
    let split = matches.is_present("split-alignment");
    let read_per_line = matches.is_present("read-per-line");
    let read_per_two_node = matches.is_present("read-per-two-range");
    let min_read_len = matches
        .value_of("min-read-length")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(0u32);
    let no_bits = matches
        .value_of("no-bits")
        .and_then(|t| t.parse::<u16>().ok())
        .unwrap_or(1796u16);
    let read_index = matches
        .value_of("read-index")
        .and_then(|a| a.parse::<usize>().ok());

    let max_coverage = matches
        .value_of("max-coverage")
        .and_then(|a| a.parse::<u32>().ok());
    let snp_frequency = matches
        .value_of("snp-frequency")
        .and_then(|a| a.parse::<f64>().ok()); // default 0.2
                                              // Calculate coverage; it won't work on sort_by_name
                                              // let mut frequency = BTreeMap::new(); // Vec::with_capacity();

    {
        for i in vis.iter() {
            let mut list = i.list.lock().unwrap();
            let mut annotation = i.annotation.lock().unwrap();
            list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
            annotation.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
            /*eprintln!(
                "{:?}",
                list.iter()
                    .map(|t| String::from_utf8_lossy(t.1.name()))
                    .collect::<Vec<_>>(),
            );*/
        }
        //let mut list = &mut *vis[0].list.get_mut().unwrap();

        //list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
    }

    if pileup {
        //let mut freq_tmp = BTreeMap::new();
        for i in vis.iter() {
            let list = i.list.lock().unwrap();
            let mut freq = i.frequency.lock().unwrap();
            //let range = i.range; //.clone();
            let prefetch_range = &i.prefetch_range; //&vis[0].prefetch_range.clone();

            list.iter().group_by(|elt| elt.0).into_iter().for_each(|t| {
                //let no_bits = matches.value_of("no-bits").and_then(|t| t.parse::<u16>().ok()).unwrap_or(1796u16);
                let mut line =
                    Vec::with_capacity((prefetch_range.end - prefetch_range.start + 1) as usize);
                for column in bam::Pileup::with_filter(&mut RecordIter::new(t.1), |record| {
                    record.flag().no_bits(1796)
                }) {
                    let column = column.unwrap();
                    /*eprintln!(
                        "Column at {}:{}, {} records",
                        column.ref_id(),
                        column.ref_pos() + 1,
                        column.entries().len()
                    );*/
                    if let Some(freq) = snp_frequency {
                        let mut seqs = Vec::with_capacity(column.entries().len()); //vec![];
                        for entry in column.entries().iter() {
                            let seq: Option<_> = entry.sequence();
                            match seq {
                                Some(a) => {
                                    let seq: Vec<_> = a.map(|nt| nt as char).take(1).collect();
                                    seqs.push(seq);
                                }
                                _ => seqs.push(vec![]),
                            }
                            //let qual: Vec<_> = entry.qualities().unwrap().iter()
                            //    .map(|q| (q + 33) as char).collect();
                            //if column.ref_pos() == 8879824 {
                            //    eprintln!("    {:?}: {:?}", entry.record(), seq);
                            //}
                        }
                        let unique_elements = seqs.iter().cloned().unique().collect_vec();
                        let mut unique_frequency = vec![];
                        for unique_elem in unique_elements.iter() {
                            unique_frequency.push((
                                seqs.iter().filter(|&elem| elem == unique_elem).count(),
                                unique_elem,
                            ));
                        }
                        unique_frequency.sort_by_key(|t| t.0);
                        unique_frequency.reverse();
                        //eprintln!("{:?} {}", unique_frequency, column.entries().len() as u32);
                        let d: usize = unique_frequency.iter().map(|t| t.0).sum();
                        let threshold = d as f64 * freq;
                        // let minor = d - seqs[0].0;
                        // let second_minor = d - unique_frequency[1].0;
                        let (car, _cdr) = unique_frequency.split_first().unwrap();
                        /*cdr.iter()
                        .filter(|t| t.0 >= threshold as usize)
                        .for_each(|t2| {
                            if !t2.1.is_empty() {
                                line.push((
                                    column.ref_pos() as u64,
                                    t2.0 as u32,
                                    t2.1[0].to_ascii_uppercase(),
                                ));
                            }
                        });*/
                        if car.0 <= d - threshold as usize && !car.1.is_empty() {
                            line.push((
                                column.ref_pos() as u64,
                                car.0 as u32,
                                car.1[0].to_ascii_uppercase(),
                            ));
                        } else if car.1.is_empty() {
                            line.push((column.ref_pos() as u64, car.0 as u32, 'N'));
                        }
                    }
                    // Should we have sparse occurrence table?
                    //eprintln!("{:?} {:?}",  range.path, lambda(column.ref_id() as usize).unwrap_or(&column.ref_id().to_string()));
                    // lambda(column.ref_id() as usize).unwrap_or(&column.ref_id().to_string())
                    // == range.path
                    // &&
                    if prefetch_range.start <= column.ref_pos() as u64
                        && column.ref_pos() as u64 <= prefetch_range.end
                    {
                        line.push((column.ref_pos() as u64, column.entries().len() as u32, '*'));
                    }
                }
                //eprintln!("{:?}", line);

                //freq_tmp.insert(t.0, line);
                freq.insert(t.0, line);
            });
        }
        //let mut freq = &mut *vis[0].frequency.get_mut().unwrap();
        //freq.extend(freq_tmp);
    }

    //eprintln!("{:?}", freq.keys());
    if sort_by_name {
        for i in vis.iter() {
            let mut list = i.list.lock().unwrap();

            list.sort_by(|a, b| {
                a.0.cmp(&b.0)
                    .then(a.1.name().cmp(&b.1.name()))
                    .then(a.1.start().cmp(&b.1.start()))
            });
        }
    }

    // Packing for each genome
    let mut prev_index = 0;
    let mut last_prev_index = 0;
    //let mut compressed_list = BTreeMap::<u64, usize>::new();
    let mut compressed_list = vec![];
    //let mut index_list = Vec::with_capacity(vis[0].list.lock().unwrap().len());
    let mut supplementary_list = vec![];
    let mut suppl_map = HashMap::new();

    // New_list is merged list to decide the order of alignments.
    let mut new_list = {
        // new_list is a tuple (sample_id, record, range_id), whose record needs to be cloned.
        let mut new_list = vec![];
        for (index, i) in vis.iter().enumerate() {
            let list = i.list.lock().unwrap();
            /*list.sort_by(|a, b| {
                a.0.cmp(&b.0)
                    .then(a.1.name().cmp(&b.1.name()))
                    .then(a.1.start().cmp(&b.1.start()))
            });*/
            new_list.append(
                &mut list
                    .iter()
                    .map(|t| (t.0, t.1.clone(), index))
                    .collect::<Vec<_>>(),
            );
        }
        // Should we mark the duplicated alignment as "this is duplicated so you shouldn't display more than once"?

        new_list
    };

    {
        //let mut list = &mut *vis[0].list.get_mut().unwrap();
        if split {
            let mut end_map = HashMap::new();

            //Avoid immutable borrow occurs here.
            let mut tmp_list = new_list.clone();
            tmp_list.sort_by(|a, b| {
                a.0.cmp(&b.0)
                    .then(a.1.name().cmp(&b.1.name()))
                    .then(a.1.start().cmp(&b.1.start()))
            });
            //eprintln!("{:#?}", tmp_list);
            tmp_list
                .iter()
                .group_by(|elt| elt.0)
                .into_iter()
                .for_each(|t| {
                    let sample_id = t.0;
                    t.1.group_by(|elt| elt.1.name()).into_iter().for_each(|s| {
                        let mut items: Vec<&(u64, Record, usize)> = s.1.into_iter().collect();
                        if items.len() > 1 {
                            let last: &(u64, Record, usize) = items
                                .iter()
                                .max_by_key(|t| t.1.calculate_end() as u64 + ((t.2 as u64) << 32))
                                .unwrap();
                            let first: &(u64, Record, usize) = items
                                .iter()
                                .min_by_key(|t| t.1.calculate_end() as u64 + ((t.2 as u64) << 32))
                                .unwrap();
                            //eprintln!("Split: {:?}", last);
                            if !(split_only && read_per_two_node) || last.2 != first.2 {
                                end_map.insert(
                                    (sample_id, s.0),
                                    (
                                        items[0].1.calculate_end(),
                                        last.1.start(),
                                        last.1.calculate_end(),
                                        items.len(),
                                        last.2,
                                    ),
                                );
                            }
                            items.sort_by(|a, b| {
                                a.0.cmp(&b.0).then(a.1.name().cmp(&b.1.name())).then(
                                    (a.1.cigar().soft_clipping(!a.1.flag().is_reverse_strand())
                                        + a.1
                                            .cigar()
                                            .hard_clipping(!a.1.flag().is_reverse_strand()))
                                    .cmp(
                                        &((b.1
                                            .cigar()
                                            .soft_clipping(!b.1.flag().is_reverse_strand()))
                                            + (b.1
                                                .cigar()
                                                .hard_clipping(!b.1.flag().is_reverse_strand()))),
                                    ),
                                )
                            });
                            //suppl_map.insert((sample_id, s.0), )
                            for (idx, item) in items.iter().enumerate() {
                                if idx > 0 {
                                    // Add head
                                    suppl_map
                                        .entry((sample_id, s.0))
                                        .or_insert_with(Vec::new)
                                        .push(item.1.start());
                                }
                                if idx < items.len() - 1 {
                                    // Add tail
                                    suppl_map
                                        .entry((sample_id, s.0))
                                        .or_insert_with(Vec::new)
                                        .push(item.1.calculate_end());
                                }
                            }
                        }
                    })
                    //group.into_iter().for_each(|t| {})
                });

            if sort_by_name {
                if false {
                    // sort_by_cigar {
                    new_list.sort_by(|a, b| {
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
                    new_list.sort_by(|a, b| {
                        /*
                        a.0.cmp(&b.0).then(a.1.name().cmp(&b.1.name())).then(
                            /*a.1.cigar()
                            .soft_clipping(!a.1.flag().is_reverse_strand())
                            .cmp(&b.1.cigar().soft_clipping(!a.1.flag().is_reverse_strand())),*/
                            a.1.aligned_query_start().cmp(&b.1.aligned_query_start()),
                        )*/
                        a.0.cmp(&b.0).then(a.1.name().cmp(&b.1.name())).then(
                            (a.1.cigar().soft_clipping(true) + a.1.cigar().hard_clipping(true))
                                .cmp(
                                    &((b.1.cigar().soft_clipping(true))
                                        + (b.1.cigar().hard_clipping(true))),
                                ),
                        )
                    });
                }
            } else {
                //new_list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
                new_list.sort_by(|a, b| {
                    a.0.cmp(&b.0).then(
                        (a.1.start() as u64 + ((a.2 as u64) << 32))
                            .cmp(&(b.1.start() as u64 + ((b.2 as u64) << 32))),
                    )
                });
            }
            if split_only {
                new_list = new_list
                    .into_iter()
                    .filter(|(sample_id, record, _range_id)| {
                        end_map.contains_key(&(*sample_id, record.name()))
                    })
                    .collect();
                for i in vis.iter() {
                    let mut list = i.list.lock().unwrap();
                    *list = list
                        .clone()
                        .into_iter()
                        .filter(|(sample_id, record)| {
                            end_map.contains_key(&(*sample_id, record.name()))
                        })
                        .collect::<Vec<_>>();
                }
            }
            if split_exclude {
                new_list = new_list
                    .into_iter()
                    .filter(|(sample_id, record, _range_id)| {
                        !end_map.contains_key(&(*sample_id, record.name()))
                    })
                    .collect();
                for i in vis.iter() {
                    let mut list = i.list.lock().unwrap();
                    *list = list
                        .clone()
                        .into_iter()
                        .filter(|(sample_id, record)| {
                            !end_map.contains_key(&(*sample_id, record.name()))
                        })
                        .collect::<Vec<_>>();
                }
            }
            //eprintln!("{:#?}", new_list);

            new_list
                .iter()
                .group_by(|elt| elt.0)
                .into_iter()
                .for_each(|t| {
                    // let mut heap = BinaryHeap::<(i64, usize)>::new();
                    let mut packing_vec = vec![0u64];
                    let mut name_index = HashMap::new();
                    prev_index += 1;
                    let sample_id = t.0;
                    (t.1).enumerate().for_each(|(e, k)| {
                        let end: u64 = if !packing {
                            (vis.len() as u64 + 1) << 32 //range.end() as i32
                        } else if let Some(end) = end_map.get(&(sample_id, k.1.name())) {
                            // eprintln!("{:?} {:?}", sample_id, k.1.name());
                            // 1 line for only split alignments.
                            if read_per_line {
                                end.2 as u64 + ((vis.len() as u64) << 32)
                            } else if read_per_two_node {
                                end.2 as u64 + ((end.4 as u64 + 1u64) << 32)
                            } else {
                                end.2 as u64 + ((end.4 as u64) << 32)
                            }
                        } else if read_per_two_node {
                            k.1.calculate_end() as u64 + ((k.2 as u64 + 1u64) << 32)
                        } else {
                            k.1.calculate_end() as u64 + ((k.2 as u64) << 32)
                        };

                        let mut index = if sort_by_name {
                            prev_index += 1;
                            e
                        } else if let Some(index) = name_index.get(k.1.name()) {
                            *index
                        } else if let Some(index) = packing_vec
                            .iter_mut()
                            .enumerate()
                            .find(|(_, item)| **item < k.1.start() as u64 + ((k.2 as u64) << 32))
                        {
                            *index.1 = end;
                            index.0
                        } else {
                            packing_vec.push(end);
                            prev_index += 1;
                            packing_vec.len() - 1
                        };
                        if let Some(end) = end_map.get(&(sample_id, k.1.name())) {
                            if None == name_index.get(k.1.name()) {
                                supplementary_list.push((
                                    k.1.name().to_vec(),
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
                        vis[k.2]
                            .index_list
                            .lock()
                            .unwrap()
                            .push(index + last_prev_index);
                        name_index.insert(k.1.name(), index);
                        //eprintln!("{:?} {:?}", index, String::from_utf8_lossy(k.1.name()));
                    });
                    /*eprintln!(
                        "{:?}\n{:?}\n{:?}\n{:?}",
                        vis[0].index_list.lock().unwrap(),
                        vis[0]
                            .list
                            .lock()
                            .unwrap()
                            .iter()
                            .map(|t| String::from_utf8_lossy(t.1.name()))
                            .collect::<Vec<_>>(),
                        vis[1].index_list.lock().unwrap(),
                        vis[1]
                            .list
                            .lock()
                            .unwrap()
                            .iter()
                            .map(|t| String::from_utf8_lossy(t.1.name()))
                            .collect::<Vec<_>>(),
                    );*/
                    if let Some(max_cov) = max_coverage {
                        prev_index = max_cov as usize + last_prev_index;
                    }
                    compressed_list.push((t.0, prev_index));
                    last_prev_index = prev_index;
                });
        } else if packing {
            new_list
                .iter()
                .group_by(|elt| elt.0)
                .into_iter()
                .for_each(|t| {
                    // let mut heap = BinaryHeap::<(i64, usize)>::new();
                    let mut packing = vec![0u64];
                    prev_index += 1;
                    (t.1).for_each(|k| {
                        let mut index =
                            if let Some(index) = packing.iter_mut().enumerate().find(|(_, item)| {
                                **item < k.1.start() as u64 + ((k.2 as u64) << 32)
                            }) {
                                //packing[index.0] = k.1.calculate_end() as u64;
                                *index.1 = k.1.calculate_end() as u64 + ((k.2 as u64) << 32);
                                index.0
                            } else {
                                packing.push(k.1.calculate_end() as u64 + ((k.2 as u64) << 32));
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
                        vis[k.2]
                            .index_list
                            .lock()
                            .unwrap()
                            .push(index + last_prev_index);
                        // eprintln!("{:?}", packing);
                        //(index, (k.0, k.1))
                    }); //.collect::<Vec<(usize, (u64, Record))>>()
                        // compressed_list.push(prev_index);
                        //compressed_list.insert(t.0, prev_index);
                        //prev_index += 1;
                    if let Some(max_cov) = max_coverage {
                        prev_index = max_cov as usize + last_prev_index;
                    }
                    //compressed_list.push((t.0, prev_index));
                    //eprintln!("{:?} {:?} {:?}", compressed_list, packing, index_list);
                    last_prev_index = prev_index;
                    //(t.0, ((t.1).0, (t.1).1))
                    // .collect::<&(u64, Record)>(). // collect::<Vec<(usize, (u64, Record))>>
                });
            let mut tmp_interval = 0;
            new_list.iter().group_by(|elt| elt.0).into_iter().for_each(
                |(sample_sequential_id, sample)| {
                    let mut count = sample.count();
                    if let Some(max_cov) = max_coverage {
                        count = max_cov as usize;
                    }
                    tmp_interval += count;
                    //prev_index += count;
                    // compressed_list.push(prev_index);
                    // compressed_list.insert(sample_sequential_id, prev_index);
                    compressed_list.push((sample_sequential_id, tmp_interval));
                },
            )
        } else {
            eprintln!("Not packing, not split; multi samples are not supported.");
            // Now does not specify the maximal length by max_coverage.
            // index_list = (0..list.len()).collect();
            //let mut prev_index = 0;
            for i in vis.iter() {
                let mut index_list = i.index_list.lock().unwrap(); // = (0..new_list.len()).collect();
                                                                   //let temp_index_list = (prev_index..prev_index + i.list.lock().unwrap().len()).collect();
                                                                   //prev_index += i.list.lock().unwrap().len();
                let temp_index_list = (0..new_list.len()).collect();
                //mem::replace(&mut *index_list, temp_index_list);
                *index_list = temp_index_list;
            }

            // list.sort_by(|a, b| a.0.cmp(&b.0));
            // eprintln!("{}", list.len());
            //new_list.sort_by_key(|elt| elt.0);
            new_list.iter().group_by(|elt| elt.0).into_iter().for_each(
                |(_sample_sequential_id, sample)| {
                    let count = sample.count();
                    prev_index += count;
                    // compressed_list.push(prev_index);
                    // compressed_list.insert(sample_sequential_id, prev_index);
                    //compressed_list.push((sample_sequential_id, count));
                },
            )
        }
        for i in vis.iter() {
            if let Some(index) = read_index {
                let mut index_list = i.index_list.lock().unwrap();
                let read_indices: Vec<_> = index_list
                    .iter()
                    .enumerate()
                    .filter(|&(_, &value)| value == index)
                    .map(|(index, _)| index)
                    .collect();

                let temp_index_list = read_indices.iter().map(|_| 0).collect();
                *index_list = temp_index_list;
                let mut list = i.list.lock().unwrap();
                let temp_list = read_indices.iter().map(|t| list[*t].clone()).collect();

                *list = temp_list;
                prev_index = 1;
                let mut freq = i.frequency.lock().unwrap();
                let tmp_freq = BTreeMap::new();
                *freq = tmp_freq;
            }
        }
    }
    eprintln!("{:?}", compressed_list);

    if matches.is_present("web") {
        let list = &*vis[0].list.lock().unwrap();
        let ann = &*vis[0].annotation.lock().unwrap();
        let freq = &*vis[0].frequency.lock().unwrap();
        let index_list = &*vis[0].index_list.lock().unwrap();
        let range = &vis[0].range;
        let prefetch_range = &vis[0].prefetch_range; //.clone();
        server(
            matches.clone(),
            range.clone(),
            prefetch_range.clone(),
            args.to_owned(),
            list.to_vec(),
            ann.to_vec(),
            freq.clone(),
            compressed_list,
            index_list.to_vec(),
            prev_index,
            supplementary_list,
            threads,
        )?;
    } else {
        let mut vis_ref = Vec::with_capacity(vis.len());
        for i in vis {
            //let list = i.list.lock().unwrap();
            //let ann = i.annotation.lock().unwrap();
            //let freq = i.frequency.lock().unwrap();
            //let index_list = i.index_list.lock().unwrap();
            let vis_item = VisOrig::new(
                i.range,
                i.list.lock().unwrap().to_vec(),
                i.annotation.lock().unwrap().to_vec(),
                i.frequency.lock().unwrap().clone(),
                &compressed_list,
                i.index_list.lock().unwrap().to_vec(),
                prev_index,
                &supplementary_list,
            );
            //let vis_ref_item = vis_item.clone();
            vis_ref.push(vis_item);
        }
        /*let vis_ref = vis
        .into_iter()
        .map(|i| {
            let val = Arc::try_unwrap(i.list).unwrap().into_inner().unwrap();
            return VisRef::new(
                range.clone(),
                &val,
                &Arc::try_unwrap(i.annotation).unwrap().into_inner().unwrap(),
                &Arc::try_unwrap(i.frequency).unwrap().into_inner().unwrap(),
                &compressed_list,
                &Arc::try_unwrap(i.index_list).unwrap().into_inner().unwrap(),
                prev_index,
                &supplementary_list,
            );
        })
        .collect::<Vec<_>>();*/
        bam_record_vis_orig(
            matches, vis_ref,
            /*
            .map(|i| {
                let list = i.list.lock().unwrap();
                let ann = i.annotation.lock().unwrap();
                let freq = i.frequency.lock().unwrap();
                let index_list = i.index_list.lock().unwrap();
                let vis_item = VisOrig::new(
                    range.clone(),
                    list.to_vec(),
                    i.annotation.lock().unwrap().to_vec(),
                    i.frequency.lock().unwrap().clone(),
                    &compressed_list,
                    i.index_list.lock().unwrap().to_vec(),
                    prev_index,
                    &supplementary_list,
                );
                //let vis_ref_item = vis_item.clone();
                return vis_item.convert().clone();
            })
            .collect::<Vec<_>>(),*/
            lambda,
        )?;
    };
    Ok(())
}

#[cfg(not(feature = "web"))]
fn server(
    matches: ArgMatches,
    range: StringRegion,
    prefetch_range: StringRegion,
    args: Vec<String>,
    mut buffer: ChromosomeBuffer,
    threads: u16,
) -> std::io::Result<()> {
    unimplemented!("Please add web as a feature.")
}
