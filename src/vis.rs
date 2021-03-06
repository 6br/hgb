use crate::{bed, VisPreset};
use bam::record::{
    tags::{StringType, TagValue},
    Cigar,
};
use bam::{Record, RecordReader};
use bio_types::strand::Strand;
use clap::ArgMatches;
use genomic_range::StringRegion;
use itertools::Itertools;
use plotters::coord::ReverseCoordTranslate;
use plotters::prelude::Palette;
use plotters::prelude::*;
use plotters::style::RGBColor;
use std::ops::Range;
use udon::{Udon, UdonPalette, UdonScaler, UdonUtils};
// use rayon::prelude::*;
use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
};

//Copied from
trait RangeUtils
where
    Self: Sized,
{
    fn has_overlap(&self, query: &Self) -> bool;
    fn clip(&self, query: &Self) -> Option<Self>;
    fn scale(&self, divisor: f64) -> (Self, f64);
}

impl RangeUtils for Range<usize> {
    fn has_overlap(&self, query: &Range<usize>) -> bool {
        if query.start > self.end {
            return false;
        }
        if self.start > query.end {
            return false;
        }
        return true;
    }

    fn clip(&self, query: &Range<usize>) -> Option<Range<usize>> {
        /* clip query range by self (window). returns local range in the window */
        if query.start > self.end {
            return None;
        }
        if self.start > query.end {
            return None;
        }

        Some(Range::<usize> {
            start: query.start.saturating_sub(self.start),
            end: query.end.min(self.end).saturating_sub(self.start),
        })
    }

    fn scale(&self, divisor: f64) -> (Range<usize>, f64) {
        let start = self.start as f64 / divisor;
        let offset = start.fract();

        let range = Range::<usize> {
            start: start as usize,
            end: (self.end as f64 / divisor).ceil() as usize,
        };

        (range, offset)
    }
}

macro_rules! predefined_color {
    ($name:ident, $r:expr, $g:expr, $b:expr, $doc:expr) => {
        #[doc = $doc]
        pub const $name: RGBColor = RGBColor($r, $g, $b);
    };
/*
    ($name:ident, $r:expr, $g:expr, $b:expr, $a: expr, $doc:expr) => {
        #[doc = $doc]
        let $name: RGBAColor = RGBColor($r, $g, $b).mix($a);
    }*/
}

/*
predefined_color!(
    NEG_COL,
    150,
    150,
    230,
    /*0.75,*/
    "The predefined negative-str color"
);
predefined_color!(
    POS_COL,
    230,
    150,
    150,
    /*0.75,*/
    "The predefined positive-str color"
);
predefined_color!(INS_COL, 138, 94, 161, "The predefined insertion color");
predefined_color!(SPL_COL, 120, 85, 43, "The predefined split-alignment color");
predefined_color!(N_COL, 80, 80, 80, "The predefined N color");
predefined_color!(A_COL, 0, 200, 0, "The predefined A color");
predefined_color!(C_COL, 0, 0, 200, "The predefined C color");
predefined_color!(T_COL, 255, 0, 40, "The predefined T color");
predefined_color!(G_COL, 209, 113, 5, "The predefined G color");
*/

// Now I adapt nordtheme. https://www.nordtheme.com/
predefined_color!(
    NEG_COL,
    0x81, //0x5e,
    0xa1, //0x81,
    0xc1, //0xac,
    /*0.75,*/
    "The predefined negative-str color"
);
predefined_color!(
    POS_COL,
    0xbf,
    0x61,
    0x6a,
    /*0.75,*/
    "The predefined positive-str color"
);
//predefined_color!(INS_COL, 0xb4, 0x8e, 0xad, "The predefined insertion color");
predefined_color!(INS_COL, 153, 0, 153, "The predefined insertion color");
predefined_color!(
    SPL_COL,
    0x8f,
    0xbc,
    0xbb,
    "The predefined split-alignment color"
);
predefined_color!(N_COL, 80, 80, 80, "The predefined N color");
predefined_color!(A_COL, 119, 217, 168, "The predefined A color"); //#00c800
predefined_color!(C_COL, 0, 90, 255, "The predefined C color"); // #0000c8
predefined_color!(T_COL, 255, 75, 0, "The predefined T color"); // #ff0000
predefined_color!(G_COL, 128, 64, 0, "The predefined G color"); // 209,113,5

//RGBColor();

//newtype RecordIter<'a> = std::slice::Iter<'a, (u64, Record)>;
struct RecordIter<'a, I: Iterator<Item = &'a (u64, Record)>>(I);

impl<'a, I> RecordReader for RecordIter<'a, I>
where
    I: Iterator<Item = &'a (u64, Record)>,
{
    fn read_into(&mut self, record: &mut Record) -> std::io::Result<bool> {
        if let Some(next_record) = self.0.next() {
            // record = next_record.1;
            let _ = std::mem::replace(record, next_record.1.clone());
            Ok(true)
        } else {
            Ok(false)
        }
    }

    fn pause(&mut self) {}
}

/// Iterator over records.
impl<'a, I> Iterator for RecordIter<'a, I>
where
    I: Iterator<Item = &'a (u64, Record)>,
{
    type Item = std::io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Record::new();
        match self.read_into(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub fn frequency_vis<'a, F>(
    matches: &ArgMatches,
    range: StringRegion,
    mut list: Vec<(u64, Record)>,
    lambda: F,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: Fn(usize) -> Option<&'a str>,
{
    let output = matches.value_of("output").unwrap();
    let max_coverage = matches
        .value_of("max-coverage")
        .and_then(|a| a.parse::<u32>().ok());
    let x = matches
        .value_of("x")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(1280u32);
    let freq_size = matches
        .value_of("freq-height")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(100u32);

    list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
    // Calculate coverage; it won't work on sort_by_name
    let mut frequency = BTreeMap::new(); // Vec::with_capacity();
    list.iter().group_by(|elt| elt.0).into_iter().for_each(|t| {
        let mut line = Vec::with_capacity((range.end - range.start + 1) as usize);
        for column in
            bam::Pileup::with_filter(&mut RecordIter(t.1), |record| record.flag().no_bits(1796))
        {
            let column = column.unwrap();
            /*println!("Column at {}:{}, {} records", column.ref_id(),
            column.ref_pos() + 1, column.entries().len());*/
            // Should we have sparse occurrence table?
            // eprintln!("{:?} {:?}",  range.path, lambda(column.ref_id() as usize).unwrap_or(&column.ref_id().to_string()));
            // lambda(column.ref_id() as usize).unwrap_or(&column.ref_id().to_string())
            // == range.path
            // &&
            if range.start <= column.ref_pos() as u64 && column.ref_pos() as u64 <= range.end {
                line.push((column.ref_pos() as u64, column.entries().len() as u32));
            }
        }
        // eprintln!("{:?}", line);
        frequency.insert(t.0, line);
    });
    let root =
        BitMapBackend::new(output, (x, frequency.len() as u32 * freq_size)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(0, 0, 0, 0);
    // After this point, we should be able to draw construct a chart context
    // let areas = root.split_by_breakpoints([], compressed_list);
    if frequency.len() > 0 {
        let coverages = root.split_evenly((frequency.len(), 1));
        for (i, (sample_sequential_id, values)) in frequency.iter().rev().enumerate() {
            let idx = *sample_sequential_id as usize;
            let y_max = match max_coverage {
                Some(a) => a,
                _ => values.iter().map(|t| t.1).max().unwrap_or(1),
            };
            let mut chart = ChartBuilder::on(&coverages[i])
                // Set the caption of the chart
                //.caption(format!("{}", range), ("sans-serif", 20).into_font())
                // Set the size of the label region
                .x_label_area_size(20)
                .y_label_area_size(40)
                // Finally attach a coordinate on the drawing area and make a chart context
                .build_ranged((range.start() - 1)..(range.end() + 1), 0..y_max)?;
            chart
                .configure_mesh()
                // We can customize the maximum number of labels allowed for each axis
                //.x_labels(5)
                // .y_labels(4)
                .y_label_style(("sans-serif", 12).into_font())
                // We can also change the format of the label text
                .x_label_formatter(&|x| format!("{:.3}", x))
                .draw()?;
            let color = Palette99::pick(idx); // BLUE
            chart
                .draw_series(
                    Histogram::vertical(&chart)
                        .style(color.filled())
                        .data(values.iter().map(|t| *t)),
                )?
                .label(format!("{}", lambda(idx).unwrap_or(&idx.to_string())))
                .legend(move |(x, y)| {
                    Rectangle::new(
                        [(x - 5, y - 5), (x + 5, y + 5)],
                        Palette99::pick(idx).filled(),
                    )
                });
            chart
                .configure_series_labels()
                .background_style(&WHITE.mix(0.8))
                .border_style(&BLACK)
                .draw()?;
        }
    }
    Ok(())
}

pub fn bam_record_vis<'a, F>(
    matches: &ArgMatches,
    range: StringRegion,
    mut list: Vec<(u64, Record)>,
    mut annotation: Vec<(u64, bed::Record)>,
    mut frequency: BTreeMap<u64, Vec<(u64, u32)>>,
    lambda: F,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: Fn(usize) -> Option<&'a str>,
{
    // let diff = range.end() - range.start();
    let preset: Option<VisPreset> = matches.value_of_t("preset").ok(); // .unwrap_or_else(|e| e.exit());
    eprintln!("Preset: {:?}", preset);
    let no_margin = matches.is_present("no-scale");
    let _prefetch_range = matches.is_present("prefetch-range");
    let output = matches.value_of("output").unwrap();
    let no_cigar = matches.is_present("no-cigar");
    let udon = matches.is_present("udon");
    let packing = !matches.is_present("no-packing");
    let quality = matches.is_present("quality");
    let legend = !matches.is_present("no-legend");
    let insertion = !matches.is_present("no-insertion");
    let split = matches.is_present("split-alignment");
    let split_only = matches.is_present("only-split-alignment");
    let sort_by_name = matches.is_present("sort-by-name");
    let sort_by_cigar = matches.is_present("sort-by-cigar");
    let pileup = matches.is_present("pileup");
    let all_bases = matches.is_present("all-bases");
    let hide_alignment = matches.is_present("hide-alignment");
    let only_translocation = matches.is_present("only-translocation");
    let output_translocation = matches.is_present("output-translocation");
    let translocation_target = matches.value_of("translocation-target");
    let square = matches.is_present("square");
    if hide_alignment {
        return frequency_vis(matches, range, list, lambda);
    }
    let max_coverage = matches
        .value_of("max-coverage")
        .and_then(|a| a.parse::<u32>().ok());
    let x = matches
        .value_of("x")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(1280u32);
    let x_scale = matches
        .value_of("x-scale")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(20u32);
    let y = matches
        .value_of("y")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(15u32);
    let freq_size = matches
        .value_of("freq-height")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(50u32);
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
    if sort_by_name {
        list.sort_by(|a, b| {
            a.0.cmp(&b.0)
                .then(a.1.name().cmp(&b.1.name()))
                .then(a.1.start().cmp(&b.1.start()))
        });
    } else {
        list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
    }

    // Calculate coverage; it won't work on sort_by_name
    // let mut frequency = BTreeMap::new(); // Vec::with_capacity();
    if pileup {
        list.iter().group_by(|elt| elt.0).into_iter().for_each(|t| {
            let mut line = Vec::with_capacity((range.end - range.start + 1) as usize);
            for column in
                bam::Pileup::with_filter(&mut RecordIter(t.1), |record| record.flag().no_bits(1796))
            {
                let column = column.unwrap();
                /*println!("Column at {}:{}, {} records", column.ref_id(),
                column.ref_pos() + 1, column.entries().len());*/
                // Should we have sparse occurrence table?
                // eprintln!("{:?} {:?}",  range.path, lambda(column.ref_id() as usize).unwrap_or(&column.ref_id().to_string()));
                // lambda(column.ref_id() as usize).unwrap_or(&column.ref_id().to_string())
                // == range.path
                // &&
                if range.start <= column.ref_pos() as u64 && column.ref_pos() as u64 <= range.end {
                    line.push((column.ref_pos() as u64, column.entries().len() as u32));
                }
            }
            // eprintln!("{:?}", line);
            frequency.insert(t.0, line);
        });
    }

    annotation.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
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
    let axis_count = axis.len(); //into_iter().unique_by(|s| s.0).count();
                                 /*
                                 let annotation = {
                                     let mut btree = BTreeMap::new();
                                     for i in annotation {
                                         //let k = i?;
                                         //btree.insert(k[0].to_string(), (k[1].parse::<u64>()?, k[2].parse::<u64>()?));
                                         //eprintln!("{:?}", k);
                                         let item = btree.entry(i.0).or_insert(vec![]);
                                         item.push((i.1.name().clone(), i.1.chrom(), i.1.start(), i.1.end()));
                                         //vec.push((k[0].to_string(), k[1].parse::<u64>()?, k[2].parse::<u64>()?))
                                     }
                                     //btree
                                     btree
                                 };*/
    let annotation_count = annotation.iter().unique_by(|s| s.0).count(); // annotation.len();

    let y_len = 40
        + (prev_index as u32 + axis_count as u32 + annotation_count as u32 * 2) * y
        + frequency.len() as u32 * freq_size;
    let x_len = if square { y_len } else { x };

    let root = BitMapBackend::new(output, (x_len, y_len)).into_drawing_area();
    let approximate_one_pixel = 1; //((range.end() - range.start()) / x as u64) as u32;
    root.fill(&WHITE)?;
    let root = root.margin(0, 0, 0, 0);
    // After this point, we should be able to draw construct a chart context
    // let areas = root.split_by_breakpoints([], compressed_list);
    let top_margin = if no_margin { 0 } else { 40 };
    let (alignment, coverage) = root.split_vertically(
        top_margin + (prev_index as u32 + axis_count as u32 + annotation_count as u32 * 2) * y,
    );
    eprintln!("{:?} {:?} {:?}", prev_index, axis_count, annotation_count);
    let y_area_size = if no_margin { 0 } else { 40 };

    let mut chart = if no_margin {
        ChartBuilder::on(&alignment)
            // Set the caption of the chart
            // Set the size of the label region
            .x_label_area_size(x_scale)
            .y_label_area_size(y_area_size)
            // Finally attach a coordinate on the drawing area and make a chart context
            .build_ranged(
                (range.start() - 1)..(range.end() + 1),
                0..(1 + prev_index + axis_count + annotation_count * 2),
            )?
    } else {
        ChartBuilder::on(&alignment)
            // Set the caption of the chart
            .caption(format!("{}", range), ("sans-serif", 20).into_font())
            // Set the size of the label region
            .x_label_area_size(20)
            .y_label_area_size(y_area_size)
            // Finally attach a coordinate on the drawing area and make a chart context
            .build_ranged(
                (range.start() - 1)..(range.end() + 1),
                0..(1 + prev_index + axis_count + annotation_count * 2),
            )?
    };
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

    // Draw annotation if there is bed-compatible annotation.
    annotation
        .into_iter()
        .group_by(|elt| elt.0)
        .into_iter()
        .enumerate()
        .for_each(|(key, value)| {
            // let sample = value.0.clone();
            value
                .1
                .into_iter()
                .enumerate()
                .for_each(|(index, (_, record))| {
                    let end = record.end();
                    let start = record.start();
                    if end > range.start() && start < range.end() {
                        let start = if start > range.start() {
                            start
                        } else {
                            range.start()
                        };
                        let end = if end < range.end() { end } else { range.end() };
                        let stroke = Palette99::pick(index as usize);
                        let outer_stroke = match record.strand() {
                            Some(Strand::Forward) => POS_COL.mix(0.8),
                            Some(Strand::Reverse) => NEG_COL.mix(0.8),
                            _ => TRANSPARENT,
                        };
                        /*let mut bar2 = Rectangle::new(
                            [(start, prev_index + key), (end, prev_index + key + 1)],
                            stroke.stroke_width(2),
                        );
                        bar2.set_margin(1, 1, 0, 0);*/
                        // prev_index += 1;
                        chart
                            .draw_series(LineSeries::new(
                                vec![
                                    (start, prev_index + key * 2 + axis_count + 1),
                                    (end, prev_index + key * 2 + axis_count + 1),
                                ],
                                outer_stroke.stroke_width(y),
                            ))
                            .unwrap();
                        chart
                            .draw_series(LineSeries::new(
                                vec![
                                    (start, prev_index + key * 2 + axis_count + 1),
                                    (end, prev_index + key * 2 + axis_count + 1),
                                ],
                                stroke.stroke_width(y / 3 * 4),
                            ))
                            .unwrap()
                            .label(format!("{}", record.name().unwrap_or(&"")))
                            .legend(move |(x, y)| {
                                Rectangle::new([(x - 5, y - 5), (x + 5, y + 5)], stroke.filled())
                            });
                        //Some(bar2)
                    }
                })
        });

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
                                vec![(start, prev_index + key), (end, prev_index + key)],
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
    if compressed_list.len() > 1 {
        if legend || no_margin {
            //list2.sort_by(|a, b| a.0.cmp(&b.0));
            // eprintln!("{}", list.len());
            // let mut prev_index = 0;

            for (sample_sequential_id, sample) in compressed_list.iter()
            // list2.into_iter().group_by(|elt| elt.0).into_iter()
            {
                // Check that the sum of each group is +/- 4.
                // assert_eq!(4, group.iter().fold(0_i32, |a, b| a + b).abs());
                let count = *sample; //.count();
                if count > 0 {
                    let idx = *sample_sequential_id as usize;
                    // let idx = sample.next().0;
                    chart
                        .draw_series(LineSeries::new(
                            vec![(range.start() - 1, count), (range.end() + 1, count)],
                            Palette99::pick(idx).stroke_width(7),
                        ))?
                        .label(format!("{}", lambda(idx).unwrap_or(&idx.to_string())))
                        .legend(move |(x, y)| {
                            Rectangle::new(
                                [(x - 5, y - 5), (x + 5, y + 5)],
                                Palette99::pick(idx).filled(),
                            )
                        });
                }
                //prev_index += count;
            }
        } else {
            // For each sample:

            let mut prev_index = 0;
            chart.draw_series(
                compressed_list
                    .into_iter()
                    .map(|(sample_sequential_id, sample)| {
                        let count = sample;
                        if count > 0 {
                            let stroke = Palette99::pick(sample_sequential_id as usize);
                            let mut bar2 = Rectangle::new(
                                [
                                    (range.start() - 1, prev_index),
                                    (range.end() + 1, prev_index + count),
                                ],
                                stroke.stroke_width(5), // filled(), //stroke_width(100),
                            );
                            bar2.set_margin(1, 0, 0, 0);
                            prev_index += count;
                            Some(bar2)
                        } else {
                            None
                        }
                    })
                    .filter_map(|t| t),
            )?;
        }
    }

    // For each supplementary alignment:

    if sort_by_cigar {
        for (idx, i) in supplementary_list.iter().enumerate() {
            let stroke = Palette99::pick(idx as usize);
            // let stroke_color = BLACK;
            chart
                .draw_series(LineSeries::new(
                    vec![(i.3 as u64, i.1), (i.4 as u64, i.2)],
                    stroke.stroke_width(y / 2),
                ))?
                .label(format!("{}", String::from_utf8_lossy(i.0)))
                .legend(move |(x, y)| {
                    Rectangle::new([(x - 5, y - 5), (x + 5, y + 5)], stroke.filled())
                });
        }
    } else {
        chart.draw_series(supplementary_list.iter().map(|i| {
            let stroke = BLACK;
            let mut bar2 = Rectangle::new(
                [(i.3 as u64, i.1), (i.4 as u64, i.2)],
                stroke.stroke_width(2), // filled(), // (y / 4), // filled(), //stroke_width(100),
            );
            bar2.set_margin(y / 4, y / 4, 0, 0);
            bar2
        }))?;
    }
    let mut split_frequency = vec![];
    // For each alignment:
    let series = {
        //list.into_iter().enumerate().map(|(index, data)| {
        let mut bars = vec![];
        index_list
            .into_iter()
            .zip(list)
            .filter(|(_, data)| {
                (data.1.start() as u64) < range.end()
                    && (data.1.calculate_end() as u64) > range.start()
            })
            .for_each(|(index, data)| {
                //chart.draw_series(index_list.into_par_iter().zip(list).map(|(index, data)| {
                //for (index, data) in list.iter().enumerate() {
                let bam = data.1;
                let color = if bam.flag().is_reverse_strand() {
                    NEG_COL
                } else {
                    POS_COL
                }
                .mix(0.8);
                let _stroke = Palette99::pick(data.0 as usize); //.unwrap(); //if data.0 % 2 == 0 { CYAN } else { GREEN };
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
                bar.set_margin(2, 2, 0, 0);

                // eprintln!("{:?}", [(start, index), (end, index + 1)]);
                bars.push(bar);
                //let mut bars =  //, bar2];
                if split {
                    match bam.tags().get(b"SA") {
                        Some(TagValue::String(array_view, StringType::String)) => {
                            // assert!(array_view.int_type() == IntegerType::U32);
                            let current_left_clip =
                                bam.cigar().soft_clipping(!bam.flag().is_reverse_strand())
                                    + bam.cigar().hard_clipping(!bam.flag().is_reverse_strand());
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
                                    /*
                                    eprintln!(
                                        "{} {:?} {:?}",
                                        String::from_utf8_lossy(bam.name()),
                                        current_left_clip,
                                        cigar.soft_clipping(strand == "+"),
                                    );*/
                                    if only_translocation && t[0] == range.path {
                                        None
                                    } else if let Some(translocation_target) = translocation_target {
                                        if translocation_target == t[0] {
                                            Some(cigar.soft_clipping(strand == "+"))
                                        } else {
                                            None
                                        }
                                    } else {
                                        Some(cigar.soft_clipping(strand == "+"))
                                    }
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

                            let color = SPL_COL;
                            if ((is_smaller && !bam.flag().is_reverse_strand())
                                || (is_larger && bam.flag().is_reverse_strand()))
                                && bam.start() as u64 > range.start()
                            {
                                // split alignment on left
                                let mut bar = Rectangle::new(
                                    [(start, index), (start + 1, index + 1)],
                                    color.stroke_width(2), //.filled(),
                                );
                                bar.set_margin(0, 0, 0, 0);
                                bars.push(bar);
                                split_frequency.push((data.0, (start, approximate_one_pixel)));
                                if output_translocation {
                                    println!("S\t{}\t{}\tL\t{}", range.path, start, end);
                                }
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
                                bars.push(bar);
                                split_frequency.push((data.0, (end, approximate_one_pixel)));
                                if output_translocation {
                                    println!("S\t{}\t{}\tR\t{}", range.path, end, start);
                                }
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
                } /*
                  if legend {
                  } else {
                      let mut bar2 =
                          Rectangle::new([(start, index), (end, index + 1)], stroke.stroke_width(2));
                      bar2.set_margin(1, 1, 0, 0);
                      //vec![bar,bar2]
                      bars.push(bar2);
                  };*/
                if !no_cigar && udon {
                    let record = bam;
                    let left_top = chart.as_coord_spec().translate(&(range.start-1, index));
                    let right_bottom = chart.as_coord_spec().translate(&(range.end+1, index + 1));

                    let opt_len = (right_bottom.0 - left_top.0) as usize;

                    let window = Range::<usize> {
                        start: range.start() as usize - 1usize,
                        end:   range.end() as usize + 1usize
                    };

                    let pixels_per_column = window.len() as f64 / opt_len as f64;
                    let scaler = UdonScaler::new(&UdonPalette::default(), pixels_per_column);
                    let base_color: [[[u8; 4]; 2]; 2] = [
                        [[255, 202, 191, 255], [255, 255, 255, 255]],
                        [[191, 228, 255, 255], [255, 255, 255, 255]]
                    ];
                    let udon = Udon::build(
                        record.cigar().raw(),
                        record.sequence().raw(),
                        if let Some(TagValue::String(s, _)) = record.tags().get(b"MD") { s } else {
                                panic!("Each BAM record must have MD string. Inspect `samtools calmd` for restoring missing MD strings.")
                        }
                    ).expect("Failed to create udon index. Would be a bug.");

                    /* compose span, skip if out of the window */
                    let range = Range::<usize> {
                            start: record.start() as usize,
                            end:   record.start() as usize + udon.reference_span()
                    };
                    //if !window.has_overlap(&range) /* || range.len() < window.len() / 8 */ { continue; }

                    /* compute local ranges */
                    let udon_range   = range.clip(&window).unwrap();
                    let window_range = window.clip(&range).unwrap();
                    //if 3 * window_range.len() < window.len() { break; }
                    let (window_range, offset_in_pixels) = window_range.scale(pixels_per_column);
                    //eprintln!("{:?}, {:?}, {:?} {:?}", range, udon_range, offset_in_pixels, window_range);

                    /* slice ribbon scaled */
                    let mut ribbon = udon.decode_scaled(
                            &udon_range,
                            offset_in_pixels,
                            &scaler
                    ).expect("Failed to decode udon ribbon. Would be a bug.");
                    ribbon.append_on_basecolor(&base_color[record.flag().is_reverse_strand() as usize]).correct_gamma();
                    let horizontal_offset = window_range.start;
                    let left_blank  = horizontal_offset;
                    let right_blank = opt_len.saturating_sub(ribbon.len() + horizontal_offset);
                    let ribbon_len  = opt_len - (left_blank + right_blank);

                    for (i, &x) in ribbon[.. ribbon_len].into_iter().enumerate() {
                        let cv = &x[0][.. 3];
                        let color = RGBColor(cv[0], cv[1], cv[2]);
                        let prev_pixel_ref = if i == 0 {
                            window.start+ (offset_in_pixels * pixels_per_column) as usize +1+ ((horizontal_offset + i) as f64 * pixels_per_column) as usize // start as usize
                        } else {
                            window.start+1+ ((horizontal_offset + i) as f64 * pixels_per_column) as usize
                        };
                        let prev_ref = if i == ribbon_len-1 {
                            end as usize
                        } else {
                            window.start + ((horizontal_offset + i + 1) as f64 * pixels_per_column) as usize
                        };
                        let mut bar = Rectangle::new(
                            [
                                (prev_pixel_ref as u64, index),
                                (prev_ref as u64, index + 1),
                            ],
                            color.filled(),
                        );
                        bar.set_margin(3, 3, 0, 0);
                        bars.push(bar);
                    }
                }
                else if !no_cigar && ! udon {
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
                                                    color = Some(INS_COL);
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
                                                if prev_ref >= range.start() as u64 && all_bases {
                                                    let record_nt =
                                                        entry.record_pos_nt().unwrap().1;
                                                    color = match record_nt as char {
                                                        'A' => Some(A_COL), //RED,
                                                        'C' => Some(C_COL), // BLUE,
                                                        'G' => Some(G_COL),
                                                        'T' => Some(T_COL), //GREEN,
                                                        _ => Some(N_COL),
                                                    };
                                                }
                                            } else {
                                                /* Mismatch */
                                                prev_ref = entry.ref_pos_nt().unwrap().0 as u64;

                                                if prev_ref > range.end() as u64 {
                                                    break;
                                                }
                                                if prev_ref >= range.start() as u64 {
                                                    let record_nt =
                                                        entry.record_pos_nt().unwrap().1;
                                                    color = match record_nt as char {
                                                        'A' => Some(A_COL), //RED,
                                                        'C' => Some(C_COL), // BLUE,
                                                        'G' => Some(G_COL),
                                                        'T' => Some(T_COL), //GREEN,
                                                        _ => Some(N_COL),
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
                                            bar.set_margin(3, 3, 0, 0);
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
                                                    RGBColor(*qual * 5, *qual * 5, *qual * 5)
                                                        .filled(),
                                                );
                                                bar.set_margin(3, 3, 0, 0);
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
                                            bar.set_margin(3, 3, 0, 0);
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
                                                INS_COL.stroke_width(1),
                                            );
                                            // eprintln!("{:?}", [(prev_ref, index), (prev_ref + 1, index + 1)]);
                                            bar.set_margin(1, 1, 0, 5);
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

    if frequency.len() > 0 {
        let coverages = coverage.split_evenly((frequency.len(), 1));
        for (i, (sample_sequential_id, values)) in frequency.iter().rev().enumerate() {
            let idx = *sample_sequential_id as usize;
            let y_max = match max_coverage {
                Some(a) => a,
                _ => values.iter().map(|t| t.1).max().unwrap_or(1),
            };
            let mut chart = ChartBuilder::on(&coverages[i])
                // Set the caption of the chart
                //.caption(format!("{}", range), ("sans-serif", 20).into_font())
                // Set the size of the label region
                .x_label_area_size(0)
                .y_label_area_size(y_area_size)
                // Finally attach a coordinate on the drawing area and make a chart context
                .build_ranged((range.start() - 1)..(range.end() + 1), 0..y_max)?;
            chart
                .configure_mesh()
                // We can customize the maximum number of labels allowed for each axis
                //.x_labels(5)
                // .y_labels(4)
                .y_label_style(("sans-serif", 12).into_font())
                // We can also change the format of the label text
                // .x_label_formatter(&|x| format!("{:.3}", x))
                .draw()?;
            let color = Palette99::pick(idx); // BLUE
            chart
                .draw_series(
                    Histogram::vertical(&chart)
                        .style(color.filled())
                        .data(values.iter().map(|t| *t)),
                )?
                .label(format!("{}", lambda(idx).unwrap_or(&idx.to_string())))
                .legend(move |(x, y)| {
                    Rectangle::new(
                        [(x - 5, y - 5), (x + 5, y + 5)],
                        Palette99::pick(idx).filled(),
                    )
                });

            chart.draw_series(
                Histogram::vertical(&chart).style(SPL_COL.filled()).data(
                    split_frequency
                        .iter()
                        .filter(|&&t| t.0 == *sample_sequential_id)
                        .map(|t| t.1),
                ),
            )?;
            /*.label(format!("{}", lambda(idx).unwrap_or(&idx.to_string())))
            .legend(move |(x, y)| {
                Rectangle::new(
                    [(x - 5, y - 5), (x + 5, y + 5)],
                    Palette99::pick(idx).filled(),
                )
            });*/
            if !no_margin {
                chart
                    .configure_series_labels()
                    .background_style(&WHITE.mix(0.8))
                    .border_style(&BLACK)
                    .draw()?;
            }
        }
    }

    if legend {
        chart
            .configure_series_labels()
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .draw()?;
    }
    Ok(())
}
