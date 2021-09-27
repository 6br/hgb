use crate::dump::*;
use crate::{color::ColorSet, color::VisColor, VisOrig, VisPreset, VisRef};
use bam::record::{
    tags::{StringType, TagValue},
    Cigar,
};
use bam::{Record, RecordReader};
use bio_types::strand::Strand;
use clap::ArgMatches;
use itertools::Itertools;
use log::{debug, info};
use num_format::{Locale, ToFormattedString};
use plotters::coord::ReverseCoordTranslate;
use plotters::prelude::Palette;
use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};
use plotters::style::RGBColor;
// use std::cmp::min;
use std::ops::Range;
use std::{collections::BTreeMap, fs::File, time::Instant};
use std::{convert::TryInto, path::PathBuf};
use twobit::TwoBitFile;
use udon::{Udon, UdonPalette, UdonScaler, UdonUtils};

const PARBASE_THRESHOLD: u64 = 5;

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
        true
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

fn nt_color(record_nt: char, preset_color: &ColorSet) -> Option<RGBColor> {
    match record_nt {
        'A' => Some(preset_color.pick(VisColor::ACol)), // RED,
        'C' => Some(preset_color.pick(VisColor::CCol)), // BLUE,
        'G' => Some(preset_color.pick(VisColor::GCol)),
        'T' => Some(preset_color.pick(VisColor::TCol)), // GREEN,
        _ => Some(preset_color.pick(VisColor::NCol)),
    }
}

fn name_to_num(name: &[u8]) -> usize {
    let mut uuid = 0usize;
    // let right_most = if name.len() > 8 { 8 } else { name.len() };
    for &i in name.iter() {
        uuid += i as usize
    }
    uuid
}

fn switch_base(c: char) -> char {
    match c {
        'a' => 't',
        'c' => 'g',
        't' => 'a',
        'g' => 'c',
        'u' => 'a',
        'A' => 'T',
        'C' => 'G',
        'T' => 'A',
        'G' => 'C',
        'U' => 'A',
        _ => 'N',
    }
}
pub struct RecordIter<'a, I: Iterator<Item = &'a (u64, Record)>>(I);

impl<'a, I> RecordIter<'a, I>
where
    I: Iterator<Item = &'a (u64, Record)>,
{
    pub fn new(item: I) -> Self {
        RecordIter(item)
    }
}

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
    //range: &StringRegion,
    //frequency: &BTreeMap<u64, Vec<(u64, u32, char)>>,
    vis: Vec<VisRef>,
    lambda: F,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: Fn(usize) -> Option<&'a str>,
{
    let no_margin = matches.is_present("no-scale");
    let y_area_size = if no_margin { 0 } else { 40 };
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
    let x_as_range = matches.is_present("x-as-range");
    let dynamic_partition = matches.is_present("dynamic-partition");
    let freq_len = vis.iter().map(|a| a.frequency.len()).max().unwrap();
    let freq_len_ids = vis
        .iter()
        .map(|a| (a.frequency.len(), a.frequency.keys()))
        .max_by_key(|t| t.0)
        .unwrap()
        .1
        .collect::<Vec<_>>();
    let x_scale = matches
        .value_of("x-scale")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(20u32)
        / if freq_len > 1 { freq_len as u32 } else { 1u32 };
    let x_len = if x_as_range {
        vis.iter().map(|t| t.range.interval() as u32).sum::<u32>()
    } else {
        x
    };
    let y_max = match max_coverage {
        Some(a) => a,
        _ => vis
            .iter()
            .map(|t| {
                t.frequency
                    .iter()
                    .map(|(_, values)| values.iter().map(|t| t.1).max().unwrap_or(1))
                    .max()
                    .unwrap_or(1)
            })
            .max()
            .unwrap_or(1),
        /*if frequency.len() > 0 {

            let coverages = coverage.split_evenly((frequency.len(), 1));
            for (i, (sample_sequential_id, values)) in frequency.iter().rev().enumerate() {
                values.iter().map(|t| t.1).max().unwrap_or(1),
            }
        }*/
    };
    // list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
    // Calculate coverage; it won't work on sort_by_name

    let root = BitMapBackend::new(output, (x, freq_size)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(0, 0, 0, 0);

    let areas = if !dynamic_partition {
        root.split_evenly((1, vis.len()))
    } else {
        let x_axis_sum: u64 = vis.iter().map(|t| t.range.interval()).sum();
        let mut x_axis = vis
            .iter()
            .map(|t| (t.range.interval() * x_len as u64 / x_axis_sum) as u32)
            .scan(0, |acc, x| {
                *acc += x;
                Some(*acc)
            })
            .collect::<Vec<u32>>();
        x_axis.pop();
        let y_axis: Vec<u32> = vec![];
        root.split_by_breakpoints(x_axis, y_axis)
    };
    //let areas_len = vis.len();
    // After this point, we should be able to draw construct a chart context
    // let areas = root.split_by_breakpoints([], compressed_list);
    if freq_len > 0 {
        for (index, area) in areas.iter().enumerate() {
            let area = area.margin(0, 0, 0, 0);
            let vis = &vis[index];
            let range = &vis.range;
            let frequency = &vis.frequency;
            let coverages = area.split_evenly((frequency.len(), 1));
            for (i, &sample_sequential_id) in freq_len_ids.iter().rev().enumerate() {
                //}
                //for (i, (sample_sequential_id, values)) in frequency.iter().rev().enumerate() {
                let idx = *sample_sequential_id as usize;

                let x_spec = if no_margin && range.end() - range.start() <= 5000 {
                    range.start()..range.end()
                } else {
                    (range.start() - 1)..(range.end() + 1)
                };
                let x_offset = if range.end() - range.start() <= PARBASE_THRESHOLD {
                    coverages[i].get_pixel_range().0.len()
                        / ((range.end() - range.start()) * 2) as usize
                } else {
                    0usize
                };
                let x_label_formatter = {
                    &|x: &u64| {
                        if *x == range.start() || range.end() - range.start() > PARBASE_THRESHOLD {
                            x.to_formatted_string(&Locale::en)
                        } else {
                            format!("")
                        }
                    }
                };
                let mut chart = ChartBuilder::on(&coverages[i])
                    // Set the caption of the chart
                    //.caption(format!("{}", range), ("sans-serif", 20).into_font())
                    // Set the size of the label region
                    .x_label_area_size(x_scale)
                    .y_label_area_size(y_area_size)
                    // Finally attach a coordinate on the drawing area and make a chart context
                    .build_cartesian_2d(x_spec.clone(), 0..y_max)?;
                chart
                    .configure_mesh()
                    // We can customize the maximum number of labels allowed for each axis
                    .x_label_offset(x_offset as u32)
                    .x_label_style(("sans-serif", x_scale / 2).into_font())
                    // .y_labels(4)
                    .y_label_style(("sans-serif", 12).into_font())
                    // We can also change the format of the label text
                    .x_label_formatter(x_label_formatter)
                    // .x_label_formatter(&|x| format!("{}", x.to_formatted_string(&Locale::en)))
                    .draw()?;
                if let Some(values) = frequency.get(sample_sequential_id) {
                    let color = Palette99::pick(idx); // BLUE
                                                      /*eprintln!("{} {:?}", y_max, values
                                                      .iter()
                                                      .filter(|t| t.0 >= range.start() && t.0 < range.end())
                                                      .map(|t| *t));*/
                    chart
                        .draw_series(
                            Histogram::vertical(&chart)
                                .style(color.filled())
                                .margin(2)
                                .data(
                                    values
                                        .iter()
                                        .filter(|t| {
                                            t.0 >= range.start() && t.0 < range.end() && t.2 == '*'
                                        })
                                        .map(|t| (t.0, t.1)),
                                ),
                        )?
                        .label(lambda(idx).unwrap_or(&idx.to_string()).to_string())
                        .legend(move |(x, y)| {
                            Rectangle::new(
                                [(x - 5, y - 5), (x + 5, y + 5)],
                                Palette99::pick(idx).filled(),
                            )
                        });
                    if !no_margin {
                        chart
                            .configure_series_labels()
                            .background_style(&WHITE.mix(0.8))
                            .border_style(&BLACK)
                            .draw()?;
                    }
                }
            }
        }
    }
    Ok(())
}

pub fn bam_record_vis_orig<'a, F>(
    matches: &ArgMatches,
    vis: Vec<VisOrig>,
    lambda: F,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: Fn(usize) -> Option<&'a str>,
{
    bam_record_vis(
        matches,
        vis.iter()
            .map(|t| *Box::new(t.convert()))
            .collect::<Vec<_>>(),
        lambda,
    )?;
    Ok(())
}

pub fn bam_record_vis<'a, F>(
    matches: &ArgMatches,
    mut vis: Vec<VisRef>,
    lambda: F,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: Fn(usize) -> Option<&'a str>,
{
    let start = Instant::now();
    let preset: Option<VisPreset> = matches.value_of_t("preset").ok(); // .unwrap_or_else(|e| e.exit());
    let preset_color: ColorSet = matches
        .value_of_t("preset-color")
        .ok()
        .unwrap_or_else(ColorSet::new);
    eprintln!("Preset: {:?}", preset);
    let show_read_id = matches.is_present("show-read-id");
    let overlapping_annotation = matches.is_present("dump-json");
    let no_bold_line = matches.is_present("no-bold-line");
    let no_margin = matches.is_present("no-scale");
    let no_ruler = matches.is_present("no-ruler");
    let _prefetch_range = matches.is_present("prefetch-range");
    let output = matches.value_of("output").unwrap();
    let no_cigar = matches.is_present("no-cigar");
    let udon = matches.is_present("udon");
    let _packing = !matches.is_present("no-packing");
    let quality = matches.is_present("quality");
    let legend = !matches.is_present("no-legend");
    let insertion = !matches.is_present("no-insertion");
    let deletion = !matches.is_present("no-deletion");
    let split = matches.is_present("split-alignment");
    let _split_only = matches.is_present("only-split-alignment");
    let _sort_by_name = matches.is_present("sort-by-name");
    let sort_by_cigar = matches.is_present("sort-by-cigar");
    let colored_by_name = matches.is_present("colored-by-name");
    let pileup = matches.is_present("pileup");
    let all_bases = matches.is_present("all-bases");
    let hide_alignment = matches.is_present("hide-alignment");
    let only_translocation = matches.is_present("only-translocation");
    let end_split = matches.is_present("end-split-callets");
    let with_caption = matches.is_present("with-caption");
    let with_caption_val = matches.value_of("with-caption").unwrap_or("");
    let output_translocation = matches.is_present("output-translocation");
    let _translocation_target = matches.value_of("translocation-target");
    let square = matches.is_present("square");
    let read_index = matches.is_present("read-index");
    let x_as_range = matches.is_present("x-as-range");
    let dump_json = matches.is_present("dump-json");
    let insertion_string = matches.is_present("insertion-string");
    let dynamic_partition = matches.is_present("dynamic-partition");
    let colored_by_motif = matches.occurrences_of("colored-by-motif") != 0;
    let colored_by_motif_vec: Option<Vec<String>> = matches
        .value_of("colored-by-motif")
        .map(|t| t.split(':').map(|t| t.to_string()).collect());
    let colored_by_tag = matches.occurrences_of("colored-by-tag") != 0;
    let colored_by_tag_vec = matches.value_of("colored-by-tag");
    let twobit = matches.value_of("ref-column");
    // let metadata = vec![];

    if hide_alignment {
        //Multi-ranged frequency vis has not yet been supported.
        /*let vis = &vis[0];
        let range = &vis.range;
        let frequency = vis.frequency;
        return frequency_vis(matches, range, frequency, lambda);*/
        return frequency_vis(matches, vis, lambda);
    }
    let vis_index = matches
        .value_of("range-index")
        .and_then(|a| a.parse::<usize>().ok());
    let max_coverage = matches
        .value_of("max-coverage")
        .and_then(|a| a.parse::<u32>().ok());
    let snp_frequency = matches
        .value_of("snp-frequency")
        .and_then(|a| a.parse::<f64>().ok());
    let x = matches
        .value_of("x")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(1280u32);
    let x_scale = matches
        .value_of("x-scale")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(40u32);
    let y = matches
        .value_of("y")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(20u32);
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

    /*
    let axis = if let Some(graph) = graph {
        graph.into_iter().collect().unique_by(|s| s[0]).count() // group_by(|elt| elt[0]).count() as usize
    } else {
        0usize
    };*/
    let axis = if let Some(mut graph) = graph {
        let mut btree = BTreeMap::new();
        for i in graph.records() {
            let k = i?;
            //btree.insert(k[0].to_string(), (k[1].parse::<u64>()?, k[2].parse::<u64>()?));
            //eprintln!("{:?}", k);
            let item = btree.entry(k[0].to_string()).or_insert_with(Vec::new);
            item.push((k[1].parse::<u64>()?, k[2].parse::<u64>()?));
        }
        btree
    //        graph.into_iter().collect().unique_by(|s| s[0]).count() // group_by(|elt| elt[0]).count() as usize
    } else {
        BTreeMap::new()
    };
    let axis_count = axis.len();

    let end0 = start.elapsed();
    eprintln!("{}.{:03} sec.", end0.as_secs(), end0.subsec_millis());

    //let annotation_count = annotation.iter().unique_by(|s| s.0).count(); // annotation.len();
    let annotation_count = vis
        .iter()
        .map(|a| a.annotation)
        .collect::<Vec<_>>()
        .into_iter()
        .flatten()
        .unique_by(|s| s.0)
        .count(); // annotation.len();
    let prev_index = vis.iter().map(|a| a.prev_index).max().unwrap();
    let freq_len = if pileup {
        vis.iter().map(|a| a.frequency.len()).max().unwrap()
    } else {
        0
    };
    let freq_len_ids = vis
        .iter()
        .map(|a| (a.frequency.len(), a.frequency.keys()))
        .max_by_key(|t| t.0)
        .unwrap()
        .1
        .collect::<Vec<_>>();

    let top_margin = if no_margin { 0 } else { 40 };
    let vis_len = vis.len();
    let y_max = match max_coverage {
        Some(a) => a,
        _ => vis
            .iter()
            .map(|t| {
                t.frequency
                    .iter()
                    .map(|(_, values)| values.iter().map(|t| t.1).max().unwrap_or(1))
                    .max()
                    .unwrap_or(1)
            })
            .max()
            .unwrap_or(1),
        /*if frequency.len() > 0 {

            let coverages = coverage.split_evenly((frequency.len(), 1));
            for (i, (sample_sequential_id, values)) in frequency.iter().rev().enumerate() {
                values.iter().map(|t| t.1).max().unwrap_or(1),
            }
        }*/
    };
    debug!("y_max: {}", y_max);
    let mut n_x_labels = if !dynamic_partition {
        vec![]
    } else {
        //let x_axis_sum: u64 = vis.iter().map(|t| t.range.interval()).sum();
        let x_axis_max: u64 = vis.iter().map(|t| t.range.interval()).max().unwrap();
        let x_axis = vis
            .iter()
            .map(|t| (t.range.interval() * 10 / x_axis_max) as usize)
            .map(|t| if t >= 1 { t } else { 1 }) // At least 1
            //.scan(0, |acc, x| {
            //    *acc = *acc + x;
            //    Some(*acc)
            //})
            .collect::<Vec<usize>>();
        //x_axis.pop();
        //eprintln!("{:?}", x_axis);
        x_axis
    };
    //eprintln!("n_x_labels: {:?}", n_x_labels);
    if let Some(val) = vis_index {
        vis = vec![vis[val].clone()];
        if dynamic_partition {
            n_x_labels = vec![n_x_labels[val]];
        }
    }
    let y_len = top_margin
        + (prev_index as u32 + axis_count as u32 + annotation_count as u32 * 2) * y
        + freq_len as u32 * freq_size;
    let x_len = if square {
        y_len
    } else if x_as_range {
        vis.iter().map(|t| t.range.interval() as u32).sum::<u32>()
    } else {
        x
    };

    let root = BitMapBackend::new(output, (x_len, y_len)).into_drawing_area();
    let approximate_one_pixel = 1; //((range.end() - range.start()) / x as u64) as u32;
    root.fill(&WHITE)?;
    let root = root.margin(0, 0, 0, 0);
    let mut json = vec![];

    let areas = if !dynamic_partition {
        root.split_evenly((1, vis.len()))
    } else {
        let x_axis_sum: u64 = vis.iter().map(|t| t.range.interval()).sum();
        let mut x_axis = vis
            .iter()
            .map(|t| (t.range.interval() * x_len as u64 / x_axis_sum) as u32)
            .scan(0, |acc, x| {
                *acc += x;
                Some(*acc)
            })
            .collect::<Vec<u32>>();
        x_axis.pop();
        eprintln!("{:?}", x_axis);
        let y_axis: Vec<u32> = vec![];
        root.split_by_breakpoints(x_axis, y_axis)
    };
    //let areas_len = vis.len();

    for (index, area) in areas.iter().enumerate() {
        let area = area.margin(0, 0, 0, 0);
        let vis = &vis[index];
        let range = &vis.range;
        let frequency = &vis.frequency;
        let list = &vis.list;
        let annotation = &vis.annotation;
        let compressed_list = &vis.compressed_list;
        let index_list = &vis.index_list;
        let prev_index = vis.prev_index;
        let supplementary_list = &vis.supplementary_list;

        // After this point, we should be able to draw construct a chart context
        // let areas = root.split_by_breakpoints([], compressed_list);
        //let (alignment, coverage) = area.split_vertically(
        //    y_len - freq_len as u32 * freq_size, //top_margin + (prev_index as u32 + axis_count as u32 + annotation_count as u32 * 2) * y,
        //);
        let (coverage, alignment) = area.split_vertically(freq_len as u32 * freq_size);
        info!(
            "{:?} {:?} {:?} <{:?}>",
            prev_index,
            axis_count,
            annotation_count,
            top_margin + (prev_index as u32 + axis_count as u32 + annotation_count as u32 * 2) * y
        );
        let y_area_size = if no_margin {
            if dynamic_partition {
                1
            } else {
                0
            }
        } else {
            25
        };
        let x_area_size = if no_ruler { 0 } else { x_scale / 2 };
        let x_spec = if no_margin && range.end() - range.start() <= 5000 {
            range.start()..range.end()
        } else {
            (range.start() - 1)..(range.end() + 1)
        };
        let x_offset = if range.end() - range.start() <= PARBASE_THRESHOLD {
            area.get_pixel_range().0.len() / ((range.end() - range.start()) * 2) as usize
        } else {
            0usize
        };
        let x_label_formatter = {
            &|x: &u64| {
                if (*x == range.start() || range.end() - range.start() > PARBASE_THRESHOLD)
                    && !no_ruler
                {
                    x.to_formatted_string(&Locale::en)
                } else {
                    format!("")
                }
            }
        };
        let y_spec_max = if read_index {
            1
        } else {
            1 + prev_index
                + axis_count
                + annotation_count * 2
                + if twobit.is_some() { 2 } else { 0 }
        };
        let y_spec = y_spec_max..0;
        let top_x_area_size = if read_index { 0 } else { x_area_size };
        let mut chart = if no_margin {
            if with_caption {
                ChartBuilder::on(&alignment)
                    // Set the caption of the chart
                    // Set the size of the label region
                    .caption(
                        if !with_caption_val.is_empty() {
                            format!("{} {}", with_caption_val, range)
                        } else {
                            format!("{}", range)
                        },
                        ("sans-serif", x_scale / 2).into_font(),
                    )
                    .y_label_area_size(y_area_size)
                    .top_x_label_area_size(x_area_size)
                    .x_label_area_size(top_x_area_size)
                    // Finally attach a coordinate on the drawing area and make a chart context
                    .build_cartesian_2d(x_spec.clone(), y_spec)?
            } else {
                ChartBuilder::on(&alignment)
                    // Set the caption of the chart
                    // Set the size of the label region
                    .y_label_area_size(y_area_size)
                    .top_x_label_area_size(x_area_size)
                    .x_label_area_size(top_x_area_size)
                    // Finally attach a coordinate on the drawing area and make a chart context
                    .build_cartesian_2d(x_spec.clone(), y_spec)?
            }
        } else if with_caption {
            ChartBuilder::on(&alignment)
                // Set the caption of the chart
                .caption(
                    if !with_caption_val.is_empty() {
                        format!("{} {}", with_caption_val, range)
                    } else {
                        format!("{}", range)
                    },
                    ("sans-serif", x_scale / 2).into_font(),
                )
                // Set the size of the label region
                .top_x_label_area_size(x_area_size)
                .x_label_area_size(top_x_area_size)
                .y_label_area_size(y_area_size)
                // Finally attach a coordinate on the drawing area and make a chart context
                .build_cartesian_2d((range.start() - 1)..(range.end() + 1), y_spec)?
        } else {
            ChartBuilder::on(&alignment)
                // Set the size of the label region
                .top_x_label_area_size(x_area_size)
                .x_label_area_size(top_x_area_size)
                .y_label_area_size(y_area_size)
                // Finally attach a coordinate on the drawing area and make a chart context
                .build_cartesian_2d((range.start() - 1)..(range.end() + 1), y_spec)?
        };
        let x_labels = n_x_labels.get(index).unwrap_or(&10); //if dynamic_partition {} else {10};
                                                             // Then we can draw a mesh
        if no_bold_line {
            chart
                .configure_mesh()
                // We can customize the maximum number of labels allowed for each axis
                .x_labels(*x_labels) // Default value is 10
                .disable_x_axis()
                .bold_line_style(RGBColor(0, 0, 0).mix(0.1))
                .light_line_style(RGBColor(0, 0, 0).mix(0.1))
                //.x_labels(1)
                .x_label_offset(x_offset as u32)
                .y_labels(1)
                .x_label_style(("sans-serif", x_scale / 4).into_font())
                // We can also change the format of the label text
                .x_label_formatter(x_label_formatter)
                .draw()?
        } else if no_margin {
            chart
                .configure_mesh()
                // We can customize the maximum number of labels allowed for each axis
                .x_labels(*x_labels) // Default value is 10
                .disable_x_axis()
                //.x_labels(1)
                .x_label_offset(x_offset as u32)
                .y_labels(1)
                .x_label_style(("sans-serif", x_scale / 4).into_font())
                // We can also change the format of the label text
                .x_label_formatter(x_label_formatter)
                .draw()?;
        } else {
            chart
                .configure_mesh()
                // We can customize the maximum number of labels allowed for each axis
                .x_labels(*x_labels) // Default value is 10
                //.x_labels(1)
                .x_label_offset(x_offset as u32)
                .y_labels(1)
                .x_label_style(("sans-serif", x_scale / 4).into_font())
                // We can also change the format of the label text
                .x_label_formatter(x_label_formatter)
                .draw()?;
        }

        let mut node_id_dict: BTreeMap<u64, (u64, u64)> = BTreeMap::new();
        let mut prev_pos = std::u64::MAX;
        let mut prev_node_id = 0u64;

        // We assume that axis[&range.path] includes all node id to be serialized.
        if let Some(axis) = axis.get(&range.path) {
            axis.iter().for_each(|(node_id, pos)| {
                if prev_pos > *pos {
                } else {
                    node_id_dict.insert(prev_node_id, (prev_pos, *pos));
                }
                prev_pos = *pos;
                prev_node_id = *node_id;
            });
            node_id_dict.insert(prev_node_id, (prev_pos, range.end()));
        }

        if let Some(twobit) = twobit {
            let left_top = chart.as_coord_spec().translate(&(range.start, index)); // range.start - 1 is better?
            let right_bottom = chart.as_coord_spec().translate(&(range.end + 1, index + 1));
            eprintln!("{:?}, {:?}", left_top, right_bottom);
            let opt_len = (right_bottom.0 - left_top.0) as usize;

            let window = Range::<usize> {
                start: range.start() as usize - 1usize,
                end: range.end() as usize + 1usize,
            };

            let pixels_per_column = opt_len as f64 / window.len() as f64;

            let enable_softmask = false;
            let tb_nosoft = TwoBitFile::open(twobit, enable_softmask).unwrap();

            if let Ok(seq) =
                tb_nosoft.sequence(&range.path, range.start as usize, range.end as usize)
            {
                let mut pos = range.start;
                let mut bars = vec![];
                let mut texts = vec![];
                for base in seq.chars() {
                    let color = nt_color(base, &preset_color)
                        .unwrap_or_else(|| preset_color.pick(VisColor::NCol));
                    let mut bar = Rectangle::new(
                        [(pos, y_spec_max - 1), (pos + 1, y_spec_max)],
                        color.filled(),
                    );
                    bar.set_margin(0, 0, 0, 0);
                    bars.push(bar);
                    if pixels_per_column as u32 >= y / 2 {
                        let text = Text::new(
                            format!("{}", base),
                            (pos, y_spec_max),
                            ("sans-serif", y / 4 * 3),
                        );
                        texts.push(text);
                    }
                    pos += 1;
                }
                chart.draw_series(bars)?;
                chart.draw_series(texts)?;
            }
        }
        let mut annotations = vec![];

        // Draw annotation if there is bed-compatible annotation.
        annotation
            .iter()
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
                                Some(Strand::Forward) => {
                                    preset_color.pick(VisColor::PosCol).mix(0.8)
                                }
                                Some(Strand::Reverse) => {
                                    preset_color.pick(VisColor::NegCol).mix(0.8)
                                }
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
                                    stroke.stroke_width(y / 2),
                                ))
                                .unwrap()
                                .label(record.name().unwrap_or(&"").to_string())
                                .legend(move |(x, y)| {
                                    Rectangle::new(
                                        [(x - 5, y - 5), (x + 5, y + 5)],
                                        stroke.filled(),
                                    )
                                });
                            if show_read_id {
                                let pos = Pos::new(HPos::Left, VPos::Bottom);
                                let style =
                                    TextStyle::from(("sans-serif", y / 3 * 2).into_font()).pos(pos);
                                let text = Text::new(
                                    record.name().unwrap_or(&"").to_string(),
                                    (start, prev_index + key * 2 + axis_count + 1),
                                    style,
                                );
                                chart.draw_series(vec![text]).unwrap();
                            }
                            if overlapping_annotation {
                                // println!("{}: {}", range, record.name().unwrap_or(&""));
                                let (lt, lb) = chart
                                    .as_coord_spec()
                                    .translate(&(start, prev_index + key * 2 + axis_count + 1)); // range.start - 1 is better?
                                let (rt, rb) = chart
                                    .as_coord_spec()
                                    .translate(&(end + 1, prev_index + key * 2 + axis_count + 1));
                                let annotation = Annotation {
                                    rectangle: (lt, lb, rt, rb),
                                    name: record.name().unwrap_or(&"").to_string(),
                                    start,
                                    end,
                                };
                                annotations.push(annotation)
                            }
                        }
                    })
            });

        // Draw graph axis if there is graph information.
        axis.iter()
            //.group_by(|elt| elt.0.as_str())
            //.into_iter()
            .enumerate()
            .for_each(|(key, value)| {
                let path_name = value.0.clone();
                value.1.iter().for_each(|(node_id, _pos)| {
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
                            let stroke = Palette99::pick(*node_id as usize);
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
                                    Rectangle::new(
                                        [(x - 5, y - 5), (x + 5, y + 5)],
                                        stroke.filled(),
                                    )
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
                for (sample_sequential_id, sample) in compressed_list.iter()
                // list2.into_iter().group_by(|elt| elt.0).into_iter()
                {
                    // Check that the sum of each group is +/- 4.
                    // assert_eq!(4, group.iter().fold(0_i32, |a, b| a + b).abs());
                    let count = *sample; //.count();
                    if count > 0 {
                        let idx = *sample_sequential_id as usize;
                        // let idx = sample.next().0;
                        let chart = chart.draw_series(LineSeries::new(
                            vec![(range.start(), count), (range.end(), count)],
                            Palette99::pick(idx).stroke_width(y / 4 * 3),
                        ))?;
                        if !no_margin {
                            chart
                                .label(lambda(idx).unwrap_or(&idx.to_string()).to_string())
                                .legend(move |(x, y)| {
                                    Rectangle::new(
                                        [(x - 5, y - 5), (x + 5, y + 5)],
                                        Palette99::pick(idx).filled(),
                                    )
                                });
                        }
                    }
                    //prev_index += count;
                }
            } else {
                // For each sample:
                let mut prev_index = 0;
                chart.draw_series(
                    compressed_list
                        .iter()
                        .map(|(sample_sequential_id, sample)| {
                            let count = *sample;
                            if count > 0 {
                                let stroke = Palette99::pick(*sample_sequential_id as usize);
                                let mut bar2 = Rectangle::new(
                                    [(range.start(), prev_index), (range.end(), count)],
                                    stroke.stroke_width(y / 2), // filled(), //stroke_width(100),
                                );
                                bar2.set_margin(1, 0, 0, 0);
                                prev_index = count;
                                Some(bar2)
                            } else {
                                None
                            }
                        })
                        .flatten(),
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
                    .label(format!("{}", String::from_utf8_lossy(&i.0)))
                    .legend(move |(x, y)| {
                        Rectangle::new([(x - 5, y - 5), (x + 5, y + 5)], stroke.filled())
                    });
            }
        } else {
            /*if no_margin {
                for (idx, i) in supplementary_list.iter().enumerate() {
                    let stroke = Palette99::pick(idx as usize);
                    // let stroke_color = BLACK;
                    chart
                        .draw_series(LineSeries::new(
                            vec![(range.start(), i.1), (range.end() as u64, i.1)],
                            stroke.stroke_width(2),
                        ))?;
                    chart
                    .draw_series(LineSeries::new(
                        vec![(i.3 as u64, i.2+1), (i.4 as u64, i.2+1)],
                        stroke.stroke_width(2),
                    ))?;
                }
            } else {*/
            if vis_len <= 1 {
                chart.draw_series(supplementary_list.iter().filter(|t| t.3 < t.4).map(|i| {
                    let stroke = BLACK;
                    let mut bar2 = Rectangle::new(
                        [(i.3 as u64, i.1), (i.4 as u64, i.2)],
                        stroke.stroke_width(1), // filled(), // (y / 4), // filled(), //stroke_width(100),
                    );
                    bar2.set_margin(y / 4, y / 4, 0, 0);
                    bar2
                }))?;
            }
            //}
        }
        let mut split_frequency = vec![];
        // let mut allele_frequency = vec![];
        // For each alignment:
        let mut reads = vec![];
        let (bars, texts) = {
            //list.into_iter().enumerate().map(|(index, data)| {
            let mut bars = vec![];
            let mut texts = vec![];
            index_list
                .iter()
                .zip(list.iter())
                .filter(|(index, data)| {
                    (data.1.start() as u64) < range.end()
                        && (data.1.calculate_end() as u64) > range.start() && **index < std::u32::MAX as usize
                })
                .for_each(|(&index, data)| {
                    //chart.draw_series(index_list.into_par_iter().zip(list).map(|(index, data)| {
                    //for (index, data) in list.iter().enumerate() {
                    let bam = &data.1;
                    let color = if colored_by_tag {
                        if let Some(colored_by_str) = colored_by_tag_vec {
                            if colored_by_str.is_empty() {
                                if bam.flag().is_reverse_strand() {
                                    preset_color.pick(VisColor::NegCol).mix(0.8)
                                } else {
                                    preset_color.pick(VisColor::PosCol).mix(0.8)
                                }
                            } else {
                                let tag: &[u8;2] = colored_by_str.as_bytes().try_into().expect("colored by tag with unexpected length: tag name must be two characters.");
                                if let Some(TagValue::Int(tag_id,_)) = bam.tags().get(tag) {
                                    Palette9999::pick(tag_id as usize).mix(0.4)
                                } else {
                                    preset_color.pick(VisColor::DefCol).mix(0.8)
                                }
                            }
                        } else {
                            preset_color.pick(VisColor::DefCol).mix(0.8)
                        }
                    } else if bam.flag().is_reverse_strand() {
                        preset_color.pick(VisColor::NegCol).mix(0.8)
                    } else {
                        preset_color.pick(VisColor::PosCol).mix(0.8)
                    };
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
                    let mut insertions = vec![];

                    /*chart
                    .draw_series(LineSeries::new(vec![(start, index), (end, index)], &color))?;*/
                    /*if let Some(val) = read_index {
                        if val == index {
                            let mut bar =
                            Rectangle::new([(start, 0), (end, 1)], color.filled());
                            bar.set_margin(2, 2, 0, 0);
                            bars.push(bar);
                        } else {
                            continue
                        }
                    } else {*/
                        let mut bar =
                            Rectangle::new([(start, index), (end, index + 1)], color.filled());
                        bar.set_margin(2, 2, 0, 0);

                        bars.push(bar);
                    //}
                    if show_read_id {
                        let pos = Pos::new(HPos::Left, VPos::Bottom);
                        let style = TextStyle::from(("sans-serif", y / 3 * 2).into_font()).pos(pos);
                        let text = Text::new(
                            format!("{}", String::from_utf8_lossy(bam.name())),
                            (start, index),
                            style,
                        );
                        texts.push(text);
                    }

                    if colored_by_name {
                        let color = Palette99::pick(name_to_num(data.1.name())).mix(0.8);
                        let mut inner_bar =
                            Rectangle::new([(start, index), (end, index + 1)], color.filled());
                        inner_bar.set_margin(3, 3, 0, 0);
                        bars.push(inner_bar);
                    }/* else if colored_by_tag {
                        if let Some(colored_by_str) = colored_by_tag_vec {
                            let tag: &[u8;2] = colored_by_str.as_bytes().try_into().expect("colored by tag with unexpected length: tag name must be two characters.");
                           
                            if let Some(TagValue::Int(tag_id,_)) = bam.tags().get(tag) {
                                //eprintln!("{:?}", tag_id);
                                let color = Palette99::pick(tag_id as usize).mix(0.3);
                                let mut inner_bar =
                                    Rectangle::new([(start, index), (end, index + 1)], color.filled());
                                inner_bar.set_margin(3, 3, 0, 0);
                                bars.push(inner_bar);
                            } else {
                                let color = DEF_COL.mix(0.3);
                                let mut inner_bar =
                                    Rectangle::new([(start, index), (end, index + 1)], color.filled());
                                inner_bar.set_margin(3, 3, 0, 0);
                                bars.push(inner_bar);
                            }
                        }
                    }*/
                    // eprintln!("{:?}", [(start, index), (end, index + 1)]);

                    //let mut bars =  //, bar2];
                    if split || end_split {
                        match bam.tags().get(b"SA") {
                            Some(TagValue::String(array_view, StringType::String)) => {
                                // assert!(array_view.int_type() == IntegerType::U32);
                                let current_left_clip = bam
                                    .cigar()
                                    .soft_clipping(!bam.flag().is_reverse_strand())
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

                                        if only_translocation && t[0] == range.path {
                                            None
                                        } else {
                                            Some(cigar.soft_clipping(strand == "+"))
                                        }
                                    })
                                    .flatten()
                                    .collect();
                                let is_smaller = sa_left_clip
                                    .iter()
                                    .any(|t| t < &current_left_clip);
                                let is_larger = sa_left_clip
                                    .iter()
                                    .any(|t| t > &current_left_clip);

                                let color = preset_color.pick(VisColor::SplCol);
                                if ((is_smaller && !bam.flag().is_reverse_strand())
                                    || (is_larger && bam.flag().is_reverse_strand()))
                                    && bam.start() as u64 > range.start()
                                {
                                    // split alignment on left
                                    let mut bar = Rectangle::new(
                                        [(start, index), (start, index + 1)],
                                        color.stroke_width(2), //.filled(),
                                    );
                                    bar.set_margin(0, 0, 0, 0);
                                    bars.push(bar);
                                    split_frequency.push((data.0, (start, approximate_one_pixel)));
                                    if output_translocation {
                                        println!("S\t{}\t{}\tR\t{}", range.path, end, start);
                                    }
                                    if dump_json {
                                        let (lt, _) = chart.as_coord_spec().translate(&(start, index));
                                        insertions.push((lt, start, "Split-alignment".to_string()));
                                    }
                                }
                                if ((is_larger && !bam.flag().is_reverse_strand())
                                    || (is_smaller && bam.flag().is_reverse_strand()))
                                    && (bam.calculate_end() as u64) < range.end()
                                {
                                    // split alignment on right
                                    let mut bar = Rectangle::new(
                                        [(end, index), (end, index + 1)],
                                        color.stroke_width(2), //.filled(),
                                    );
                                    bar.set_margin(0, 0, 0, 0);
                                    bars.push(bar);
                                    split_frequency.push((data.0, (end, approximate_one_pixel)));
                                    if output_translocation {
                                        println!("S\t{}\t{}\tR\t{}", range.path, end, start);
                                    }
                                    if dump_json {
                                        let (lt, _) = chart.as_coord_spec().translate(&(end, index));
                                        insertions.push((lt, end, "Split-alignment".to_string()));
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
                    let left_top = chart.as_coord_spec().translate(&(range.start-1, index)); // range.start - 1 is better?
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
                    ).unwrap_or_else(|| panic!("Failed to decode udon ribbon. Would be a bug. Read id: {}, start: {}", String::from_utf8_lossy(record.name()), record.start()));

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

                    /* slice ribbon scaled */
                    let mut ribbon = udon.decode_scaled(
                            &udon_range,
                            offset_in_pixels,
                            &scaler
                    ).unwrap_or_else(|| panic!("Failed to decode udon ribbon. Would be a bug. Read id: {}, start: {}", String::from_utf8_lossy(record.name()), record.start()));
                    ribbon.append_on_basecolor(&base_color[record.flag().is_reverse_strand() as usize]).correct_gamma();
                    let horizontal_offset = window_range.start;
                    let left_blank  = horizontal_offset;
                    let right_blank = opt_len.saturating_sub(ribbon.len() + horizontal_offset);
                    let ribbon_len  = opt_len - (left_blank + right_blank);

                    for (i, &x) in ribbon[.. ribbon_len].iter().enumerate() {
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
                else if !no_cigar && !udon {
                    let mut prev_ref = bam.start() as u64 - 1;
                    let mut prev_nt = ' ';
                    let mut prev_pixel_ref = start-1;
                    let left_top = chart.as_coord_spec().translate(&(start, index));
                    let right_bottom = chart.as_coord_spec().translate(&(end, index + 1));
                    let offset = 1; /* if range.end() - range.start() == 1 {
                        1 // 0
                    } else {
                        1
                    };*/

                        //if let Ok(mut a) = bam.alignment_entries() {
                        match (quality, bam.alignment_entries()) {
                            (false, Ok(mut a)) => {
                                for i in left_top.0 + 1..right_bottom.0 + offset {
                                    let k =
                                        chart.as_coord_spec().reverse_translate((i, left_top.1));
                                    //eprintln!("{:?} {:?}", i,k);

                                    let mut color = None;
                                    let mut insertion_flag = false;
                                    let mut insertion_str = vec![];
                                    if let Some(k) = k {
                                        while k.0 > prev_ref
                                            && prev_ref != bam.calculate_end() as u64
                                        {
                                            /*if prev_min < entry.ref_pos_nt().unwrap().0 as u64 {
                                                prev_min = entry.ref_pos_nt().unwrap().0 as u64;
                                            }*/
                                            let entry = a.next();
                                            if let Some(entry) = entry {
                                                /*eprintln!(
                                                    "{:?} {:?} {:?} {:?}",
                                                    entry.ref_pos_nt(),
                                                    k.0,
                                                    prev_ref,
                                                    i
                                                );*/
                                                if entry.is_insertion() {
                                                    if prev_ref >= range.start() as u64 && insertion
                                                    {
                                                        //let mut bar = Rectangle::new([(prev_ref, index), (prev_ref+1, index + 1)], MAGENTA.stroke_width(1));
                                                        //eprintln!("{:?}", [(prev_ref, index), (prev_ref + 1, index + 1)]);
                                                        //bar.set_margin(0, 0, 0, 3);
                                                        //bars.push(bar);
                                                        //prev_ref = 0;
                                                        // color = Some(preset_color.pick(VisColor::INS_COL));
                                                        insertion_flag = true;
                                                        insertion_str.push((prev_ref, entry.record_nt().unwrap() as char));
                                                    }
                                                } else if entry.is_deletion() {
                                                    prev_ref = entry.ref_pos_nt().unwrap().0 as u64;
                                                    if prev_ref > range.end() as u64 {
                                                        break;
                                                    }
                                                    if prev_ref >= range.start() as u64 && deletion {
                                                        //let mut bar = Rectangle::new([(prev_ref , index), (prev_ref + 1, index + 1)], WHITE.filled());
                                                        //bar.set_margin(2, 2, 0, 0);
                                                        //eprintln!("White: {:?}", [(prev_pixel_ref, index), (prev_ref, index + 1)]);
                                                        //bars.push(bar);
                                                        color = Some(WHITE);
                                                    }

                                                } else if entry.is_seq_match() {
                                                    prev_ref = entry.ref_pos_nt().unwrap().0 as u64;
                                                    if prev_ref > range.end() as u64 {
                                                        break;
                                                    }
                                                    // If all bases shows the SNPs
                                                    if prev_ref >= range.start() as u64 && all_bases
                                                    {
                                                        let record_nt =
                                                            entry.record_pos_nt().unwrap().1;
                                                        color = nt_color(record_nt as char, &preset_color);
                                                    }
                                                } else {
                                                    /* Mismatch */
                                                    prev_ref = entry.ref_pos_nt().unwrap().0 as u64;

                                                    if prev_ref > range.end() as u64 {
                                                        break;
                                                    }
                                                    if prev_ref >= range.start() as u64
                                                        && !colored_by_motif
                                                    {
                                                        let record_nt =
                                                            entry.record_pos_nt().unwrap().1;
                                                        color = nt_color(record_nt as char, &preset_color);
                                                        //eprintln!("Mismatch: {:?}", [(prev_pixel_ref, index), (prev_ref, index + 1)]);
                                                        /*if let Some(_f) = snp_frequency {
                                                            allele_frequency.push((data.0, (record_nt, start, approximate_one_pixel)));
                                                        }*/

                                                        /*let mut bar =
                                                        Rectangle::new([(prev_ref as u64, index), (prev_ref as u64 + 1, index + 1)], color.filled());
                                                        bar.set_margin(2, 2, 0, 0);
                                                        //eprintln!("{:?}", [(prev_ref, index), (prev_ref + 1, index + 1)]);
                                                        bars.push(bar);*/
                                                    }
                                                }
                                                if !entry.is_deletion()
                                                    && colored_by_motif && prev_ref >= range.start() as u64
                                                {
                                                    if let Some(colored_by_motif_vec) =
                                                        &colored_by_motif_vec
                                                    {
                                                        let record_nt =
                                                            entry.record_pos_nt().unwrap().1;
                                                        let ref_nt = entry
                                                            .ref_nt().map(|t| t as char)
                                                            .unwrap_or(' ');
                                                        let next_ref_nt = a.clone()
                                                            .next()
                                                            .and_then(|t| t.ref_nt()).map(|t| t as char)
                                                            .unwrap_or(' ');
                                                        let string = format!(
                                                            "{}{}",
                                                            ref_nt, next_ref_nt
                                                        );
                                                        let revcomp_string = format!("{}{}", prev_nt, ref_nt);
                                                        debug!(
                                                            "Read: {} {} {} {} {}",
                                                            prev_ref, string, revcomp_string, bam.flag().is_reverse_strand(), record_nt as char
                                                        );
                                                        if bam.flag().is_reverse_strand() && revcomp_string == colored_by_motif_vec[2] {
                                                            /*let next_record_nt = a
                                                                .next()
                                                                .and_then(|t| t.record_pos_nt())
                                                                .and_then(|t| Some(t.1 as char))
                                                                .unwrap_or(' ');*/
                                                            eprintln!(
                                                                "REV: {} {} {} {}",
                                                                string, revcomp_string, colored_by_motif_vec[2], record_nt as char
                                                            );
                                                            if colored_by_motif_vec[0].starts_with(switch_base(record_nt as char))
                                                            {
                                                                color = Some(RED);
                                                            } else if colored_by_motif_vec[1].starts_with(switch_base(record_nt as char))
                                                            {
                                                                color = Some(BLUE);
                                                            }
                                                        } else if !bam.flag().is_reverse_strand() && string == colored_by_motif_vec[2] {
                                                            eprintln!(
                                                                "FWD: {} {} {}",
                                                                string, colored_by_motif_vec[2], record_nt as char
                                                            );
                                                            if colored_by_motif_vec[0].starts_with(record_nt as char)
                                                            {
                                                                color = Some(RED);
                                                            } else if colored_by_motif_vec[1].starts_with(record_nt as char)
                                                            {
                                                                color = Some(BLUE);
                                                            }
                                                        }
                                                    }
                                                }
                                                if let Some(ref_pos_nt) = entry.ref_pos_nt() {
                                                    prev_nt = ref_pos_nt.1 as char;
                                                }
                                            } else {
                                                break;
                                            }
                                        }
                                        if prev_ref >= range.start() as u64 {
                                            if let Some(color) = color {
                                                    let mut bar = Rectangle::new(
                                                        [
                                                            (prev_pixel_ref as u64+1, index),
                                                            (prev_ref as u64+1, index + 1),
                                                        ],
                                                        color.filled(),
                                                    );
                                                    bar.set_margin(3, 3, 0, 0);
                                                    bars.push(bar);
                                                /*eprintln!(
                                                "{:?}",
                                                [
                                                    (prev_pixel_ref, index),
                                                    (prev_ref, index + 1)
                                                ]);*/
                                            }
                                            if insertion_flag {
                                                // If insertion
                                                let mut bar = Rectangle::new(
                                                    [
                                                        (prev_pixel_ref as u64+1, index),
                                                        (prev_pixel_ref as u64+1, index + 1),
                                                    ],
                                                    preset_color.pick(VisColor::InsCol).stroke_width(1),
                                                );
                                                bar.set_margin(1, 1, 0, 0);
                                                bars.push(bar);
                                                if insertion_string {
                                                    //Note that 
                                                    let text = Text::new(
                                                        insertion_str.iter().map(|t| t.1).join("").to_string(),
                                                        (prev_pixel_ref + 1, index),
                                                        ("sans-serif", y / 4 * 3),
                                                    );
                                                    texts.push(text);
                                                }
                                                if dump_json {
                                                    insertion_str.iter().group_by(|elt| elt.0).into_iter().for_each(|(ge0, group)| {
                                                        let (lt, _) = chart.as_coord_spec().translate(&(ge0+1, index));
                                                        insertions.push((lt, ge0, group.map(|t| t.1).join("").to_string()));
                                                    }
                                                    );
                                                    // push(lt, prev_ref, ins);
                                                }
                                            }
                                            prev_pixel_ref = k.0;
                                        }
                                    }
                                }
                            }
                            _ => {
                                let mut insertion_str = vec![];
                                for entry in bam.aligned_pairs() {
                                    match entry {
                                        // (Seq_idx, ref_idx)
                                        (Some(_record), Some(reference)) => {
                                            if prev_ref > range.start() as u64 && quality {
                                                if let Some(qual) =
                                                    bam.qualities().raw().get(_record as usize)
                                                {
                                                    let qual_color = if *qual < 51 {*qual * 5 } else { 255 };

                                                    let mut bar = Rectangle::new(
                                                        [
                                                            (reference as u64, index),
                                                            (reference as u64 + 1, index + 1),
                                                        ],
                                                        RGBColor(qual_color, qual_color, qual_color)
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
                                        (Some(record), None) => {
                                            //Insertion
                                            if prev_ref > range.start() as u64 && insertion {
                                                let mut bar = Rectangle::new(
                                                    [(prev_ref, index), (prev_ref + 1, index + 1)],
                                                    preset_color.pick(VisColor::InsCol).stroke_width(1),
                                                );
                                                // eprintln!("{:?}", [(prev_ref, index), (prev_ref + 1, index + 1)]);
                                                bar.set_margin(1, 1, 0, 5);
                                                bars.push(bar);
                                                prev_ref = 0;
                                                //insertion_str.push((prev_ref, entry.record_nt().unwrap() as char));
                                                if insertion_string {
                                                    //Note that 
                                                    let text = Text::new(
                                                        (bam.sequence().at(record as usize) as char).to_string(),
                                                        (prev_ref, index),
                                                        ("sans-serif", y / 4 * 3),
                                                    );
                                                    texts.push(text);
                                                    insertion_str.push((prev_ref, bam.sequence().at(record as usize) as char));
                                                }
                                            }
                                            // eprintln!("{}", prev_ref)
                                        }
                                        _ => {}
                                    }
                                }
                                if dump_json {
                                    insertion_str.iter().group_by(|elt| elt.0).into_iter().for_each(|(ge0, group)| {
                                        let (lt, _) = chart.as_coord_spec().translate(&(ge0, index));
                                        insertions.push((lt, ge0, group.map(|t| t.1).join("").to_string()));
                                    }
                                    );
                                    // push(lt, prev_ref, ins);
                                }
                            }
                        }
                    }
                    if dump_json {
                        //println!("{}", String::from_utf8_lossy(bam.name()));
                        let (lt, lb) = chart.as_coord_spec().translate(&(range.start, index));
                        let (rt, rb) = chart.as_coord_spec().translate(&(range.end+1, index + 1));
                        let sastr = match bam.tags().get(b"SA") {
                            Some(TagValue::String(array_view, StringType::String)) => {
                                 String::from_utf8_lossy(array_view).to_string()
                            }
                        _ => "".to_string()
                        };
                        let mut readable: Vec<u8> = Vec::new();
                        bam.cigar().write_readable(&mut readable).unwrap();
                        let readable_string = String::from_utf8_lossy(&readable);
                        let read = Read{rectangle: (lt, lb, rt ,rb), read_id: String::from_utf8_lossy(&bam.name()).to_string(), start: bam.start(), end: bam.calculate_end(), insertions: insertions, strand: bam.flag().is_reverse_strand(), flag: bam.flag().0, track: data.0, mapq: bam.mapq(), query_len: bam.query_len(), sa: sastr, cigar: readable_string.to_string()};//bam.tags().get(b"SA").}; //, tags: tags  }
                        reads.push(read);
                    }
                });
            (bars, texts)
        };
        if dump_json {
            json.push(Area {
                pileups: reads,
                annotations,
            });
        }
        chart.draw_series(bars)?;
        chart.draw_series(texts)?;
        //let dump = {reads: [], annotation: };
        let end1 = start.elapsed();
        eprintln!("{}.{:03} sec.", end1.as_secs(), end1.subsec_millis());

        if freq_len > 0 {
            let coverages = coverage.split_evenly((freq_len, 1));
            for (i, &sample_sequential_id) in freq_len_ids.iter().rev().enumerate() {
                //}
                //for (i, (sample_sequential_id, values)) in frequency.iter().rev().enumerate() {
                let idx = *sample_sequential_id as usize;
                /*let y_max = match max_coverage {
                    Some(a) => a,
                    _ => values.iter().map(|t| t.1).max().unwrap_or(1),
                };*/
                /*let absolute_index = freq_len_ids
                .iter()
                .position(|x| x == &sample_sequential_id)
                .unwrap();*/
                let mut chart = ChartBuilder::on(&coverages[i])
                    // Set the caption of the chart
                    //.caption(format!("{}", range), ("sans-serif", 20).into_font())
                    // Set the size of the label region
                    .x_label_area_size(0)
                    .top_x_label_area_size(5)
                    .y_label_area_size(y_area_size)
                    // Finally attach a coordinate on the drawing area and make a chart context
                    .build_cartesian_2d(x_spec.clone(), 0..y_max)?;
                chart
                    .configure_mesh()
                    // We can customize the maximum number of labels allowed for each axis
                    //.x_labels(5)
                    //.y_labels(4)
                    //.x_label_style(("sans-serif", x_scale / 2).into_font())
                    .y_label_style(("sans-serif", 12).into_font())
                    // We can also change the format of the label text
                    // .x_label_formatter(&|x| format!("{:.3}", x))
                    .draw()?;
                if let Some(values) = frequency.get(sample_sequential_id) {
                    let color = if snp_frequency.is_some() {
                        Palette99::pick(idx).stroke_width(1)
                    } else {
                        Palette99::pick(idx).filled()
                    }; // BLUE
                    chart
                        .draw_series(
                            Histogram::vertical(&chart).style(color).margin(2).data(
                                values
                                    .iter()
                                    .filter(|t| {
                                        t.0 >= range.start() && t.0 < range.end() && t.2 == '*'
                                    })
                                    .map(|t| (t.0, t.1)),
                            ),
                        )?
                        .label(lambda(idx).unwrap_or(&idx.to_string()).to_string())
                        .legend(move |(x, y)| {
                            Rectangle::new(
                                [(x - 5, y - 5), (x + 5, y + 5)],
                                Palette99::pick(idx).filled(),
                            )
                        });
                    if let Some(_f) = snp_frequency {
                        /*eprintln!("{:?}", values.iter()
                        .filter(|t| {
                            t.0 >= range.start()
                                && t.0 < range.end()
                        }).collect::<Vec<_>>());*/
                        for i in ['A', 'C', 'G', 'T'].iter() {
                            let color = nt_color(*i, &preset_color).unwrap();
                            chart.draw_series(
                                Histogram::vertical(&chart)
                                    .style(color.filled())
                                    .margin(1)
                                    .data(
                                        values
                                            .iter()
                                            .filter(|t| {
                                                t.0 >= range.start()
                                                    && t.0 < range.end()
                                                    && t.2 == *i
                                            })
                                            .map(|t| (t.0, t.1)),
                                    ),
                            )?;
                        }
                    }

                    chart.draw_series(
                        Histogram::vertical(&chart)
                            .style(preset_color.pick(VisColor::SplCol).filled())
                            .margin(1)
                            .data(
                                split_frequency
                                    .iter()
                                    .filter(|&&t| {
                                        t.0 == *sample_sequential_id
                                            && (t.1).0 >= range.start()
                                            && (t.1).0 < range.end()
                                    })
                                    .map(|t| t.1),
                            ),
                    )?;

                    /*if snp_frequency {
                        [('A', A_COL), ]
                    }*/
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
                            .label_font(("sans-serif", 20))
                            .draw()?;
                    }
                }
            }
        }

        if legend && (!no_margin || annotation_count + axis_count > 0) {
            chart
                .configure_series_labels()
                .background_style(&WHITE.mix(0.8))
                .border_style(&BLACK)
                .draw()?;
        }
        let end2 = start.elapsed();
        eprintln!("{}.{:03} sec.", end2.as_secs(), end2.subsec_millis());
    }

    if dump_json {
        let mut json_path = PathBuf::from(output);
        json_path.set_extension("json");
        serde_json::to_writer(&File::create(json_path)?, &json)?
    }
    Ok(())
}
