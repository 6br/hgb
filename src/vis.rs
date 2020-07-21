use crate::bed;
use bam::record::{
    tags::{StringType, TagValue},
    Cigar,
};
use bam::Record;
use clap::ArgMatches;
use genomic_range::StringRegion;
use itertools::Itertools;
use plotters::coord::ReverseCoordTranslate;
use plotters::prelude::Palette;
//use plotters::prelude::RGBAColor;
use plotters::prelude::*;
use plotters::style::{RGBAColor, RGBColor};
use rayon::prelude::*;
use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    mem,
    sync::{Arc, Mutex},
};

pub struct Vis<'a> {
    range: StringRegion,
    list: Mutex<Vec<(u64, Record)>>,
    annotation: Mutex<Vec<(u64, bed::Record)>>,
    index_list: Mutex<Vec<usize>>,
    suppl_list: Mutex<Vec<(&'a [u8], usize, usize, usize, usize)>>,
}

impl<'a> Vis<'a> {
    pub fn new(
        range: StringRegion,
        list: Vec<(u64, Record)>,
        annotation: Vec<(u64, bed::Record)>,
    ) -> Self {
        Vis {
            range,
            list: Mutex::new(list),
            annotation: Mutex::new(annotation),
            index_list: Mutex::new(vec![]),
            suppl_list: Mutex::new(vec![]),
        }
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
//RGBColor();

pub fn bam_record_vis<'a, F>(
    matches: &ArgMatches,
    mut vis: Vec<Vis>,
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
    let split_only = matches.is_present("only-split-alignment");
    let sort_by_name = matches.is_present("sort-by-name");
    let sort_by_cigar = matches.is_present("sort-by-cigar");
    let x = matches
        .value_of("x")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(1280u32);
    let y = matches
        .value_of("y")
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
    for i in vis.iter() {
        let mut list = i.list.lock().unwrap();
        let mut annotation = i.annotation.lock().unwrap();
        if sort_by_name {
            list.sort_by(|a, b| {
                a.0.cmp(&b.0)
                    .then(a.1.name().cmp(&b.1.name()))
                    .then(a.1.start().cmp(&b.1.start()))
            });
        } else {
            list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
        }

        annotation.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
    }

    // New_list is merged list to decide the order of alignments.
    let mut new_list = {
        // new_list is a tuple (sample_id, record, range_id), whose record needs to be cloned.
        let mut new_list = vec![];
        for (index, i) in vis.iter().enumerate() {
            let mut list = i.list.lock().unwrap();
            list.sort_by(|a, b| {
                a.0.cmp(&b.0)
                    .then(a.1.name().cmp(&b.1.name()))
                    .then(a.1.start().cmp(&b.1.start()))
            });
            new_list.append(
                &mut list
                    .iter()
                    .map(|t| (t.0, t.1.clone(), index))
                    .collect::<Vec<_>>(),
            );
        }
        // Should we mark the duplicated alignment as "this is duplicated so you sholdn't display more than once"?

        new_list
    };
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
    let mut compressed_list = vec![]; //This is shared among samples.
                                      // let mut index_list = Vec::with_capacity(list.len());
    let mut supplementary_list = vec![];

    if split {
        let mut end_map = HashMap::new();

        let tmp_list = new_list.clone();
        tmp_list
            .iter()
            .group_by(|elt| elt.0)
            .into_iter()
            .for_each(|t| {
                let sample_id = t.0.clone();
                t.1.group_by(|&elt| elt.1.name()).into_iter().for_each(|s| {
                    let items: Vec<&(u64, Record, usize)> = s.1.into_iter().collect();
                    if items.len() > 1 {
                        let last: &(u64, Record, usize) =
                            items.iter().max_by_key(|t| t.1.calculate_end()).unwrap();
                        end_map.insert(
                            (sample_id, s.0),
                            (
                                items[0].1.calculate_end(), // The end of first item
                                last.1.start(),             // The start of last item
                                last.1.calculate_end(),     // The end of last item
                                items.len(), // How many alignments correspond to one read
                            ),
                        );
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
                        (a.1.cigar().soft_clipping(true) + a.1.cigar().hard_clipping(true)).cmp(
                            &((b.1.cigar().soft_clipping(true))
                                + (b.1.cigar().hard_clipping(true))),
                        ),
                    )
                });
            }
        } else {
            new_list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
        }

        if split_only {
            new_list = new_list
                .into_iter()
                .filter(|(sample_id, record, _range_id)| {
                    end_map.contains_key(&(*sample_id, record.name()))
                })
                .collect();
        }

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
                        end.2 as u64 + ((vis.len() as u64) << 32)
                    } else {
                        k.1.calculate_end() as u64 + ((k.2 as u64) << 32)
                    };

                    let index = if sort_by_name {
                        prev_index += 1;
                        e
                    } else if let Some(index) = name_index.get(k.1.name()) {
                        *index
                    } else if let Some(index) = packing_vec
                        .iter_mut()
                        .enumerate()
                        .find(|(_, item)| **item < k.1.start() as u64)
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
                                k.1.name(),
                                index + last_prev_index,
                                index + last_prev_index + end.3,
                                end.0,
                                end.1,
                            ));
                        }
                    }
                    vis[k.2]
                        .index_list
                        .lock()
                        .unwrap()
                        .push(index + last_prev_index);
                    name_index.insert(k.1.name(), index);
                });
                compressed_list.push((t.0, prev_index));
                last_prev_index = prev_index;
            });
    } else {
        if packing {
            new_list
                .iter()
                .group_by(|elt| elt.0)
                .into_iter()
                .for_each(|t| {
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
                    compressed_list.push((t.0, prev_index));
                    //eprintln!("{:?} {:?} {:?}", compressed_list, packing, index_list);
                    last_prev_index = prev_index;
                    //(t.0, ((t.1).0, (t.1).1))
                    // .collect::<&(u64, Record)>(). // collect::<Vec<(usize, (u64, Record))>>
                });
        } else {
            //index_list = (0..list.len()).collect();

            // Now we asssume that we'd like to vertically stack all of reads.
            for i in vis.iter() {
                let mut index_list = i.index_list.lock().unwrap(); // = (0..new_list.len()).collect();
                let temp_index_list = (0..new_list.len()).collect();
                mem::replace(&mut *index_list, temp_index_list);
            }

            // list.sort_by(|a, b| a.0.cmp(&b.0));
            // eprintln!("{}", list.len());
            new_list.iter().group_by(|elt| elt.0).into_iter().for_each(
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
        Arc::new(btree)
    //        graph.into_iter().collect().unique_by(|s| s[0]).count() // group_by(|elt| elt[0]).count() as usize
    } else {
        Arc::new(BTreeMap::new())
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

    // We asssume that the length of annotation is uniform.
    let annotation_count = vis[0]
        .annotation
        .lock()
        .unwrap()
        .iter()
        .unique_by(|s| s.0)
        .count(); // annotation.len();

    let root = BitMapBackend::new(
        output,
        (
            x * vis.len() as u32,
            40 + (prev_index as u32 + axis_count as u32 + annotation_count as u32 * 2) * y,
        ),
    )
    .into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(10, 10, 10, 10);
    let areas = root.split_evenly((1, vis.len()));
    // After this point, we should be able to draw construct a chart context
    // let areas = root.split_by_breakpoints([], compressed_list);
    eprintln!("{:?} {:?} {:?}", prev_index, axis_count, annotation_count);

    for (index, area) in areas.iter().enumerate() {
        let range = &vis[index].range;
        let index_list = vis[index].index_list.lock().unwrap();
        let list = &*vis[index].list.lock().unwrap();
        let annotation = &*vis[index].annotation.lock().unwrap();

        // let
        let mut chart = ChartBuilder::on(&area)
            // Set the caption of the chart
            .caption(format!("{}", range), ("sans-serif", 20).into_font())
            // Set the size of the label region
            .x_label_area_size(20)
            .y_label_area_size(40)
            // Finally attach a coordinate on the drawing area and make a chart context
            .build_ranged(
                (range.start() - 1)..(range.end() + 1),
                0..(1 + prev_index + axis_count + annotation_count * 2),
            )?;
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

        // Draw annotation axis if there is graph information.
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
                                    stroke.stroke_width(y / 2),
                                ))
                                .unwrap()
                                .label(format!("{}", record.name().unwrap_or(&"")))
                                .legend(move |(x, y)| {
                                    Rectangle::new(
                                        [(x - 5, y - 5), (x + 5, y + 5)],
                                        stroke.filled(),
                                    )
                                });
                            //Some(bar2)
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
                    .iter()
                    .map(|(sample_sequential_id, sample)| {
                        let count = *sample;
                        if count > 0 {
                            let stroke = Palette99::pick(*sample_sequential_id as usize);
                            let mut bar2 = Rectangle::new(
                                [
                                    (range.start() - 1, prev_index),
                                    (range.end() + 1, prev_index + count),
                                ],
                                stroke.stroke_width(7), // filled(), //stroke_width(100),
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
        // For each alignment:
        let series = {
            //list.into_iter().enumerate().map(|(index, data)| {
            let mut bars = vec![];
            (*index_list).iter().zip(list).for_each(|(&index, data)| {
                //chart.draw_series(index_list.into_par_iter().zip(list).map(|(index, data)| {
                //for (index, data) in list.iter().enumerate() {
                let bam = &data.1;
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
                            let color = SPL_COL;
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
                } /*
                  if legend {
                  } else {
                      let mut bar2 =
                          Rectangle::new([(start, index), (end, index + 1)], stroke.stroke_width(2));
                      bar2.set_margin(1, 1, 0, 0);
                      //vec![bar,bar2]
                      bars.push(bar2);
                  };*/
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

        if legend {
            chart
                .configure_series_labels()
                .background_style(&WHITE.mix(0.8))
                .border_style(&BLACK)
                .draw()?;
        }
    }
    Ok(())
}
