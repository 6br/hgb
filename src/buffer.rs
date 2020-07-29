use crate::index::Region;
use crate::range::Default;
use crate::{
    bed, checker_index::Reference, header::Header, range::Format, reader::IndexedReader,
    twopass_alignment::Alignment, vis::RecordIter, Vis,
};
use clap::ArgMatches;
use genomic_range::StringRegion;
use itertools::Itertools;
use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    io::BufReader,
};

pub struct ChromosomeBuffer<'a> {
    chr_id: usize,
    matches: &'a ArgMatches,
    bins: BTreeMap<usize, (Vec<(u64, bam::Record)>, Vec<(u64, bed::Record)>)>,
    freq: BTreeMap<u64, Vec<(u64, u32)>>,
    header: Header,
    reader: &'a IndexedReader<BufReader<File>>,
}

impl<'a> ChromosomeBuffer<'a> {
    fn new(
        header: Header,
        reader: &'a IndexedReader<BufReader<File>>,
        matches: &'a ArgMatches,
    ) -> Self {
        ChromosomeBuffer {
            chr_id: 0,
            matches,
            bins: BTreeMap::new(),
            freq: BTreeMap::new(),
            header,
            reader,
        }
    }

    // This range is completely included or not?
    fn included(&mut self, range: Region) -> bool {
        let bins = self.bins.keys();
        let bins_iter =
            self.reader.index().references()[range.ref_id() as usize].region_to_bins(range);
        for i in bins_iter {
            if bins.any(|t| !i.bin_disp_range().contains(t)) {
                return false;
            }
        }
        true
    }

    //This range is at least once overlaps with bins?
    fn overlap(&mut self, range: Region) -> bool {
        let bins = self.bins.keys();
        let bins_iter =
            self.reader.index().references()[range.ref_id() as usize].region_to_bins(range);
        for i in bins_iter {
            if bins.any(|t| i.bin_disp_range().contains(t)) {
                return true;
            }
        }
        false
    }

    fn add(&mut self, range: Region) {
        let bins_iter =
            self.reader.index().references()[range.ref_id() as usize].region_to_bins(range);
        let matches = self.matches;
        let min_read_len = matches
            .value_of("min-read-length")
            .and_then(|a| a.parse::<u32>().ok())
            .unwrap_or(0u32);
        for i in bins_iter {
            for bins in i.slice {
                for bin in bins {
                    let viewer = self.reader.chunk(bin.chunks_mut()).unwrap();

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
                                || std::mem::discriminant(&format_type)
                                    == std::mem::discriminant(&data)
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

                    // Before append into BTreeMap, we need to calculate pileup.
                    list.iter().group_by(|elt| elt.0).into_iter().for_each(|t| {
                /*let mut line =
                    Vec::with_capacity((prefetch_range.end - prefetch_range.start + 1) as usize);*/
                    let line = self.freq.entry(bin.bin_id as u64s);
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
                    // Should we have sparse occurrence table?
                    //eprintln!("{:?} {:?}",  range.path, lambda(column.ref_id() as usize).unwrap_or(&column.ref_id().to_string()));
                    // lambda(column.ref_id() as usize).unwrap_or(&column.ref_id().to_string())
                    // == range.path
                    // &&
                    if prefetch_range.start <= column.ref_pos() as u64
                        && column.ref_pos() as u64 <= prefetch_range.end
                    {
                        line.push((column.ref_pos() as u64, column.entries().len() as u32));
                    }
                }
                //eprintln!("{:?}", line);
                freq.insert(t.0, line);
            });
                }
            }
        }
    }

    fn vis(&mut self, range: StringRegion) -> Vis {
        let closure = |x: &str| self.reader.reference_id(x);
        let _reference_name = &range.path;
        let range = Region::convert(&range, closure).unwrap();

        let list: Vec<(u64, bam::Record)> = self
            .bins
            .values()
            .into_iter()
            .map(|t| t.0)
            .flatten()
            .collect();
        let ann: Vec<(u64, bed::Record)> = self
            .bins
            .values()
            .into_iter()
            .map(|t| t.1)
            .flatten()
            .collect();
        let matches = self.matches;
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

        ann.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
        list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));

        //TODO(FIX: pileup changed to append )

        eprintln!("{:?}", freq.keys());
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
                            (a.1.cigar().soft_clipping(true) + a.1.cigar().hard_clipping(true))
                                .cmp(
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
                    .filter(|(sample_id, record)| {
                        end_map.contains_key(&(*sample_id, record.name()))
                    })
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

        Vis {
            range,
            list: list,
            annotation: ann,
            freq: self.freq,
            compressed_list: compressed_list,
            index_list: index_list,
            prev_index: prev_index,
            supplementary_list: suppementary_list,
            prefetch_max: self.reader.reference_len(range.path),
        }
    }

    // If the region is completely not overlapped,
    fn drop(&mut self) {}
}
