use crate::index::Region;
use crate::range::Default;
use crate::{
    bed, range::Format, reader::IndexedReader, twopass_alignment::Alignment, vis::RecordIter, Vis,
};
use bam::{record::tags::TagValue, Record};
use clap::ArgMatches;
use genomic_range::StringRegion;
use itertools::Itertools;
use log::debug;
use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fs::File,
    io::BufReader,
};

pub struct ChromosomeBuffer {
    ref_id: u64,
    matches: ArgMatches,
    bins: BTreeMap<usize, Vec<(u64, bed::Record)>>,
    freq: BTreeMap<u64, Vec<(u64, u32, char)>>,
    reader: IndexedReader<BufReader<File>>,
}

impl ChromosomeBuffer {
    pub fn new(reader: IndexedReader<BufReader<File>>, matches: ArgMatches) -> Self {
        ChromosomeBuffer {
            ref_id: 0,
            matches,
            bins: BTreeMap::new(),
            freq: BTreeMap::new(),
            reader,
        }
    }

    // If the region is completely not overlapped,
    fn drop(&mut self) {
        eprintln!("Buffer has dropped.");
        self.ref_id = 0;
        self.bins = BTreeMap::new();
        self.freq = BTreeMap::new();
    }

    fn set_ref_id(&mut self, ref_id: u64) {
        self.ref_id = ref_id;
    }

    // This range is completely included or not?
    fn included(&self, range: Region) -> bool {
        if range.ref_id() != self.ref_id {
            return false;
        }
        let bins = self.bins.keys().collect::<Vec<_>>();
        let bins_iter =
            self.reader.index().references()[range.ref_id() as usize].region_to_bins(range);
        for i in bins_iter {
            if i.bin_disp_range().into_iter().any(|t| !bins.contains(&&t)) {
                return false;
            }
        }
        true
    }

    //This range is at least once overlaps with bins?
    /*
    fn overlap(&mut self, range: Region) -> bool {
        if range.ref_id() != self.ref_id {
            return false;
        }
        let bins = self.bins.keys().collect::<Vec<_>>();
        let bins_iter =
            self.reader.index().references()[range.ref_id() as usize].region_to_bins(range);
        for i in bins_iter {
            if bins.iter().any(|t| i.bin_disp_range().contains(t)) {
                return true;
            }
        }
        false
    }
    */

    pub fn bins(&self) -> Vec<&usize> {
        self.bins.keys().collect::<Vec<&usize>>()
    }

    fn size_limit(&self) -> bool {
        //std::mem::size_of(self)
        //This is a very heuristic way.
        debug!("Current bin size: {}", self.bins.keys().len());
        self.bins.keys().len() > 5000
    }

    fn size_limit_local(&self, bins: &BTreeSet<u64>) -> bool {
        //std::mem::size_of(self)
        //This is a very heuristic way.
        debug!("Current local bin size: {}", bins.iter().len());
        bins.iter().len() > 2000
    }

    pub fn add(&mut self, range: &StringRegion) -> (bool, Vec<(u64, Record)>, BTreeSet<u64>) {
        debug!("Add: {}", range);
        let closure = |x: &str| self.reader.reference_id(x);
        let reference_name = &range.path;
        let range = Region::convert(&range, closure).unwrap();
        let mut reload_flag = false;
        // Check if overlap, and drop them if the chrom_id is different.
        if range.ref_id() != self.ref_id || self.size_limit() {
            self.drop();
            self.set_ref_id(range.ref_id());
            reload_flag = true;
        }
        let matches = self.matches.clone();
        let min_read_len = matches
            .value_of("min-read-length")
            .and_then(|a| a.parse::<u32>().ok())
            .unwrap_or(0u32);
        let snp_frequency = matches
            .value_of("snp-frequency")
            .and_then(|a| a.parse::<f64>().ok());
        let no_bits = matches
            .value_of("no-bits")
            .and_then(|t| t.parse::<u16>().ok())
            .unwrap_or(1796u16);

        let mut chunks = BTreeMap::new();
        let mut bin_ids = BTreeSet::new();

        for i in
            self.reader.index().references()[range.ref_id() as usize].region_to_bins(range.clone())
        {
            if let Some(bins) = i.slice {
                for bin in bins {
                    if !self.bins.contains_key(&(bin.bin_id() as usize)) {
                        chunks.insert(bin.bin_id(), bin.clone().chunks_mut());
                        bin_ids.insert(bin.bin_id() as u64);
                    }
                }
            }
        }

        let mut merged_list = vec![];

        for (bin_id, chunks) in chunks {
            let viewer = self.reader.chunk(chunks).unwrap();
            let mut ann = vec![];
            let mut list = vec![];
            let sample_ids_opt: Option<Vec<u64>> = matches
                .values_of("id")
                .map(|a| a.map(|t| t.parse::<u64>().unwrap()).collect());
            let sample_id_cond = sample_ids_opt.is_some();
            let sample_ids = sample_ids_opt.unwrap_or_default();
            let filter = matches.is_present("filter");

            let format_type_opt = matches.value_of_t::<Format>("type");
            let format_type_cond = format_type_opt.is_ok();
            let format_type = format_type_opt.unwrap_or(Format::Default(Default {}));

            let _ = viewer.into_iter().for_each(|t| {
                let f = t.unwrap();
                if !sample_id_cond || sample_ids.iter().any(|&i| i == f.sample_id()) {
                    let sample_id = f.sample_id();
                    // eprintln!("{:?}", sample_id);
                    let data = f.data();
                    if !format_type_cond
                        || std::mem::discriminant(&format_type) == std::mem::discriminant(&data)
                    {
                        match data {
                            Format::Range(rec) => {
                                for i in rec.to_record(reference_name) {
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
                                    // && i.flag().no_bits(no_bits)
                                    // && i.query_len() > min_read_len
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

            // Before append into BTreeMap, we need to calculate pileup.
            list.iter().group_by(|elt| elt.0).into_iter().for_each(|t| {
                /*let mut line =
                Vec::with_capacity((prefetch_range.end - prefetch_range.start + 1) as usize);*/

                let line = self.freq.entry(t.0).or_insert_with(Vec::new);
                for column in bam::Pileup::with_filter(&mut RecordIter::new(t.1), |record| {
                    record.flag().no_bits(1796) //&& record.query_len() > min_read_len
                }) {
                    let column = column.unwrap();
                    if let Some(freq) = snp_frequency {
                        let mut seqs: Vec<String> = Vec::with_capacity(column.entries().len()); //vec![];
                        for entry in column.entries().iter() {
                            let seq: Option<_> = entry.sequence();
                            match seq {
                                Some(a) => {
                                    let seq: Vec<_> = a.map(|nt| nt as char).take(1).collect();
                                    seqs.push(seq.into_iter().collect());
                                }
                                _ => seqs.push(String::new()),
                            }
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
                                (car.1.as_bytes()[0] as char).to_ascii_uppercase(),
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
                    line.push((column.ref_pos() as u64, column.entries().len() as u32, '*'));
                }
            });
            merged_list.extend(list);
            self.bins.insert(bin_id as usize, ann);
        }
        (reload_flag, merged_list, bin_ids)
    }

    pub fn included_string(&self, string_range: &StringRegion) -> bool {
        let closure = |x: &str| self.reader.reference_id(x);
        let _reference_name = &string_range.path;
        let range = Region::convert(string_range, closure).unwrap();
        self.included(range)
    }

    // This range is completely included or not?
    fn included_local(&self, range: Region, bins: &BTreeSet<u64>) -> bool {
        if range.ref_id() != self.ref_id {
            return false;
        }
        let bins_iter =
            self.reader.index().references()[range.ref_id() as usize].region_to_bins(range);
        for i in bins_iter {
            if i.bin_disp_range()
                .into_iter()
                .any(|t| !bins.contains(&(t as u64)))
            {
                //if bins.iter().any(|t| !i.bin_disp_range().contains(t)) {
                return false;
            }
        }
        true
    }

    pub fn included_string_local(&self, string_range: &StringRegion, bins: &BTreeSet<u64>) -> bool {
        let closure = |x: &str| self.reader.reference_id(x);
        let _reference_name = &string_range.path;
        let range = Region::convert(string_range, closure).unwrap();
        self.included_local(range, bins)
    }

    pub fn add_local(
        &mut self,
        range: &StringRegion,
        local_bins: &mut BTreeSet<u64>,
    ) -> (bool, Vec<(u64, Record)>) {
        debug!("Add: {}", range);
        let closure = |x: &str| self.reader.reference_id(x);
        let reference_name = &range.path;
        let range = Region::convert(&range, closure).unwrap();
        let mut reload_flag = false;
        // Check if overlap, and drop them if the chrom_id is different.
        if self.size_limit_local(local_bins) {
            //self.drop();
            *local_bins = BTreeSet::new();
            reload_flag = true;
        }
        let matches = self.matches.clone();
        let min_read_len = matches
            .value_of("min-read-length")
            .and_then(|a| a.parse::<u32>().ok())
            .unwrap_or(0u32);
        let no_bits = matches
            .value_of("no-bits")
            .and_then(|t| t.parse::<u16>().ok())
            .unwrap_or(1796u16);
        //let mut chunks = BTreeMap::new();
        let mut chunks = vec![];

        for i in
            self.reader.index().references()[range.ref_id() as usize].region_to_bins(range.clone())
        {
            if let Some(bins) = i.slice {
                for bin in bins {
                    if !local_bins.contains(&(bin.bin_id() as u64)) && !reload_flag {
                        chunks.extend(bin.clone().chunks_mut());
                        local_bins.insert(bin.bin_id() as u64);
                    }
                }
            }
        }
        //eprintln!("3");

        let mut merged_list = vec![];

        let viewer = self.reader.chunk(chunks).unwrap();
        let mut ann = vec![];
        let sample_ids_opt: Option<Vec<u64>> = matches
            .values_of("id")
            .map(|a| a.map(|t| t.parse::<u64>().unwrap()).collect());
        let sample_id_cond = sample_ids_opt.is_some();
        let sample_ids = sample_ids_opt.unwrap_or_default();
        let filter = matches.is_present("filter");

        let format_type_opt = matches.value_of_t::<Format>("type");
        let format_type_cond = format_type_opt.is_ok();
        let format_type = format_type_opt.unwrap_or(Format::Default(Default {}));

        let _ = viewer.into_iter().for_each(|t| {
            let f = t.unwrap();
            if !sample_id_cond || sample_ids.iter().any(|&i| i == f.sample_id()) {
                let sample_id = f.sample_id();
                let data = f.data();
                if !format_type_cond
                    || std::mem::discriminant(&format_type) == std::mem::discriminant(&data)
                {
                    match data {
                        Format::Range(rec) => {
                            for i in rec.to_record(reference_name) {
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
                                //&& i.flag().no_bits(no_bits)
                                //&& i.query_len() > min_read_len
                                {
                                    merged_list.push((sample_id, i));
                                }
                            }
                        }
                        _ => {}
                    }
                }
            }
        });
        (reload_flag, merged_list)
    }

    pub fn retrieve(
        &mut self,
        string_range: &StringRegion,
        list: &mut Vec<(u64, bam::Record)>,
        list_btree: &mut BTreeSet<u64>,
    ) {
        let closure = |x: &str| self.reader.reference_id(x);
        let _reference_name = &string_range.path;
        let range = Region::convert(string_range, closure).unwrap();
        debug!("List len: {}", list.len());

        if !self.included(range.clone()) {
            let (reset_flag, new_list, bins) = self.add(string_range);
            debug!("Loaded len: {}", new_list.len());
            if reset_flag {
                *list = new_list;
                *list_btree = bins;
            } else {
                list.extend(new_list);
                list_btree.extend(bins)
            }
            debug!("After load list len: {}", list.len());
        }

        if !self.included_local(range, list_btree) {
            let (reset_flag, new_list) = self.add_local(string_range, list_btree);
            if reset_flag {
                *list = new_list
            } else {
                list.extend(new_list);
            }
            debug!("After load list len: {}", list.len());
        }
    }
    pub fn vis(
        &self,
        string_range: &StringRegion,
        list: &mut Vec<(u64, bam::Record)>,
        _list_btree: &mut BTreeSet<u64>,
    ) -> Option<Vis> {
        let closure = |x: &str| self.reader.reference_id(x);
        let _reference_name = &string_range.path;
        let range = Region::convert(string_range, closure).unwrap();
        let _freq = &self.freq;

        let mut ann: Vec<(u64, bed::Record)> =
            self.bins.values().into_iter().cloned().flatten().collect();
        let matches = self.matches.clone();
        let _split_only = matches.is_present("only-split-alignment");
        let sort_by_name = matches.is_present("sort-by-name");
        let packing = !matches.is_present("no-packing");
        let split = matches.is_present("split-alignment");
        let max_coverage = matches
            .value_of("max-coverage")
            .and_then(|a| a.parse::<u32>().ok());
        let min_read_len = matches
            .value_of("min-read-length")
            .and_then(|a| a.parse::<u32>().ok())
            .unwrap_or(0u32);
        let no_bits = matches
            .value_of("no-bits")
            .and_then(|t| t.parse::<u16>().ok())
            .unwrap_or(1796u16);
        // Calculate coverage; it won't work on sort_by_name
        // let mut frequency = BTreeMap::new(); // Vec::with_capacity();

        ann.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));
        list.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.start().cmp(&b.1.start())));

        //TODO(FIX: pileup changed to append )

        //eprintln!("{:?}", freq.keys());
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
                    let sample_id = t.0;
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
            /*
            if split_only {
                list = list
                    .into_iter()
                    .filter(|(sample_id, record)| {
                        end_map.contains_key(&(*sample_id, record.name()))
                    })
                    .collect();
            }
            */
            list.iter()
                .filter(|(_, i)| i.flag().no_bits(no_bits) && i.query_len() > min_read_len)
                .group_by(|elt| elt.0)
                .into_iter()
                .for_each(|t| {
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
        } else if packing {
            list.iter()
                .filter(|(_, i)| i.flag().no_bits(no_bits) && i.query_len() > min_read_len)
                .group_by(|elt| elt.0)
                .into_iter()
                .for_each(|t| {
                    // let mut heap = BinaryHeap::<(i64, usize)>::new();
                    let mut packing = vec![0u64];
                    prev_index += 1;
                    (t.1).for_each(|k| {
                        let mut index =
                            if let Some(TagValue::Int(array_view, _)) = k.1.tags().get(b"YY") {
                                array_view as usize
                            } else if let Some(index) = packing
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
                    });
                    if let Some(max_cov) = max_coverage {
                        prev_index = max_cov as usize + last_prev_index;
                    }
                    compressed_list.push((t.0, prev_index));
                    //eprintln!("{:?} {:?} {:?}", compressed_list, packing, index_list);
                    last_prev_index = prev_index;
                });
        } else {
            // Now does not specify the maximal length by max_coverage.
            index_list = (0..list.len()).collect();

            // eprintln!("{}", list.len());
            list.iter().group_by(|elt| elt.0).into_iter().for_each(
                |(sample_sequential_id, sample)| {
                    let count = sample.count();
                    prev_index += count;
                    compressed_list.push((sample_sequential_id, prev_index));
                },
            )
        }

        return Some(Vis {
            range: string_range.clone(),
            //list: list,
            annotation: ann,
            freq: self.freq.clone(),
            compressed_list,
            index_list,
            prev_index,
            supplementary_list,
            prefetch_max: self.reader.header().reference_len(0).unwrap(), // The max should be the same as the longest ?
        });
    }
}
