use crate::gff;
use bio_types::strand::*;
use genomic_range::StringRegion;
//use lib::{Config, ConfigFeature};
//use lib::{Database, GeneNameEachReference, GeneNameTree, Region};
//use rocks::rocksdb::*;
use log::{debug, info};
use std::collections::{BTreeMap, HashMap};
use std::path::Path;
//use vg::GraphDB;
//use vg::GraphDB::VG;

pub type GeneNameTree = BTreeMap<String, StringRegion>;
pub type GeneNameEachReference = BTreeMap<String, GeneNameTree>;
// NodeId to corresponding feature items.
type Features = HashMap<u64, Vec<Feature>>;
pub type FeatureDB = Vec<Features>;

#[derive(Debug, PartialEq)]
pub struct Database {
    pub gene_name_tree: GeneNameEachReference,
}

impl Database {
    pub fn new() -> Database {
        let mut gene_name_tree = BTreeMap::new();
        gene_name_tree.insert("default".to_string(), BTreeMap::new());
        Database { gene_name_tree }
    }

    pub fn convert(self: &Database, str: &String) -> Option<&StringRegion> {
        self.gene_name_tree["default"].get(str)
    }
}

#[derive(Debug, PartialEq)]
pub struct Feature {
    pub start_offset: u64,
    pub stop_offset: u64,
    pub id: u64,
    pub name: String,
    pub is_reverse: Option<bool>,
    pub attributes: Vec<String>,
    pub value: Option<f32>,
}

#[derive(Debug, PartialEq)]
struct FeatureSet {
    feature_set_id: u64,
    dataset_id: Vec<u64>,
    attributes: Option<String>,
}

fn opt_strand_to_opt_bool(strand: Option<Strand>) -> Option<bool> {
    strand.and_then(|strand| match strand {
        Strand::Forward => Some(false),
        Strand::Reverse => Some(true),
        Strand::Unknown => None,
    })
}
/*
fn record_to_nodes(
    record: bed::Record,
    coord_map: &CoordToNodeId,
    bed_id: u64,
    chr_prefix: &Option<String>,
) -> HashMap<u64, Feature> {
    let mut hash_map: HashMap<u64, Feature> = HashMap::new();
    let chr = match *chr_prefix {
        Some(ref k) => record.chrom().replace(k, ""),
        None => record.chrom().to_string(),
    };
    let ref vec = match coord_map.get(&chr) {
        Some(k) => k,
        None => return hash_map,
    };
    let lower_bound_index = match vec.binary_search_by_key(&record.start(), |b| b.coord) {
        Ok(x) => x,
        Err(x) => x,
    };

    hash_map.insert(
        vec[lower_bound_index].id,
        Feature {
            start_offset: vec[lower_bound_index].coord - record.start(),
            stop_offset: 0,
            id: bed_id,
            name: record.name().unwrap_or_default().to_string(),
            is_reverse: opt_strand_to_opt_bool(record.strand()),
            attributes: vec![],
            value: None,
        },
    );

    let mut index = lower_bound_index;
    while vec.len() > index + 1 && vec[index + 1].coord < record.end() {
        index += 1;
        hash_map.insert(
            vec[index].id,
            Feature {
                start_offset: 0,
                stop_offset: 0,
                id: bed_id,
                name: record.name().unwrap_or_default().to_string(),
                is_reverse: opt_strand_to_opt_bool(record.strand()),
                attributes: vec![],
                value: None,
            },
        );
    }
    return hash_map;
}
*/
// tmpNew should be replicated with a novel implementation.
// Required input list is sorted by coordinates.
//pub fn tmp_new(graph: Arc<Graph>, config: &Config) -> Database {
pub fn tmp_new(features: Vec<String>) -> Database {
    let mut gene_per_ref = GeneNameEachReference::new();
    //for data in config.reference.data.iter() {
    let mut gene: GeneNameTree = GeneNameTree::new();
    for feature in features.into_iter() {
        // It limits only "config,reference Items."
        let path = Path::new(&feature);
        info!("Parsing:  {:?}", path);
        match path.extension().unwrap_or_default().to_str() {
            /*Some("bed") => {
                vec.push(tmp_new_internal(feature, &graph, &hashmap));
            }*/
            Some("gff") => {
                tmp_new_gene_internal(feature, &mut gene, gff::GffType::GFF3);
            }
            Some("gff3") => {
                tmp_new_gene_internal(feature, &mut gene, gff::GffType::GFF3);
            }
            Some("gtf") => {
                tmp_new_gene_internal(feature, &mut gene, gff::GffType::GTF2);
            }
            _ => println!("Unsupported format {:?}", path),
        }
    }
    gene_per_ref.insert("default".to_string(), gene);
    //}
    return Database {
        gene_name_tree: gene_per_ref,
    };
}

// It includes only "gene" row.
fn tmp_new_gene_internal(url: String, gene: &mut GeneNameTree, gff_type: gff::GffType) {
    let gff3 = url;
    let path = Path::new(&gff3);
    let mut reader = match gff::Reader::from_file(path, gff_type) {
        Ok(f) => f,
        Err(e) => {
            debug!("could not open; skipping.");
            //return result?;
            return;
        }
    };
    let mut index = 0;
    for record in reader.records() {
        index += 1;
        match record {
            Ok(rec) => match rec.feature_type() {
                "gene" | "transcript" => {
                    let reg = match opt_strand_to_opt_bool(rec.strand()) {
                        Some(false) => StringRegion {
                            path: rec.seqname().to_string(),
                            end: *rec.start(),
                            start: *rec.end(),
                        },
                        _ => StringRegion {
                            path: rec.seqname().to_string(),
                            start: *rec.start(),
                            end: *rec.end(),
                        },
                    };
                    match rec.attributes().get("gene_name") {
                        Some(name) => gene.insert(name.clone().to_string(), reg),
                        None => continue,
                    };
                }
                _ => continue,
            },
            Err(_) => continue,
        }
    }
    debug!("{} lines processed. end.", index);
}
/*
fn tmp_new_internal(
    feature: &ConfigFeature,
    _graph: &GraphDB,
    hashmap: &CoordToNodeId,
) -> Features {
    let bed = &feature.url;
    let path = Path::new(&bed);
    let mut features: Features = Features::new();

    let mut reader = match bed::Reader::from_file(path) {
        Ok(f) => f,
        Err(e) => {
            debug!("could not open {}; skipping.", e.description());
            return features;
        }
    };
    let mut index: u64 = 0;

    for record in reader.records() {
        let rec = record.ok().expect("Error reading record.");
        let nodes = record_to_nodes(rec, &hashmap, index, &feature.chr_prefix);
        for (key, value) in nodes.into_iter() {
            features.entry(key).or_insert(Vec::new()).push(value);
        }
        index += 1;
    }
    return features;
}
*/
