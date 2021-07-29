use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Area {
    pub annotations: Vec<Annotation>,
    pub pileups: Vec<Read>,
}
/*
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Pileup {
    pub reads: Vec<Read>,
}
*/
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Read {
    pub rectangle: (i32, i32, i32, i32),
    pub read_id: String,
    pub start: i32,
    pub end: i32,
    pub insertions: BTreeMap<u64, String>,
    // pub tags
}
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Annotation {
    pub rectangle: (i32, i32, i32, i32),
    pub start: u64,
    pub end: u64,
    pub name: String,
}
