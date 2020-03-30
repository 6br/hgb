use std::collections::HashMap;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct VirtualOffset(u64);

/// Chunk `[start-end)`, where `start` and `end` are [virtual offsets](struct.VirtualOffset.html).
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Chunk {
    sample_id: u32,
    file_id: u32, // When we split files of inverted data structure
    // format_type: u32, // Enum
    start: VirtualOffset,
    end: VirtualOffset, // Or, "length"?
}

/// Single bin that stores chunks of 
#[derive(Clone)]
pub struct Bin {
    bin_id: u32,
    chunks: Vec<Chunk>, // chunks for each sample, each format type
}

/// Index for a single reference sequence. Contains [bins](struct.Bin.html).
#[derive(Clone)]
pub struct Reference {
    bins: HashMap<u32, Bin>, // bin_id is not continuous.
    // linear_index: LinearIndex,
}

#[derive(Clone)]
pub struct Index {
    samples: Vec<String>, // sample names, we expect sample id is continuous.
    references: Vec<Reference>, // the length is inherited from n_ref
}


/// Returns a BAI bin for the record with alignment `[beg-end)`.
pub fn region_to_bin(beg: i32, end: i32) -> u32 {
    let end = end - 1;
    let mut res = 0_i32;
    for i in (14..27).step_by(3) {
        if beg >> i == end >> i {
            res = ((1 << 29 - i) - 1) / 7 + (beg >> i);
            break;
        }
    }
    res as u32
}

/// Returns all possible BAI bins for the region `[beg-end)`.
pub fn region_to_bins(start: i32, end: i32) -> BinsIter {
    BinsIter {
        i: -1,
        t: 0,
        start,
        end,
        curr_bin: 0,
        bins_end: 0,
    }
}

/// Iterator over bins.
pub struct BinsIter {
    i: i32,
    t: i32,
    start: i32,
    end: i32,
    curr_bin: u32,
    bins_end: u32,
}

impl Iterator for BinsIter {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr_bin == self.bins_end {
            if self.i >= 4 {
                return None;
            }
            self.i += 1;
            self.t += 1 << (self.i * 3);
            self.curr_bin = (self.t + (self.start >> 26 - 3 * self.i)) as u32 - 1;
            self.bins_end = (self.t + (self.end >> 26 - 3 * self.i)) as u32;

            if self.i == 0 {
                return Some(0);
            }
        }
        self.curr_bin += 1;
        Some(self.curr_bin)
    }
}

impl std::iter::FusedIterator for BinsIter {}

/// Maximal possible bin value.
pub const MAX_BIN: u16 = 37448;

/// Returns a maximal region for a given bin.
pub fn bin_to_region(bin: u16) -> (i32, i32) {
    if bin == 0 {
        return (std::i32::MIN, std::i32::MAX);
    }
    let mut left = 1;
    for i in 1..6 {
        let right = left + (1 << 3 * i);
        if bin >= left && bin < right {
            let beg = (bin - left) as i32;
            return (beg << 29 - 3 * i, beg + 1 << 29 - 3 * i);
        }
        left = right;
    }
    panic!("Bin id should be not bigger than MAX_BIN ({} > {})", bin, MAX_BIN);
}