use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use genomic_range::StringRegion;
use regex::Regex;
use std::collections::BTreeMap;
use std::fmt::{self, Debug, Display, Formatter};
use std::fs::File;
use std::io::ErrorKind::InvalidData;
use std::io::{BufReader, BufWriter, Error, Read, Result, Write};
use std::path::Path;
use std::result;

/// Per BAM specification, bin with `bin_id == SUMMARY_BIN` contains summary over the reference.
// const SUMMARY_BIN: u32 = 37450; // TODO(FIXME)

/// Genomic coordinates, used in [struct.IndexedReader.html#method.fetch] and [struct.IndexedReader.html#method.pileup].
/// `ref_id` is 0-based, `start-end` is 0-based half-open interval.
#[derive(Clone, PartialEq, Debug)]
pub struct Region {
    ref_id: u64,
    start: u64,
    end: u64,
}

impl Region {
    /// Creates new region. `ref_id` is 0-based, `start-end` is 0-based half-open interval.
    pub fn new(ref_id: u64, start: u64, end: u64) -> Region {
        assert!(
            start <= end,
            "Region: start should not be greater than end ({} > {})",
            start,
            end
        );
        Region { ref_id, start, end }
    }

    pub fn convert<F>(
        path: &StringRegion,
        to_id: F,
    ) -> std::result::Result<Self, Box<dyn std::error::Error>>
    where
        F: Fn(&str) -> Option<u64>,
    {
        Ok(Region {
            ref_id: to_id(&path.path).ok_or("Error: the reference id is not recognized.")?,
            start: path.start,
            end: path.end,
        })
    }

    pub fn parse<F>(path: &str, to_id: F) -> std::result::Result<Self, Box<dyn std::error::Error>>
    where
        F: Fn(&str) -> Option<u64>,
    {
        let re = Regex::new(r"^(.+):(\d*)-?(\d*)$").unwrap();
        let caps = re.captures(path).ok_or("Parse Error")?;
        let path = caps
            .get(1)
            .and_then(|t| Some(t.as_str()))
            .ok_or("Parse Path Error")?;
        let start = caps
            .get(2)
            .and_then(|t| t.as_str().parse::<u64>().ok())
            .ok_or("Error: the reference start is not recognized.")?;
        let stop = caps
            .get(3)
            .and_then(|t| t.as_str().parse::<u64>().ok())
            .ok_or("Error: the reference end is not recognized.")?;

        return Ok(Region {
            ref_id: to_id(path).ok_or("Error: the reference id is not recognized.")?,
            start: start,
            end: stop,
        });
    }

    pub fn ref_id(&self) -> u64 {
        self.ref_id
    }

    pub fn start(&self) -> u64 {
        self.start
    }

    pub fn end(&self) -> u64 {
        self.end
    }

    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    pub fn set_ref_id(&mut self, ref_id: u64) {
        self.ref_id = ref_id;
    }

    pub fn set_start(&mut self, start: u64) {
        assert!(
            start <= self.end,
            "Region: start should not be greater than end ({} > {})",
            start,
            self.end
        );
        self.start = start;
    }

    pub fn set_end(&mut self, end: u64) {
        assert!(
            self.start <= end,
            "Region: start should not be greater than end ({} > {})",
            self.start,
            end
        );
        self.end = end;
    }

    pub fn contains(&self, ref_id: u64, pos: u64) -> bool {
        self.ref_id == ref_id && self.start <= pos && pos < self.end
    }

    pub fn include(&self, range: &Region) -> bool {
        self.ref_id == range.ref_id && self.start <= range.start && range.end < self.end
    }
}

/// Virtual offset. Represents `block_offset << 16 | contents_offset`, where
/// `block_offset` is `u48` and represents the offset in the bgzip file to the beginning of the
/// block (also known as `coffset` or `compressed_offset`).
///
/// `contents_offset` is `u16` and represents offset in the uncompressed data in a single block
/// (also known as `uoffset` or `uncompressed_offset`).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct VirtualOffset(u64);

impl VirtualOffset {
    /// Construct virtual offset from raw value.
    pub fn from_raw(raw: u64) -> VirtualOffset {
        VirtualOffset(raw)
    }

    fn from_stream<R: Read>(stream: &mut R) -> Result<Self> {
        Ok(VirtualOffset(stream.read_u64::<LittleEndian>()?))
    }

    fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        stream.write_u64::<LittleEndian>(self.raw())?;
        Ok(())
    }

    /// Constructs virtual offset from `block_offset` and `contents_offset`.
    pub fn new(block_offset: u64, contents_offset: u16) -> Self {
        VirtualOffset(block_offset << 16 | contents_offset as u64)
    }

    /// Returns the raw value.
    pub fn raw(&self) -> u64 {
        self.0
    }

    /// Returns the block offset. Represents the offset in the Bgzip file to the beginning of the block.
    pub fn block_offset(&self) -> u64 {
        self.0 >> 16
    }

    /// Represents the contents offset. Represents the offset into the uncompressed contents of the block.
    pub fn contents_offset(&self) -> u16 {
        self.0 as u16
    }

    /// Checks if the `self` is the same as `block_offset << 16 | contents_offset`.
    pub fn equal(&self, block_offset: u64, contents_offset: u16) -> bool {
        self.0 == (block_offset << 16 | contents_offset as u64)
    }

    /// Minimal possible offset, same as `VirtualOffset::from_raw(0)`.
    pub const MIN: VirtualOffset = VirtualOffset(0);
    /// Maximal possible offset, same as `VirtualOffset::from_raw(std::u64::MAX)`.
    pub const MAX: VirtualOffset = VirtualOffset(std::u64::MAX);
}

impl Display for VirtualOffset {
    fn fmt(&self, f: &mut Formatter) -> result::Result<(), fmt::Error> {
        write!(f, "c={},u={}", self.block_offset(), self.contents_offset())
    }
}

/// Chunk `[start-end)`, where `start` and `end` are [virtual offsets](struct.VirtualOffset.html).
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Copy)]
pub struct Chunk {
    start: VirtualOffset,
    end: VirtualOffset, // We don't need "end"
    sample_id: u64,
    file_id: u32, // When we split files of inverted data structure
                  // format_type: u32, // Enum
}

impl Chunk {
    /// Constructs a `Chunk` from two virtual offsets.
    pub fn new(sample_id: u64, file_id: u32, start: VirtualOffset, end: VirtualOffset) -> Self {
        Chunk {
            sample_id,
            file_id,
            start,
            end,
        }
    }

    pub fn from_stream<R: Read>(stream: &mut R, check: bool) -> Result<Self> {
        let sample_id = stream.read_u64::<LittleEndian>()?;
        let file_id = stream.read_u32::<LittleEndian>()?;
        let start = VirtualOffset::from_stream(stream)?;
        let end = VirtualOffset::from_stream(stream)?;
        if check && end <= start {
            Err(Error::new(
                InvalidData,
                format!("GHI chunk end < start ({}  <  {})", end, start),
            ))
        } else {
            Ok(Chunk {
                sample_id,
                file_id,
                start,
                end,
            })
        }
    }

    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        stream.write_u64::<LittleEndian>(self.sample_id)?;
        stream.write_u32::<LittleEndian>(self.file_id)?;
        self.start.to_stream(stream)?;
        self.end.to_stream(stream)?;
        Ok(())
    }

    /// Checks if two chunks intersect.
    pub fn intersect(&self, other: &Chunk) -> bool {
        self.start < other.end && other.start < self.end
    }

    /// Checks if two chunks intersect or one of the chunks goes right after another.
    pub fn can_combine(&self, other: &Chunk) -> bool {
        self.start <= other.end && other.start <= self.end
    }

    ///Checks if two chunks is adjacent
    pub fn is_adjacent(&self, next: &Chunk) -> bool {
        self.can_combine(next) && self.end == next.start
    }

    pub fn concat(&mut self, next: &Chunk) {
        self.end = next.end;
    }

    /// Returns the start of the chunk.
    pub fn start(&self) -> VirtualOffset {
        self.start
    }

    /// Returns the end of the chunk.
    pub fn end(&self) -> VirtualOffset {
        self.end
    }
}

impl Debug for Chunk {
    fn fmt(&self, f: &mut Formatter) -> result::Result<(), fmt::Error> {
        write!(f, "{{{}__{}}}", self.start, self.end)
    }
}

impl Display for Chunk {
    fn fmt(&self, f: &mut Formatter) -> result::Result<(), fmt::Error> {
        write!(f, "{{{}__{}}}", self.start, self.end)
    }
}

/// Single bin that stores chunks of
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Bin {
    bin_id: u32,
    chunks: Vec<Chunk>, // chunks for each sample, each format type
}

impl Bin {
    pub fn new(bin_id: u32, chunks: Vec<Chunk>) -> Bin {
        Bin {
            bin_id: bin_id,
            chunks: chunks,
        }
    }
    pub fn dummy() -> Bin {
        Bin {
            bin_id: 0,
            chunks: vec![],
        }
    }
    pub fn from_stream<R: Read>(stream: &mut R) -> Result<Self> {
        let bin_id = stream.read_u32::<LittleEndian>()?;
        let n_chunks = stream.read_i32::<LittleEndian>()? as usize;
        // let check_chunks = bin_id != SUMMARY_BIN;
        let check_chunks = false;
        let mut chunks = Vec::with_capacity(n_chunks);
        for i in 0..n_chunks {
            chunks.push(Chunk::from_stream(stream, check_chunks)?);
            if check_chunks && i > 0 && chunks[i].start() < chunks[i - 1].end() {
                return Err(Error::new(
                    InvalidData,
                    format!("Invalid index: chunks are not sorted for bin {}", bin_id),
                ));
            }
        }
        Ok(Bin { bin_id, chunks })
    }

    pub fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        stream.write_u32::<LittleEndian>(self.bin_id)?;
        let chunks = self.merge_chunks();
        /*if chunks.len() > 1 {
            eprintln!("{:?} {:?}", self.chunks(), chunks);
        }*/
        let n_chunks = chunks.len() as i32;
        stream.write_i32::<LittleEndian>(n_chunks)?;

        for chunk in chunks {
            chunk.to_stream(stream)?;
        }
        Ok(())
    }

    /// Returns the bin ID.
    pub fn bin_id(&self) -> u32 {
        self.bin_id
    }

    /// Returns all the chunks in the bin.
    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }

    fn merge_chunks(&self) -> Vec<Chunk> {
        self.chunks.iter().fold(vec![], |mut acc: Vec<Chunk>, x| {
            if let Some(last) = acc.last_mut() {
                if last.is_adjacent(x) {
                    last.concat(x);
                    acc
                } else {
                    acc.push(*x);
                    acc
                }
            } else {
                vec![*x]
            }
        })
    }
}

impl Display for Bin {
    fn fmt(&self, f: &mut Formatter) -> result::Result<(), fmt::Error> {
        write!(f, "Bin {}:  ", self.bin_id)?;
        for (i, chunk) in self.chunks.iter().enumerate() {
            write!(f, "{}{}", if i > 0 { ",  " } else { "" }, chunk)?;
        }
        Ok(())
    }
}

impl Debug for Bin {
    fn fmt(&self, f: &mut Formatter) -> result::Result<(), fmt::Error> {
        write!(f, "Bin {}:  ", self.bin_id)?;
        for (i, chunk) in self.chunks.iter().enumerate() {
            write!(f, "{}{}", if i > 0 { ",  " } else { "" }, chunk)?;
        }
        Ok(())
    }
}

/// Index for a single reference sequence. Contains [bins](struct.Bin.html).
#[derive(Clone, PartialEq, Eq)]
pub struct Reference {
    bins: BTreeMap<u32, Bin>, // bin_id is not continuous.
                              // Should we set `bins` as Btree-map?
                              // linear_index: LinearIndex,
}

impl Reference {
    pub fn new(bins: BTreeMap<u32, Bin>) -> Reference {
        return Reference { bins: bins };
    }

    fn from_stream<R: Read>(stream: &mut R) -> Result<Self> {
        let n_bins = stream.read_i32::<LittleEndian>()? as usize;
        let mut bins = BTreeMap::new(); //with_capacity(n_bins);
        for _ in 0..n_bins {
            let bin = Bin::from_stream(stream)?;
            bins.insert(bin.bin_id, bin);
        }

        //let linear_index = LinearIndex::from_stream(stream)?;
        //Ok(Reference { bins, linear_index })
        Ok(Reference { bins })
    }

    /// Returns all bins for the reference.
    pub fn bins(&self) -> &BTreeMap<u32, Bin> {
        &self.bins
    }

    fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        let n_ref = self.bins.len() as i32;
        stream.write_i32::<LittleEndian>(n_ref)?;

        for (_id, bin) in &self.bins {
            bin.to_stream(stream)?;
        }
        Ok(())
    }
}

impl Display for Reference {
    fn fmt(&self, f: &mut Formatter) -> result::Result<(), fmt::Error> {
        if !self.bins.is_empty() {
            writeln!(f, "    Bins:")?;
            for bin in self.bins.values() {
                writeln!(f, "        {}", bin)?;
            }
        }
        Ok(())
    }
}

#[derive(Clone, PartialEq, Eq)]
pub struct Index {
    // samples: Vec<String>, // sample names, we expect sample id is continuous.
    references: Vec<Reference>, // the length is inherited from n_ref
}

impl Index {
    pub fn new(references: Vec<Reference>) -> Index {
        Index {
            references: references,
        }
    }
    /// Loads index from stream.
    pub fn from_stream<R: Read>(mut stream: R) -> Result<Index> {
        let mut magic = [0_u8; 4];
        stream.read_exact(&mut magic)?;
        if magic != [b'G', b'H', b'I', 1] {
            return Err(Error::new(InvalidData, "Input is not in GHI format"));
        }

        let n_ref = stream.read_i32::<LittleEndian>()? as usize;
        let mut references = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            references.push(Reference::from_stream(&mut stream)?);
        }
        //        let n_unmapped = stream.read_u64::<LittleEndian>().ok();
        //        Ok(Index { references, n_unmapped })
        Ok(Index { references })
    }

    /// Loads index from path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Index> {
        let f = File::open(&path)?;
        let r = BufReader::with_capacity(262144, f);
        Index::from_stream(r)
    }

    pub fn to_stream<W: Write>(&self, mut stream: W) -> Result<()> {
        stream.write_all(&[b'G', b'H', b'I', 1])?;
        let n_ref = self.references.len() as i32;
        stream.write_i32::<LittleEndian>(n_ref)?;

        for reference in &self.references {
            reference.to_stream(&mut stream)?;
        }
        Ok(())
    }

    pub fn to_path<P: AsRef<Path> + Write>(&self, path: P) -> Result<()> {
        let f = File::open(&path)?;
        let b = BufWriter::with_capacity(262144, f);
        Index::to_stream(&self, b)
    }

    /// Returns all [references](struct.Reference.html) present in the BAI index.
    pub fn references(&self) -> &[Reference] {
        &self.references
    }

    /// Fetches [chunks](struct.Chunk.html) of the BAM file that contain all records for a given region.
    pub fn fetch_chunks(&self, ref_id: u64, start: u64, end: u64) -> Vec<Chunk> {
        let mut chunks = Vec::new();
        let ref_id = ref_id as usize;

        for bin_id in region_to_bins(start, end) {
            if let Some(bin) = self.references[ref_id].bins.get(&bin_id) {
                chunks.extend(bin.chunks.iter().copied());
            }
        }
        chunks.sort();
        chunks
        /*
        let mut res = Vec::new();
        for i in 0..chunks.len() {
            res.push(chunks[i].clone());
        }
        res*/
        //chunks
        /*
        let mut res = Vec::new();
        if chunks.is_empty() {
            return res;
        }
        chunks

        chunks.sort();
        let mut curr = chunks[0].clone();
        for i in 1..chunks.len() {
            if !curr.can_combine(&chunks[i]) {
                res.push(curr);
                curr = chunks[i].clone();
            } else {
                curr = curr.combine(&chunks[i]);
            }
        }
        res.push(curr);
        res*/
    }
}

impl Display for Index {
    fn fmt(&self, f: &mut Formatter) -> result::Result<(), fmt::Error> {
        for (i, reference) in self.references.iter().enumerate() {
            writeln!(f, "Reference {}:", i)?;
            reference.fmt(f)?;
        }
        Ok(())
    }
}

/// Returns a bin for the record with alignment `[beg-end)` for multiplex lower-level
/// Fixed version
pub fn region_to_bin_3(beg: u64, end: u64) -> u32 {
    let end = end - 1;
    let mut res = 0_i32;

    if beg >> 14 == end >> 14 {
        return (((1 << 29 - 14) - 1) / 7 + ((beg >> 14) << 1)) as u32;
    }
    if (beg + (1 << 13)) >> 14 == (end + (1 << 13)) >> 14 {
        return (((1 << 29 - 14) - 1) / 7 + ((beg >> 14) << 1) + 1) as u32;
    }

    for i in (17..27).step_by(3) {
        if beg >> i == end >> i {
            res = ((1 << 29 - i) - 1) / 7 + (beg >> i) as i32;
            break;
        }
    }
    res as u32
}

/// Returns all possible BAI bins for the region `[beg-end)`.
pub fn region_to_bins(start: u64, end: u64) -> BinsIter {
    BinsIter::new(-1, 0, start, end, 0, 0)
    /*    BinsIter {
        i: -1,
        t: 0,
        start,
        end,
        curr_bin: 0,
        bins_end: 0,
    }*/
}

/// Iterator over bins.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct BinsIter {
    i: i32,
    t: i32,
    start: i32,
    end: i32,
    curr_bin: u32,
    bins_end: u32,
}

impl BinsIter {
    fn new(i: i32, t: i32, start: u64, end: u64, curr_bin: u32, bins_end: u32) -> BinsIter {
        return {
            let start = start as i32;
            let end = end as i32;
            BinsIter {
                i,
                t,
                start,
                end,
                curr_bin,
                bins_end,
            }
        };
    }
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
            if self.i == 4 {
                // when the last line
                self.curr_bin = (self.t + (self.start >> 26 - 3 * self.i - 1)) as u32 - 2;
                if self.curr_bin < 4680 {
                    self.curr_bin = 4680;
                }
                self.bins_end = (self.t + (self.end >> 26 - 3 * self.i - 1)) as u32;
            } else {
                self.curr_bin = (self.t + (self.start >> 26 - 3 * self.i)) as u32 - 1;
                self.bins_end = (self.t + (self.end >> 26 - 3 * self.i)) as u32;
            }
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
pub const MAX_BIN: u32 = 70217;

/// Returns a maximal region for a given bin.
pub fn bin_to_region(bin: u32) -> (i32, i32) {
    if bin == 0 {
        return (std::i32::MIN, std::i32::MAX);
    }
    let mut left = 1;
    for i in 1..5 {
        let right = left + (1 << 3 * i);
        if bin >= left && bin < right {
            let beg = (bin - left) as i32;
            return (beg << 29 - 3 * i, beg + 1 << 29 - 3 * i);
        }
        left = right;
    }
    let i = 5;
    let right = left + (1 << 3 * i);
    if bin >= left && bin < right {
        let beg = ((bin - left) / 2u32) as i32;
        if bin % 2 == 0 {
            return (
                (1 << 13) + (beg << 29 - 3 * i),
                (1 << 13) + (beg + 1 << 29 - 3 * i),
            );
        } else {
            return (beg << 29 - 3 * i, beg + 1 << 29 - 3 * i);
        }
    }
    panic!(
        "Bin id should be not bigger than MAX_BIN ({} > {})",
        bin, MAX_BIN
    );
}

#[cfg(test)]
mod tests {
    use crate::index::bin_to_region;
    use crate::index::region_to_bins;
    use crate::index::BinsIter;

    #[test]
    fn bin_to_region_works() {
        assert_eq!(bin_to_region(1), (0, 1 << 26));
        assert_eq!(bin_to_region(9), (0, 1 << 23));
        assert_eq!(bin_to_region(73), (0, 1 << 20));
        assert_eq!(bin_to_region(585), (0, 1 << 17));
        assert_eq!(bin_to_region(4681), (0, 16384));
        assert_eq!(bin_to_region(4683), (16384, 16384 * 2));
        assert_eq!(bin_to_region(4682), (8192, 16384 + 8192));
    }

    #[test]
    fn bin_iter_4681() {
        let mut bin_iter = region_to_bins(0, 8191);
        assert_eq!(bin_iter, BinsIter::new(-1, 0, 0, 8191, 0, 0));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 0, 8191, 0, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 0, 8191, 1, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(1, 9, 0, 8191, 9, 9));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(2, 73, 0, 8191, 73, 73));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(3, 585, 0, 8191, 585, 585));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(4, 4681, 0, 8191, 4681, 4681));
        assert_eq!(bin_iter.next(), None);
    }

    #[test]
    fn bin_iter_4681_4683() {
        let mut bin_iter = region_to_bins(0, 16384);
        assert_eq!(bin_iter, BinsIter::new(-1, 0, 0, 16384, 0, 0));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 0, 16384, 0, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 0, 16384, 1, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(1, 9, 0, 16384, 9, 9));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(2, 73, 0, 16384, 73, 73));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(3, 585, 0, 16384, 585, 585));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(4, 4681, 0, 16384, 4681, 4683));
        assert_eq!(bin_iter.next(), Some(4682));
    }

    #[test]
    fn bin_iter_4682_4683() {
        let mut bin_iter = region_to_bins(16485, 16486);
        assert_eq!(bin_iter, BinsIter::new(-1, 0, 16485, 16486, 0, 0));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 16485, 16486, 0, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 16485, 16486, 1, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(1, 9, 16485, 16486, 9, 9));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(2, 73, 16485, 16486, 73, 73));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(3, 585, 16485, 16486, 585, 585));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(4, 4681, 16485, 16486, 4682, 4683));
        assert_eq!(bin_iter.next(), Some(4683));
    }

    #[test]
    fn bin_iter_4681_4683_2() {
        let mut bin_iter = region_to_bins(0, 16486);
        assert_eq!(bin_iter, BinsIter::new(-1, 0, 0, 16486, 0, 0));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 0, 16486, 0, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 0, 16486, 1, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(1, 9, 0, 16486, 9, 9));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(2, 73, 0, 16486, 73, 73));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(3, 585, 0, 16486, 585, 585));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(4, 4681, 0, 16486, 4681, 4683));
        assert_eq!(bin_iter.next(), Some(4682));
    }

    #[test]
    fn bin_iter_4681_4683_3() {
        let mut bin_iter = region_to_bins(8199, 16384);
        assert_eq!(bin_iter, BinsIter::new(-1, 0, 8199, 16384, 0, 0));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 8199, 16384, 0, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 8199, 16384, 1, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(1, 9, 8199, 16384, 9, 9));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(2, 73, 8199, 16384, 73, 73));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(3, 585, 8199, 16384, 585, 585));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(4, 4681, 8199, 16384, 4681, 4683));
        assert_eq!(bin_iter.next(), Some(4682));
    }

    #[test]
    fn bin_iter_4683_4684() {
        let mut bin_iter = region_to_bins(24578, 24589);
        assert_eq!(bin_iter, BinsIter::new(-1, 0, 24578, 24589, 0, 0));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 24578, 24589, 0, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(0, 1, 24578, 24589, 1, 1));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(1, 9, 24578, 24589, 9, 9));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(2, 73, 24578, 24589, 73, 73));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(3, 585, 24578, 24589, 585, 585));
        bin_iter.next();
        assert_eq!(bin_iter, BinsIter::new(4, 4681, 24578, 24589, 4683, 4684));
        assert_eq!(bin_iter.next(), Some(4684));
    }
}
