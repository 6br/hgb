use std::fs::File;
use crate::index::{Bin, Chunk, Region};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use std::io::{Write, Read, Result, Error};
use std::io::ErrorKind::InvalidData;
use std::mem::MaybeUninit;
use std::fmt;
use std::{path::Path, fmt::{Display, Formatter}, result};
/*
fn get_file_as_byte_vec(filename: &String) -> Vec<u8> {
    let mut f = File::open(&filename).expect("no file found");
    let metadata = fs::metadata(&filename).expect("unable to read metadata");
    let mut buffer = vec![0; metadata.len() as usize];
    f.read(&mut buffer).expect("buffer overflow");

    buffer
}
*/
#[derive(Clone)]
pub struct Reference {
    bin_count_mask: u64,
    bin_pitch_indices: [u8;64],
    bins: Vec<Bin> // or bins: Vec<Chunk>,
}

impl Reference {
    pub fn new(bin_count_mask: u64, bin_pitch_indices: [u8;64], bins: Vec<Bin>) -> Reference {
        Reference {
            bin_count_mask,
            bin_pitch_indices,
            bins
        }
    }

    /// The reference length is at most 512Mbp.
    pub fn new_with_bai_half_overlapping() -> Reference {
        let mut bin_pitch_indices = [0;64];
        bin_pitch_indices[63] = 29;
        bin_pitch_indices[62] = 26;
        bin_pitch_indices[61] = 23;
        bin_pitch_indices[60] = 20;
        bin_pitch_indices[59] = 17;
        bin_pitch_indices[58] = 14;
        Reference {
            bin_count_mask: 0b10010010010010001,
            bin_pitch_indices, //: [0,..,14,17,20,23,26,29],
            bins: vec![]
        }
    }

    fn from_stream<R: Read>(stream: &mut R) -> Result<Self> {
        // https://tyfkda.github.io/blog/2020/03/19/rust-init-array.html
        let bin_count_mask = stream.read_u64::<LittleEndian>()?;
        const LEN: usize = 64;

        let mut array: [MaybeUninit<u8>; LEN] = unsafe { MaybeUninit::uninit().assume_init() };
        for (i, slot) in array.iter_mut().enumerate() {
            *slot = MaybeUninit::new((i * 2) as u8);
        }
        let mut array: [u8; LEN] = unsafe { std::mem::transmute::<_, [u8; LEN]>(array) };

        for i in 0..63 {
            array[i] = stream.read_u8()?;
        }

        // Just copy and paste from old Bins
        let n_chunks = stream.read_u64::<LittleEndian>()? as usize;
        // let check_chunks = true; // bin_id != SUMMARY_BIN;
        let mut chunks = Vec::with_capacity(n_chunks);
        for _i in 0..n_chunks {
            chunks.push(Bin::from_stream(stream)?);
            /*if check_chunks && i > 0 && chunks[i].start() < chunks[i - 1].end() {
                return Err(Error::new(InvalidData, format!("Invalid index: chunks are not sorted")));
            }*/
        }

        Ok(Reference{bin_count_mask, bin_pitch_indices: array, bins: chunks})
    }

    fn to_stream<W: Write>(&self, stream: &mut W) -> Result<()> {
        stream.write_u64::<LittleEndian>(self.bin_count_mask)?;
        for i in self.bin_pitch_indices.iter() {
            stream.write_u8(*i)?;
        }
        //stream.write_all(self.bins.as_ref());
        let n_chunks = self.bins.len() as i32;
        stream.write_i32::<LittleEndian>(n_chunks)?;

        for chunk in &self.bins {
            chunk.to_stream(stream)?;
        }

        Ok(())
    }

    pub fn region_to_bin(&self, range: Region) -> usize {
        let mut iter = self.region_to_bins(range);
        let mut prev = iter.next().unwrap();
        loop {
            let next = iter.next();
            match next {
                None => { return prev.bin_disp_start }
                Some(slice) => {
                    if slice.bin_size > 1 {
                        return slice.bin_disp_start
                    } else {
                        prev = slice;
                    }
                }
            }
        }
    }

    pub fn region_to_bins(&self, range: Region) -> BinsIter {
        BinsIter{
            range,
            finished: 0,
            header: self
        }
    }
}

impl Display for Reference {
    fn fmt(&self, f: &mut Formatter) -> result::Result<(), fmt::Error> {
        if !self.bins.is_empty() {
            writeln!(f, "    Bins:")?;
            for bin in &self.bins {
                writeln!(f, "        {}", bin)?;
            }
        }
        Ok(())
    }
}


#[derive(Clone)]
pub struct Index {
    // samples: Vec<String>, // sample names, we expect sample id is continuous.
    references: Vec<Reference>, // the length is inherited from n_ref
}

impl Index {
    pub fn new(references: Vec<Reference>) -> Index {
        Index{references: references}
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
        Index::from_stream(f)
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
        let _f = File::open(&path)?;
        Index::to_stream(&self, path)
    }

    /// Returns all [references](struct.Reference.html) present in the BAI index.
    pub fn references(&self) -> &[Reference] {
        &self.references
    }
    /// Fetches [chunks](struct.Chunk.html) of the BAM file that contain all records for a given region.
    ///pub fn fetch_chunks(&self, ref_id: u64, start: u64, end: u64) -> Vec<Chunk> {
        pub fn fetch_chunks(&self, region: Region) -> Vec<Chunk> {
        let mut chunks = Vec::new();
        let ref_id = region.ref_id() as usize;
        let ref_container = &self.references[ref_id];

        for slice in ref_container.region_to_bins(region) {
            for bin in slice.slice {
//            if let Some(bin) = self.references[ref_id].bins.get(&bin_id) {
                chunks.extend(bin.chunks().iter());
            }
        }
        chunks.sort();

        let mut res = Vec::new();
        for  i in 0..chunks.len() {
            res.push(chunks[i].clone());
        }
        res
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

pub struct BinsIter<'a> {
    range: Region,
    finished: u64,
    header: &'a Reference,
}

pub struct Slice<'a> {
    slice: &'a [Bin], // or Chunk
    bin_size: usize,
    bin_disp_start: usize, // For debugging and 
    range: Region,
}

impl<'a> Iterator for BinsIter<'a> {
    type Item = Slice<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        /* load constant and previous state */
        let finished = self.finished;
        let everything = self.header.bin_count_mask;

        /* mask finished bits out to compute remaining bits */
        let remaining = everything & (!finished);

        /* if there is no remaining bits, iteration is done */
        if remaining == 0 {
            return None
        }
        
        /* compute bin offset; offset base is the sum of the number of finished bins */
        let bin_ofs_base = (everything & finished) as usize;
        
        /* determine bin pitch and size; note that bin span (size) is twice as large as bin pitch, because bins are half-overlapping */
        let iterations_done;
        unsafe {
            iterations_done = core::arch::x86_64::_popcnt64(bin_ofs_base as i64) as usize;
        }
        let bin_pitch_index = self.header.bin_pitch_indices[iterations_done];
        let bin_pitch = (0x01u64<<bin_pitch_index) as usize;
        let bin_size = 2 * bin_pitch;

        let bin_ofs_disp_start = {
            let mut t = self.range.start() >> bin_pitch_index;
            if t > 0 { t -= 1; }
            t
        };

        let bin_ofs_disp_end = {
            let mut t = self.range.end() >> bin_pitch_index;
            if t > finished { t = finished; }
            /* previous finished (not new one) is <bin count for this depth> - 1 */
            t + 1									/* make the range margined (tail margin) */
        };

        /* compute bin range */
        let bin_range_start = bin_ofs_disp_start << bin_pitch_index;
        let bin_range_end = (bin_ofs_disp_end + 1) << bin_pitch_index;
        let bin_disp_start: usize = bin_ofs_base + bin_ofs_disp_start as usize;
        let bin_disp_end: usize = bin_ofs_disp_end as usize + bin_disp_start;

        /* update finished mask; find least significant set bit to locate the next unfinished bit, then xor to squish further remaining bits out */
        self.finished = remaining ^ (remaining - 1);

        /* done; compute bin pointer */
        Some(Slice{
            slice: &self.header.bins[bin_disp_start..bin_disp_end],
            bin_size,
            range: Region::new(self.range.ref_id(), bin_range_start, bin_range_end),
            bin_disp_start, //Region{start: bin_range_start, end: bin_range_end, ref_id: self.ref_id},
        })
    }
}

mod tests {
    use super::Reference;
    use crate::index::Region;

    #[test]
    fn iterator_works() {
        let bai = Reference::new_with_bai_half_overlapping();
        assert_eq!(bai.bin_pitch_indices[63], 29);
        let mut iter = bai.region_to_bins(Region::new(0,0,8192));
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2);
        assert_eq!(slice.range, Region::new(0,0,0));
        assert_eq!(slice.bin_disp_start, 0);

    }

    #[test]
    fn region_to_bin_works() {
        let bai = Reference::new_with_bai_half_overlapping();
        let bin = bai.region_to_bin(Region::new(0, 0, 100_000_000));
        assert_eq!(bin, 0);
        let bin = bai.region_to_bin(Region::new(0, 0, 58_000_000));
        assert_eq!(bin, 1);
        let bin = bai.region_to_bin(Region::new(0, 58_000_000, 112_000_000));
        assert_eq!(bin, 2);
    }
}