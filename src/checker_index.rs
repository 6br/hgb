use std::fs::File;
use crate::index::{Bin, Chunk, Region};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use std::io::{Write, Read, Result, Error};
use std::io::ErrorKind::InvalidData;
use std::mem::MaybeUninit;
use std::fmt;
use std::{path::Path, fmt::{Display, Formatter}, result, ops::Range};
use bam::header::HeaderEntry;
use fmt::Debug;
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

impl Debug for Reference {
    fn fmt(&self, f: &mut Formatter) -> result::Result<(), fmt::Error> {
        write!(f, "Mask: {}", self.bin_count_mask)
    }
}

impl Reference {
    pub fn new(bin_count_mask: u64, bin_pitch_indices: [u8;64], bins: Vec<Bin>) -> Reference {
        Reference {
            bin_count_mask,
            bin_pitch_indices,
            bins
        }
    }

    pub fn new_from_len(len: u64) -> Reference {
        if len > 2_u64.pow(32) {
            panic!("Unsupported reference size");
        } else {
            return Reference::new_with_bai_half_overlapping()
        }
    }

    pub fn new_from_reference(reference: &HeaderEntry) -> Reference {
        if let Some(len) = reference.get(b"LN") {
            if let Some(len) = len.parse::<u64>().ok() {
                return Reference::new_from_len(len)
            }
        }
        return Reference::new_with_bai_half_overlapping()
    }

    /// The reference length is at most 512Mbp.
    pub fn new_with_bai_half_overlapping() -> Reference {
        let mut bin_pitch_indices = [0;64];
        bin_pitch_indices[0] = 28;
        bin_pitch_indices[1] = 25;
        bin_pitch_indices[2] = 22;
        bin_pitch_indices[3] = 19;
        bin_pitch_indices[4] = 16;
        bin_pitch_indices[5] = 13;
        Reference {
            bin_count_mask: 0b10010010010010001,
            bin_pitch_indices, //: [0,..,14,17,20,23,26,29],
            bins: vec![]
        }
    }

    pub fn bins(&mut self) -> &Vec<Bin> {
        return &self.bins
    }

    pub fn update(&mut self, bin_id: usize, bin: Bin) {
        self.bins[bin_id] = bin;
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
        let n_chunks = self.bins.len() as u64;
        stream.write_u64::<LittleEndian>(n_chunks)?;

        for chunk in &self.bins {
            chunk.to_stream(stream)?;
        }

        Ok(())
    }
    /*
    pub fn region_to_bin(&self, range: Region) -> usize {
        let range_middle = (range.start() + range.end()) / 2;
        let middle = Region::new(range.ref_id(), range_middle, range_middle + 1);
        println!("Range: {:?} {:?}", range, range_middle);
        let mut iter = self.region_to_bins(middle);
        let mut prev = 0_usize; //iter.next().unwrap();
        loop {
            let next = iter.next();
            match next {
                None => { return prev }
                Some(slice) => {
                    println!("{:?} {:?} {}", slice.bin_disp_range, slice.range,  (slice.bin_size as u64));
                    println!("{} {}", !slice.range.include(&range), (slice.bin_size as u64) < range.len() );
                    if !slice.range.include(&range) || (slice.bin_size as u64) < range.len() {
                        return prev
                    } else {
                        prev = slice.bin_disp_range.start;
                    }
                }
            }
        }
    }
*/
    pub fn region_to_bin(&self, range: Region) -> usize {
        // println!("Range: {:?}", range);
        let len = range.len();
        let start = range.start();
        let mut iter = self.region_to_bins(range);
        // let mut prev = iter.next().unwrap();
        let mut prev = 0;
        loop {
            let next = iter.next();
            match next {
                None => { return prev } //v.bin_disp_range.start }
                Some(slice) => {
                    // println!("{:?} {:?} {}", slice.bin_disp_range, slice.range,  (slice.bin_size as u64));
                    if (slice.bin_size as u64) < len {
                        return prev //.bin_disp_range.start
                    } else {
                        prev = slice.bin_disp_range.start;
                        let mut pos = slice.range.start();
                        while pos < start {
                            prev += 1;
                            pos += slice.bin_size as u64;
                        }
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
    pub fn fetch_chunks(&self, ref_id: u64, start: u64, end: u64) -> Vec<Chunk> {
    //pub fn fetch_chunks(&self, region: Region) -> Vec<Chunk> {
        let mut chunks: Vec<&Chunk> = Vec::new();
        //let ref_id = region.ref_id() as usize;
        let region = Region::new(ref_id, start, end);
        let ref_id = ref_id as usize;
        let ref_container = &self.references[ref_id];
        // println!("{:?}", region);
        for slice in ref_container.region_to_bins(region) {
            if let Some(bin) = slice.slice {
                for chunk in bin {
                    chunks.extend(chunk.chunks());
                }
            }
        }
        chunks.sort();

        let mut res = Vec::new();
        for i in 0..chunks.len() {
            res.push(chunks[i].clone());
        }
        res
    }
}


impl Display for Index {
    fn fmt(&self, f: &mut Formatter) -> result::Result<(), fmt::Error> {
        for (i, reference) in self.references.iter().enumerate() {
            writeln!(f, "Reference {}:", i)?;
            //reference.fmt::<dyn Display>(f)?;
            Display::fmt(&reference, f)?
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
    slice: Option<&'a [Bin]>, // or Chunk
    bin_size: usize,
    //bin_disp_start: usize, // For debugging and 
    bin_disp_range: Range<usize>,
    range: Region, // The entire range the slice covers.
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
            //if t > finished { t = finished; }
            if t > remaining ^ (remaining - 1) { t = remaining ^ (remaining - 1); }
            /* previous finished (not new one) is <bin count for this depth> - 1 */
            t + 1									/* make the range margined (tail margin) */
        };

        /* compute bin range */
        let bin_range_start = bin_ofs_disp_start << bin_pitch_index;
        //let bin_range_end = (bin_ofs_disp_start + bin_ofs_disp_end + 1) << bin_pitch_index;
        let bin_range_end = (bin_ofs_disp_end + 1) << bin_pitch_index;
        let bin_disp_start: usize = bin_ofs_base + bin_ofs_disp_start as usize;
        let bin_disp_end: usize = bin_ofs_disp_end as usize + bin_ofs_base;

        /* update finished mask; find least significant set bit to locate the next unfinished bit, then xor to squish further remaining bits out */
        self.finished = remaining ^ (remaining - 1);

        println!("{:?} {} {} {} {} {:?} {:?} {} {}", self.range, bin_range_start, bin_range_end, bin_size, finished, bin_disp_start..bin_disp_end, bin_ofs_disp_start..bin_ofs_disp_end, iterations_done, bin_pitch_index);

        /* done; compute bin pointer */
        Some(Slice{
            //slice: &self.header.bins[bin_disp_start..bin_disp_end],
            slice: self.header.bins.get(bin_disp_start..bin_disp_end),
            bin_size,
            range: Region::new(self.range.ref_id(), bin_range_start, bin_range_end),
            bin_disp_range: bin_disp_start..bin_disp_end, //Region{start: bin_range_start, end: bin_range_end, ref_id: self.ref_id},
        })
    }
}

#[cfg(test)]
mod tests {
    use super::Reference;
    use crate::index::Region;

    #[test]
    fn iterator_works() {
        let bai = Reference::new_with_bai_half_overlapping();
        //assert_eq!(bai.bin_pitch_indices[0], 28);
        let mut iter = bai.region_to_bins(Region::new(0,0,8191));
        let slice = iter.next().unwrap();
        /* Because we adopt half-overlapping, the most  */
        assert_eq!(slice.bin_size, 2_usize.pow(29));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(29)));
        assert_eq!(slice.bin_disp_range, 0..1); //1 is not included.
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(26));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(26)));
        assert_eq!(slice.bin_disp_range, 1..2); //2 is not included.
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(23));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(23)));
        assert_eq!(slice.bin_disp_range, 17..18); //18 is not included.
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(20));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(20)));
        assert_eq!(slice.bin_disp_range, 145..146); //2 is not included
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(17));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(17)));
        assert_eq!(slice.bin_disp_range, 1169..1170); //2 is not included
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(14));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(14)));
        assert_eq!(slice.bin_disp_range, 9361..9362); //2 is not included
    }
    #[test]
    fn iterator_2_works() {
        let bai = Reference::new_with_bai_half_overlapping();
        //assert_eq!(bai.bin_pitch_indices[0], 28);
        let mut iter = bai.region_to_bins(Region::new(0,0,8192));
        let slice = iter.next().unwrap();
        /* Because we adopt half-overlapping, the most  */
        assert_eq!(slice.bin_size, 2_usize.pow(29));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(29)));
        assert_eq!(slice.bin_disp_range, 0..1); //1 is not included.
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(26));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(26)));
        assert_eq!(slice.bin_disp_range, 1..2); //2 is not included.
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(23));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(23)));
        assert_eq!(slice.bin_disp_range, 17..18); //18 is not included.
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(20));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(20)));
        assert_eq!(slice.bin_disp_range, 145..146); //2 is not included
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(17));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(17)));
        assert_eq!(slice.bin_disp_range, 1169..1170); //2 is not included
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(14));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(14) + 2_u64.pow(13)));
        assert_eq!(slice.bin_disp_range, 9361..9363); //2 is not included
    }
    #[test]
    fn iterator3_works() {
        let bai = Reference::new_with_bai_half_overlapping();
        //assert_eq!(bai.bin_pitch_indices[0], 28);
        let mut iter = bai.region_to_bins(Region::new(0,17000,17500));
        let slice = iter.next().unwrap();
        /* Because we adopt half-overlapping, the most  */
        assert_eq!(slice.bin_size, 2_usize.pow(29));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(29)));
        assert_eq!(slice.bin_disp_range, 0..1); //1 is not included.
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(26));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(26)));
        assert_eq!(slice.bin_disp_range, 1..2); //2 is not included.
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(23));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(23)));
        assert_eq!(slice.bin_disp_range, 17..18); //18 is not included.
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(20));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(20)));
        assert_eq!(slice.bin_disp_range, 145..146); //2 is not included
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(17));
        assert_eq!(slice.range, Region::new(0,0,2_u64.pow(17)));
        assert_eq!(slice.bin_disp_range, 1169..1170); //2 is not included
        let slice = iter.next().unwrap();
        assert_eq!(slice.bin_size, 2_usize.pow(14));
        assert_eq!(slice.range, Region::new(0,2_u64.pow(13),2_u64.pow(14)+2_u64.pow(14)));
        assert_eq!(slice.bin_disp_range, 9362..9364); //2 is not included
    }
    #[test]
    fn region_to_bin_works() {
        let bai = Reference::new_with_bai_half_overlapping();
        let bin = bai.region_to_bin(Region::new(0, 0, 100_000_000));
        assert_eq!(bin, 0);
        let bin = bai.region_to_bin(Region::new(0, 0, 2_u64.pow(25)-1));
        assert_eq!(bin, 1);
        let bin = bai.region_to_bin(Region::new(0, 66_000_000, 112_000_000));
        assert_eq!(bin, 2);
        let bin = bai.region_to_bin(Region::new(0, 0, 1<<17));
        assert_eq!(bin, 1169);
        let bin = bai.region_to_bin(Region::new(0, 0, 16384));
        assert_eq!(bin, 9361);
        let bin = bai.region_to_bin(Region::new(0, 0, 8191));
        assert_eq!(bin, 9361);
        let bin = bai.region_to_bin(Region::new(0, 8192, 8192+16384));
        assert_eq!(bin, 9362);
        let bin = bai.region_to_bin(Region::new(0, 16384, 16384*2));
        assert_eq!(bin, 9363);
        let bin = bai.region_to_bin(Region::new(0, 16382, 16385));
        assert_eq!(bin, 9362);
    }
}