use rstar::{Envelope, Point, PointDistance, PointExt, RTree, RTreeObject, AABB};
use serde::{Deserialize, Serialize};
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
    pub track: u64,
    pub read_id: String,
    pub start: i32,
    pub end: i32,
    pub query_len: u32,
    pub strand: bool,
    pub flag: u16,
    pub mapq: u8,
    pub sa: String,
    pub cigar: String,
    pub insertions: Vec<(i32, u64, String)>, // Pixel, Real_pos, Ins_seq
}
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Annotation {
    pub rectangle: (i32, i32, i32, i32),
    pub start: u64,
    pub end: u64,
    pub name: String,
}

impl RTreeObject for Read {
    type Envelope = AABB<[i32; 2]>;

    fn envelope(&self) -> Self::Envelope {
        AABB::from_corners(
            [self.rectangle.0, self.rectangle.1],
            [self.rectangle.2, self.rectangle.3],
        )
    }
}

impl Read {
    /// Returns the nearest point within this rectangle to a given point.
    ///
    /// If `query_point` is contained within this rectangle, `query_point` is returned.
    pub fn nearest_point(&self, query_point: [i32; 2]) -> [i32; 2] {
        AABB::from_corners(
            [self.rectangle.0, self.rectangle.1],
            [self.rectangle.2, self.rectangle.3],
        )
        .min_point(&query_point)
    }
}

impl PointDistance for Read {
    fn distance_2(
        &self,
        point: &<Self::Envelope as Envelope>::Point,
    ) -> <<Self::Envelope as Envelope>::Point as Point>::Scalar {
        self.nearest_point(*point).sub(point).length_2()
    }

    fn contains_point(&self, point: &<Self::Envelope as Envelope>::Point) -> bool {
        AABB::from_corners(
            [self.rectangle.0, self.rectangle.1],
            [self.rectangle.2, self.rectangle.3],
        )
        .contains_point(point)
    }

    fn distance_2_if_less_or_equal(
        &self,
        point: &<Self::Envelope as Envelope>::Point,
        max_distance_2: <<Self::Envelope as Envelope>::Point as Point>::Scalar,
    ) -> Option<<<Self::Envelope as Envelope>::Point as Point>::Scalar> {
        let distance_2 = self.distance_2(point);
        if distance_2 <= max_distance_2 {
            Some(distance_2)
        } else {
            None
        }
    }
}

pub struct ReadTree(RTree<Read>);

impl ReadTree {
    pub fn new(reads: Area) -> Self {
        return ReadTree(RTree::bulk_load(reads.pileups));
    }
    pub fn find(self, x: i32, y: i32) -> Option<(Read, (u64, String))> {
        let read = self.0.locate_at_point(&[x, y]);
        match read {
            Some(read) => {
                let insertion = read.insertions.binary_search_by_key(&x, |(a, _, _)| *a);
                let insertions_str = match insertion {
                    Ok(index) => (read.insertions[index].1, read.insertions[index].2.clone()),
                    Err(index) => {
                        if index > 0 && (read.insertions[index - 1].0 - x).abs() <= 5 {
                            (
                                read.insertions[index - 1].1,
                                read.insertions[index - 1].2.clone(),
                            )
                        } else {
                            (0, "".to_string())
                        }
                    }
                };
                Some((read.clone(), insertions_str))
            }
            _ => None,
        }
    }
}
