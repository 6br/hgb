use enum_map::Enum;
use plotters::style::Color;
use plotters::{prelude::RGBColor, style::Palette, style::PaletteColor};
use std::str::FromStr;

#[derive(Debug, Enum, Copy, Clone)]
pub enum VisColor {
    ACol,
    CCol,
    GCol,
    TCol,
    NCol,
    DefCol,
    InsCol,
    DelCol,
    PosCol,
    NegCol,
    SplCol,
}

fn f(foo: &VisColor) -> usize {
    *foo as usize
}

pub enum ColorSet {
    HgbColor(HgbColor),
    IgvColor(IgvColor),
    //    UCSC_COLOR,
}

impl FromStr for ColorSet {
    type Err = &'static str;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "hgb" => Ok(ColorSet::HgbColor(HgbColor)),
            "igv" => Ok(ColorSet::IgvColor(IgvColor)),
            "" => Ok(ColorSet::HgbColor(HgbColor)),
            _ => Ok(ColorSet::HgbColor(HgbColor)),
        }
    }
}

impl ColorSet {
    pub fn new() -> Self {
        ColorSet::HgbColor(HgbColor)
    }
    pub fn pick(&self, idx: VisColor) -> RGBColor {
        let (r, g, b) = match self {
            ColorSet::HgbColor(HgbColor) => HgbColor::p(idx).rgb(),
            ColorSet::IgvColor(IgvColor) => IgvColor::p(idx).rgb(),
        };
        RGBColor(r, g, b)
    }
}

pub trait Palette2: Palette {
    // const COLORS: &'static [(u8, u8, u8)];
    fn p(idx: VisColor) -> PaletteColor<Self>
    where
        Self: Sized,
    {
        PaletteColor::<Self>::pick(f(&idx))
    }
}

pub struct HgbColor;
pub struct IgvColor;

// Inspired by nordtheme. https://www.nordtheme.com
impl Palette for HgbColor {
    const COLORS: &'static [(u8, u8, u8)] = &[
        (119, 217, 168),
        (0, 90, 255),
        (128, 64, 0),
        (255, 75, 0),
        (80, 80, 80),
        (150, 150, 150),
        (153, 0, 153),
        (120, 85, 43),
        (0xbf, 0x61, 0x6a),
        (0x81, 0xa1, 0xc1),
        (0x8f, 0xbc, 0xbb),
    ];
}

impl Palette2 for HgbColor {}

impl Palette for IgvColor {
    const COLORS: &'static [(u8, u8, u8)] = &[
        (0, 200, 0),
        (0, 0, 200),
        (209, 113, 5),
        (255, 0, 40),
        (80, 80, 80),
        (150, 150, 150),
        (153, 0, 153),
        (120, 85, 43),
        (230, 150, 150),
        (150, 150, 230),
        (120, 85, 43),
    ];
}

impl Palette2 for IgvColor {}
