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

fn f(color: &VisColor) -> usize {
    *color as usize
}

pub enum ColorSet {
    HgbColor(HgbColor),
    IgvColor(IgvColor),
    JBrowseColor(JBrowseColor),
    //    UCSC_COLOR,
}

impl FromStr for ColorSet {
    type Err = &'static str;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "hgb" => Ok(ColorSet::HgbColor(HgbColor)),
            "igv" => Ok(ColorSet::IgvColor(IgvColor)),
            "jbrowse" => Ok(ColorSet::JBrowseColor(JBrowseColor)),
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
            ColorSet::JBrowseColor(_) => JBrowseColor::p(idx).rgb(),
        };
        RGBColor(r, g, b)
    }
}

impl Default for ColorSet {
    fn default() -> Self {
        Self::new()
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
pub struct JBrowseColor;

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

impl Palette for JBrowseColor {
    const COLORS: &'static [(u8, u8, u8)] = &[
        (76, 175, 80), //#4caf50
        (33, 150, 243), //#2196f3
        (255, 193, 7), //#ffc107
        (238, 34, 34), //#ee2222
        (170, 170, 170), //#aaaaaa
        (153, 153, 153), //#999999
        (153, 0, 153), //#in
        (120, 85, 43), //#del
        (236, 139, 139), //#EC8B8B 
        (143, 143, 216), //#8F8FD8
        (120, 85, 43), //#Spl
    ];
}

impl Palette2 for JBrowseColor {}

//    ACol,
//CCol,
//GCol,
//TCol,
//NCol,
//DefCol,
//InsCol,
//DelCol,
//PosCol,
//NegCol,
//SplCol,