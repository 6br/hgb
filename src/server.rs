
use actix_files::NamedFile;
use rand::Rng;
use actix_web::http::header::{ContentDisposition, DispositionType};
use actix_web::{HttpRequest, Result, web, Responder, error, middleware::Logger};
use std::{sync::{RwLock, Arc},  collections::BTreeMap, fs};
use clap::{App,  ArgMatches, Arg, AppSettings};
use ghi::{bed, vis::bam_record_vis};
// use crate::subcommands::bam_vis;
use genomic_range::StringRegion;
use bam::Record;

fn id_to_range<'a>(range: &StringRegion, args: &Vec<String>, zoom: u64, path: u64, param: &Param, path_string: String) -> (ArgMatches, StringRegion) {

    let app = App::new("vis")
            .setting(AppSettings::AllArgsOverrideSelf)
            .about("Visualize GHB and other genomic data")
            .arg(
                Arg::new("range")
                    .short('r')
                    .takes_value(true)
                    .multiple(true)
                    .about("Genomic range to visualize. Format is chr:from-to."),
            )
            .arg(
                Arg::new("prefetch-range")
                    .short('R')
                    .takes_value(true)
                    .multiple(true)
                    .about("Genomic range to pre-fetch. Optional. The same order as range option. Format is chr:from-to."),
            )
            .arg(
                Arg::new("type")
                    .short('t')
                    .takes_value(true)
                    // .default_value("default")
                    .possible_values(&["alignment", "range", "default"])
                    .about("The type of GHB data structure to display (used for debugging.)"),
            )
            .arg(
                Arg::new("preset")
                    .short('z')
                    //.takes_value(true)
                    // .empty_values(true)
                    // .default_value("auto")
                    //.setting(ArgSettings::AllowEmptyValues)
                    .possible_values(&["","auto", "base", "gene", "chrom", "sv", "qual"])
                    .about("Preset (always overwrites other options) ['auto', 'base', 'gene', 'chrom', 'sv', 'qual']"),
            )
            .arg(
                Arg::new("id")
                    .short('i')
                    .takes_value(true)
                    .multiple(true)
                    .about("Indices on GHB data structure to display (used for debugging.)"),
            )
            .arg(
                Arg::new("frequency")
                    .short('F')
                    .takes_value(true)
                    .multiple(true)
                    .about("A subset of sorted bed for coverage plot (start and score fields are used)"),
            )
            .arg(Arg::new("no-filter").short('f').about("Disable pre-filtering on loading BAM index (used for debugging)"))
            .arg(Arg::new("no-cigar").short('n').about("Do not show cigar string"))
            .arg(Arg::new("no-scale").short('S').about("Do not show y-axis scale"))
            .arg(Arg::new("no-packing").short('p').about("Disable read packing"))
            .arg(Arg::new("no-legend").short('l').about("Hide legend"))
            .arg(Arg::new("square").short('Q').about("Set square width (overwritten)"))
            .arg(Arg::new("quality").short('q').about("Display reads by quality value"))
            .arg(Arg::new("x").short('x').takes_value(true).about("The width of image"))
            .arg(Arg::new("y").short('y').takes_value(true).about("The height of each read alignment"))
            .arg(Arg::new("web").short('w').takes_value(true).about("Serve the web server"))
            .arg(
                Arg::new("freq-height")
                    .short('Y')
                    .takes_value(true)
                    .about("The height of each coverage track"),
            )
            .arg(
                Arg::new("max-coverage")
                    .short('m')
                    .takes_value(true)
                    .about("Max coverage value on coverage track"),
            )
            .arg(
                Arg::new("x-scale")
                    .short('X')
                    .takes_value(true)
                    .about("Size of x-scale legend"),
            )
            .arg(Arg::new("split-alignment").short('s').about("Display split alignments in the same "))
            .arg(Arg::new("only-split-alignment").short('u').about("Display only split alignments or mate-paired reads on alignment track"))
            .arg(
                Arg::new("sort-by-name")
                    .short('N')
                    .about("Sort alignments by read id (for split-alignment visualization)"),
            )
            .arg(
                Arg::new("sort-by-cigar")
                    .short('C')
                    .about("Display each split-alignment connection as a colored line"),
            )
            .arg(Arg::new("no-insertion").short('I').about("Hide insertion callets on read alignments"))
            .arg(
                Arg::new("hide-alignment")
                    .short('A')
                    .about("Hide read alignments (display only coverage; for chromosome scale view)"),
            )
            .arg(Arg::new("all-bases").short('B').about("Show all nucleotides by color"))
            .arg(
                Arg::new("pileup")
                    .short('P')
                    .about("Show pileup as coverage plot"),
            )
            .arg(Arg::new("only-translocation").short('T').about("Show callets on ends of read alignments if the read contains translocation split-alignment"))
            .arg(
                Arg::new("graph")
                    .short('c')
                    .takes_value(true)
                    .about("[Input] Graph genome coordinates tsv (tab-separated, generated by vg view -N command)"),
            )
            // .arg(Arg::new("insertion").short('').about("No-md"))
            .arg(
                Arg::new("output")
                    .short('o')
                    .takes_value(true)
                    .about("[Output] image file (prefixed as .bmp / .png)"),
            )
            .arg(
                Arg::new("bed")
                    .short('L')
                    .multiple(true)
                    .takes_value(true)
                    .about("[Input] A subset of sorted bed to display as annotation track"),
            )
            .arg(
                Arg::new("bam")
                    .short('a')
                    .takes_value(true)
                    .multiple(true)
                    .about("[Input] Sorted BAM to display read alignment track"),
            )
            .arg(
                Arg::new("gff3")
                    .short('g')
                    .takes_value(true)
                    .multiple(true)
                    .about("[Input] A subset of sorted gff3 to display as annotation track"),
            )
            .arg(
                Arg::new("INPUT")
                    .about("(Optional) GHB format to display both alignment and annotation tracks")
                    .index(1),
            );
    let prefetch_max = 25000000; // prefetch_range.end() - prefetch_range.start();
    let criteria = param.criteria; // range.end() - range.start();
    // let x_width = prefetch_max / criteria * (740); //max_x
    // Here 2**17 < x_width < 2**18 so maxZoom should be 18+ 1;
    let max_zoom = param.max_zoom as u64; // 18;
    let max_y = param.max_y as u64;
    let freq_y = param.y_freq as u64;
    let y = param.y as u64;
    let scalex_default = param.x_scale as u64;
    // let path = path << (max_zoom - zoom);
    let b: Vec<String> = if criteria << (max_zoom - zoom) >= 3000000 { // i.e. 2**22s
        vec!["-A".to_string(), "-Y".to_string(), (max_y >> (max_zoom-zoom)).to_string(), "-y".to_string(), (y >> (max_zoom-zoom)).to_string(), "-n".to_string(), "-I".to_string(), "-o".to_string(), path_string, "-X".to_string(), ((scalex_default >> (max_zoom-zoom) )* (max_y / freq_y)).to_string()]
    } else if criteria << (max_zoom - zoom ) <= 200000 { // Base-pair level
        vec!["-Y".to_string(), (freq_y >> (max_zoom-zoom)).to_string(), "-y".to_string(), (y >> (max_zoom-zoom)).to_string(), "-X".to_string(), (scalex_default >> (max_zoom-zoom)).to_string(), "-o".to_string(), path_string]
    } else {
        vec!["-Y".to_string(), (freq_y >> (max_zoom-zoom)).to_string(), "-y".to_string(), (y >> (max_zoom-zoom)).to_string(), "-X".to_string(), (scalex_default >> (max_zoom-zoom)).to_string(), "-n".to_string(), "-I".to_string(), "-o".to_string(), path_string]
    };
    let mut args = args.clone();
    args.extend(b); //b.extend(args.clone());
    args.remove(0);
    args = args.into_iter().skip_while(|t| t != "vis").collect();
    //b.insert(0, "vis".to_string());
    let start =  (criteria << (max_zoom - zoom)) * path + range.start;
    let range = StringRegion::new(&format!("{}:{}-{}", range.path, start, start + (criteria << (max_zoom - zoom)) - 1).to_string()).unwrap();
    eprintln!("{:?} {:?}", args, range);
    let matches = app.get_matches_from(args);
    // let args: Vec<String> = args.into_iter().chain(b.into_iter()).collect(); 
    (matches, range)
}

async fn get_dzi<'a>(data: web::Data<RwLock<Vis<'a>>>) -> impl Responder {
    return web::Json(data.read().unwrap().dzi.clone());
}

async fn index<'a>(data: web::Data<RwLock<Vis<'a>>>, list: web::Data<Vec<(u64, Record)>>, req: HttpRequest) -> Result<NamedFile> {
    let zoom: u64 = req.match_info().query("zoom").parse().unwrap();
    let path: u64 = req.match_info().query("filename").parse().unwrap();// .parse().unwrap();
    let data = data.read().unwrap();
    // printlnlet list = &data.list;
    //let data = data.read().unwrap().vis; //(*data).lock().unwrap(); //get_mut();
    let cache_dir = &data.params.cache_dir;
    match NamedFile::open(format!("{}/{}/{}_0.png", cache_dir, zoom, path)) {
        Ok(file) => Ok(file
            .use_last_modified(true)
            .set_content_disposition(ContentDisposition {
                disposition: DispositionType::Attachment,
                parameters: vec![],
            })),
        _ => {
            fs::create_dir( format!("{}/{}", cache_dir, zoom)); //error is permitted.
            let ann = &data.annotation;
            let params = &data.params;
            let freq = &data.freq;
            let compressed_list = &data.compressed_list;
            let index_list = &data.index_list;
            let supplementary_list = &data.supplementary_list;
            let prev_index = data.prev_index;

            let min_zoom = 13;
            
            let path_string = format!("{}/{}/{}_0.png", cache_dir, zoom, path);
            let max_zoom = (&data.params).max_zoom as u64;
            if zoom < min_zoom || zoom > max_zoom as u64 {
                return Err(error::ErrorBadRequest("zoom level is not appropriate"));
            }
            let (matches, string_range) = id_to_range(&data.range, &data.args, zoom, path, params, path_string.clone());

            // If the end is exceeds the prefetch region, raise error.
            // let arg_vec = vec!["ghb", "vis", "-t", "1", "-r",  "parse"];
            bam_record_vis(&matches, string_range, list.to_vec(), ann, freq, compressed_list, index_list, prev_index, supplementary_list,|_| None).unwrap();
            // bam_vis(matches, 1);
            Ok(NamedFile::open(path_string)?
                .use_last_modified(true)
                .set_content_disposition(ContentDisposition {
                    disposition: DispositionType::Attachment,
                    parameters: vec![],
                }))
        }
    }
}

/*
pub struct Vis<'a> {
    range: StringRegion,
    list: Mutex<Vec<(u64, Record)>>,
    annotation: Mutex<Vec<(u64, bed::Record)>>,
    index_list: Mutex<Vec<usize>>,
    suppl_list: Mutex<Vec<(&'a [u8], usize, usize, usize, usize)>>,
}*/



#[derive(Clone)]
pub struct Vis<'a> {
    range: StringRegion, args: Vec<String>, annotation: Vec<(u64, bed::Record)>, freq: BTreeMap<u64, Vec<(u64, u32)>>,     compressed_list: Vec<(u64, usize)>,
    index_list: Vec<usize>,
    prev_index: usize,
    supplementary_list: Vec<(&'a[u8], usize, usize, i32, i32)>,dzi: DZI, params: Param
}

impl<'a> Vis<'a> {
    fn new(range: StringRegion, args: Vec<String>, annotation: Vec<(u64, bed::Record)>, freq: BTreeMap<u64, Vec<(u64, u32)>>,     compressed_list: Vec<(u64, usize)>,
    index_list: Vec<usize>,
    prev_index: usize,
    supplementary_list: Vec<(&[u8], usize, usize, i32, i32)>,dzi: DZI, params: Param) -> Vis {
        Vis{range, args, annotation, freq, compressed_list, index_list, prev_index, supplementary_list, dzi, params}
    }
}

pub struct Item {
    list: Vec<(u64, Record)>,
    vis: Vis,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct DZI {
    Image: Image
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct Image {
    xmlns: String,
    Format: String,
    Overlap: u64,
    TileSize: u32,
    Size: Size
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct Size {
    height: u32,
    width: u32
}
#[derive(Debug, Clone)]
pub struct Param {
    x_scale: u32,
    max_y: u32,
    prefetch_max: u64,
    criteria: u64,
    max_zoom: u32,
    y_freq: u32,
    y: u32,
    cache_dir: String
}

const fn num_bits<T>() -> usize { std::mem::size_of::<T>() * 8 }

fn log_2(x: i32) -> u32 {
    assert!(x > 0);
    num_bits::<i32>() as u32 - x.leading_zeros() - 1
}

#[actix_rt::main]
pub async fn server(matches: ArgMatches, range: StringRegion, prefetch_range: StringRegion, args: Vec<String>, list: Vec<(u64, Record)>, annotation: Vec<(u64, bed::Record)>, freq: BTreeMap<u64, Vec<(u64, u32)>>,     compressed_list: Vec<(u64, usize)>,
index_list: Vec<usize>,
prev_index: usize,
supplementary_list: Vec<(&[u8], usize, usize, i32, i32)>,threads: u16) -> std::io::Result<()> {

    use actix_web::{web, HttpServer};
    let bind = matches.value_of("web").unwrap_or(&"0.0.0.0:4000");
    let x = matches
        .value_of("x")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(1280u32);
    let diff = range.end - range.start;
    let all = prefetch_range.end - prefetch_range.start;
    let size = Size{height: x, width: (all as u32 / diff as u32 + 1) * x};
    let image = Image{xmlns: "http://schemas.microsoft.com/deepzoom/2008".to_string(), Format: "png".to_string(), Overlap: 0, TileSize: x, Size: size};
    let dzi = DZI{Image: image};
    let x_scale = matches
        .value_of("x-scale")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(20u32);
    let freq_size = matches
        .value_of("freq-height")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(50u32);
    let y = matches
        .value_of("y")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(15u32);
    let mut rng = rand::thread_rng();
    let cache_dir = matches
        .value_of("cache-dir")
        .and_then(|a| Some(a.to_string()))
        .unwrap_or(rng.gen::<u32>().to_string());
    
    let x_width = all as u32 / diff as u32 * x;
    let max_zoom = 18; //log_2(x_width as i32) + 1;

    match fs::create_dir(&cache_dir) {
        Err(e) => panic!("{}: {}", &cache_dir, e),
        Ok(_) => {},
    };
    let params = Param{x_scale, max_y:x, prefetch_max: all, max_zoom, criteria:diff, y_freq: freq_size, y, cache_dir};
    eprintln!("{:?}, threads:{}", params, threads);
    // Create some global state prior to building the server
    //#[allow(clippy::mutex_atomic)] // it's intentional.
    //let counter1 = web::Data::new(Mutex::new((matches.clone(), range, list, annotation, freq)));
    // let counter = RwLock::new(Vis::new(range.clone(), args.clone(), list.clone(), annotation.clone(), freq.clone(), dzi.clone(), params.clone()));
    //let counter = Arc::new(RwLock::new(Vis::new(range, args, annotation, freq, dzi, params)));
    //#[allow(clippy::mutex_atomic)] 
    let counter = web::Data::new(RwLock::new(Vis::new(range, args, annotation, freq, compressed_list, index_list, prev_index, supplementary_list,dzi, params)));

    //https://github.com/actix/examples/blob/master/state/src/main.rs
    HttpServer::new(move|| {
        let list = list.clone();
        //let counter = Arc::new(RwLock::new(Item{list: list, vis: Vis::new(range, args, annotation, freq, dzi, params)}));
        //let counter = Cell::new(Vis::new( range.clone(),args.clone(), list.clone(), annotation.clone(), freq.clone()));
        actix_web::App::new().data(list).app_data(counter.clone()).route("genome.dzi", web::get().to(get_dzi)).route("/{zoom:.*}/{filename:.*}_0.png", web::get().to(index)).wrap(Logger::default())
    })
    .bind(bind)?
    .workers(threads as usize)
    .run()
    .await
}
