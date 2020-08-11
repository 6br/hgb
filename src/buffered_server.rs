
use actix_files::NamedFile;
use rand::Rng;
use std::time::Instant;
use actix_cors::Cors;
use actix_web::http::header::{ContentDisposition, DispositionType};
use actix_web::{HttpRequest, Result, web, Responder, error, middleware::Logger};
use std::{sync::{RwLock},  collections::{BTreeSet}, fs};
use clap::{App,  ArgMatches, Arg, AppSettings};
use ghi::{vis::bam_record_vis, Vis, simple_buffer::ChromosomeBuffer};
use genomic_range::StringRegion;
use bam::Record;
use itertools::Itertools;

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
            )                .arg(Arg::new("whole-chromosome").short('W').about("Pretend as if the prefetch range is the whole chromosome"))
            .arg(Arg::new("no-filter").short('f').about("Disable pre-filtering on loading BAM index (used for debugging)"))
            .arg(Arg::new("no-cigar").short('c').about("Do not show cigar string"))
            .arg(Arg::new("no-scale").short('S').about("Do not show y-axis scale"))
            .arg(Arg::new("no-packing").short('p').about("Disable read packing"))
            .arg(Arg::new("no-legend").short('l').about("Hide legend"))
            .arg(Arg::new("colored-by-name").short('n').about("Set read colors by read name"))
            .arg(Arg::new("meaningless").short('Q').about("Set square width (overwritten)"))
            .arg(Arg::new("quality").short('q').about("Display reads by quality value"))
            .arg(Arg::new("x").short('x').takes_value(true).about("The width of image"))
            .arg(Arg::new("y").short('y').takes_value(true).about("The height of each read alignment"))
            .arg(Arg::new("web").short('w').takes_value(true).about("Serve the web server"))
            .arg(Arg::new("end-split-callets").short('e').about("Show callets on ends of read alignments if the read contains split-alignment"))
            .arg(
                Arg::new("freq-height")
                    .short('Y')
                    .takes_value(true)
                    .about("The height of each coverage track"),
            )
            .arg(
                Arg::new("snp-frequency")
                    .short('Z')
                    .takes_value(true)
                    .about("The portion of SNP frequency to display on each coverage track"),
            )
            .arg(
                Arg::new("max-coverage")
                    .short('m')
                    .takes_value(true)
                    .about("Max coverage value on coverage track"),
            )
            .arg(
                Arg::new("cache-dir")
                    .short('d')
                    .takes_value(true)
                    .about("Cache directory for server (generated randomly if not specified)"),
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
                    .short('G')
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
    // let prefetch_max = 25000000; // prefetch_range.end() - prefetch_range.start();
    let criteria = param.criteria; // range.end() - range.start();
    // let x_width = prefetch_max / criteria * (740); //max_x
    // Here 2**17 < x_width < 2**18 so maxZoom should be 18+ 1;
    let max_zoom = param.max_zoom as u64; // 18;
    let max_y = param.max_y as u64; // It is the same as x;
    let freq_y = param.y_freq as u64;
    let y = param.y as u64;
    let x = param.x;
    let scalex_default = param.x_scale as u64;
    // let path = path << (max_zoom - zoom);
    let b: Vec<String> = if (y >> (max_zoom-zoom)) <= 2 { // i.e. 2**22s
        vec!["-A".to_string(), "-Y".to_string(), (max_y >> (max_zoom-zoom)).to_string(), "-y".to_string(), (y >> (max_zoom-zoom)).to_string(), "-c".to_string(), "-I".to_string(), "-o".to_string(), path_string, "-X".to_string(), ((max_y >> (max_zoom-zoom)) / 2).to_string(), "-x".to_string(), x.to_string(), "-l".to_string(), "-e".to_string()]
    } else if criteria << (max_zoom - zoom) <= 10000 { // Base-pair level with legend and insertion
        vec!["-Y".to_string(), (freq_y >> (max_zoom-zoom)).to_string(), "-y".to_string(), (y >> (max_zoom-zoom)).to_string(), "-X".to_string(), (scalex_default >> (max_zoom-zoom)).to_string(), "-o".to_string(), path_string, "-x".to_string(), x.to_string(), "-e".to_string()]
    } else if criteria << (max_zoom - zoom) <= 25000 && (y >> (max_zoom-zoom)) >= 8 { // Base-pair level
        vec!["-Y".to_string(), (freq_y >> (max_zoom-zoom)).to_string(), "-y".to_string(), (y >> (max_zoom-zoom)).to_string(), "-X".to_string(), (scalex_default >> (max_zoom-zoom)).to_string(), "-I".to_string(), "-o".to_string(), path_string, "-x".to_string(), x.to_string(), "-l".to_string(), "-e".to_string()]
    } else { // No alignment
        vec!["-Y".to_string(), (freq_y >> (max_zoom-zoom)).to_string(), "-y".to_string(), (y >> (max_zoom-zoom)).to_string(), "-X".to_string(), (scalex_default >> (max_zoom-zoom)).to_string(), "-c".to_string(), "-I".to_string(), "-o".to_string(), path_string, "-x".to_string(), x.to_string(), "-l".to_string(), "-e".to_string()]
    };
    let mut args = args.clone();
    args.extend(b);
    args.remove(0);
    args = args.into_iter().skip_while(|t| t != "vis").collect();
    let start =  (criteria << (max_zoom - zoom)) * path + range.start;
    let range = StringRegion::new(&format!("{}:{}-{}", range.path, start, start + (criteria << (max_zoom - zoom)) - 1).to_string()).unwrap();
    eprintln!("{:?} {:?}", args.join(" "), range);
    let matches = app.get_matches_from(args);
    (matches, range)
}

async fn get_dzi(data: web::Data<RwLock<Item>>) -> impl Responder {
    return web::Json(data.read().unwrap().dzi.clone());
}

async fn get_index(data: web::Data<RwLock<Item>>) -> Result<NamedFile> {
    return Ok(NamedFile::open(format!("static/index.html"))?);
}

async fn get_js(data: web::Data<RwLock<Item>>) -> Result<NamedFile> {
    return Ok(NamedFile::open(format!("static/openseadragon.min.js"))?);
}

async fn get_js_map(data: web::Data<RwLock<Item>>) -> Result<NamedFile> {
    return Ok(NamedFile::open(format!("static/openseadragon.min.js.map"))?);
}


async fn index(item: web::Data<RwLock<Item>>, vis: web::Data<RwLock<Vis>>, list: web::Data<RwLock<Vec<(u64, Record)>>>, list_btree: web::Data<RwLock<BTreeSet<u64>>>, buffer: web::Data<RwLock<ChromosomeBuffer>>, req: HttpRequest) -> Result<NamedFile> {
    let start = Instant::now();
    let zoom: u64 = req.match_info().query("zoom").parse().unwrap();
    let path: u64 = req.match_info().query("filename").parse().unwrap();
    let data = item.read().unwrap();

    let cache_dir = &data.params.cache_dir;
    match NamedFile::open(format!("{}/{}/{}_0.png", cache_dir, zoom, path)) {
        Ok(file) => Ok(file
            .set_content_disposition(ContentDisposition {
                disposition: DispositionType::Attachment,
                parameters: vec![],
            })),
        _ => {
            let end0 = start.elapsed();
            eprintln!(
                "match named file: {}.{:03} sec.",
                end0.as_secs(),
                end0.subsec_nanos() / 1_000_000
            );
            let params = &data.params;
            let args = &data.args;



            //let min_zoom = 13;
            
            let path_string = format!("{}/{}/{}_0.png", cache_dir, zoom, path);
            let max_zoom = params.max_zoom as u64;
            //let min_zoom = ((&data.params).criteria << (max_zoom - zoom)) >= 10000000;
            let min_zoom = zoom < params.min_zoom as u64;

            if min_zoom || zoom > max_zoom as u64 {
                return Err(error::ErrorBadRequest("zoom level is not appropriate"));
            }
            fs::create_dir( format!("{}/{}", cache_dir, zoom)); //error is permitted.
            let end1 = start.elapsed();
            eprintln!(
                "create dir: {}.{:03} sec.",
                end1.as_secs(),
                end1.subsec_nanos() / 1_000_000
            );

            let (matches, string_range) = id_to_range(&data.range, args, zoom, path, params, path_string.clone());
            if !buffer.read().unwrap().included_string_local(&string_range, &list_btree.read().unwrap()) {
                // TODO() This ignores 
                let endx = start.elapsed();
                eprintln!("Fallback to reload: {}.{:03} sec.",
                    endx.as_secs(),
                    endx.subsec_nanos() / 1_000_000
                );
                //let (mut list, mut list_btree) = &*;
                buffer.write().unwrap().retrieve(&string_range, &mut list.write().unwrap(), &mut list_btree.write().unwrap());
                let new_vis = buffer.read().unwrap().vis(&string_range, &mut list.write().unwrap(), &mut list_btree.write().unwrap());
                let endy = start.elapsed();
                eprintln!("Fallback to reload2: {}.{:03} sec.",
                    endy.as_secs(),
                    endy.subsec_nanos() / 1_000_000
                );
                let mut old_vis = vis.write().unwrap();
                *old_vis = new_vis.unwrap();
                let endz = start.elapsed();
                eprintln!("Fallback to reload3: {}.{:03} sec.",
                    endz.as_secs(),
                    endz.subsec_nanos() / 1_000_000
                );
            }


            let data = vis.read().unwrap();
            let ann = &data.annotation;
            let freq = &data.freq;
            let compressed_list = &data.compressed_list;
            let index_list = &data.index_list;
            let supplementary_list = &data.supplementary_list;
            let prev_index = data.prev_index;
            let end2 = start.elapsed();
            eprintln!(
                "id_to_range: {}.{:03} sec.",
                end2.as_secs(),
                end2.subsec_nanos() / 1_000_000
            );
            // If the end is exceeds the prefetch region, raise error.
            // let arg_vec = vec!["ghb", "vis", "-t", "1", "-r",  "parse"];
            bam_record_vis(&matches, string_range, &list.read().unwrap(), ann, freq, compressed_list, index_list, prev_index, supplementary_list,|_| None).unwrap();
            let end3 = start.elapsed();
            eprintln!(
                "img_saved: {}.{:03} sec.",
                end3.as_secs(),
                end3.subsec_nanos() / 1_000_000
            );
            // bam_vis(matches, 1);
            Ok(NamedFile::open(path_string)?
                .set_content_disposition(ContentDisposition {
                    disposition: DispositionType::Attachment,
                    parameters: vec![],
                }))
        }
    }
}

pub struct Item {
    //vis: Vis,
    range: StringRegion,
    args: Vec<String>,
    params: Param,
    dzi: DZI,
}

impl Item {
    fn new(range: StringRegion, args: Vec<String>, params: Param, dzi: DZI) -> Self {
        Item{range, args:args,params:params,dzi:dzi}
    }
}

/*
impl Vis {
    fn new(range: StringRegion, args: Vec<String>, annotation: Vec<(u64, bed::Record)>, freq: BTreeMap<u64, Vec<(u64, u32, char)>>, compressed_list: Vec<(u64, usize)>,
    index_list: Vec<usize>,
    prev_index: usize,
    supplementary_list: Vec<(Vec<u8>, usize, usize, i32, i32)>) -> Vis {
        Vis{range, args, annotation, freq, compressed_list, index_list, prev_index, supplementary_list}
    }
}*/
#[derive(Debug, Serialize, Deserialize, Clone)]
struct DZI {
    Image: Image,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct Image {
    xmlns: String,
    Url: String,
    Format: String,
    Overlap: String,
    TileSize: String,
    Size: Size,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct Size {
    Height: String,
    Width: String
}
#[derive(Debug, Clone)]
pub struct Param {
    x_scale: u32,
    max_y: u32,
    prefetch_max: u64,
    criteria: u64,
    max_zoom: u32,
    min_zoom: u32,
    y_freq: u32,
    x: u32,
    y: u32,
    cache_dir: String
}


const fn num_bits<T>() -> usize { std::mem::size_of::<T>() * 8 }

fn log_2(x: i64) -> u32 {
    assert!(x > 0);
    num_bits::<i64>() as u32 - x.leading_zeros() - 1
}

#[actix_rt::main]
pub async fn server(matches: ArgMatches, range: StringRegion, prefetch_range: StringRegion, args: Vec<String>, mut buffer:  ChromosomeBuffer, threads: u16) -> std::io::Result<()> {

    use actix_web::{web, HttpServer};
    // let list = buffer.add(&prefetch_range);
    let mut list = vec![];
    let mut list_btree = BTreeSet::new();
    let bind = matches.value_of("web").unwrap_or(&"0.0.0.0:4000");
    let no_margin = matches.is_present("no-scale");
    buffer.retrieve(&prefetch_range, &mut list, &mut list_btree);
    let vis = buffer.vis(&prefetch_range, &mut list, &mut list_btree).unwrap();
    //eprintln!("{:#?}", vis);
    let annotation = &vis.annotation;
    let prev_index = vis.prev_index;
    let freq = &vis.freq;
    let annotation_count = annotation.iter().unique_by(|s| s.0).count(); // annotation.len();
    let top_margin = if no_margin { 0 } else { 40 };
    let axis_count = 0;
    let y = matches
        .value_of("y")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(15u32);
    let freq_size = matches
        .value_of("freq-height")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(50u32);
    let zoom_range = matches
        .value_of("zoom-range")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(if matches.is_present("whole-chromosome") {6u32} else{8u32});
    let square = matches.is_present("square");
    let x = if square {
        top_margin
        + (prev_index as u32 + axis_count as u32 + annotation_count as u32 * 2) * y
        + freq.len() as u32 * freq_size
    } else {matches
        .value_of("x")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(1280u32)
    };
    let diff = range.end - range.start;
    let all = if matches.is_present("whole-chromosome") {vis.prefetch_max} else {prefetch_range.end - prefetch_range.start};
    let view_range = if matches.is_present("whole-chromosome") {StringRegion{path: prefetch_range.path, start: 1, end: vis.prefetch_max}} else {prefetch_range}; 
    let size = Size{Height: x.to_string(), Width: ((all as u32 / diff as u32 + 1) * x).to_string()};
    let image = Image{xmlns: "http://schemas.microsoft.com/deepzoom/2008".to_string(), Url: format!("http://{}/", bind).to_string(), Format: "png".to_string(), Overlap: "0".to_string(), TileSize: x.to_string(), Size: size};
    let dzi = DZI{Image: image};
    let x_scale = matches
        .value_of("x-scale")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(20u32);
    let mut rng = rand::thread_rng();
    let cache_dir = matches
        .value_of("cache-dir")
        .and_then(|a| Some(a.to_string()))
        .unwrap_or(rng.gen::<u32>().to_string());
    
    let x_width = all as u32 / diff as u32 * x;
    eprintln!("Total len: {}, partial len: {}, x_width: {}({}), x: {}", all,diff,x_width, x_width as i64, x);
    let max_zoom = log_2(x_width as i64) + 1;
    let min_zoom = max_zoom - zoom_range;

    match fs::create_dir(&cache_dir) {
        Err(e) => panic!("{}: {}", &cache_dir, e),
        Ok(_) => {},
    };
    let params = Param{x_scale, max_y:x, prefetch_max: all, max_zoom, min_zoom, criteria:diff, y_freq: freq_size, x, y, cache_dir};
    eprintln!("{:?}, threads: {}, zoom: {}", params, threads, log_2(x_width as i64) + 1);
    println!("Server is running on {}", bind);
    // Create some global state prior to building the server
    //#[allow(clippy::mutex_atomic)] // it's intentional.
    //let counter1 = web::Data::new(Mutex::new((matches.clone(), range, list, annotation, freq)));
    // let counter = RwLock::new(Vis::new(range.clone(), args.clone(), list.clone(), annotation.clone(), freq.clone(), dzi.clone(), params.clone()));
    //let counter = Arc::new(RwLock::new(Vis::new(range, args, annotation, freq, dzi, params)));
    //#[allow(clippy::mutex_atomic)] 
    // let vis = Vis::new(range, annotation, freq, compressed_list, index_list, prev_index, supplementary_list, 250000000);
    let counter = web::Data::new(RwLock::new(Item::new(view_range, args, params, dzi)));
    let vis = web::Data::new(RwLock::new(vis));
    let buffer = web::Data::new(RwLock::new(buffer));

    //https://github.com/actix/examples/blob/master/state/src/main.rs
    HttpServer::new(move|| {
        let list = RwLock::new(list.clone());
        let list_btree = RwLock::new(list_btree.clone());//buffer = 
        //let buffer = buffer.clone();
        actix_web::App::new().data(list).data(list_btree)/*.app_data(web::Data::new(RwLock::new(list.clone()))).*/.app_data(counter.clone()).app_data(vis.clone())
        .app_data(buffer.clone()).route("/", web::get().to(get_index))
        .route("openseadragon.min.js", web::get().to(get_js))
        .route("openseadragon.min.js.map", web::get().to(get_js_map))
        .route("genome.dzi", web::get().to(get_dzi))
        .route("/{zoom:.*}/{filename:.*}_0.png", web::get().to(index)).service(actix_files::Files::new("/images", "static/images").show_files_listing()).wrap(Logger::default()).wrap(
            Cors::new().supports_credentials() /*allowed_origin("*").allowed_methods(vec!["GET", "POST"])
            .allowed_headers(vec![http::header::AUTHORIZATION, http::header::ACCEPT])
            .allowed_header(http::header::CONTENT_TYPE)
            .max_age(3600)*/
            .finish()
        )
    })
    .bind(bind)?
    .workers(threads as usize)
    .run()
    .await
}
