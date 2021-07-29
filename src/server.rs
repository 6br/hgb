use actix_cors::Cors;
use actix_files::NamedFile;
use actix_web::http::header::{ContentDisposition, DispositionType};
use actix_web::{error, middleware::Logger, web, HttpRequest, Responder, Result};
use bam::Record;
use clap::{App, AppSettings, Arg, ArgMatches};
use genomic_range::StringRegion;
use ghi::{bed, vis::bam_record_vis, Vis, VisRef};
use itertools::Itertools;
use rand::Rng;
use std::time::Instant;
use std::{collections::BTreeMap, fs, sync::RwLock};

fn id_to_range(
    range: &StringRegion,
    args: &[String],
    zoom: u64,
    path: u64,
    param: &Param,
    path_string: String,
) -> (ArgMatches, StringRegion) {
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
            .arg(Arg::new("bed-range").short('J').takes_value(true).about("Visualize multiple region where BED file specifies"))
            .arg(Arg::new("neighbor").short('K').takes_value(true).about("Visualize specified base-pair neighbor of BED region"))
            .arg(Arg::new("whole-chromosome").short('W').about("Pretend as if the prefetch range is the whole chromosome"))
            .arg(Arg::new("adjust-y").short('&').about("Do not adjust y on server mode"))
            .arg(Arg::new("no-filter").short('f').about("Disable pre-filtering on loading BAM index (used for debugging)"))
            .arg(Arg::new("no-cigar").short('c').about("Do not show cigar string"))
            .arg(Arg::new("no-scale").short('S').about("Do not show y-axis scale"))
            .arg(Arg::new("no-packing").short('p').about("Disable read packing"))
            .arg(Arg::new("no-legend").short('l').about("Hide legend"))
            .arg(Arg::new("colored-by-name").short('n').about("Set read colors by read name"))
            .arg(
                Arg::new("colored-by-motif")
                    .short('E')
                    .takes_value(true)
                    .default_value("C:T:CG")
                    .about("Colored by specified motif (for bisulfite sequencing)"),
            )
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
                    .short('V')
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
                Arg::new("min-read-length")
                    .short('M')
                    .takes_value(true)
                    .about("Minimum read length on coverage/alignment track"),
            )
            .arg(
                Arg::new("x-scale")
                    .short('X')
                    .takes_value(true)
                    .about("Size of x-scale legend"),
            )
            .arg(Arg::new("split-alignment").short('s').about("Display split alignments in the same "))
            .arg(Arg::new("only-split-alignment").short('u').about("Display only split alignments or mate-paired reads on alignment track"))
            .arg(Arg::new("udon").short('U').about("Colored by udon library"))
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
                Arg::new("format")
                    .short('O')
                    .takes_value(true)
                    .about("Output format (automatically detected; optionally used for server mode)"),
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
    let scalex_default = if param.y_adjust {
        param.x_scale as u64
    } else {
        param.x_scale as u64 >> (max_zoom - zoom)
    };
    let adjusted_y = if param.y_adjust {
        y
    } else {
        y >> (max_zoom - zoom)
    };
    let adjusted_large_y = if param.y_adjust {
        freq_y
    } else if (y >> (max_zoom - zoom)) <= 2 {
        max_y >> (max_zoom - zoom)
    } else {
        freq_y >> (max_zoom - zoom)
    };
    let adjusted_scale_x = if param.y_adjust {
        param.x_scale as u64
    } else {
        adjusted_large_y / 2
    };
    let only_coverage_condition = if param.y_adjust {
        criteria << (max_zoom - zoom) >= 500000
    } else {
        y >> (max_zoom - zoom) <= 2
    };

    let b: Vec<String> = if only_coverage_condition {
        // Only coverage
        vec![
            "-A".to_string(),
            "-Y".to_string(),
            adjusted_large_y.to_string(),
            "-y".to_string(),
            adjusted_y.to_string(),
            "-c".to_string(),
            "-I".to_string(),
            "-o".to_string(),
            path_string,
            "-X".to_string(),
            adjusted_scale_x.to_string(),
            "-x".to_string(),
            x.to_string(),
            "-l".to_string(),
            "-e".to_string(),
        ]
    } else if criteria << (max_zoom - zoom) <= 10000 {
        // Base-pair level with legend and insertion
        vec![
            "-Y".to_string(),
            adjusted_large_y.to_string(),
            "-y".to_string(),
            adjusted_y.to_string(),
            "-X".to_string(),
            scalex_default.to_string(),
            "-o".to_string(),
            path_string,
            "-x".to_string(),
            x.to_string(),
            "-e".to_string(),
        ]
    } else if criteria << (max_zoom - zoom) <= 25000 && (y >> (max_zoom - zoom)) >= 8 {
        // Base-pair level
        vec![
            "-Y".to_string(),
            adjusted_large_y.to_string(),
            "-y".to_string(),
            adjusted_y.to_string(),
            "-X".to_string(),
            scalex_default.to_string(),
            "-I".to_string(),
            "-o".to_string(),
            path_string,
            "-x".to_string(),
            x.to_string(),
            "-l".to_string(),
            "-e".to_string(),
        ]
    } else {
        // No alignment
        vec![
            "-Y".to_string(),
            adjusted_large_y.to_string(),
            "-y".to_string(),
            adjusted_y.to_string(),
            "-X".to_string(),
            scalex_default.to_string(),
            "-c".to_string(),
            "-I".to_string(),
            "-o".to_string(),
            path_string,
            "-x".to_string(),
            x.to_string(),
            "-l".to_string(),
            "-e".to_string(),
        ]
    };
    let mut args = args.to_owned();
    args.extend(b);
    args.remove(0);
    args = args.into_iter().skip_while(|t| t != "vis").collect();
    let start = (criteria << (max_zoom - zoom)) * path + range.start;
    let range = StringRegion::new(&format!(
        "{}:{}-{}",
        range.path,
        start,
        start + (criteria << (max_zoom - zoom)) - 1
    ))
    .unwrap();
    eprintln!("{:?} {:?}", args.join(" "), range);
    let matches = app.get_matches_from(args);
    (matches, range)
}

async fn get_dzi(data: web::Data<RwLock<Item>>) -> impl Responder {
    return web::Json(data.read().unwrap().dzi.clone());
}

async fn get_index(_data: web::Data<RwLock<Item>>) -> Result<NamedFile> {
    return Ok(NamedFile::open("static/index.html".to_string())?);
}

async fn get_js(_data: web::Data<RwLock<Item>>) -> Result<NamedFile> {
    return Ok(NamedFile::open("static/openseadragon.min.js".to_string())?);
}

async fn get_js_map(_data: web::Data<RwLock<Item>>) -> Result<NamedFile> {
    return Ok(NamedFile::open(
        "static/openseadragon.min.js.map".to_string(),
    )?);
}

async fn index(
    data: web::Data<RwLock<Item>>,
    list: web::Data<RwLock<Vec<(u64, Record)>>>,
    req: HttpRequest,
) -> Result<NamedFile> {
    let start = Instant::now();
    let zoom: u64 = req.match_info().query("zoom").parse().unwrap();
    let path: u64 = req.match_info().query("filename").parse().unwrap();
    let data = data.read().unwrap();

    let cache_dir = &data.params.cache_dir;
    match NamedFile::open(format!("{}/{}/{}_0.png", cache_dir, zoom, path)) {
        Ok(file) => Ok(file.set_content_disposition(ContentDisposition {
            disposition: DispositionType::Attachment,
            parameters: vec![],
        })),
        _ => {
            let end0 = start.elapsed();
            eprintln!(
                "match named file: {}.{:03} sec.",
                end0.as_secs(),
                end0.subsec_millis()
            );
            let params = &data.params;
            let args = &data.args;
            let data = &data.vis;
            let ann = &data.annotation;
            let freq = &data.freq;
            let compressed_list = &data.compressed_list;
            let index_list = &data.index_list;
            let supplementary_list = &data.supplementary_list;
            let prev_index = data.prev_index;

            //let min_zoom = 13;

            let path_string = format!("{}/{}/{}_0.png", cache_dir, zoom, path);
            let max_zoom = params.max_zoom as u64;
            //let min_zoom = ((&data.params).criteria << (max_zoom - zoom)) >= 10000000;
            let min_zoom = params.min_zoom as u64;
            if zoom < min_zoom || zoom > max_zoom as u64 {
                return Err(error::ErrorBadRequest(format!(
                    "zoom level {} should be between {} and {}",
                    zoom, min_zoom, max_zoom
                )));
            }
            fs::create_dir_all(format!("{}/{}", cache_dir, zoom))?; //error is permitted.
            let end1 = start.elapsed();
            eprintln!(
                "create dir: {}.{:03} sec.",
                end1.as_secs(),
                end1.subsec_millis()
            );
            let (matches, string_range) =
                id_to_range(&data.range, args, zoom, path, params, path_string.clone());
            let end2 = start.elapsed();
            eprintln!(
                "id_to_range: {}.{:03} sec.",
                end2.as_secs(),
                end2.subsec_millis()
            );
            // If the end is exceeds the prefetch region, raise error.
            // let arg_vec = vec!["ghb", "vis", "-t", "1", "-r",  "parse"];
            //bam_record_vis(&matches, vec![VisOrig::new(string_range, list.read().unwrap().to_vec(), ann.to_vec(), *freq, compressed_list, index_list.to_vec(), prev_index, supplementary_list)],|_| None).unwrap();
            bam_record_vis(
                &matches,
                vec![VisRef::new(
                    string_range,
                    &list.read().unwrap(),
                    ann,
                    freq,
                    compressed_list,
                    index_list,
                    prev_index,
                    supplementary_list,
                )],
                |_| None,
            )
            .unwrap();
            let end3 = start.elapsed();
            eprintln!(
                "img_saved: {}.{:03} sec.",
                end3.as_secs(),
                end3.subsec_millis()
            );
            // bam_vis(matches, 1);
            Ok(
                NamedFile::open(path_string)?.set_content_disposition(ContentDisposition {
                    disposition: DispositionType::Attachment,
                    parameters: vec![],
                }),
            )
        }
    }
}

pub struct Item {
    vis: Vis,
    args: Vec<String>,
    params: Param,
    dzi: DZI,
}

impl Item {
    fn new(vis: Vis, args: Vec<String>, params: Param, dzi: DZI) -> Self {
        Item {
            vis,
            args,
            params,
            dzi,
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DZI {
    #[serde(rename = "Image")]
    pub image: Image,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Image {
    pub xmlns: String,
    #[serde(rename = "Url")]
    pub url: String,
    #[serde(rename = "Format")]
    pub format: String,
    #[serde(rename = "Overlap")]
    pub overlap: String,
    #[serde(rename = "TileSize")]
    pub tile_size: String,
    #[serde(rename = "Size")]
    pub size: Size,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Size {
    #[serde(rename = "Height")]
    pub height: String,
    #[serde(rename = "Width")]
    pub width: String,
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
    cache_dir: String,
    y_adjust: bool,
}

const fn num_bits<T>() -> usize {
    std::mem::size_of::<T>() * 8
}

fn log_2(x: i32) -> u32 {
    assert!(x > 0);
    num_bits::<i32>() as u32 - x.leading_zeros() - 1
}

#[actix_rt::main]
pub async fn server(
    matches: ArgMatches,
    range: StringRegion,
    prefetch_range: StringRegion,
    args: Vec<String>,
    list: Vec<(u64, Record)>,
    annotation: Vec<(u64, bed::Record)>,
    freq: BTreeMap<u64, Vec<(u64, u32, char)>>,
    compressed_list: Vec<(u64, usize)>,
    index_list: Vec<usize>,
    prev_index: usize,
    supplementary_list: Vec<(Vec<u8>, usize, usize, i32, i32)>,
    threads: u16,
) -> std::io::Result<()> {
    use actix_web::{web, HttpServer};

    let bind = matches.value_of("web").unwrap_or(&"0.0.0.0:4000");
    let no_margin = matches.is_present("no-scale");
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
    let square = matches.is_present("square");
    let x = if square {
        top_margin
            + (prev_index as u32 + axis_count as u32 + annotation_count as u32 * 2) * y
            + freq.len() as u32 * freq_size
    } else {
        matches
            .value_of("x")
            .and_then(|a| a.parse::<u32>().ok())
            .unwrap_or(1280u32)
    };
    let diff = range.end - range.start;
    let all = if matches.is_present("whole-chromosome") {
        250000000
    } else {
        prefetch_range.end - prefetch_range.start
    };

    let size = Size {
        height: x.to_string(),
        width: ((all as u32 / diff as u32 + 1) * x).to_string(),
    };
    let image = Image {
        xmlns: "http://schemas.microsoft.com/deepzoom/2008".to_string(),
        url: format!("http://{}/", bind),
        format: "png".to_string(),
        overlap: "0".to_string(),
        tile_size: x.to_string(),
        size,
    };
    let dzi = DZI { image };
    let x_scale = matches
        .value_of("x-scale")
        .and_then(|a| a.parse::<u32>().ok())
        .unwrap_or(20u32);
    let mut rng = rand::thread_rng();
    let cache_dir = matches
        .value_of("cache-dir")
        .map(|a| a.to_string())
        .unwrap_or_else(|| rng.gen::<u32>().to_string());
    let y_adjust = matches.is_present("adjust-y");

    let x_width = all as u32 / diff as u32 * x;
    let max_zoom = log_2(x_width as i32) + 1;
    let min_zoom = if y_adjust { 2 } else { max_zoom - 8 }; // + (prev_index );

    if let Err(e) = fs::create_dir(&cache_dir) {
        panic!("{}: {}", &cache_dir, e)
    }
    let params = Param {
        x_scale,
        max_y: x,
        prefetch_max: all,
        max_zoom,
        min_zoom,
        criteria: diff,
        y_freq: freq_size,
        x,
        y,
        cache_dir,
        y_adjust,
    };
    eprintln!(
        "{:?}, threads: {}, zoom: {}",
        params,
        threads,
        log_2(x_width as i32) + 1
    );
    println!("Server is running on {}", bind);
    // Create some global state prior to building the server
    let vis = Vis::new(
        range,
        annotation,
        freq,
        compressed_list,
        index_list,
        prev_index,
        supplementary_list,
        250000000,
    );
    let counter = web::Data::new(RwLock::new(Item::new(vis, args, params, dzi)));
    //let buffer = web::Data::new(RwLock::new(ChromosomeBuffer::new()));
    let cross_origin_bool = matches.is_present("production");

    //https://github.com/actix/examples/blob/master/state/src/main.rs
    HttpServer::new(move || {
        let cross_origin = if cross_origin_bool {
            Cors::default()
        } else {
            Cors::permissive()
        };

        actix_web::App::new()
            .data(())
            .app_data(web::Data::new(RwLock::new(list.clone())))
            .app_data(counter.clone())
            .route("/", web::get().to(get_index))
            .route("openseadragon.min.js", web::get().to(get_js))
            .route("openseadragon.min.js.map", web::get().to(get_js_map))
            .route("genome.dzi", web::get().to(get_dzi))
            .route("/{zoom:.*}/{filename:.*}_0.png", web::get().to(index))
            .service(actix_files::Files::new("/images", "static/images").show_files_listing())
            .wrap(Logger::default())
            .wrap(cross_origin)
    })
    .bind(bind)?
    .workers(threads as usize)
    .run()
    .await
}
