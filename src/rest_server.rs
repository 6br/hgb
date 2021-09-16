use actix_cors::Cors;
use actix_files::NamedFile;
use actix_web::{
    http::header::{ContentDisposition, DispositionType},
    HttpResponse,
};
use qstring::QString;
use actix_web::HttpRequest;
use actix_web::{middleware::Logger, web, Result};
use bam::Record;
use clap::{App, AppSettings, Arg, ArgMatches, Error};
use genomic_range::StringRegion;
use ghi::{simple_buffer::ChromosomeBuffer, vis::bam_record_vis, Vis, VisRef};
use rand::Rng;
use serde::Deserialize;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::time::Instant;
use std::{collections::BTreeSet, fs, sync::RwLock};

use crate::subcommands::{bam_vis, vis_query};

struct Item {
    //vis: Vis,
    range: StringRegion,
    args: Vec<String>,
    cache_dir: String,
}

impl Item {
    fn new(range: StringRegion, args: Vec<String>, cache_dir: String) -> Self {
        Item {
            range,
            args,
            cache_dir,
        }
    }
}

#[derive(Deserialize, Hash)]
pub struct RequestBody {
    params: String,
    format: String,
    prefetch: bool,
}

#[derive(Serialize)]
pub struct ResponseBody {
    message: String,
}

fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}

fn get_matches_from(args: Vec<String>) -> Result<ArgMatches, Error> {
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
    .arg(Arg::new("bed-range").short('J').takes_value(true).about("Visualize multiple region where BED file specifies"))
    .arg(Arg::new("neighbor").short('K').takes_value(true).about("Visualize specified base-pair neighbor of BED region"))
    .arg(
        Arg::new("format")
            .short('O')
            .takes_value(true)
            .about("Output format (automatically detected; optionally used for server mode)"),
    )
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
    .arg(Arg::new("meaningless2").short('>').about("Serve the web server with accepting any parameter"))
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
        Arg::new("x-scale")
            .short('X')
            .takes_value(true)
            .about("Size of x-scale legend"),
    )
    .arg(Arg::new("split-alignment").short('s').about("Display split alignments in the same "))
    .arg(Arg::new("only-split-alignment").short('u').about("Display only split alignments or mate-paired reads on alignment track"))
    .arg(Arg::new("exclude-split-alignment").short('3').about("Display only NOT split alignments or mate-paired reads on alignment track"))
    .arg(Arg::new("full-length").short('4').about("Display only full-length match reads against prefetch range."))
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
    app.try_get_matches_from(args)
}

fn id_to_range(
    _range: &StringRegion,
    args: &[String],
    params: String,
    path_string: String,
) -> Result<(ArgMatches, StringRegion), Error> {
    let a: Vec<String> = params.split(' ').map(|t| t.to_string()).collect();
    let b: Vec<String> = vec!["-o".to_string(), path_string];
    let mut args = args.to_owned();
    args.extend(b);
    args.extend(a);
    args.remove(0);
    args = args.into_iter().skip_while(|t| t != "vis").collect();
    eprintln!("{:?}", args.join(" "));
    let matches = get_matches_from(args)?;
    let ranges_str = matches.values_of("range").unwrap();
    let range = {
        let ranges_tmp: Vec<String> = ranges_str.into_iter().map(|t| t.to_string()).collect();
        StringRegion::new(&ranges_tmp[ranges_tmp.len() - 1]).unwrap()
    };
    eprintln!("{:?}", range);
    Ok((matches, range))
}

fn id_to_range_ab_initio(
    params: String,
    path_string: String,
) -> Result<(ArgMatches, Vec<String>), Error> {
    let a: Vec<String> = params.split(' ').map(|t| t.to_string()).collect();
    let b: Vec<String> = vec!["-o".to_string(), path_string];
    let mut args = vec!["vis".to_string()];
    args.extend(a);
    args.extend(b);
    eprintln!("{:?}", args);
    let matches = get_matches_from(args.clone())?;
    eprintln!("{:?}", matches.value_of("INPUT"));
    Ok((matches, args))
}

//     
async fn get_index(
    req: HttpRequest,
    item: web::Data<RwLock<Item>>,
    vis: web::Data<RwLock<Vis>>,
    list: web::Data<RwLock<Vec<(u64, Record)>>>,
    list_btree: web::Data<RwLock<BTreeSet<u64>>>,
    buffer: web::Data<RwLock<ChromosomeBuffer>>,
    //query: web::Query<RequestBody>
) -> Result<NamedFile> {
    let qs = QString::from(req.query_string());
    let format = qs.get("format").unwrap().clone(); // &query.format.clone();
    let params = qs.get("params").unwrap().clone();
    let prefetch = qs.get("params").is_some();
    let hash: u64 = calculate_hash(&RequestBody{format: format.to_string(), params: params.to_string(), prefetch});
    return index2(item, vis, list, list_btree, buffer, format.to_string(), params.to_string(), prefetch, hash);
}

async fn index(
    item: web::Data<RwLock<Item>>,
    vis: web::Data<RwLock<Vis>>,
    list: web::Data<RwLock<Vec<(u64, Record)>>>,
    list_btree: web::Data<RwLock<BTreeSet<u64>>>,
    buffer: web::Data<RwLock<ChromosomeBuffer>>,
    request_body: web::Json<RequestBody>,
) -> Result<NamedFile> {
    let format = &request_body.format.clone();
    let params = &request_body.params.clone();
    let prefetch = &request_body.prefetch.clone();
    let hash: u64 = calculate_hash(&request_body.into_inner());
    return index2(item, vis, list, list_btree, buffer, format.to_string(), params.to_string(), *prefetch, hash);
}

fn index2(
        item: web::Data<RwLock<Item>>,
        vis: web::Data<RwLock<Vis>>,
        list: web::Data<RwLock<Vec<(u64, Record)>>>,
        list_btree: web::Data<RwLock<BTreeSet<u64>>>,
        buffer: web::Data<RwLock<ChromosomeBuffer>>,
        format: String,
        params: String,
        prefetch: bool,
        hash: u64
    ) -> Result<NamedFile> {   
    let data = item.read().unwrap();
    let cache_dir = &data.cache_dir;
    let args = &data.args;
    let start = Instant::now();
    let path_string = format!("{}/{}.{}", cache_dir, hash, format);
    eprintln!("{} {} {:?}", format, params, path_string);

    match NamedFile::open(path_string.clone()) {
        Ok(file) => Ok(file.set_content_disposition(ContentDisposition {
            disposition: DispositionType::Attachment,
            parameters: vec![],
        })),
        _ => {
            if prefetch {
                let end0 = start.elapsed();
                eprintln!(
                    "match named file: {}.{:03} sec.",
                    end0.as_secs(),
                    end0.subsec_millis()
                );

                let end1 = start.elapsed();
                eprintln!(
                    "create dir: {}.{:03} sec.",
                    end1.as_secs(),
                    end1.subsec_millis()
                );

                let (matches, string_range) =
                    id_to_range(&data.range, args, params.to_string(), path_string.clone())
                        .map_err(|t| {
                            HttpResponse::BadRequest().json(ResponseBody {
                                message: format!("parameter error: {}", t), //String::from("parameter error: " + t.description()),
                            })
                        })?;
                let end2 = start.elapsed();
                eprintln!(
                    "id_to_range: {}.{:03} sec.",
                    end2.as_secs(),
                    end2.subsec_millis()
                );
                if !buffer
                    .read()
                    .unwrap()
                    .included_string_local(&string_range, &list_btree.read().unwrap())
                {
                    // TODO() This ignores
                    let endx = start.elapsed();
                    eprintln!(
                        "Fallback to reload: {}.{:03} sec.",
                        endx.as_secs(),
                        endx.subsec_millis()
                    );
                    //let (mut list, mut list_btree) = &*;
                    {
                        buffer.write().unwrap().retrieve(
                            &string_range,
                            &mut list.write().unwrap(),
                            &mut list_btree.write().unwrap(),
                        );
                    }
                    let new_vis = buffer.read().unwrap().vis(
                        &string_range,
                        &mut list.write().unwrap(),
                        &mut list_btree.write().unwrap(),
                    );
                    let endy = start.elapsed();
                    eprintln!(
                        "Fallback to reload2: {}.{:03} sec.",
                        endy.as_secs(),
                        endy.subsec_millis()
                    );
                    let mut old_vis = vis.write().unwrap();
                    *old_vis = new_vis.unwrap();
                    let endz = start.elapsed();
                    eprintln!(
                        "Fallback to reload3: {}.{:03} sec.",
                        endz.as_secs(),
                        endz.subsec_millis()
                    );
                }

                let data = vis.read().unwrap();
                let ann = &data.annotation;
                let freq = &data.freq;
                let compressed_list = &data.compressed_list;
                let index_list = &data.index_list;
                let supplementary_list = &data.supplementary_list;
                let prev_index = data.prev_index;

                // If the end is exceeds the prefetch region, raise error.
                // let arg_vec = vec!["ghb", "vis", "-t", "1", "-r",  "parse"];
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
            } else {
                //Visualization for unprefetch data.
                let (matches, args) =
                    id_to_range_ab_initio(params.to_string(), path_string.clone()).map_err(
                        |t| {
                            HttpResponse::BadRequest().json(ResponseBody {
                                message: format!("parameter error: {}", t), //String::from("parameter error: " + t.description()),
                            })
                        },
                    )?;
                let threads = matches
                    .value_of("threads")
                    .and_then(|t| t.parse::<u16>().ok())
                    .unwrap_or(1u16);
                match matches.is_present("INPUT") {
                    true => vis_query(&matches, args, threads).unwrap(),
                    false => bam_vis(&matches, args, threads).unwrap(),
                }
            }
            Ok(
                NamedFile::open(path_string)?.set_content_disposition(ContentDisposition {
                    disposition: DispositionType::Attachment,
                    parameters: vec![],
                }),
            )
        }
    }
}

#[actix_rt::main]
pub async fn server(
    matches: ArgMatches,
    _range: StringRegion,
    prefetch_range: StringRegion,
    args: Vec<String>,
    mut buffer: ChromosomeBuffer,
    threads: u16,
) -> std::io::Result<()> {
    use actix_web::{web, HttpServer};
    // let list = buffer.add(&prefetch_range);
    let mut list = vec![];
    let mut list_btree = BTreeSet::new();
    let bind = matches.value_of("web").unwrap_or(&"0.0.0.0:4000");
    buffer.retrieve(&prefetch_range, &mut list, &mut list_btree);
    let vis = buffer
        .vis(&prefetch_range, &mut list, &mut list_btree)
        .unwrap();
    let view_range = if matches.is_present("whole-chromosome") {
        StringRegion {
            path: prefetch_range.path,
            start: 1,
            end: vis.prefetch_max,
        }
    } else {
        prefetch_range
    };

    let mut rng = rand::thread_rng();
    let cache_dir = matches
        .value_of("cache-dir")
        .map(|a| a.to_string())
        .unwrap_or_else(|| rng.gen::<u32>().to_string());

    if let Err(e) = fs::create_dir(&cache_dir) {
        panic!("{}: {}", e, &cache_dir)
    }
    println!("REST Server is running on {}", bind);
    // Create some global state prior to building the server

    let counter = web::Data::new(RwLock::new(Item::new(view_range, args, cache_dir)));
    let buffer = web::Data::new(RwLock::new(buffer));

    let cross_origin_bool = matches.is_present("production");

    // https://github.com/actix/examples/blob/master/state/src/main.rs
    HttpServer::new(move || {
        let cross_origin = if cross_origin_bool {
            Cors::default()
        } else {
            Cors::permissive()
        };

        let list = RwLock::new(list.clone());
        let list_btree = RwLock::new(list_btree.clone()); //buffer =
        let vis = RwLock::new(vis.clone());
        //let buffer = buffer.clone();
        actix_web::App::new()
            .data(list)
            .data(list_btree)
            .app_data(counter.clone())
            .data(vis) //.app_data(vis.clone())
            .app_data(buffer.clone())
            .route("/", web::post().to(index))
            .route("/", web::get().to(get_index))
            .wrap(Logger::default())
            .wrap(cross_origin)
    })
    .bind(bind)?
    .workers(threads as usize)
    .run()
    .await
}
