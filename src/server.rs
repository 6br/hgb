
use actix_files::NamedFile;

use actix_web::http::header::{ContentDisposition, DispositionType};
use actix_web::{HttpRequest, Result, web};
use std::{sync::{Arc, Mutex}, path::PathBuf, cell::Cell, collections::BTreeMap};
use clap::{App,  ArgMatches, Arg};
use ghi::{bed, vis::bam_record_vis};
// use crate::subcommands::bam_vis;
use genomic_range::StringRegion;
use bam::Record;

fn id_to_range<'a>(range: &StringRegion, args: &Vec<String>, zoom: i64, path: i64) -> (ArgMatches, StringRegion) {

    let app = App::new("GHB/GHI genomic data visualization tool")
    // .setting(AppSettings::ArgsNegateSubcommands)
    .version("0.1")
    .author("6br. ")
    .about(
        "Command-line visualization tool for read alignment, ranged annotation and frequency.",
    )
    .arg(
        Arg::new("verbose")
            .short('v')
            .multiple(true)
            .about("Sets the level of verbosity"),
    )
    .arg(
        Arg::new("threads")
            .short('t')
            .default_value("1")
            .about("Sets the number of threads"),
    )
    .subcommand(
        App::new("build")
            .about("construct GHB and GHI index from bam/")
            .arg(
                Arg::new("bam")
                    .short('a')
                    .takes_value(true)
                    .multiple(true)
                    .about("sorted bam"),
            )
            .arg(
                Arg::new("bed")
                    .short('b')
                    .takes_value(true)
                    .multiple(true)
                    .about("sorted bed"),
            )
            .arg(
                Arg::new("gff3")
                    .short('g')
                    .takes_value(true)
                    .multiple(true)
                    .about("sorted gff3"),
            )
            .arg(
                Arg::new("chrom")
                    .short('c')
                    .takes_value(true)
                    .about("chroms.sizes"),
            )
            .arg(Arg::new("force").short('f').about("Outputs only header"))
            .arg(
                Arg::new("OUTPUT")
                    .about("Sets the output file to use")
                    .takes_value(true)
                    .required(true)
                    .index(1),
            ), /*.arg(Arg::new("gff")
               .short('g')
               .multiple(true)
               .about("sorted gff"))*/
    )
    .subcommand(
        App::new("query")
            .about("construct GHB and GHI index from bam/")
            .arg(
                Arg::new("range")
                    .short('r')
                    .takes_value(true)
                    .multiple(true)
                    .about("sorted bam"),
            )
            .arg(
                Arg::new("bed")
                    .short('L')
                    .takes_value(true)
                    .about("sorted bam"),
            )
            .arg(
                Arg::new("compression")
                    .short('c')
                    .takes_value(true)
                    .about("compression level"),
            )
            .arg(
                Arg::new("type")
                    .short('t')
                    .takes_value(true)
                    // .default_value("default")
                    .possible_values(&["alignment", "range", "default"])
                    .about("annotation type to fetch"),
            )
            .arg(
                Arg::new("id")
                    .short('i')
                    .takes_value(true)
                    .multiple(true)
                    .about("annotation sample to fetch"),
            )
            .arg(Arg::new("filter").short('f').about("Pre-filter"))
            .arg(Arg::new("binary").short('b').about("Binary"))
            /*                .arg(Arg::new("vis").short('v').about("Binary"))
            .arg(Arg::new("no-cigar").short('n').about("Binary"))
            .arg(Arg::new("packing").short('p').about("Binary"))
            .arg(Arg::new("legend").short('l').about("Legend"))
            .arg(Arg::new("quality").short('q').about("Quality"))
            .arg(Arg::new("x").short('x').takes_value(true).about("x"))
            .arg(Arg::new("y").short('y').takes_value(true).about("y"))
            .arg(
                Arg::new("freq-height")
                    .short('Y')
                    .takes_value(true)
                    .about("yf"),
            )
            .arg(
                Arg::new("max-coverage")
                    .short('m')
                    .takes_value(true)
                    .about("Max coverage"),
            )
            .arg(Arg::new("split-alignment").short('s').about("No-md"))
            .arg(Arg::new("no-insertion").short('z').about("No-md"))
            .arg(Arg::new("only-split-alignment").short('u').about("No-md"))
            .arg(
                Arg::new("hide-alignment")
                    .short('A')
                    .about("Hide alignment"),
            )
            .arg(Arg::new("all-bases").short('B').about("Show all bases"))
            .arg(Arg::new("pileup").short('P').about("Show coverage plot"))
            .arg(
                Arg::new("sort-by-name")
                    .short('N')
                    .about("Sort-by-name (read-id)"),
            )
            .arg(
                Arg::new("sort-by-cigar")
                    .short('C')
                    .about("Display split alignment name"),
            )
            // .arg(Arg::new("insertion").short('').about("No-md"))
            .arg(
                Arg::new("graph")
                    .short('g')
                    .takes_value(true)
                    .about("graph"),
            )*/
            .arg(
                Arg::new("output")
                    .short('o')
                    .takes_value(true)
                    .about("Output format"),
            )
            .arg(
                Arg::new("INPUT")
                    .about("Sets the input file to use")
                    .required(true)
                    .index(1),
            ),
    )
    .subcommand(
        App::new("decompose")
            .about("extract bam/bed from GHB and GHI index: (*) Now it needs to ")
            .arg(
                Arg::new("id")
                    .short('i')
                    .takes_value(true)
                    //            .multiple(true)
                    .about("annotation sample to fetch (alignment | annotation)"),
            )
            .arg(Arg::new("header").short('z').about("Outputs only header"))
            .arg(
                Arg::new("INPUT")
                    .about("Sets the input file to use")
                    .required(true)
                    .index(1),
            )
            .arg(
                Arg::new("output")
                    .short('o')
                    .takes_value(true)
                    .about("Output format"),
            ),
    )
    .subcommand(
        App::new("bin")
            .about("extract bam/bed from GHB and GHI index: (*) Now it needs to ")
            .arg(
                Arg::new("id")
                    .short('i')
                    .takes_value(true)
                    //            .multiple(true)
                    .about("annotation sample to fetch (alignment | annotation)"),
            )
            .arg(
                Arg::new("INPUT")
                    .about("Sets the input file to use")
                    .required(true)
                    .index(1),
            )
            .arg(
                Arg::new("bin")
                    .short('b')
                    .takes_value(true)
                    //            .multiple(true)
                    .about("annotation sample to fetch (alignment | annotation)"),
            )
            .arg(
                Arg::new("start")
                    .short('s')
                    .takes_value(true)
                    //            .multiple(true)
                    .about("annotation sample to fetch (alignment | annotation)"),
            )
            .arg(
                Arg::new("end")
                    .short('e')
                    .takes_value(true)
                    //            .multiple(true)
                    .about("annotation sample to fetch (alignment | annotation)"),
            )
            .arg(
                Arg::new("ref_id")
                    .short('c')
                    .takes_value(true)
                    //            .multiple(true)
                    .about("annotation sample to fetch (alignment | annotation)"),
            )
            .arg(
                Arg::new("range")
                    .short('r')
                    .takes_value(true)
                    .multiple(true)
                    .about("sorted bam"),
            )
            .arg(
                Arg::new("layer")
                    .short('l')
                    .takes_value(true)
                    .about("sorted bam"),
            )
            .arg(
                Arg::new("output")
                    .short('o')
                    .takes_value(true)
                    .about("Output format"),
            ),
    )
    .subcommand(
        App::new("server")
            .about("start the web server.")
            .arg(Arg::new("web").short('w').about("bind host and port"))
            .arg(
                Arg::new("range")
                    .short('r')
                    .takes_value(true)
                    .multiple(true)
                    .about("sorted bam"),
            )
            .arg(
                Arg::new("max-coverage")
                    .short('m')
                    .takes_value(true)
                    .about("Max coverage value on coverage track"),
            )
            .arg(
                Arg::new("bam")
                    .short('a')
                    .takes_value(true)
                    .multiple(true)
                    .about("sorted bam"),
            )
            .arg(
                Arg::new("INPUT")
                    .about("Sets the input file to use")
                    .required(true)
                    .index(1),
            ),
    )
    .subcommand(
        App::new("vis")
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
            ),
    )
    .subcommand(
        App::new("split")
            .about("start the web server.")
            .arg(
                Arg::new("output")
                    .short('o')
                    .takes_value(true)
                    .about("Output format"),
            )
            .arg(
                Arg::new("INPUT")
                    .about("Sets the input file to use")
                    .required(true)
                    .index(1),
            ),
    );
    let prefetch_max = 25000000;
    let criteria = 1000000;
    let x_width = prefetch_max / criteria * (740);
    // Here 2**17 < x_width < 2**18 so maxZoom should be 18+ 1;
    let max_zoom = 19;
    let freq_y = 100;
    let y = 40;
    let path = path * (max_zoom - zoom);
    if criteria * (19 - zoom) > 4000000 {

    }
    let mut b: Vec<String> = vec!["-Y".to_string(), (freq_y / (max_zoom-zoom)).to_string(), "-y".to_string(), (y / (max_zoom-zoom)).to_string()];
    b.extend(args.clone());
    let matches = app.get_matches_from(b);
    // let args: Vec<String> = args.into_iter().chain(b.into_iter()).collect(); 
    (matches, StringRegion::new(&format!("{}:{}-{}", range.path, criteria * path + 1, criteria * (path + 1)).to_string()).unwrap())
}

async fn index(data: web::Data<Arc<Vis>>, req: HttpRequest) -> Result<NamedFile> {
    let zoom: i64 = req.match_info().query("zoom").parse().unwrap();
    let path: i64 = req.match_info().query("filename").parse().unwrap();// .parse().unwrap();
    match NamedFile::open(format!("{}/{}_0.png", zoom, path)) {
        Ok(file) => Ok(file
            .use_last_modified(true)
            .set_content_disposition(ContentDisposition {
                disposition: DispositionType::Attachment,
                parameters: vec![],
            })),
        _ => {
            let data = &*data; //(*data).lock().unwrap(); //get_mut();
            let (matches, string_range) = id_to_range(&data.range, &data.args, zoom, path);
            let list = &data.list;
            let ann = &data.annotation;
            let freq = &data.freq;

            // let arg_vec = vec!["ghb", "vis", "-t", "1", "-r",  "parse"];
            bam_record_vis(&matches, string_range, list.to_vec(), ann.to_vec(), BTreeMap::new(), |_| None).unwrap();
            // bam_vis(matches, 1);
            Ok(NamedFile::open(format!("{}/{}_0.png", zoom, path))?
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

// #[derive(Copy)]
pub struct Vis {
    range: StringRegion, args: Vec<String>, list: Vec<(u64, Record)>, annotation: Vec<(u64, bed::Record)>, freq: BTreeMap<u64, Vec<(u64, u32)>>
}

impl Vis {
    fn new(range: StringRegion, args: Vec<String>, list: Vec<(u64, Record)>, annotation: Vec<(u64, bed::Record)>, freq: BTreeMap<u64, Vec<(u64, u32)>>) -> Vis {
        Vis{range, args, list, annotation, freq}
    }
}

#[actix_rt::main]
pub async fn server(matches: &'static ArgMatches, range: StringRegion, args: Vec<String>, list: Vec<(u64, Record)>, annotation: Vec<(u64, bed::Record)>, freq: BTreeMap<u64, Vec<(u64, u32)>>, threads: u16) -> std::io::Result<()> {

    use actix_web::{web, HttpServer};

    // Create some global state prior to building the server
    //#[allow(clippy::mutex_atomic)] // it's intentional.
    //let counter1 = web::Data::new(Mutex::new((matches.clone(), range, list, annotation, freq)));
    HttpServer::new(move || {
        //let counter = Cell::new(Vis::new( range.clone(),args.clone(), list.clone(), annotation.clone(), freq.clone()));
        let counter = Arc::new(Vis::new(range.clone(), args.clone(), list.clone(), annotation.clone(), freq.clone()));

        actix_web::App::new().data(counter).route("/{zoom:.*}/{filename:.*}.png", web::get().to(index))
    })
        .bind("0.0.0.0:4000")?
        .workers(threads as usize)
        .run()
        .await
}
