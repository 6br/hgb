extern crate log;
pub mod subcommands;

use clap::{App, Arg};
use env_logger;
use subcommands::*;

fn main() {
    env_logger::init();
    let matches = App::new("GHB/GHI genomic data visualization tool")
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
                .arg(Arg::new("host").short('i').about("host"))
                .arg(Arg::new("port").short('p').about("port"))
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
                    Arg::new("bed")
                        .short('L')
                        .multiple(true)
                        .takes_value(true)
                        .about("A subset of sorted bed to display as annotation track"),
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
                .arg(Arg::new("no-packing").short('p').about("Disable read packing"))
                .arg(Arg::new("legend").short('l').about("Show legend"))
                .arg(Arg::new("quality").short('q').about("Display reads by quality value"))
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
                .arg(Arg::new("only-split-alignment").short('u').about("No-md"))
                .arg(
                    Arg::new("sort-by-name")
                        .short('N')
                        .about("Sort-by-name (read-id)"),
                )
                .arg(
                    Arg::new("sort-by-cigar")
                        .short('C')
                        .about("Sort-by-name (read-id)"),
                )
                .arg(Arg::new("no-insertion").short('z').about("No insertion"))
                .arg(
                    Arg::new("hide-alignment")
                        .short('A')
                        .about("Hide alignment"),
                )
                .arg(Arg::new("all-bases").short('B').about("Show all bases"))
                .arg(
                    Arg::new("pileup")
                        .short('P')
                        .about("Show pileup as coverage plot"),
                )
                .arg(Arg::new("only-translocation").short('T').about("Show all bases"))
                .arg(
                    Arg::new("graph")
                        .short('c')
                        .takes_value(true)
                        .about("graph csv (tab-separated)"),
                )
                // .arg(Arg::new("insertion").short('').about("No-md"))
                .arg(
                    Arg::new("output")
                        .short('o')
                        .takes_value(true)
                        .about("Output format"),
                )
                .arg(
                    Arg::new("bam")
                        .short('a')
                        .takes_value(true)
                        .multiple(true)
                        .about("sorted bam"),
                )
                .arg(
                    Arg::new("gff3")
                        .short('g')
                        .takes_value(true)
                        .multiple(true)
                        .about("sorted bed"),
                )
                .arg(
                    Arg::new("INPUT")
                        .about("Sets the input file to use")
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
        )
        .get_matches();

    let threads = matches
        .value_of("threads")
        .and_then(|t| t.parse::<u16>().ok())
        .unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads as usize - 1)
        .build_global()
        .unwrap();

    if let Some(ref matches) = matches.subcommand_matches("build") {
        build(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("query") {
        match matches.is_present("binary") {
            true => bam_query(matches, threads),
            false => match matches.is_present("vis") {
                true => vis_query(matches, threads).unwrap(),
                false => query(matches, threads),
            },
        }
    } else if let Some(ref matches) = matches.subcommand_matches("decompose") {
        decompose(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("split") {
        split(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("bin") {
        bin(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("vis") {
        eprintln!("{:?}", matches.is_present("input"));
        match matches.is_present("input") {
            true => vis_query(matches, threads).unwrap(),
            false => bam_vis(matches, threads).unwrap(),
        }
    }
}
