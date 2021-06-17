extern crate log;
#[macro_use]
extern crate serde_derive;
pub mod subcommands;
#[cfg(feature = "web")]
pub mod server;
#[cfg(feature = "web")]
pub mod buffered_server;
#[cfg(feature = "web")]
pub mod rest_server;

use clap::{App, AppSettings, Arg, ArgSettings};
use env_logger;
use subcommands::*;
use std::env;

fn main() {
//    env::set_var("RUST_LOG", "actix_web=info");
    env_logger::init();
    let app = App::new("A hybrid genomic data visualization tool")
        // .setting(AppSettings::ArgsNegateSubcommands)
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .setting(AppSettings::UnifiedHelpMessage)
        .setting(AppSettings::ArgRequiredElseHelp)
        .setting(AppSettings::ColoredHelp)
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
                .setting(AppSettings::ArgRequiredElseHelp)
                .setting(AppSettings::ColoredHelp)
                .about("Constructs hybrid genome index from bam/bed files")
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
                .arg(Arg::new("header").short('z').about("Outputs only header"))
                .arg(Arg::new("formatted-header").short('f').about("Outputs formatted header"))
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
                .setting(AppSettings::ArgRequiredElseHelp)
                .setting(AppSettings::ColoredHelp)
                .about("Extracts bam/bed from hybrid genome index with ranged query")
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
                .arg(Arg::new("bench").short('B').about("For benchmarking"))
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
                .about("Extracts an entire bam/bed from a hybrid genome index. (*) used for debugging")
                .arg(
                    Arg::new("id")
                        .short('i')
                        .takes_value(true)
                        //            .multiple(true)
                        .about("annotation sample to fetch (alignment | annotation)"),
                )
                .arg(Arg::new("header").short('z').about("Outputs only header"))
                .arg(Arg::new("formatted-header").short('f').about("Outputs formatted header"))
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
                .about("Extracts bam/bed from a hybrid genome index. (*) used for debugging")
                .arg(
                    Arg::new("id")
                        .short('i')
                        .takes_value(true)
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
                .about("Starts the web server. (*) used for debugging")
                .arg(Arg::new("host").short('i').about("host"))
                .arg(Arg::new("port").short('p').about("port"))
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
                .setting(AppSettings::ArgRequiredElseHelp)
                .setting(AppSettings::ColoredHelp)
                .about("Visualizes genomic data e.g. alignments and annotations")
                .arg(
                    Arg::new("range")
                        .short('r')
                        .takes_value(true)
                        .multiple(true)
                        .required(true)
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
                    Arg::new("preset-color")
                        .short('#')
                        .possible_values(&["", "hgb", "igv"])
                        .about("Preset color scheme ['igv', 'hgb'] (default is hgb"),
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
                .arg(Arg::new("no-filter").short('f').about("Disable pre-filtering on loading BAM index (used for debugging)"))
                .arg(Arg::new("no-cigar").short('c').about("Do not show cigar string"))
                .arg(Arg::new("no-scale").short('S').about("Do not show y-axis scale and legends"))
                .arg(Arg::new("no-ruler").short('*').setting(ArgSettings::MultipleOccurrences).about("Do not show x-axis ruler"))
                .arg(Arg::new("no-packing").short('p').about("Disable read packing"))
                .arg(Arg::new("no-legend").short('l').about("Hide legend"))
                .arg(Arg::new("square").short('Q').about("Set square width (overwritten)"))
                .arg(Arg::new("quality").short('q').about("Display reads by quality value"))
                .arg(Arg::new("colored-by-name").short('n').about("Set read colors by read name"))
                .arg(Arg::new("read-per-line").short('1').about("Show one line one read on split-alignment mode"))
                .arg(Arg::new("read-per-two-range").short('2').about("Show one read per one line across two range on split-alignment mode"))
                .arg(Arg::new("x-as-range").short('9').about("Set X as the interval of the input range"))
                .arg(
                    Arg::new("range-index")
                        .short('8')
                        .takes_value(true)
                        .about("Choose a range index to visualize when multiple ranges are specified."),
                )
                .arg(
                    Arg::new("read-index")
                        .short('_')
                        .takes_value(true)
                        .about("Choose a read index to visualize of a single input file."),
                )
                .arg(
                    Arg::new("colored-by-motif")
                        .short('E')
                        .takes_value(true)
                        .default_value("C:T:CG")
                        .about("Colored by specified motif (for bisulfite sequencing)"),
                )
                .arg(
                    Arg::new("colored-by-tag")
                        .short('0')
                        .takes_value(true)
                        .default_value("")
                        .setting(ArgSettings::AllowEmptyValues)
                        .about("Colored by specified tags on read alignments"),
                )
                .arg(
                    Arg::new("separated-by-tag")
                        .short('j')
                        .takes_value(true)
                        .default_value("")
                        .setting(ArgSettings::AllowEmptyValues)
                        .about("Separated tracks by specified tags on read alignments"),
                )
                .arg(
                    Arg::new("separated-by-tag-offset")
                        .short('k')
                        .takes_value(true)
                        .default_value("3")
                        .about("The maximal number of tracks by specified tags on read alignments"),
                )
                .arg(Arg::new("x").short('x').takes_value(true).about("The width of image"))
                .arg(Arg::new("y").short('y').takes_value(true).about("The height of each read alignment"))
                .arg(Arg::new("web").short('w').takes_value(true).about("Serve the web server"))
                .arg(Arg::new("whole-chromosome").short('W').about("Pretend as if the prefetch range is the whole chromosome"))
                .arg(Arg::new("rest").short('>').about("Serve the web server with accepting any parameter"))
                .arg(Arg::new("production").short('$').about("Serve the web server on production mode (no cross-origin request is allowed)"))
                .arg(Arg::new("dump-json").short('%').about("Dump JSON of read metadata"))
                .arg(Arg::new("adjust-y").short('&').about("Do not adjust y on server mode"))
                .arg(Arg::new("ref-column").short('!').takes_value(true).about("Show the base colors of reference genome (input must be 2bit format)"))
                .arg(Arg::new("insertion-string").short('{').about("Show the insertion sequence along with insertion callets"))
                .arg(
                    Arg::new("labels")
                        .short('}')
                        .takes_value(true)
                        .about("Labels displayed on legend"),
                )
                .arg(Arg::new("dynamic-partition").short('D').about("Divide multiple genomic range with dynamic width"))
                .arg(
                    Arg::new("format")
                        .short('O')
                        .takes_value(true)
                        .about("Output format (automatically detected; optionally used for server mode)"),
                )
                .arg(
                    Arg::new("show-read-id")
                        .short('H')
                        .about("Write a read id on the beginning of each read"),
                )
                .arg(
                    Arg::new("overlapping-reads")
                        .short('7')
                        .about("Output a list of read id overlapping the range"),
                )
                .arg(
                    Arg::new("snp-frequency")
                        .short('V')
                        .takes_value(true)
                        .about("The portion of allele frequency to display on each coverage track"),
                )
                .arg(
                    Arg::new("zoom-range")
                        .short('Z')
                        .takes_value(true)
                        .about("Permitted zoom range on web server mode"),
                )
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
                        .about("Maximum coverage value on coverage track"),
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
                        .about("Minimum read mapping length on coverage/alignment track"),
                )
                .arg(
                    Arg::new("x-scale")
                        .short('X')
                        .takes_value(true)
                        .about("Size of x-scale legend"),
                )
                .arg(Arg::new("split-alignment").short('s').about("Display split alignments in the same line"))
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
                .arg(Arg::new("no-deletion").short('+').about("Hide deletion callets on read alignments"))
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
                .arg(Arg::new("with-caption").short('<').about("Show caption on the top of chart"))
                .arg(Arg::new("end-split-callets").short('e').about("Show callets on ends of read alignments if the read contains split-alignment"))
                .arg(Arg::new("output-translocation").short('5').about("Write a list of translocation split-alignment to stdout"))
                .arg(Arg::new("translocation-target").short('6').takes_value(true).about("Set a translocation target chromosome"))
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
                ),
        )
        .subcommand(
            App::new("split")
                .about("Splits a hybrid genome index (not implemented).")
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
                )
                .arg(Arg::new("header").short('z').about("Outputs only header"))
                .arg(Arg::new("formatted-header").short('f').about("Outputs formatted header")),
        );
    let matches = app.get_matches();
    let args: Vec<String> = env::args().collect();
    let threads = matches
        .value_of("threads")
        .and_then(|t| t.parse::<u16>().ok())
        .unwrap_or(1u16);

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads as usize - 1)
        .build_global()
        .unwrap();

    if let Some(ref matches) = matches.subcommand_matches("build") {
        build(matches, threads);
    } else if let Some(ref matches) = matches.subcommand_matches("query") {
        match matches.is_present("binary") {
            true => bam_query(matches, threads),
            false => match matches.is_present("bench") {
                true => bench_query(matches, args, threads),
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
        // eprintln!("{:?}", matches.is_present("INPUT"));
        match matches.is_present("INPUT") {
            true => vis_query(matches, args, threads).unwrap(),
            false => bam_vis(matches, args, threads).unwrap(),
        }
    } else if let Some(ref _matches) = matches.subcommand_matches("server") {
    }
}
