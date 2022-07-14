extern crate log;
#[macro_use]
#[cfg(feature = "web")]
extern crate serde_derive;
#[cfg(feature = "web")]
pub mod buffered_server;
#[cfg(feature = "web")]
pub mod rest_server;
#[cfg(feature = "web")]
pub mod server;

//#[cfg(feature = "web")]
//pub mod bam_server;

pub mod subcommands;

use clap::{App, AppSettings, Arg, ArgSettings};

use std::env;
use subcommands::*;

fn main() {
    //    env::set_var("RUST_LOG", "actix_web=info");
    env_logger::init();
    let app = App::new("A hybrid genomic data visualization tool")
        // .setting(AppSettings::ArgsNegateSubcommands)
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .setting(AppSettings::UnifiedHelpMessage)
        .setting(AppSettings::ArgRequiredElseHelp)
        .setting(AppSettings::ColoredHelp)
        .version("0.3")
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
                .arg(Arg::new("yy").short('y').about("Calculate y coordinate and attach as a YY tag"))
                .arg(Arg::new("formatted-header").short('f').about("Outputs formatted header"))
                .arg(
                    Arg::new("max-coverage")
                        .short('m')
                        .takes_value(true)
                        .about("Max coverage value on coverage track"),
                )
                .arg(
                    Arg::new("no-bits")
                        .short('(')
                        .long("flag-exclude")
                        .takes_value(true)
                        .default_value("1796")
                        .about("Read must have NONE of these flags"),
                )
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
                .about("Visualizes genomic data, e.g. read alignments and annotations")
                .arg(
                    Arg::new("range")
                        .short('r')
                        .long("range")
                        .takes_value(true)
                        .multiple(true)
                        .required(true)
                        .about("Genomic range to visualize. Format is chr:from-to"),
                )
                .arg(
                    Arg::new("prefetch-range")
                        .long("prefetch-range")
                        .short('R')
                        .takes_value(true)
                        .multiple(true)
                        .about("(Optional) Genomic range to pre-fetch. The same order as range option. Format is chr:from-to"),
                )
                .arg(
                    Arg::new("type")
                        .long("ghb-type")
                        .short('t')
                        .takes_value(true)
                        // .default_value("default")
                        .possible_values(&["alignment", "range", "default"])
                        .about("The type of GHB data structure to display (used for debugging)"),
                )
                .arg(
                    Arg::new("preset")
                        .short('z')
                        .long("preset")
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
                        .long("preset-color")
                        .possible_values(&["", "hgb", "igv", "jbrowse"])
                        .about("Preset color scheme ['igv', 'hgb', 'jbrowse'] (default is hgb"),
                )
                .arg(
                    Arg::new("id")
                        .short('i')
                        .long("id")
                        .takes_value(true)
                        .multiple(true)
                        .about("Indices on GHB data structure to display (used for debugging)"),
                )
                .arg(
                    Arg::new("frequency")
                        .short('F')
                        .long("frequency")
                        .takes_value(true)
                        .multiple(true)
                        .about("[Input] A subset of sorted bed for coverage plot (start and score fields are used)"),
                )
                .arg(Arg::new("bed-range").short('J').long("bed-range").takes_value(true).about("BED file to specify multiple regions to display"))
                .arg(Arg::new("neighbor").short('K').long("bed-neighbor-bases").takes_value(true).about("Visualizes specified base-pair neighbor of BED region"))
                .arg(Arg::new("no-filter").short('f').long("disable-read-prefilter").about("Disables pre-filtering on loading BAM index (used for debugging)"))
                .arg(Arg::new("no-cigar").short('c').long("hide-cigar").about("Do not show cigar string"))
                .arg(Arg::new("no-scale").short('S').long("hide-y-scale").about("Do not show y-axis scale and legends"))
                .arg(Arg::new("no-ruler").short('*').long("hide-x-scale").setting(ArgSettings::MultipleOccurrences).about("Do not show x-axis ruler"))
                .arg(Arg::new("no-packing").short('p').long("disable-read-packing").about("Disables read packing"))
                .arg(Arg::new("no-legend").short('l').long("hide-legend").about("Hides legend"))
                .arg(Arg::new("square").short('Q').long("square").about("Sets width as the same as height to generate s square image (overwritten)"))
                .arg(Arg::new("quality").short('q').long("read-quality").about("Displays reads colored by quality value"))
                .arg(Arg::new("colored-by-name").short('n').long("colored-by-name").about("Sets read colors by read name"))
                .arg(Arg::new("colored-by-track").short('?').long("colored-by-track").about("Sets read colors by track on hgb"))
                .arg(Arg::new("read-per-line").short('1').long("read-per-line").about("Shows a read per line on split-alignment mode"))
                .arg(Arg::new("read-per-two-range").short('2').long("read-per-line-two-ranges").about("Shows a read per line across two ranges on split-alignment mode"))
                .arg(Arg::new("x-as-range").short('9').long("input-range-is-width").about("Sets the genomic interval of the input range as width pixels instead of x parameter"))
                .arg(
                    Arg::new("range-index")
                        .short('8')
                        .long("range-index")
                        .takes_value(true)
                        .about("Range index to visualize when multiple ranges are specified"),
                )
                .arg(
                    Arg::new("read-index")
                        .short('_')
                        .long("read-index")
                        .takes_value(true)
                        .about("Read index to visualize of a single input file"),
                )
                .arg(
                    Arg::new("colored-by-motif")
                        .short('E')
                        .long("colored-by-motif")
                        .takes_value(true)
                        .default_value("C:T:CG")
                        .about("Sequence motif to be colored (for bisulfite sequencing)"),
                )
                .arg(
                    Arg::new("colored-by-tag")
                        .short('0')
                        .long("colored-by-tag")
                        .takes_value(true)
                        .default_value("")
                        .setting(ArgSettings::AllowEmptyValues)
                        .about("Tags on read alignments to be colored"),
                )
                .arg(
                    Arg::new("separated-by-tag")
                        .short('j')
                        .long("grouped-by-tag")
                        .takes_value(true)
                        .default_value("")
                        .setting(ArgSettings::AllowEmptyValues)
                        .about("Tracks are grouped by specified tags on read alignments"),
                )
                .arg(
                    Arg::new("separated-by-tag-offset")
                        .short('k')
                        .long("grouped-by-tag-offset")
                        .takes_value(true)
                        .default_value("3")
                        .about("The maximal number of tracks by a specified tag on read alignments"),
                )
                .arg(
                    Arg::new("filtered-by-tag")
                        .short('~')
                        .long("filtered-by-tag")
                        .takes_value(true)
                        .about("Tag on read alignments as <tag>:<value> to be filtered by (e.g. HP:0)"),
                )
                .arg(Arg::new("border-height").short('^').takes_value(true).about("The height of border between samples"))
                .arg(Arg::new("x").short('x').takes_value(true).about("The width of image"))
                .arg(Arg::new("y").short('y').takes_value(true).about("The height of each read alignment"))
                .arg(Arg::new("web").short('w').takes_value(true).long("web-server").about("Web server mode running on TCP socket (host:port)"))
                .arg(Arg::new("whole-chromosome").short('W').long("whole-chromosome").about("Pretends as if the prefetch range is the whole chromosome"))
                .arg(Arg::new("rest").short('>').long("rest-server").about("Serves a web server that accepts any parameters"))
                .arg(Arg::new("production").short('$').long("serve-as-production").about("Serves a web server on production mode (no cross-origin request is allowed)"))
                .arg(Arg::new("dump-json").short('%').long("write-json").about("Dumps JSON of read metadata"))
                .arg(Arg::new("adjust-y").short('&').long("not-adjust-y").about("Do not adjust y parameter on server mode"))
                .arg(Arg::new("ref-column").short('!').long("2bit").takes_value(true).about("[Input] 2 bit file of reference genome to display the base colors of reference genome"))
                .arg(Arg::new("insertion-string").short('{').long("show-insertion-sequence").about("Shows insertion sequences along with insertion symbols"))
                .arg(
                    Arg::new("labels")
                        .short('}')
                        .long("labels")
                        .takes_value(true)
                        .multiple(true)
                        .about("Labels displayed on legend"),
                )
                .arg(Arg::new("dynamic-partition").short('D').long("dynamic-range-partition").about("Divides multiple genomic ranges of fields with dynamic width"))
                .arg(
                    Arg::new("format")
                        .short('O')
                        .long("output-format")
                        .takes_value(true)
                        .about("Output format (automatically detected; optionally used for server mode)"),
                )
                .arg(
                    Arg::new("show-read-id")
                        .short('H')
                        .long("show-read-ids")
                        .about("Writes a read id on the beginning of each read"),
                )
                .arg(
                    Arg::new("no-bold-line")
                        .short('7')
                        .long("no-bold-line")
                        .about("Draws a mesh without bold lines"),
                )
                .arg(
                    Arg::new("snp-frequency")
                        .short('V')
                        .long("heterozygous-frequency")
                        .takes_value(true)
                        .about("The portion of heterozygous allele frequency to display on each coverage track"),
                )
                .arg(
                    Arg::new("zoom-range")
                        .short('Z')
                        .long("zoom-range")
                        .takes_value(true)
                        .about("Maximum zoom level of DZI format on web server mode"),
                )
                .arg(
                    Arg::new("freq-height")
                        .short('Y')
                        .long("y-scale")
                        .takes_value(true)
                        .about("The height of each coverage track"),
                )
                .arg(
                    Arg::new("max-coverage")
                        .short('m')
                        .long("max-coverage")
                        .takes_value(true)
                        .about("Maximum coverage value on coverage track"),
                )
                .arg(
                    Arg::new("cache-dir")
                        .short('d')
                        .long("cache-directory")
                        .takes_value(true)
                        .about("Cache directory for web server (generated randomly if not specified)"),
                )
                .arg(
                    Arg::new("static-dir")
                        .short('.')
                        .long("static-directory")
                        .default_value("./static")
                        .takes_value(true)
                        .about("Static serve directory for web server")
                )
                .arg(
                    Arg::new("min-read-length")
                        .short('M')
                        .long("min-read-length")
                        .takes_value(true)
                        .about("Minimum read mapping length on coverage/alignment track"),
                )
                .arg(
                    Arg::new("x-scale")
                        .short('X')
                        .long("x-scale")
                        .takes_value(true)
                        .about("Size of legends on x-scale"),
                )
                .arg(Arg::new("split-alignment").short('s').long("align-reads-horizontally").about("Displays split alignments in the same line"))
                .arg(Arg::new("only-split-alignment").short('u').long("only-split-alignments").about("Displays only split alignments or mate-paired reads on alignment track"))
                .arg(Arg::new("exclude-split-alignment").short('3').long("only-full-length-alignments").about("Displays only NOT split alignments or mate-paired reads on alignment track"))
                .arg(Arg::new("full-length").short('4').long("full-length-reads").about("Displays only full-length match reads against prefetch range"))
                .arg(Arg::new("udon").short('U').long("udon").about("Colors by udon library"))
                .arg(
                    Arg::new("sort-by-name")
                        .short('N')
                        .long("sort-by-name")
                        .about("Sorts alignments by read id (for split-alignment visualization)"),
                )
                .arg(
                    Arg::new("sort-by-cigar")
                        .short('C')
                        .long("show-split-alignment-line")
                        .about("Displays each split-alignment connection as a colored line"),
                )
                .arg(Arg::new("no-insertion").short('I').long("hide-insertion").about("Hides insertion symbols on read alignments"))
                .arg(Arg::new("no-deletion").short('+').long("hide-deletion").about("Hides deletion symbols on read alignments"))
                .arg(
                    Arg::new("hide-alignment")
                        .short('A')
                        .long("only-coverage-plot")
                        .about("Hides read alignments (display only coverage; for chromosome scale view)"),
                )
                .arg(Arg::new("all-bases").short('B').long("all-bases").about("Shows all nucleotides with color"))
                .arg(
                    Arg::new("pileup")
                        .short('P')
                        .long("coverage-plot")
                        .about("Shows both pileup and coverage plot"),
                )
                .arg(Arg::new("only-translocation").short('T').long("show-translocation-callets").about("Shows symbols on ends of read alignments if the read contains translocational split-alignment"))
                .arg(Arg::new("with-caption").short('<').long("caption").takes_value(true).setting(ArgSettings::AllowEmptyValues)
                .default_value("").about("Caption on the top of chart"))
                .arg(Arg::new("end-split-callets").short('e').long("show-split-alignment-callets").about("Shows symbols on ends of read alignments if the read contains split-alignment"))
                .arg(Arg::new("output-translocation").short('5').long("write-split-alignment").about("Writes a list of translocation split-alignment to stdout"))
                .arg(Arg::new("translocation-target").short('6').long("translocation-target-chromosome").takes_value(true).about("Target chromosome on writing a list of translocation"))
                .arg(
                    Arg::new("graph")
                        .short('G')
                        .long("graph")
                        .takes_value(true)
                        .about("[Input] Graph genome coordinates tsv (tab-separated, generated by vg view -N command)"),
                )
                .arg(
                    Arg::new("basic-auth")
                        .short('[')
                        .long("basic-auth")
                        .takes_value(true)
                        .about("Basic authentication on web server (username:password)"),
                )
                .arg(
                    Arg::new("unix-socket")
                        .short(']')
                        .long("unix-socket")
                        .takes_value(true)
                        .about("Web server mode running on unix socket"),
                )
                .arg(
                    Arg::new("no-bits")
                        .short('(')
                        .long("flag-exclude")
                        .takes_value(true)
                        .default_value("1796")
                        .about("Read must have NONE of these flags"),
                )
                .arg(
                    Arg::new("read-name")
                        .short(')')
                        .long("filtered-by-read-name")
                        .takes_value(true)
                        .about("Read name to filter by"),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .takes_value(true)
                        .about("[Output] Image file (prefixed as .bmp / .png)"),
                )
                .arg(
                    Arg::new("bed")
                        .short('L')
                        .long("bed")
                        .multiple(true)
                        .takes_value(true)
                        .about("[Input] A subset of sorted bed to display as annotation track"),
                )
                .arg(
                    Arg::new("bam")
                        .short('a')
                        .long("bam")
                        .takes_value(true)
                        .multiple(true)
                        .about("[Input] Sorted BAM to display read alignment track"),
                )
                .arg(
                    Arg::new("gff3")
                        .short('g')
                        .long("gff3")
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
            App::new("precomp")
                .about("Attach y positions on each read by default view options.")
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
    } else if let Some(ref matches) = matches.subcommand_matches("precomp") {
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
