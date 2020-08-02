# hgb

A hybrid genome browser for zooming out long-read alignments

You can get this image in seconds with command line option.

![alignments](fig/alignments.png)

![large](fig/large.png)

## Feature

* A light-weight binary to visualize read alignments.
* Visualize read alignments more than 100 samples at once.
* Much more options to visualize read alignments.

## Install

If you didn't install Rust/cargo, you need to install rustup.

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

```bash
git clone https://github.com/6br/hgb
cargo build --release
cargo run -- vis --help
```

## Usage

* A simple example:

```bash
cargo run -- vis -a input.bam -o output.png -r chr1:91234-92334
```

* A simple examples with two bam files:

```bash
cargo run -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334
```

* A simple examples with two bam files of split alignments:

```bash
cargo run -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -s
```

* A simple examples with two bam files of only split alignments:

```bash
cargo run -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -p -s -u
```

* A simple examples with two bam files with read quality:

```bash
cargo run -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -q
```

* A simple examples with two bam files with hidden legends:

```bash
cargo run -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -l
```

* A simple examples with two bam files without cigars:

```bash
cargo run -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -n -I
```

## Author

6br
