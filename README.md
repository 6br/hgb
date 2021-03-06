# hgb

[![Build Status](https://travis-ci.org/6br/hgb.svg?branch=master)](https://travis-ci.org/6br/hgb)

A hybrid genome browser for zooming out long-read alignments

Screenshots in seconds.

![alignments](fig/alignments.png)

Integrated with [Udon](https://github.com/ocxtal/udon):

![udon](fig/udon.png)

## Features

* A light-weight binary to visualize read alignments as PNG/JPG/BMP files.
* Visualize read alignments more than 100 samples at once.
* Much more options to visualize read alignments.

## Install

If you didn't install Rust/cargo, you need to install rustup.

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

```bash
git clone https://github.com/6br/hgb
cd hgb
cargo build --release
cargo run --release -- vis --help
```

## Usage

The input BAM file **must** be indexed using `samtools index`. The input BAM file *needed to* be calculated MD tag using `samtools calmd` if mismatches are to be visualized.

```bash
samtools calmd -b aln.bam ref.fasta > input.bam
samtools index input.bam
```

* A simple example:

```bash
cargo run --release -- vis -a input.bam -o output.png -r chr1:91234-92334
```

* A simple example with [Udon](https://github.com/ocxtal/udon) mode:

```bash
cargo run --release -- vis -a input.bam -o output.png -l -U -r chr1:91234-92334
```

* A simple examples with two bam files:

```bash
cargo run --release -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334
```

* A simple examples with two bam files of split alignments:

```bash
cargo run --release -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -s
```

* A simple examples with two bam files of only split alignments:

```bash
cargo run --release -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -s -u
```

* A simple examples with two bam files with read quality:

```bash
cargo run --release -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -q
```

* A simple examples with two bam files with hidden legends:

```bash
cargo run --release -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -l
```

* A simple examples with two bam files within 30x coverages at most:

```bash
cargo run --release -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -m 30
```

* A simple examples with two bam files without insertion cigars:

```bash
cargo run --release -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -I
```

## Gallery

![large](fig/large.png)

## Author

6br
