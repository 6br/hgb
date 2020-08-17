# hgb

A hybrid genome browser for zooming out long-read alignments

You can get these images in seconds.

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
cd hgb
cargo build --release
cargo run --release -- vis --help
```

## Usage

The input BAM file **must** be indexed using `samtools index`. The input BAM file *needed to* be calculated MD tag using `samtools calmd` if mismatches are to be visualized.

* A simple example:

```bash
cargo run --release -- vis -a input.bam -o output.png -r chr1:91234-92334
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

* A simple examples with two bam files without cigars:

```bash
cargo run --release -- vis -a input1.bam input2.bam -o output.png -r chr1:91234-92334 -n -I
```

## Author

6br
