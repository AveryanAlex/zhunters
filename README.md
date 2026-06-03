# zhunters

`zhunters` is a Rust port of the Z-HUNT 3 Z-DNA scanner. It provides the `zhunt` command-line tool for scanning DNA sequences and writing legacy-compatible `.Z-SCORE` output files.

This project is based on the original Z-HUNT implementation by Ho Lab, Colorado State University:

https://github.com/Ho-Lab-Colostate/zhunt

It also incorporates correctness and optimization ideas from Carlos Bederián's Fast Z-Hunt implementation:

https://github.com/zzzoom/fast-zhunt

## Installation

Install Rust, then build the release binary:

```bash
cargo build --release
```

The binary will be available at:

```bash
target/release/zhunt
```

## Usage

```bash
zhunt [--threads <threads>] <windowsize> <minsize> <maxsize> <datafile>
```

Example:

```bash
target/release/zhunt 12 8 12 input.fa
```

Use a fixed number of worker threads:

```bash
target/release/zhunt --threads 8 12 8 12 input.fa
```

The output file is written next to the input file as:

```text
input.fa.Z-SCORE
```

## Arguments

- `windowsize`: maximum window size in dinucleotides
- `minsize`: minimum region size in dinucleotides
- `maxsize`: maximum region size in dinucleotides
- `datafile`: input DNA file
- `--threads <threads>`: optional number of scoring threads; defaults to available parallelism

Input parsing follows the legacy behavior: the scanner reads all `A/T/G/C/N` bases from the file and ignores other bytes.

## Development

Run formatting, linting, and tests:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --lib
cargo test --test cli
```

Generate a release build:

```bash
cargo build --release --bin zhunt
```

## Profiling

If `cargo-flamegraph` is installed, use the helper script:

```bash
./profile.sh 12 8 12 input.fa --threads 16
```

This writes `flamegraph.svg`.
