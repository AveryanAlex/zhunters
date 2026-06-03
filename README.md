# ZHunters

`zhunters` is a Rust port of the Z-HUNT 3 Z-DNA scanner. It provides the `zhunt` command-line tool for scanning DNA sequences and writing legacy-compatible `.Z-SCORE` output files.

This project is based on the original Z-HUNT implementation by Ho Lab, Colorado State University:

https://github.com/Ho-Lab-Colostate/zhunt

It also incorporates correctness and optimization ideas from Carlos Bederián's Fast Z-Hunt implementation:

https://github.com/zzzoom/fast-zhunt

## Changes from the original Z-HUNT

The command-line interface, input handling, and `.Z-SCORE` output format are kept compatible with the original scanner, but this port includes a few correctness and usability updates:

- **Physical constants:** `RT` is computed from the CODATA molar gas constant using the `dimensional_quantity` crate and an explicit 25 °C temperature (`298.15 K`). This replaces the older rounded `298 K`/hand-written constant calculation while keeping the Z-HUNT energy model and DBZED table intact.
- **Anti/syn scoring fix:** anti/syn path selection uses exact integer centi-kcal energy sums instead of accumulating recursive `float` additions/subtractions. This avoids the numerical drift described in Ho-Lab-Colostate/zhunt issue [#9](https://github.com/Ho-Lab-Colostate/zhunt/issues/9).
- **Anti/syn dynamic programming optimization:** exact anti/syn scoring is evaluated with a compact two-state dynamic program that keeps only the best prefix ending in `AS` and the best prefix ending in `SA`. This removes per-position hash-map/frontier bookkeeping while preserving the corrected integer scoring behavior.
- **Slope fix:** the reported slope is computed from the log coefficients for the actual best dinucleotide length, not whatever length happened to run last. This fixes the wrong-`logcoef` bug described in Ho-Lab-Colostate/zhunt issue [#10](https://github.com/Ho-Lab-Colostate/zhunt/issues/10).
- **Progress reporting:** long scans show a progress bar while positions are scored and streamed to disk.
- **Rayon multithreading:** scoring is parallelized with Rayon work stealing. Results are streamed as ordered score blocks so the writer can consume completed work without waiting for a whole large compute chunk. By default, systems with 8 or more logical CPUs reserve one CPU for writing/progress work and use `cores - 1` scoring workers; smaller systems use all available CPUs. `--threads` overrides this default.

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
zhunt [--threads <threads>] [-o <output>] <windowsize> <minsize> <maxsize> <datafile>
```

Example:

```bash
target/release/zhunt 12 8 12 input.fa
```

Use a fixed number of worker threads:

```bash
target/release/zhunt --threads 8 12 8 12 input.fa
```

Write results to a custom path:

```bash
target/release/zhunt -o results.Z-SCORE 12 8 12 input.fa
```

By default, the output file is written next to the input file as:

```text
input.fa.Z-SCORE
```

## Arguments

- `windowsize`: maximum window size in dinucleotides
- `minsize`: minimum region size in dinucleotides
- `maxsize`: maximum region size in dinucleotides
- `datafile`: input DNA file
- `--threads <threads>`: optional number of Rayon scoring workers; defaults to available parallelism, except systems with 8 or more logical CPUs default to `cores - 1` to leave room for writing/progress work
- `-o, --output <output>`: optional output `.Z-SCORE` path; defaults to `<datafile>.Z-SCORE`

Input parsing follows the legacy behavior: the scanner reads all `A/T/G/C/N` bases from the file and ignores other bytes.

Note: the current physical-constant implementation uses `dimensional_quantity`, which requires nightly Rust.

## Development

Run formatting, linting, and tests:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --lib
cargo test --test cli
```

## Profiling

If `cargo-flamegraph` is installed, use the helper script:

```bash
./profile.sh 12 8 12 input.fa --threads 16
```

This writes `flamegraph.svg`.
