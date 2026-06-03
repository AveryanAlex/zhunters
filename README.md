# ZHunters

`zhunters` is a Rust port of the Z-HUNT 3 Z-DNA scanner. It provides the `zhunt` command-line tool for scanning DNA sequences and writing legacy-compatible `.Z-SCORE` output files.

This project is based on the original Z-HUNT implementation by Ho Lab, Colorado State University:

https://github.com/Ho-Lab-Colostate/zhunt

It also incorporates correctness and optimization ideas from Carlos Bederián's Fast Z-Hunt implementation:

https://github.com/zzzoom/fast-zhunt

## Benchmark

Local benchmark (Ryzen 7 5800X) on `NC_043715.1[1..59306649].fa` (59.3 Mbases) with `zhunt 12 8 12`:

| Scanner | Real time | User time | System time | CPU |
| --- | ---: | ---: | ---: | ---: |
| `zhunters` | `20.568s` | `251.41s` | `5.11s` | `1247%` |
| `fast-zhunt` | `1:36.23` | `710.58s` | `8.53s` | `747%` |

## Changes from the original Z-HUNT

The command-line interface, input handling, and `.Z-SCORE` output format are kept compatible with the original scanner. The main changes are grouped by whether they were added in this Rust port or adapted from Fast Z-Hunt.

### Rust-port improvements

- **Delta-linking candidate pruning:** candidate lengths that cannot beat the current best delta-linking value are skipped using the monotonicity of the delta-linking equation. This preserves the legacy result while avoiding unnecessary root searches.
- **Delta-linking root search acceleration:** root search uses a safeguarded Newton predictor, snaps to the legacy bisection grid, and verifies the final point with the original equation. The predictor also uses a factorized form of the delta-linking equation to reduce runtime `exp()` calls while preserving the verified legacy-grid result.
- **Rayon multithreading:** scoring is parallelized with Rayon work stealing. Results are streamed as ordered score blocks so the writer can consume completed work without waiting for a whole large compute chunk. By default, systems with 8 or more logical CPUs reserve one CPU for writing/progress work and use `cores - 1` scoring workers; smaller systems use all available CPUs. `--threads` overrides this default.
- **Progress reporting:** long scans show a progress bar while positions are scored and streamed to disk.
- **Worker-side output formatting:** score blocks are formatted into byte buffers by Rayon workers, then written in order by the writer thread. This keeps record formatting parallel and reduces writer-side bottlenecks.
- **Physical constants:** `RT` is computed from the CODATA molar gas constant using the `dimensional_quantity` crate and an explicit 25 °C temperature (`298.15 K`). This replaces the older rounded `298 K`/hand-written constant calculation while keeping the Z-HUNT energy model and DBZED table intact.

### Fast Z-Hunt-derived fixes and optimizations

- **Anti/syn scoring fix:** anti/syn path selection uses exact integer centi-kcal energy sums instead of accumulating recursive `float` additions/subtractions. This avoids the numerical drift described in Ho-Lab-Colostate/zhunt issue [#9](https://github.com/Ho-Lab-Colostate/zhunt/issues/9).
- **Slope fix:** the reported slope is computed from the log coefficients for the actual best dinucleotide length, not whatever length happened to run last. This fixes the wrong-`logcoef` bug described in Ho-Lab-Colostate/zhunt issue [#10](https://github.com/Ho-Lab-Colostate/zhunt/issues/10).
- **Anti/syn dynamic programming optimization:** exact anti/syn scoring is evaluated with a compact two-state dynamic program that keeps only the best prefix ending in `AS` and the best prefix ending in `SA`. This removes per-position hash-map/frontier bookkeeping while preserving the corrected integer scoring behavior. This optimization is based on the compact anti/syn DP used in Carlos Bederián's Fast Z-Hunt implementation.

## Installation

Install Rust, then build the release binary:

```bash
cargo build --release
```

The binary will be available at:

```bash
target/release/zhunt
```

### Jupyter/Colab

In a Jupyter or Google Colab notebook, run the following commands in a cell to install Rust, build `zhunt`, and make it available on `PATH`:

```bash
!curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
!rm -rf /tmp/zhunters && GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 https://github.com/AveryanAlex/zhunters.git /tmp/zhunters
!cd /tmp/zhunters && PATH="$PATH:/root/.cargo/bin" cargo build --release && cd -
!cp /tmp/zhunters/target/release/zhunt /usr/local/bin/zhunt && rm -rf /tmp/zhunters
!zhunt --help
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
