# Scoring optimization report

## Benchmark setup

- Input: `/home/alex/Documents/ФКН/bioinf/hw4-2026/NC_043715.1[1..59306649].fa`
- Arguments: `12 8 12`
- Output validation: full `.Z-SCORE` SHA-256 against the exact reference output
- Reference hash: `9107c5945925b3cd229fc088f8115819d7c7e00cfc7da43369b6467fb3324850`
- Reference output size: `5348384067` bytes

## Summary

The best accepted changes are exact-output optimizations. They preserve the full output hash and reduce the full benchmark from about `20.6s` after constant tuning to about `14.7s`.

| Version | Time | Output |
| --- | ---: | --- |
| Initial post-tuning exact best | ~`20.6s` | byte-identical |
| Final optimized exact path | `14.66s`–`14.71s` | byte-identical |

This is roughly a **29% time reduction** versus the previous exact best.

## Accepted exact optimizations

### 1. Delta-linking coefficient math

Previous delta-linking code stored `ln(sum(products))` coefficients, rebuilt exponent arrays, and used an exponent-offset path for equation evaluation.

The optimized version stores raw coefficient sums directly and evaluates the same equations as:

```text
weight = coefficient * exp(K_RT * z * z)
```

This removes repeated `ln()` work, removes the `exponents` scratch buffer, and removes the unused `EXP_LIMIT` path for the tested workload. The full benchmark remained byte-identical.

### 2. Remove endpoint checks in normal Newton search

The previous search evaluated `f_min` and `f_max` before using the Newton predictor. For this workload the normal path is monotone and valid, so the optimized code starts from the midpoint and lets the verified grid snap/fallback handle safety.

Accepted setting:

```rust
const DELTA_LINKING_NEWTON_STEPS: usize = 4;
const DELTA_LINKING_GRID_ADJUST_LIMIT: usize = 0;
```

This kept the output byte-identical and was faster than the previous `6` Newton steps.

### 3. Manual record formatting

`format_records` previously used `writeln!(&mut bytes, "{record}")`, which showed up clearly in profiling through generic float formatting.

The optimized version writes directly into the output byte buffer:

- integer fields via small decimal helpers,
- fixed-width 3-decimal floats via direct scaling,
- scientific probability format via a custom legacy-compatible path,
- anti/syn strings via direct byte appends.

An initial direct formatter changed `0.05045%` of lines due to float tie rounding (`x.xxx5` cases). Adding ties-to-even rounding restored the exact reference hash.

## Experiments that did not make it

### Reduced delta-linking grid size

Tried reducing `DELTA_LINKING_GRID_STEPS` below `65_536`.

Results:

- `8192` and `4096` grids still ran around `20.5s` before later algorithmic changes.
- Output hash changed.
- Output size changed.

Conclusion: grid size is not a hot loop count in the current code; reducing it mostly changes quantization/precision, not runtime.

### Lower Newton counts with verified snapping

With exact verified snapping:

- `NEWTON_STEPS=4` was good and exact, around `19.1s` before formatter/coefficient optimizations.
- `NEWTON_STEPS=3` stayed exact but was slower, around `21.3s`.
- Very low Newton counts in the older search setup were much slower (`~36s`–`40s`).

Conclusion: `4` is a good exact setting for the current path.

### Fast snap without verification

Skipping grid verification gave large speedups but changed output.

Examples:

| Variant | Time | Output |
| --- | ---: | --- |
| no-endpoint `n4` fast snap | ~`15.7s` before formatter | changed |
| manual formatter + `n4` fast snap | ~`13.1s` | changed |
| coefficient math + `n0` fast snap | ~`9.1s` | changed substantially |

Detailed comparison for one approximate `n4` fast-snap output:

- changed lines: `18.305%`
- changed winning length/path/sequence: `2.910%`
- max printed delta-linking difference: `10.0`

Conclusion: this can reach ~10s, but the precision/output loss is far above the requested approximate tolerance.

### Fast exponential approximations

Tried replacing predictor exponentials with approximate exponentials.

- Bit-level approximation: much slower (`~31s`) despite exact final hash in the tested path.
- Polynomial approximation: around `21s`, not an improvement.

Conclusion: not useful here.

### Cache best coefficients for slope

Tried storing the best candidate's coefficients instead of recomputing them for the final slope.

- Output stayed exact.
- Runtime was not improved (`~20.9s` in that experiment).

Conclusion: extra copying outweighed the saved recomputation.

## Profiling notes

Before formatting/coefficient optimizations, `perf` showed major time in:

- `exp` / `__ieee754_exp_fma`
- `log` / `__ieee754_log_fma`
- `delta_linking_equation_into`
- `format_records` / generic float formatting

After the accepted optimizations, formatting dropped to a small fraction of runtime and the remaining bottleneck is mostly scoring math, especially exponentials in delta-linking.

## Current changed files

- `src/scoring.rs`
  - constants tuned,
  - direct coefficient delta-linking math,
  - no endpoint checks in the normal Newton predictor path,
  - no `exponents` scratch buffer.
- `src/output.rs`
  - custom byte-level record formatter.
- `src/constants.rs`
  - removed unused `EXP_LIMIT`.
- `scripts/tune_scoring_constants.sh`
  - sed-based constants tuner.
- `scripts/tune_delta_linking_constants.sh`
  - sed-based delta-linking tuner with SHA checking.

## Validation

- `cargo fmt`
- `cargo test`: passed all tests
  - 11 lib tests
  - 8 CLI tests
  - 0 doc tests
- Final full benchmark after cleanup:
  - `14.713s`, byte-identical
  - repeat `14.664s`, byte-identical
