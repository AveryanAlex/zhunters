#!/usr/bin/env bash
set -Eeuo pipefail

# Sed-based tuner for src/scoring.rs constants.
#
# Defaults benchmark the requested full FASTA with zhunt args `12 8 12`, write
# each multi-GB result file under OUT_DIR, and delete it after recording size.
# The current CLI intentionally deletes pre-existing output paths, so /dev/null
# is not a valid output target for this benchmark.
#
# Useful overrides:
#   INPUT=/path/to.fa scripts/tune_scoring_constants.sh
#   THREADS=16 scripts/tune_scoring_constants.sh
#   REPEATS=3 scripts/tune_scoring_constants.sh
#   RUN_LIMIT=5 scripts/tune_scoring_constants.sh
#   KEEP_OUTPUTS=1 HASH_OUTPUT=1 scripts/tune_scoring_constants.sh

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCORING_RS="$REPO_DIR/src/scoring.rs"
INPUT="${INPUT:-$HOME/Documents/ФКН/bioinf/hw4-2026/NC_043715.1[1..59306649].fa}"
WINDOW_SIZE="${WINDOW_SIZE:-12}"
MIN_SIZE="${MIN_SIZE:-8}"
MAX_SIZE="${MAX_SIZE:-12}"
THREADS="${THREADS:-}"
REPEATS="${REPEATS:-1}"
RUN_LIMIT="${RUN_LIMIT:-0}"
OUT_DIR="${OUT_DIR:-/tmp/opencode/zhuntrs-scoring-tune}"
KEEP_OUTPUTS="${KEEP_OUTPUTS:-0}"
HASH_OUTPUT="${HASH_OUTPUT:-0}"

mkdir -p "$OUT_DIR"

if [[ ! -f "$INPUT" ]]; then
  printf 'Input FASTA not found: %s\n' "$INPUT" >&2
  exit 1
fi

ORIGINAL_SCORING="$(mktemp "$OUT_DIR/scoring.rs.original.XXXXXX")"
cp "$SCORING_RS" "$ORIGINAL_SCORING"

cleanup() {
  cp "$ORIGINAL_SCORING" "$SCORING_RS"
  rm -f "$ORIGINAL_SCORING"
}
trap cleanup EXIT

TARGET_DIR="$(cargo metadata --format-version 1 --no-deps 2>/dev/null \
  | sed -n 's/.*"target_directory":"\([^"]*\)".*/\1/p')"
if [[ -z "$TARGET_DIR" ]]; then
  TARGET_DIR="$REPO_DIR/target"
fi
ZHUNT_BIN="$TARGET_DIR/release/zhunt"

CSV="$OUT_DIR/results.csv"
printf 'variant,score_block_positions,progress_batch_size,newton_steps,grid_adjust_limit,repeat,seconds,output_bytes,sha256,status\n' > "$CSV"

variants=(
  'baseline|8192|1000|6|8'

  'block_1024|1024|1000|6|8'
  'block_2048|2048|1000|6|8'
  'block_4096|4096|1000|6|8'
  'block_16384|16384|1000|6|8'
  'block_32768|32768|1000|6|8'
  'block_65536|65536|1000|6|8'
  'block_131072|131072|1000|6|8'
  'block_262144|262144|1000|6|8'

  'progress_512|8192|512|6|8'
  'progress_2000|8192|2000|6|8'
  'progress_5000|8192|5000|6|8'
  'progress_10000|8192|10000|6|8'
  'progress_50000|8192|50000|6|8'

  'newton_0|8192|10000|0|8'
  'newton_2|8192|10000|2|8'
  'newton_4|8192|10000|4|8'
  'newton_8|8192|10000|8|8'
  'newton_10|8192|10000|10|8'
  'newton_12|8192|10000|12|8'

  'adjust_0|8192|10000|6|0'
  'adjust_2|8192|10000|6|2'
  'adjust_4|8192|10000|6|4'
  'adjust_16|8192|10000|6|16'
  'adjust_32|8192|10000|6|32'

  'combo_b16k_n4|16384|10000|4|8'
  'combo_b32k_n4|32768|10000|4|8'
  'combo_b64k_n4|65536|10000|4|8'
  'combo_b16k_n6_a16|16384|10000|6|16'
  'combo_b32k_n6_a16|32768|10000|6|16'
)

apply_variant() {
  local block_positions="$1"
  local progress_batch="$2"
  local newton_steps="$3"
  local grid_adjust="$4"

  cp "$ORIGINAL_SCORING" "$SCORING_RS"
  sed -i -E \
    -e "s/^(const PROGRESS_BATCH_SIZE: usize = )[0-9_]+;/\1${progress_batch};/" \
    -e "s/^(pub\(crate\) const SCORE_WORK_BLOCK_POSITIONS: usize = )[0-9_]+;/\1${block_positions};/" \
    -e "s/^(const DELTA_LINKING_NEWTON_STEPS: usize = )[0-9_]+;/\1${newton_steps};/" \
    -e "s/^(const DELTA_LINKING_GRID_ADJUST_LIMIT: usize = )[0-9_]+;/\1${grid_adjust};/" \
    "$SCORING_RS"

  local actual_progress actual_block actual_newton actual_adjust
  actual_progress="$(sed -n -E 's/^const PROGRESS_BATCH_SIZE: usize = ([0-9_]+);/\1/p' "$SCORING_RS")"
  actual_block="$(sed -n -E 's/^pub\(crate\) const SCORE_WORK_BLOCK_POSITIONS: usize = ([0-9_]+);/\1/p' "$SCORING_RS")"
  actual_newton="$(sed -n -E 's/^const DELTA_LINKING_NEWTON_STEPS: usize = ([0-9_]+);/\1/p' "$SCORING_RS")"
  actual_adjust="$(sed -n -E 's/^const DELTA_LINKING_GRID_ADJUST_LIMIT: usize = ([0-9_]+);/\1/p' "$SCORING_RS")"

  if [[ "$actual_progress" != "$progress_batch" \
     || "$actual_block" != "$block_positions" \
     || "$actual_newton" != "$newton_steps" \
     || "$actual_adjust" != "$grid_adjust" ]]; then
    printf 'sed replacement failed for block=%s progress=%s newton=%s adjust=%s\n' \
      "$block_positions" "$progress_batch" "$newton_steps" "$grid_adjust" >&2
    exit 2
  fi
}

printf 'Benchmark input : %s\n' "$INPUT"
printf 'zhunt args      : %s %s %s\n' "$WINDOW_SIZE" "$MIN_SIZE" "$MAX_SIZE"
if [[ -n "$THREADS" ]]; then
  printf 'threads         : %s\n' "$THREADS"
else
  printf 'threads         : default zhunt policy\n'
fi
printf 'results CSV     : %s\n' "$CSV"
printf 'variant count   : %s\n\n' "${#variants[@]}"

variant_index=0
for variant in "${variants[@]}"; do
  variant_index=$((variant_index + 1))
  if (( RUN_LIMIT > 0 && variant_index > RUN_LIMIT )); then
    break
  fi

  IFS='|' read -r name block_positions progress_batch newton_steps grid_adjust <<< "$variant"
  printf '[%02d/%02d] %s: block=%s progress=%s newton=%s adjust=%s\n' \
    "$variant_index" "${#variants[@]}" "$name" "$block_positions" "$progress_batch" "$newton_steps" "$grid_adjust"

  apply_variant "$block_positions" "$progress_batch" "$newton_steps" "$grid_adjust"

  build_log="$OUT_DIR/${name}.build.log"
  cargo build --release --bin zhunt > "$build_log" 2>&1

  for ((repeat = 1; repeat <= REPEATS; repeat++)); do
    run_log="$OUT_DIR/${name}.${repeat}.run.log"
    status='ok'
    output_bytes='0'
    output_hash=''

    output_file="$OUT_DIR/${name}.${repeat}.Z-SCORE"
    rm -f "$output_file"

    thread_args=()
    if [[ -n "$THREADS" ]]; then
      thread_args=(--threads "$THREADS")
    fi

    start_ns="$(date +%s%N)"
    if ! "$ZHUNT_BIN" "${thread_args[@]}" -o "$output_file" \
      "$WINDOW_SIZE" "$MIN_SIZE" "$MAX_SIZE" "$INPUT" \
      > "$run_log" 2>&1; then
      status='fail'
    fi
    end_ns="$(date +%s%N)"
    elapsed_ns=$((end_ns - start_ns))
    seconds="$(printf '%d.%03d' \
      "$((elapsed_ns / 1000000000))" \
      "$(((elapsed_ns / 1000000) % 1000))")"

    if [[ -f "$output_file" ]]; then
      output_bytes="$(wc -c < "$output_file")"
      if [[ "$HASH_OUTPUT" == '1' ]]; then
        output_hash="$(sha256sum "$output_file" | sed -n 's/ .*//p')"
      fi
      if [[ "$KEEP_OUTPUTS" != '1' ]]; then
        rm -f "$output_file"
      fi
    fi

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
      "$name" "$block_positions" "$progress_batch" "$newton_steps" "$grid_adjust" \
      "$repeat" "$seconds" "$output_bytes" "$output_hash" "$status" >> "$CSV"
    printf '  repeat %s: %ss (%s)\n' "$repeat" "$seconds" "$status"
  done
done

printf '\nFastest successful runs:\n'
sed '1d' "$CSV" | sed -n '/,ok$/p' | sort -t, -k7,7n | sed -n '1,10p'
printf '\nRestored original src/scoring.rs. Full CSV: %s\n' "$CSV"
