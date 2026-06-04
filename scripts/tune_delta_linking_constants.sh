#!/usr/bin/env bash
set -Eeuo pipefail

# Sed-based tuner for the delta-linking root/snap constants in src/scoring.rs.
#
# This tuner is intentionally stricter than tune_scoring_constants.sh: it hashes
# every generated Z-SCORE file by default and compares each variant with the
# first reference run. Treat sha_match=yes as the conservative "no precision
# change" filter.
#
# Useful overrides:
#   INPUT=/path/to.fa scripts/tune_delta_linking_constants.sh
#   THREADS=16 scripts/tune_delta_linking_constants.sh
#   REPEATS=3 scripts/tune_delta_linking_constants.sh
#   RUN_LIMIT=5 scripts/tune_delta_linking_constants.sh
#   KEEP_OUTPUTS=1 scripts/tune_delta_linking_constants.sh
#   HASH_OUTPUT=0 scripts/tune_delta_linking_constants.sh

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCORING_RS="$REPO_DIR/src/scoring.rs"
INPUT="${INPUT:-$HOME/Documents/ФКН/bioinf/hw4-2026/NC_043715.1[1..59306649].fa}"
WINDOW_SIZE="${WINDOW_SIZE:-12}"
MIN_SIZE="${MIN_SIZE:-8}"
MAX_SIZE="${MAX_SIZE:-12}"
THREADS="${THREADS:-}"
REPEATS="${REPEATS:-1}"
RUN_LIMIT="${RUN_LIMIT:-0}"
OUT_DIR="${OUT_DIR:-/tmp/opencode/zhuntrs-delta-linking-tune}"
KEEP_OUTPUTS="${KEEP_OUTPUTS:-0}"
HASH_OUTPUT="${HASH_OUTPUT:-1}"

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
printf 'variant,grid_steps,newton_steps,grid_adjust_limit,repeat,seconds,output_bytes,sha256,sha_match,status\n' > "$CSV"

# Reference is the legacy/full snap configuration. The current performance
# tuning changed PROGRESS_BATCH_SIZE, which should not affect output; this
# script only edits the delta-linking constants below.
variants=(
  'reference_g65536_n6_a8|65536|6|8'

  'current_g65536_n6_a0|65536|6|0'
  'grid64512_n6_a0|64512|6|0'
  'grid63488_n6_a0|63488|6|0'
  'grid61440_n6_a0|61440|6|0'
  'grid57344_n6_a0|57344|6|0'
  'grid49152_n6_a0|49152|6|0'
  'grid40960_n6_a0|40960|6|0'
  'grid32768_n6_a0|32768|6|0'
  'grid24576_n6_a0|24576|6|0'
  'grid16384_n6_a0|16384|6|0'
  'grid12288_n6_a0|12288|6|0'
  'grid8192_n6_a0|8192|6|0'
  'grid4096_n6_a0|4096|6|0'

  'adjust1_g65536_n6|65536|6|1'
  'adjust2_g65536_n6|65536|6|2'
  'adjust4_g65536_n6|65536|6|4'

  'newton0_g65536_a8|65536|0|8'
  'newton1_g65536_a8|65536|1|8'
  'newton2_g65536_a8|65536|2|8'
  'newton3_g65536_a8|65536|3|8'
  'newton4_g65536_a8|65536|4|8'
  'newton5_g65536_a8|65536|5|8'
  'newton7_g65536_a8|65536|7|8'
  'newton8_g65536_a8|65536|8|8'
  'newton10_g65536_a8|65536|10|8'
  'newton12_g65536_a8|65536|12|8'

  'grid8192_n6_a8|8192|6|8'
  'grid16384_n6_a8|16384|6|8'
  'grid32768_n6_a8|32768|6|8'
  'grid49152_n6_a8|49152|6|8'
  'grid98304_n6_a8|98304|6|8'
  'grid131072_n6_a8|131072|6|8'

  'combo_g32768_n5_a8|32768|5|8'
  'combo_g32768_n7_a8|32768|7|8'
  'combo_g98304_n5_a8|98304|5|8'
  'combo_g98304_n7_a8|98304|7|8'
  'combo_g131072_n5_a8|131072|5|8'
  'combo_g131072_n7_a8|131072|7|8'
)

apply_variant() {
  local grid_steps="$1"
  local newton_steps="$2"
  local grid_adjust="$3"

  cp "$ORIGINAL_SCORING" "$SCORING_RS"
  sed -i -E \
    -e "s/^(const DELTA_LINKING_GRID_STEPS: usize = )[0-9_]+;/\1${grid_steps};/" \
    -e "s/^(const DELTA_LINKING_NEWTON_STEPS: usize = )[0-9_]+;/\1${newton_steps};/" \
    -e "s/^(const DELTA_LINKING_GRID_ADJUST_LIMIT: usize = )[0-9_]+;/\1${grid_adjust};/" \
    "$SCORING_RS"

  local actual_grid actual_newton actual_adjust
  actual_grid="$(sed -n -E 's/^const DELTA_LINKING_GRID_STEPS: usize = ([0-9_]+);/\1/p' "$SCORING_RS")"
  actual_newton="$(sed -n -E 's/^const DELTA_LINKING_NEWTON_STEPS: usize = ([0-9_]+);/\1/p' "$SCORING_RS")"
  actual_adjust="$(sed -n -E 's/^const DELTA_LINKING_GRID_ADJUST_LIMIT: usize = ([0-9_]+);/\1/p' "$SCORING_RS")"

  if [[ "$actual_grid" != "$grid_steps" \
     || "$actual_newton" != "$newton_steps" \
     || "$actual_adjust" != "$grid_adjust" ]]; then
    printf 'sed replacement failed for grid=%s newton=%s adjust=%s\n' \
      "$grid_steps" "$newton_steps" "$grid_adjust" >&2
    exit 2
  fi
}

original_progress="$(sed -n -E 's/^const PROGRESS_BATCH_SIZE: usize = ([0-9_]+);/\1/p' "$ORIGINAL_SCORING")"
original_block="$(sed -n -E 's/^pub\(crate\) const SCORE_WORK_BLOCK_POSITIONS: usize = ([0-9_]+);/\1/p' "$ORIGINAL_SCORING")"

printf 'Benchmark input : %s\n' "$INPUT"
printf 'zhunt args      : %s %s %s\n' "$WINDOW_SIZE" "$MIN_SIZE" "$MAX_SIZE"
printf 'fixed progress  : %s\n' "$original_progress"
printf 'fixed block     : %s\n' "$original_block"
if [[ -n "$THREADS" ]]; then
  printf 'threads         : %s\n' "$THREADS"
else
  printf 'threads         : default zhunt policy\n'
fi
printf 'hash output     : %s\n' "$HASH_OUTPUT"
printf 'results CSV     : %s\n' "$CSV"
printf 'variant count   : %s\n\n' "${#variants[@]}"

reference_hash=''
variant_index=0
for variant in "${variants[@]}"; do
  variant_index=$((variant_index + 1))
  if (( RUN_LIMIT > 0 && variant_index > RUN_LIMIT )); then
    break
  fi

  IFS='|' read -r name grid_steps newton_steps grid_adjust <<< "$variant"
  printf '[%02d/%02d] %s: grid=%s newton=%s adjust=%s\n' \
    "$variant_index" "${#variants[@]}" "$name" "$grid_steps" "$newton_steps" "$grid_adjust"

  apply_variant "$grid_steps" "$newton_steps" "$grid_adjust"

  build_log="$OUT_DIR/${name}.build.log"
  cargo build --release --bin zhunt > "$build_log" 2>&1

  for ((repeat = 1; repeat <= REPEATS; repeat++)); do
    run_log="$OUT_DIR/${name}.${repeat}.run.log"
    status='ok'
    output_bytes='0'
    output_hash=''
    sha_match=''

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
        if [[ -z "$reference_hash" && "$status" == 'ok' ]]; then
          reference_hash="$output_hash"
        fi
        if [[ -n "$reference_hash" ]]; then
          if [[ "$output_hash" == "$reference_hash" ]]; then
            sha_match='yes'
          else
            sha_match='no'
          fi
        fi
      fi
      if [[ "$KEEP_OUTPUTS" != '1' ]]; then
        rm -f "$output_file"
      fi
    fi

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
      "$name" "$grid_steps" "$newton_steps" "$grid_adjust" "$repeat" \
      "$seconds" "$output_bytes" "$output_hash" "$sha_match" "$status" >> "$CSV"
    printf '  repeat %s: %ss (%s, sha_match=%s)\n' "$repeat" "$seconds" "$status" "$sha_match"
  done
done

printf '\nFastest successful byte-identical runs:\n'
sed '1d' "$CSV" | sed -n '/,yes,ok$/p' | sort -t, -k6,6n | sed -n '1,10p'

printf '\nFastest successful changed-output runs:\n'
sed '1d' "$CSV" | sed -n '/,no,ok$/p' | sort -t, -k6,6n | sed -n '1,10p'

printf '\nRestored original src/scoring.rs. Full CSV: %s\n' "$CSV"
