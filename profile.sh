#!/bin/sh
set -eu

CARGO_PROFILE_RELEASE_DEBUG=true \
RUSTFLAGS="-C force-frame-pointers=yes" \
cargo flamegraph --release --bin zhunt \
  --no-inline \
  -c "record -e cycles:u -g -F 199" \
  -o flamegraph.svg \
  -- "$@"
