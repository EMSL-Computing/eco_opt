#!/usr/bin/env bash
# Copy example PNGs into docs/assets/examples/<name>/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
for d in "$ROOT"/examples/*/figs; do
  name="$(basename "$(dirname "$d")")"
  dest="$ROOT/docs/assets/examples/$name"
  mkdir -p "$dest"
  cp -a "$d"/*.png "$dest/"
done
