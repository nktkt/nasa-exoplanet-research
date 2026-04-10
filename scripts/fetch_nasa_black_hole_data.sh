#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
RAW_DIR="$ROOT_DIR/data/raw"
OUT_FILE="$RAW_DIR/intbsc_black_hole_candidates.psv"

mkdir -p "$RAW_DIR"

QUERY_URL="https://heasarc.gsfc.nasa.gov/xamin/query?table=intbsc&fields=name,source_type,ra,dec,lii,bii,class&constraint=source_type+like+%27%25BHC%25%27+or+source_type+like+%27%25Microquasar%25%27+or+source_type+like+%27%25Sgr+A*%25%27&sortvar=name&format=stream&messages=none"

echo "Fetching HEASARC INTEGRAL black-hole candidate snapshot..."
curl -sSL "$QUERY_URL" | awk -F'|' 'NF > 1 { print }' > "$OUT_FILE"

echo "Wrote $OUT_FILE"
