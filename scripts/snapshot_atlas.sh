#!/usr/bin/env bash
# Copy the AFDB pocket atlas (monomers, tiers 1-2) off purgeable scratch into the workspace and checksum it.
set -euo pipefail
SRC=${ATLAS_SRC:-/cluster/scratch/gmagnusson/pocketcompare_afdb}
DST=${ATLAS_DST:-$(dirname "$0")/../data/raw/pocketcompare_afdb_snapshot}
mkdir -p "$DST"
cp -rn "$SRC/pockets" "$SRC/structures" "$SRC/POCKETOME.md" "$DST/"
( cd "$DST" && find . -type f -name '*.parquet' -o -name '*.md' | sort | xargs sha256sum > MANIFEST.sha256 )
echo "snapshot done $(date -Iseconds)" >> "$DST/SNAPSHOT_INFO.txt"
