#!/usr/bin/env bash
# Count + sha256 of sorted directory listings (AF3 MSA folders, AFDB v6 tars).
OUT=$(cd "$(dirname "$0")/.." && pwd)/results/checksums
mkdir -p "$OUT"
h() { # name dir
  ( cd "$2" && LC_ALL=C ls | LC_ALL=C sort > "$OUT/$1.listing.tmp" )
  wc -l < "$OUT/$1.listing.tmp" > "$OUT/$1.count"
  sha256sum < "$OUT/$1.listing.tmp" | cut -d' ' -f1 > "$OUT/$1.listing.sha256"
  rm -f "$OUT/$1.listing.tmp"
}
h af3_msas_human "${AF3_MSA_HUMAN:-/cluster/project/beltrao/shared/25.05_alphafold3_msas_human}"
h af3_msas_yeast "${AF3_MSA_YEAST:-/cluster/project/beltrao/shared/25.06_alphafold3_msas_yeast}"
h afdb_v6_tars "${AFDB_V6_ROOT:-/cluster/work/beltrao/shared/alphafold/v6}"
echo "done $(date -Iseconds)" > "$OUT/LISTINGS_DONE"
