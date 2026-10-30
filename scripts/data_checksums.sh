#!/usr/bin/env bash
# Computes sha256 of small index/manifest files of the read-only data mounts.
# Output: data/raw/CHECKSUMS.sha256 (rerun to regenerate). msa_depths.csv (~GB) is skipped; see DATA.md.
set -euo pipefail
OUT="$(dirname "$0")/../data/raw/CHECKSUMS.sha256"
D=${DIMERS_ROOT:-/cluster/scratch/gmagnusson/dimers}
P=${ATLAS_ROOT:-/cluster/scratch/gmagnusson/pocketcompare_afdb}
{
  sha256sum $D/*/*/candidates.csv $D/*/*/download_manifest.csv
  sha256sum $P/POCKETOME.md $P/worklist.json $P/gcs_listing.json $P/gcs_listing.tsv.gz \
            $P/skip_accessions.txt $P/worklist.tsv $P/worklist.parquet $P/ncbi/*.dmp
} > "$OUT"
wc -l "$OUT"
