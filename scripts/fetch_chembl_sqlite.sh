#!/usr/bin/env bash
# Download the ChEMBL SQLite release (used for max_phase and pChEMBL>=6 flags per CCD-mapped compound).
set -euo pipefail
REL=${CHEMBL_RELEASE:-chembl_37}
URL=https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/${REL}/${REL}_sqlite.tar.gz
DST=$(dirname "$0")/../data/raw/chembl_sqlite
mkdir -p "$DST" && cd "$DST"
[ -f ${REL}_sqlite.tar.gz ] || curl -sS -o ${REL}_sqlite.tar.gz "$URL"
sha256sum ${REL}_sqlite.tar.gz > ${REL}_sqlite.tar.gz.sha256
tar -xzf ${REL}_sqlite.tar.gz
echo "$(date -Iseconds) $URL $(cut -d' ' -f1 ${REL}_sqlite.tar.gz.sha256)" >> FETCH_LOG.txt
