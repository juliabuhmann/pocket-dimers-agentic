#!/usr/bin/env bash
# All-vs-all human paralog search with mmseqs2 (conda_envs/mmseqs, installed from bioconda).
# Paralog = mmseqs2 hit with sequence identity >= 0.30 (DESIGN.md §3) and e-value <= 1e-3 (our choice); no coverage
# requirement (conservative: a shared domain is enough to exclude). Self-hits are dropped downstream.
# Input: data/interim/q22/human_seqs.fasta (scripts/q22_fetch_uniprot.py). Output: data/interim/q22/paralogs_mmseqs.tsv
# Status: first attempt, not yet verified.
set -euo pipefail
WS="$(cd "$(dirname "$0")/.." && pwd)"
MM="$WS/conda_envs/mmseqs/bin/mmseqs"
D="$WS/data/interim/q22"
TMP="$D/mmseqs_tmp"
mkdir -p "$TMP"
# restrict to accessions used here (dimer set + human atlas); the proteome stream also contains unreviewed entries
"$WS/.venvs/q22/bin/python" - "$D" <<'PY'
import sys, pandas as pd
D = sys.argv[1]
d = pd.read_parquet(f"{D}/../../part2/human_dimer_pockets_annotated.parquet", columns=["uniprot_A", "uniprot_B"])
a = pd.read_parquet(f"{D}/atlas_human_meta.parquet", columns=["accession"])
need = set(d.uniprot_A) | set(d.uniprot_B) | set(a.accession)
out, keep = [], False
for ln in open(f"{D}/human_seqs.fasta"):
    if ln.startswith(">"):
        keep = ln[1:].strip() in need
    if keep:
        out.append(ln)
open(f"{D}/human_seqs_needed.fasta", "w").writelines(out)
print("sequences kept", sum(l.startswith(">") for l in out), "of needed", len(need))
PY
"$MM" version
"$MM" easy-search "$D/human_seqs_needed.fasta" "$D/human_seqs_needed.fasta" "$D/paralogs_mmseqs.tsv" "$TMP" \
  --min-seq-id 0.30 -e 1e-3 -s 7.5 --max-seqs 2000 --threads 7 \
  --format-output "query,target,fident,alnlen,evalue,qcov,tcov"
wc -l "$D/paralogs_mmseqs.tsv"
rm -rf "$TMP"
