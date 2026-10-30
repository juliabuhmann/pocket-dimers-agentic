#!/usr/bin/env bash
# Q2.4: reciprocal best hits (RBH) between human dimer proteins and each species' dimer proteins with mmseqs2 (conda_envs/mmseqs).
# Sequences come from the dimer model files (scripts/q24_prepare.py -> data/interim/q24/seqs_<species>.fasta).
# Both directions: easy-search, e-value <= 1e-3 (same cut-off as q22_paralogs.sh), -s 7.5, no identity floor (distant orthologs),
# alignment strings kept (qaln/taln) so pocket residues can be mapped downstream (scripts/q24_conservation.py picks the RBH by bitscore).
# Output: data/interim/q24/hits_<species>_{fwd,rev}.tsv. Status: first attempt, not yet verified.
set -euo pipefail
WS="$(cd "$(dirname "$0")/.." && pwd)"
source "$WS/scripts/env_threads.sh"
MM="$WS/conda_envs/mmseqs/bin/mmseqs"
D="$WS/data/interim/q24"
THREADS=${THREADS:-12}
FMT="query,target,fident,alnlen,evalue,bits,qstart,qend,tstart,tend,qlen,tlen,qaln,taln"
"$MM" version
for SP in mouse rat zebrafish fly worm yeast pombe arabidopsis dicty ecoli mtb paeruginosa; do
  for DIR in fwd rev; do
    OUT="$D/hits_${SP}_${DIR}.tsv"
    [ -s "$OUT" ] && continue
    [ -s "$D/seqs_${SP}.fasta" ] || { echo "$SP: no sequences yet, skipped"; continue; }
    if [ "$DIR" = fwd ]; then Q="$D/seqs_human.fasta"; T="$D/seqs_${SP}.fasta"; else Q="$D/seqs_${SP}.fasta"; T="$D/seqs_human.fasta"; fi
    TMP="$D/mmseqs_tmp_${SP}_${DIR}"
    "$MM" easy-search "$Q" "$T" "$OUT.part" "$TMP" -e 1e-3 -s 7.5 --max-seqs 300 --threads "$THREADS" \
      --format-output "$FMT" > "$D/mmseqs_${SP}_${DIR}.log" 2>&1
    mv "$OUT.part" "$OUT"; rm -rf "$TMP"
    echo "$SP $DIR $(wc -l < "$OUT")"
  done
done
