#!/bin/bash
# Sequential driver (16 workers). Resumable: finished chunks are skipped. Priority: human (A, then B isolated), then others.
cd "$(dirname "$0")/.." && source scripts/env_tools.sh pockets
W=${WORKERS:-16}
for s in human; do for k in homo hetero isolated; do python scripts/run_dimer_pockets.py --species $s --set $k --workers $W; done; done
echo "HUMAN_DONE $(date -Is)"
for s in mouse yeast ecoli pombe dicty fly worm zebrafish rat arabidopsis mtb paeruginosa; do
  for k in homo hetero; do python scripts/run_dimer_pockets.py --species $s --set $k --workers $W; done
  echo "SPECIES_DONE $s $(date -Is)"
done
echo "ALL_DONE $(date -Is)"
