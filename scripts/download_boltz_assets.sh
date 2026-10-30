#!/usr/bin/env bash
# Pre-fetch Boltz-2 weights/CCD from HuggingFace into the workspace cache (boltz's own downloader hung via model-gateway.boltz.bio).
set -euo pipefail
WS=${WS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
C=$WS/data/boltz/cache; mkdir -p $C
B2=https://huggingface.co/boltz-community/boltz-2/resolve/main
for f in boltz2_conf.ckpt mols.tar; do [ -s $C/$f ] || curl -sSfL --retry 3 -o $C/$f $B2/$f; done
[ -s $C/ccd.pkl ] || curl -sSfL --retry 3 -o $C/ccd.pkl https://huggingface.co/boltz-community/boltz-1/resolve/main/ccd.pkl
[ -f $C/.mols_done ] || { tar -xf $C/mols.tar -C $C && touch $C/.mols_done; }
ls -la $C
