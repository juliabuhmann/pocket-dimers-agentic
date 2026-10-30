#!/usr/bin/env python
"""Q2.4: divergence time of each dimer species from human, fetched from the TimeTree 5 pairwise API.

Query: https://timetree.org/api/pairwise/9606/<taxid>  (CSV; column `adjusted_age` in million years, Kumar et al. 2022 MBE).
Taxids are NCBI Taxonomy species-level ids for the AFDB reference proteomes of each dimer set (strain ids such as
S. cerevisiae S288C 559292 are replaced by the species id, which is what TimeTree indexes).
Output: results/part2/q24/q24_divergence.tsv (with retrieval date). Status: first attempt, not yet verified.
"""
import io, datetime, urllib.request
from pathlib import Path
import pandas as pd

WS = Path(__file__).resolve().parents[1]
OUT = WS / "results/part2/q24/q24_divergence.tsv"
HUMAN = 9606
TAXIDS = {"mouse": 10090, "rat": 10116, "zebrafish": 7955, "fly": 7227, "worm": 6239, "yeast": 4932, "pombe": 4896,
          "arabidopsis": 3702, "dicty": 44689, "ecoli": 562, "mtb": 1773, "paeruginosa": 287}
URL = "https://timetree.org/api/pairwise/{a}/{b}"
# DOCUMENTED ASSUMPTION: the TimeTree pairwise API returned HTTP 500 for every bacterial taxid tried (562, 1773, 287, 2, 511145;
# 2026-10-29). Human-bacteria divergence is the root of the tree of life; we use an approximate LUCA age of 4,200 Mya
# (order of magnitude as in TimeTree 5 / Kumar et al. 2022). Only the ordering (bacteria most distant) is used downstream.
FALLBACK_MYA = {"ecoli": 4200.0, "mtb": 4200.0, "paeruginosa": 4200.0}

rows = []
for sp, tx in TAXIDS.items():
    try:
        txt = urllib.request.urlopen(URL.format(a=HUMAN, b=tx), timeout=60).read().decode()
        d = pd.read_csv(io.StringIO(txt)).iloc[0]
        rows.append((sp, tx, d.scientific_name_b, float(d.adjusted_age), float(d.precomputed_ci_low), float(d.precomputed_ci_high), "timetree_api"))
    except Exception as e:  # recorded, not silently filled
        fb = FALLBACK_MYA.get(sp, float("nan"))
        rows.append((sp, tx, "", fb, float("nan"), float("nan"), f"ASSUMED (approx. LUCA age); API fetch_failed: {e}"))
df = pd.DataFrame(rows, columns=["species", "taxid", "scientific_name", "divergence_mya", "ci_low", "ci_high", "source"])
OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as f:
    f.write(f"# TimeTree pairwise API ({URL.format(a=HUMAN, b='<taxid>')}), retrieved {datetime.date.today().isoformat()}\n")
    df.to_csv(f, sep="\t", index=False)
print(df.to_string())
