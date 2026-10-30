#!/usr/bin/env python3
"""Optional extra AF3 seed (modelSeeds [1]) for the DIMER runs of the on-axis stratum only (budget-limited; task:
"add seed 1 only if time allows"). Writes inputs named <case>_dimer_seed1 and splits them across two GPU job dirs.
usage: cofold_seed1_onaxis.py GPU_A GPU_B   then run AF3 on data/cofold/onaxis/jobs_gpu{A,B} (scripts/cofold_run_onaxis.sh --no-split)
Status: first attempt, UNVERIFIED."""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
W = ROOT / "data/cofold/onaxis"
SEED = 1
STRATA = ["on_axis_druglike"]

ga, gb = sys.argv[1], sys.argv[2]
c = pd.read_csv(ROOT / "results/cofold/onaxis/cases.tsv", sep="\t")
c = c[c.stratum.isin(STRATA)].sort_values("tokens_dimer")
for g in (ga, gb):
    d = W / f"jobs_gpu{g}"
    d.mkdir(parents=True, exist_ok=True)
    for f in d.glob("*.json"):
        f.unlink()
load = {ga: 0.0, gb: 0.0}
for i, (_, r) in enumerate(c.iterrows()):
    name = f"{r.case_id}_dimer_seed{SEED}"
    if (W / "af3" / name / f"{name}_ranking_scores.csv").exists():
        continue
    js = json.loads((ROOT / r.input_dimer).read_text())
    js["name"], js["modelSeeds"] = name, [SEED]
    p = W / "inputs" / f"{name}.json"
    p.write_text(json.dumps(js))
    g = min(load, key=load.get)
    load[g] += (max(r.tokens_dimer, 256) / 512) ** 2
    (W / f"jobs_gpu{g}" / f"{i:03d}_{name}.json").symlink_to(p)
print("seed-1 jobs split, est. load", load)
