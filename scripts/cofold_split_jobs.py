#!/usr/bin/env python3
"""Split AF3 inputs of results/cofold/onaxis/cases.tsv across two GPU job directories (symlinks).
Order = priority (dimer runs of on-axis cases first, then controls' dimers, then all monomers), so a budget cut
loses monomer/control runs first. Greedy balance on estimated runtime ~ tokens^2 (XLA attention, docs/COFOLDING.md).
Jobs whose output already exists (ranking_scores.csv) are skipped.
usage: cofold_split_jobs.py GPU_A GPU_B"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
W = ROOT / "data/cofold/onaxis"
PRIO = {"on_axis_druglike": 0, "offaxis_cofactor": 1, "sbs_druglike": 2}

ga, gb = sys.argv[1], sys.argv[2]
c = pd.read_csv(ROOT / "results/cofold/onaxis/cases.tsv", sep="\t")
jobs = []
for _, r in c.iterrows():
    for cond, tok in (("dimer", r.tokens_dimer), ("monomer", r.tokens_monomer)):
        name = f"{r.case_id}_{cond}"
        if (W / "af3" / name / f"{name}_ranking_scores.csv").exists():
            continue
        jobs.append((0 if cond == "dimer" else 1, PRIO[r.stratum], name, ROOT / r[f"input_{cond}"], tok))
jobs.sort(key=lambda j: (j[0], j[1], j[4]))
load = {ga: 0.0, gb: 0.0}
for g in (ga, gb):
    d = W / f"jobs_gpu{g}"
    d.mkdir(parents=True, exist_ok=True)
    for f in d.glob("*.json"):
        f.unlink()
for i, (_, _, name, path, tok) in enumerate(jobs):
    g = min(load, key=load.get)
    load[g] += (max(tok, 256) / 512) ** 2
    (W / f"jobs_gpu{g}" / f"{i:03d}_{name}.json").symlink_to(path)
print("jobs", len(jobs), "est. load (units of a 512-token job)", load)
