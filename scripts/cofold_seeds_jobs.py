#!/usr/bin/env python3
"""Replicate-seed AF3 jobs for the on-axis (docs/COFOLD_ONAXIS.md) and dimer-only (docs/COFOLD_DIMERONLY.md)
co-folding experiments. Inputs are the EXISTING seed-0 input JSONs with only `modelSeeds` (and `name`) changed.
Seeds follow the project's fixed list [0,1,2]; seed 0 is the finished original run, so --seeds defaults to 1 2.

Scope: on-axis experiment = on_axis_druglike + offaxis_cofactor cases (SBS controls skipped), dimer + monomer;
dimer-only experiment = all cases, dimer + monomer. Jobs whose outputs exist are skipped (so the seed-1 on-axis
dimers from scripts/cofold_seed1_onaxis.py, named <case>_dimer_seed1 in data/cofold/onaxis/af3, are not redone).
Order: seed-major (every case/condition of seed 1 before any seed 2); inside a seed, cases by dimer token count
ascending (both experiments pooled, dimer and monomer of a case adjacent) so that a time cut loses the largest
cases of the last seed. Static greedy split over the GPU job dirs (est. runtime ~ tokens^2); each AF3 process
runs its dir in sorted-name order, so global order is preserved on each GPU.
usage: cofold_seeds_jobs.py GPU [GPU ...] [--seeds 1 2]      Status: first attempt, UNVERIFIED."""
import argparse
import json
from pathlib import Path

import csv
from types import SimpleNamespace  # stdlib only: pandas import from the Lustre venv stalled

ROOT = Path(__file__).resolve().parents[1]
S = ROOT / "data/cofold/seeds"
EXPERIMENTS = {  # name -> (cases.tsv, strata to include (None = all), dir with existing AF3 outputs)
    "onaxis": (ROOT / "results/cofold/onaxis/cases.tsv", {"on_axis_druglike", "offaxis_cofactor"}, ROOT / "data/cofold/onaxis/af3"),
    "dimeronly": (ROOT / "results/cofold/dimeronly/cases.tsv", None, ROOT / "data/cofold/dimeronly/af3"),
}
CONDS = ("dimer", "monomer")
PROJECT_SEEDS = [0, 1, 2]

ap = argparse.ArgumentParser()
ap.add_argument("gpus", nargs="+")
ap.add_argument("--seeds", type=int, nargs="+", default=[s for s in PROJECT_SEEDS if s != 0])
a = ap.parse_args()

SKIP = S / "skip.txt"   # jobs that failed as first job of a fresh process (written by cofold_run_seeds.sh)
skip = set(SKIP.read_text().split()) if SKIP.exists() else set()
jobs = []
for exp, (cases_tsv, strata, af3_old) in EXPERIMENTS.items():
    with open(cases_tsv) as fh:
        c = [SimpleNamespace(**row) for row in csv.DictReader(fh, delimiter="\t")]
    if strata is not None:
        c = [r for r in c if r.stratum in strata]
    for r in c:
        for seed in a.seeds:
            for cond in CONDS:
                name = f"{r.case_id}_{cond}_seed{seed}"
                if name in skip or any((d / name / f"{name}_ranking_scores.csv").exists() for d in (S / "af3", af3_old)):
                    continue
                tok = int(getattr(r, f"tokens_{cond}"))
                jobs.append((seed, int(r.tokens_dimer), r.case_id, 0 if cond == "dimer" else 1, name, tok, getattr(r, f"input_{cond}"), exp))
jobs.sort(key=lambda j: j[:4])
load = {g: 0.0 for g in a.gpus}
for g in a.gpus:
    d = S / f"jobs_gpu{g}"
    d.mkdir(parents=True, exist_ok=True)
    for f in d.glob("*.json"):
        f.unlink()
for i, (seed, _, _, _, name, tok, inp, exp) in enumerate(jobs):
    js = json.loads((ROOT / inp).read_text())
    js["name"], js["modelSeeds"] = name, [seed]
    p = S / "inputs" / f"{name}.json"
    p.write_text(json.dumps(js))
    g = min(load, key=load.get)
    load[g] += (max(tok, 256) / 512) ** 2
    (S / f"jobs_gpu{g}" / f"{i:03d}_{name}.json").symlink_to(p)
print("jobs", len(jobs), "est. load (units of a 512-token job)", load)
