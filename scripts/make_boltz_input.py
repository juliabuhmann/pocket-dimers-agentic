#!/usr/bin/env python3
"""Build a Boltz-2 YAML (+ a3m files) for a dimer(+CCD ligand) from the cached AF3 MSAs.
Uses the cached unpairedMsa as a3m per entity; identical accessions share one entity (id [A,B]).
No cross-chain MSA pairing is provided (Boltz would otherwise pair via ColabFold server).
First attempt, not verified. usage: make_boltz_input.py --acc P60174 P60174 --ligand-ccd PGA --outdir DIR"""
import os
import argparse, gzip, json, string, sys
from pathlib import Path
MSA_DIR = os.environ.get("AF3_MSA_HUMAN", "/cluster/project/beltrao/shared/25.05_alphafold3_msas_human")
ap = argparse.ArgumentParser()
ap.add_argument("--acc", nargs=2, required=True)
ap.add_argument("--ligand-ccd", nargs="*", default=[])
ap.add_argument("--outdir", required=True)
ap.add_argument("--msa-dir", default=MSA_DIR)
ap.add_argument("--name", default=None)
a = ap.parse_args()
out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
ents = []
for acc in a.acc:
    if ents and ents[-1][0] == acc: ents[-1][2] += 1; continue
    d = json.load(gzip.open(Path(a.msa_dir) / f"{acc.lower()}_data.json.gz", "rt"))
    p = [s["protein"] for s in d["sequences"] if "protein" in s][0]
    (out / f"{acc}.a3m").write_text(p["unpairedMsa"])
    ents.append([acc, p["sequence"], 1])
letters = iter(string.ascii_uppercase)
lines = ["version: 1", "sequences:"]
for acc, seq, n in ents:
    ids = [next(letters) for _ in range(n)]
    lines += ["  - protein:", f"      id: [{', '.join(ids)}]", f"      sequence: {seq}", f"      msa: {(out / (acc + '.a3m')).resolve()}"]
for c in a.ligand_ccd:
    lines += ["  - ligand:", f"      id: {next(letters)}", f"      ccd: {c}"]
name = a.name or "_".join(x[0].lower() for x in ents)
(out / f"{name}.yaml").write_text("\n".join(lines) + "\n")
print(out / f"{name}.yaml")
