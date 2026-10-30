#!/usr/bin/env python3
"""Build an AF3 dimer(+ligand) input JSON from cached per-protein AF3 data JSONs
(unpairedMsa, pairedMsa, templates), so run_alphafold.py can use --norun_data_pipeline.

Status: first attempt; see docs/COFOLDING.md for what has been verified.

Examples
  homodimer + CCD ligand:
    make_af3_input.py --acc P60174 P60174 --ligand-ccd PGA --out tpi1_pga.json
  heterodimer + SMILES ligand, 2 seeds:
    make_af3_input.py --acc P04637 P38936 --ligand-smiles 'CCO' --seeds 1 2 --out x.json

Identical accessions -> one AF3 entity with ids [A, B] (copy number 2).
The cached sequence is checked against the current UniProt canonical sequence
(fetched from rest.uniprot.org); mismatch aborts unless --skip-uniprot-check.
Note: cached MSAs may cover only a sub-sequence ("subseq from"), the cached query is used as is.
"""
import os
import argparse, gzip, json, string, sys, urllib.request
from pathlib import Path

MSA_DIR_DEFAULT = os.environ.get("AF3_MSA_HUMAN", "/cluster/project/beltrao/shared/25.05_alphafold3_msas_human")


def load_cached(acc, msa_dir):
    p = Path(msa_dir) / f"{acc.lower()}_data.json.gz"
    if not p.exists():
        sys.exit(f"No cached MSA for {acc}: {p}")
    with gzip.open(p, "rt") as f:
        d = json.load(f)
    prots = [s["protein"] for s in d["sequences"] if "protein" in s]
    if len(prots) != 1:
        sys.exit(f"{p}: expected exactly 1 protein entry, found {len(prots)}")
    return prots[0]


def uniprot_seq(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc.upper()}.fasta"
    with urllib.request.urlopen(url, timeout=60) as r:
        txt = r.read().decode()
    return "".join(txt.splitlines()[1:])


def chain_ids():
    return iter(string.ascii_uppercase)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--acc", nargs=2, required=True, metavar=("ACC_A", "ACC_B"))
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--ligand-ccd", nargs="*", default=[], help="CCD codes, e.g. ATP HEM")
    g.add_argument("--ligand-smiles", nargs="*", default=[])
    ap.add_argument("--ligand-copies", type=int, default=1)
    ap.add_argument("--msa-dir", default=MSA_DIR_DEFAULT)
    ap.add_argument("--seeds", nargs="+", type=int, default=[1])
    ap.add_argument("--no-templates", action="store_true")
    ap.add_argument("--name", default=None)
    ap.add_argument("--skip-uniprot-check", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    entities = []  # (acc, protein dict, copies)
    for acc in a.acc:
        if entities and entities[-1][0].upper() == acc.upper():
            entities[-1][2] += 1
        else:
            entities.append([acc, load_cached(acc, a.msa_dir), 1])

    ids = chain_ids()
    seqs = []
    for acc, p, n in entities:
        if not a.skip_uniprot_check:
            ref = uniprot_seq(acc)
            if p["sequence"] != ref:
                sys.exit(f"{acc}: cached sequence (len {len(p['sequence'])}) != UniProt sequence (len {len(ref)})")
        if not p.get("pairedMsa"):
            print(f"WARNING: {acc} has empty pairedMsa", file=sys.stderr)
        prot = {
            "id": [next(ids) for _ in range(n)],
            "sequence": p["sequence"],
            "modifications": p.get("modifications", []),
            "unpairedMsa": p["unpairedMsa"],
            "pairedMsa": p.get("pairedMsa", ""),
            "templates": [] if a.no_templates else p.get("templates", []),
        }
        seqs.append({"protein": prot})
    for code in a.ligand_ccd:
        seqs.append({"ligand": {"id": [next(ids) for _ in range(a.ligand_copies)], "ccdCodes": [code]}})
    for smi in a.ligand_smiles:
        seqs.append({"ligand": {"id": [next(ids) for _ in range(a.ligand_copies)], "smiles": smi}})

    name = a.name or "_".join(x[0].lower() for x in entities) + "_x" + str(len(a.acc)) + "".join("_" + l for l in a.ligand_ccd)
    out = {"name": name, "modelSeeds": a.seeds, "sequences": seqs, "dialect": "alphafold3", "version": 2}
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(out, f)
    ntok = sum(len(e[1]["sequence"]) * e[2] for e in entities)
    print(f"wrote {a.out}: protein tokens={ntok} (+ligand atoms as tokens)")


if __name__ == "__main__":
    main()
