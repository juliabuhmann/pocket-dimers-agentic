#!/usr/bin/env python
"""Fetch / link the small set of structures used to smoke-test burrow, pocketcompare and adock.

Writes to data/interim/tool_tests/:
  afdb/      AFDB monomer models downloaded from the AlphaFold DB file server
  dimers/    symlinks to predicted human dimers (AFDB quaternary models, local scratch copy)
  pdb/       RCSB mmCIF files (as deposited) + a ligand/water-stripped copy (gemmi)
  dimer_batch.txt, dimer_batch_large.txt  first N_BATCH / N_BATCH_LARGE dimers of each class (sorted) for timing
"""
from __future__ import annotations
import os

import argparse
import urllib.request
from pathlib import Path

import gemmi

WS = Path(__file__).resolve().parents[1]
OUT = WS / "data/interim/tool_tests"

# 14-3-3 beta (P31946): chosen because the same protein exists as AFDB monomer and as a predicted homodimer.
AFDB_ACCESSIONS = ["P31946"]
AFDB_VERSION = 6
AFDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{v}.cif"

DIMER_ROOT = Path(os.environ.get("DIMERS_ROOT", "/cluster/scratch/gmagnusson/dimers")) / "afdb_human_dimers"
# homodimer of P31946 (DBREF A,B = P31946) and a heterodimer Q04917/P27348 (14-3-3 eta/theta)
DIMER_FILES = ["homodimers/AF-0000000065760335.pdb", "heterodimers/AF-0000000203452075.pdb"]

# Progesterone receptor LBD crystal structure: 2 chains (A, B) with progesterone (STR) bound to each.
PDB_IDS = ["1A28"]
RCSB_URL = "https://files.rcsb.org/download/{pid}.cif"

N_BATCH = 50  # dimers per class for the quick batch test
N_BATCH_LARGE = 500  # dimers per class for the GPU-vs-CPU throughput test (GPU warm-up ~12 s needs a big batch)


def fetch(url: str, dest: Path) -> None:
    if dest.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=60) as r:
        dest.write_bytes(r.read())
    print("fetched", url)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-batch", type=int, default=N_BATCH)
    a = ap.parse_args()

    for acc in AFDB_ACCESSIONS:
        fetch(AFDB_URL.format(acc=acc, v=AFDB_VERSION), OUT / f"afdb/AF-{acc}-F1-model_v{AFDB_VERSION}.cif")

    (OUT / "dimers").mkdir(parents=True, exist_ok=True)
    for rel in DIMER_FILES:
        link = OUT / "dimers" / Path(rel).name
        if not link.exists():
            link.symlink_to(DIMER_ROOT / rel)

    for pid in PDB_IDS:
        raw = OUT / f"pdb/{pid}.cif"
        fetch(RCSB_URL.format(pid=pid), raw)
        st = gemmi.read_structure(str(raw))
        st.setup_entities()
        st.remove_ligands_and_waters()   # drops HETATM non-polymer residues and waters
        st.remove_empty_chains()
        doc = st.make_mmcif_document()
        doc.write_file(str(OUT / f"pdb/{pid}_protein.cif"))

    for n, name in ((a.n_batch, "dimer_batch.txt"), (N_BATCH_LARGE, "dimer_batch_large.txt")):
        batch = []
        for cls in ("homodimers", "heterodimers"):
            files = sorted((DIMER_ROOT / cls).glob("*.pdb"))[:n]
            batch += [str(f) for f in files]
        (OUT / name).write_text("\n".join(batch) + "\n")
        print(f"{len(batch)} dimers listed in {OUT/name}")


if __name__ == "__main__":
    main()
