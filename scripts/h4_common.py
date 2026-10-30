"""Shared constants and helpers for the H4 / Q2.1 recovery analysis (scripts/h4_*.py).

Status: first attempt, NOT verified. Every threshold is a named constant with its source.
"""
import os
import hashlib
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

WS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WS))
sys.path.insert(0, str(WS / "scripts"))

DIMER_ROOT = Path(os.environ.get("DIMERS_ROOT", "/path/to/afdb_dimers")) / "afdb_human_dimers"  # docs/DATA_RECON.md, task brief
KINDS = {"homo": "homodimers", "hetero": "heterodimers"}
PDB_MIRROR = Path(os.environ.get("PDB_MIRROR", "/path/to/pdb_mirror"))  # <id[1:3]>/<id>.cif.gz (docs/DATA_RECON.md)
POCKETS_DIMER = WS / "data/part2/human_dimer_pockets_annotated.parquet"  # scripts/annotate_pocket_chains.py
POCKETS_ISO = WS / "data/part2/pockets/human_isolated"                     # scripts/run_dimer_pockets.py (isolated chains)
ATLAS = WS / "data/raw/pocketcompare_afdb_snapshot"                         # scripts/snapshot_atlas.sh
OUT_DATA = WS / "data/part2/h4"
OUT_RES = WS / "results/part2/h4"
OUT_FIG = WS / "results/part2/figures"

# ---- pre-registered thresholds (DESIGN.md v3 §3 Q2.1, §1 H4)
DCC_HIT = 4.0              # Å, ligand centroid vs pocket centre
TOPK_EXTRA = 2             # top-(n+2)
IFACE_DIST = 5.0           # Å heavy-atom distance for interface residues (task brief step 2)
IFACE_JACCARD_MIN = 0.3    # interface-correspondence filter
RESNAME_IDENTITY_MIN = 0.9 # per-site residue-name identity PDB vs UniProt (SIFTS)
SEEDS = [0, 1, 2]          # DESIGN §5; seed 0 primary
N_BOOT = 1000              # DESIGN §5
N_PERM = 10000             # permutation resamples (not fixed by DESIGN; chosen for p resolution 1e-4)
POWER_GATE_CLUSTERS = 30   # DESIGN §1 power gate
HETERO_MIN_MBS = 100       # DESIGN §3: hetero <100 non-redundant MBS -> descriptive only
AF_TRAINING_CUTOFF = "2021-09-30"  # LITERATURE.md (post-cutoff split, Runs N' Poses / AF3); AF-Multimer v2.3 PDB cutoff
                                   # is believed to be the same date [UNVERIFIED for the AFDB-complex models]

# ---- analysis choices (not fixed by DESIGN; documented in docs/H4_RESULTS.md)
SUPERPOSE_RADIUS = 15.0    # Å: CA atoms (both chains) within this distance of the ligand centroid define the local frame
SUPERPOSE_MIN_CA = 10      # fewer mapped CA pairs -> fall back to the whole mapped chain pair (flagged)
SURFACE_NB_RADIUS = 10.0   # Å, CB neighbour count radius for the surface definition of null patches
SURFACE_NB_MAX = 20        # residues with < this many CB neighbours are 'surface' seeds (half-sphere-exposure-like rule)
MAX_ENTRIES_PER_UNIT = 3   # best-resolution entries tried per non-redundant unit (fallbacks if the first fails mapping)

ORGANIC_EXCLUDE = {"metal_ion", "metal_cluster", "inorganic"}   # organic = carbon-containing classes
METAL_CLASSES = {"metal_ion", "metal_cluster"}                   # reported separately (DESIGN §0)


def census_paths(version):
    """v1 -> data/part1/v1/; v2 -> data/part1/ (must carry census_version == 'v2')."""
    base = WS / "data/part1" / ("v1" if version == "v1" else "")
    return base / "ligand_sites.parquet", base / "entries.parquet"


def min_atoms_primary(version):
    # binding residue rule of the census primary: v1 >=1 atom pair, v2 >=2 (docs/PART1_PIPELINE.md v2 changes 1)
    return 1 if version == "v1" else 2


def site_group_col(version):
    return "sym_group" if version == "v1" else "site_group"


def provenance():
    d = hashlib.sha256((WS / "DESIGN.md").read_bytes()).hexdigest()
    try:
        c = subprocess.run(["git", "-C", str(WS), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:
        c = "unknown"
    return {"design_sha256": d, "git_commit": c}


def write_tsv(df, path, extra=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    prov = provenance()
    if extra:
        prov.update(extra)
    with open(path, "w") as f:
        f.write("# " + " ".join(f"{k}={v}" for k, v in prov.items()) + " status=UNVERIFIED\n")
        df.to_csv(f, sep="\t", index=False)


def kabsch(P, Q):
    """Rotation R, translation t minimising |(P @ R.T + t) - Q|. Returns R, t, rmsd."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = qc - pc @ R.T
    rmsd = float(np.sqrt(((P @ R.T + t - Q) ** 2).sum(1).mean()))
    return R, t, rmsd


def parse_res(s):
    """'A:12 B:40' -> set of (chain, int pos); insertion codes dropped (AF models have none)."""
    out = set()
    for x in (s or "").split():
        c, r = x.split(":")
        try:
            out.add((c, int(r)))
        except ValueError:
            out.add((c, int("".join(ch for ch in r if ch.isdigit() or ch == "-"))))
    return out


def load_candidates():
    rows = []
    for kind, d in KINDS.items():
        c = pd.read_csv(DIMER_ROOT / d / "candidates.csv")
        c["kind"] = kind
        rows.append(c)
    return pd.concat(rows, ignore_index=True)
