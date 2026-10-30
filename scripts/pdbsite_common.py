"""Shared constants for the PDB-site pocket descriptor run (H2 secondary descriptor + EXPLORATORY physicochemistry).

Status: first attempt, not yet verified.
"""
import os
import hashlib, subprocess
from pathlib import Path

WS = Path(__file__).resolve().parents[1]
PDB_MIRROR = Path(os.environ.get("PDB_MIRROR", "/path/to/pdb_mirror"))          # same mirror as src/part1/config.py (layout <id[1:3]>/<id>.cif.gz)
UNITS = WS / "data/part1/analysis/units_primary_v2.parquet"  # scripts/part1_analysis.py output (census v2, primary)
SITES = WS / "data/part1/ligand_sites.parquet"
INTERIM = WS / "data/interim/pdbsite"
OUT = WS / "results/part1/pocket_descriptors"
FIG = WS / "results/part1/figures"
FP_OUT = WS / "data/part2/pdb_site_fingerprints.parquet"
Q22_PROPS = WS / "results/part2/q22/q22_props.tsv"          # predicted-dimer comparison (docs/Q22_RESULTS.md)
for p in (INTERIM, OUT, FIG):
    p.mkdir(parents=True, exist_ok=True)

SEED = 0                       # DESIGN.md §5: seed 0 primary
SEEDS = (0, 1, 2)
N_BOOT = 1000                  # DESIGN.md §5
MAX_SITES = 12000              # task cap
POWER_GATE_MIN_CLUSTERS = 30   # DESIGN.md §1 power gate
SPLINE_DF = 4                  # same B-spline df for heavy atoms as part1_analysis.py (our choice there)
CROP_A = 30.0                  # SENSITIVITY only (first run): keep chains with a heavy atom <= 30 Å of the ligand; primary = uncropped
BIOLIP_MIN_PAIRS = 2           # src/part1/config.py MIN_ATOM_CONTACTS_PRIMARY: vdW+0.5 Å atom pairs per residue (BioLiP2 rule)
LIG_MATCH_A = 1.5              # our choice: max distance between census centroid and re-built ligand centroid
INTERFACE_MIN_RES = 3          # DESIGN.md §0 interface pocket: >=3 lining residues on each of >=2 chains

# residue classes for lining chemistry: identical to scripts/q22_common.py (apples-to-apples with Q2.2)
HYDROPHOBIC = set("AVLIMFWC")
AROMATIC = set("FWYH")
CHARGED = set("DEKR")

# ligand class strata (task): drug-like (primary ChEMBL rule, census `cls`), cofactor (Marsh cofactor class:
# PDBe cofactor + heme), other (all remaining organic classes)
def lig_stratum(cls, marsh_class):
    if cls == "drug_like":
        return "drug_like"
    if marsh_class == "cofactor":
        return "cofactor"
    return "other"


def provenance():
    sha = hashlib.sha256((WS / "DESIGN.md").read_bytes()).hexdigest()
    try:
        commit = subprocess.check_output(["git", "-C", str(WS), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        commit = "unknown"
    return {"design_sha256": sha, "git_commit": commit}


def write_tsv(df, path, note=""):
    pv = provenance()
    with open(path, "w") as fh:
        fh.write(f"# design_sha256={pv['design_sha256']} git_commit={pv['git_commit']} status=UNVERIFIED "
                 f"H2-secondary(burrow_score)+EXPLORATORY {note}\n")
        df.to_csv(fh, sep="\t", index=False, float_format="%.6g")
    print("wrote", path, len(df), flush=True)
