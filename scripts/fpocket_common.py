"""Shared helpers for the fpocket druggability analysis (EXPLORATORY, not pre-registered; docs/FPOCKET_DRUGGABILITY.md).

fpocket is installed in conda_envs/fpocket (bioconda); the version is read from the conda metadata at run time.
Status: first attempt, not yet verified.
"""
import glob, json, re, shutil, subprocess, tempfile
from pathlib import Path
import numpy as np

WS = Path(__file__).resolve().parents[1]
FPOCKET_ENV = WS / "conda_envs/fpocket"
FPOCKET_BIN = FPOCKET_ENV / "bin/fpocket"
OUT = WS / "results/part1/fpocket"
FIG = OUT / "figures"
INTERIM = WS / "data/interim/fpocket"
for p in (OUT, FIG, INTERIM):
    p.mkdir(parents=True, exist_ok=True)

OVERLAP_A = 3.0          # task spec: ligand heavy atom within 3 Å of a pocket alpha-sphere centre
FPOCKET_TIMEOUT_S = 900   # our choice (pilot max 331 s); a timeout is logged as an error, never dropped silently
DRUG_THR = 0.5           # fpocket's conventional druggability threshold (Schmidtke & Barril 2010, J Med Chem 53:5858)

# fields of <name>_info.txt we keep (fpocket 4 labels -> our column names)
INFO_FIELDS = {"Score": "fp_score", "Druggability Score": "fp_drug", "Number of Alpha Spheres": "fp_n_spheres",
               "Total SASA": "fp_sasa", "Polar SASA": "fp_polar_sasa", "Apolar SASA": "fp_apolar_sasa",
               "Volume": "fp_volume", "Mean local hydrophobic density": "fp_hyd_density",
               "Hydrophobicity score": "fp_hydrophobicity", "Volume score": "fp_volume_score",
               "Polarity score": "fp_polarity", "Charge score": "fp_charge",
               "Proportion of polar atoms": "fp_prop_polar_atoms", "Apolar alpha sphere proportion": "fp_apolar_prop"}


def fpocket_version():
    metas = sorted(glob.glob(str(FPOCKET_ENV / "conda-meta/fpocket-*.json")))
    if not metas:
        return "unknown"
    m = json.load(open(metas[-1]))
    return f"fpocket {m['version']} ({m.get('channel', '?')}, build {m.get('build', '?')})"


def parse_info(path):
    pockets, cur = {}, None
    for line in open(path):
        m = re.match(r"Pocket (\d+) :", line)
        if m:
            cur = int(m.group(1)); pockets[cur] = {"fp_pocket": cur}
            continue
        if cur is not None and ":" in line:
            k, v = [s.strip() for s in line.split(":", 1)]
            if k in INFO_FIELDS:
                pockets[cur][INFO_FIELDS[k]] = float(v)
    return pockets


def parse_spheres(out_pdb):
    """Alpha-sphere centres per pocket from <name>_out.pdb (HETATM STP, residue number = pocket number)."""
    cen = {}
    for line in open(out_pdb):
        if line.startswith("HETATM") and line[17:20] == "STP":
            k = int(line[22:26])
            cen.setdefault(k, []).append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return {k: np.asarray(v) for k, v in cen.items()}


def parse_residues(pockets_dir):
    """Lining residues per pocket from pockets/pocketN_atm.pdb as 'chain:resnum' (author numbering)."""
    res = {}
    for f in glob.glob(str(Path(pockets_dir) / "pocket*_atm.pdb")):
        k = int(re.search(r"pocket(\d+)_atm", f).group(1))
        s = set()
        for line in open(f):
            if line.startswith(("ATOM", "HETATM")):
                s.add(f"{line[21].strip()}:{int(line[22:26])}")
        res[k] = s
    return res


def run_fpocket(pdb_path, keep_residues=False, extra_args=()):
    """Run fpocket on a PDB file inside a private temp copy; return (pockets dict, sphere centres, residues|None)."""
    pdb_path = Path(pdb_path)
    with tempfile.TemporaryDirectory(prefix="fpk_") as td:   # honours $TMPDIR (node-local is fastest)
        p = Path(td) / ("s" + pdb_path.suffix)
        shutil.copy(pdb_path, p)
        r = subprocess.run([str(FPOCKET_BIN), "-f", str(p), *extra_args], capture_output=True, text=True,
                           timeout=FPOCKET_TIMEOUT_S)
        od = Path(td) / "s_out"
        if not (od / "s_info.txt").exists():
            raise RuntimeError(f"fpocket failed rc={r.returncode}: {(r.stderr or r.stdout)[-300:]}")
        pk = parse_info(od / "s_info.txt")
        cen = parse_spheres(od / "s_out.pdb")
        res = parse_residues(od / "pockets") if keep_residues else None
    return pk, cen, res

