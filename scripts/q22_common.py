"""Shared constants and helpers for the Q2.2 / Q2.3a analyses (EXPLORATORY, DESIGN.md v3 §3).

Status: first attempt, not yet verified.
"""
import os
import hashlib, subprocess
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
DIMER_ROOT = Path(os.environ.get("DIMERS_ROOT", "/cluster/scratch/gmagnusson/dimers")) / "afdb_human_dimers"
KINDS = {"homo": "homodimers", "hetero": "heterodimers"}
ANNOT = WS / "data/part2/human_dimer_pockets_annotated.parquet"
POCKETS_DIMER = WS / "data/part2/pockets/human"
POCKETS_ISO = WS / "data/part2/pockets/human_isolated"
ATLAS = WS / "data/raw/pocketcompare_afdb_snapshot/pockets"
INTERIM = WS / "data/interim/q22"
OUT = WS / "results/part2/q22"
FIG = WS / "results/part2/figures"
for p in (INTERIM, OUT, FIG):
    p.mkdir(parents=True, exist_ok=True)

HUMAN_TAXID = 9606                 # NCBI taxonomy id (atlas column `taxid`)
HUMAN_PROTEOME = "UP000005640"     # UniProt reference proteome (AFDB tier-1 human tar name)
ATLAS_SCORE, ATLAS_PLDDT = 0.03, 70.0          # DESIGN.md §0 atlas filter
DIMER_ONLY_THRESHOLDS = (0.2, 0.5, 0.8)       # DESIGN.md §0 (0.5 primary)
DIMER_ONLY_PRIMARY = 0.5
CONTACT_A = 5.0                    # task spec: interface residue = heavy atom within 5 Å of partner chain
REL_SASA_SURFACE = 0.25            # surface residue: relative SASA >= 25% in isolated chain (Levy 2010 convention)
DRUGLIKE_VOL = (300.0, 1500.0)     # task spec fallback (no PDB drug-like pocket volumes available yet), on volume_A3
PARALOG_MIN_ID = 0.30              # DESIGN.md §3: exclude >=30% identity paralogs
PARALOG_MAX_EVALUE = 1e-3          # mmseqs2 hit e-value cut-off (our choice, documented)
SEEDS = (0, 1, 2)                  # DESIGN.md §5; 0 primary
N_BOOT = 1000                      # DESIGN.md §5

# residue classes for lining chemistry (side-chain based; our choice, documented in docs/Q22_RESULTS.md)
HYDROPHOBIC = set("AVLIMFWC")
AROMATIC = set("FWYH")
CHARGED = set("DEKR")


def provenance():
    """design sha256 + git commit, written into every output."""
    sha = hashlib.sha256((WS / "DESIGN.md").read_bytes()).hexdigest()
    try:
        commit = subprocess.check_output(["git", "-C", str(WS), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        commit = "unknown"
    return {"design_sha256": sha, "git_commit": commit}


def write_tsv(df, path, note=""):
    pv = provenance()
    path = Path(path)
    with open(path, "w") as fh:
        fh.write(f"# design_sha256={pv['design_sha256']} git_commit={pv['git_commit']} status=UNVERIFIED EXPLORATORY {note}\n")
        df.to_csv(fh, sep="\t", index=False)


def oc(a, b):
    return len(a & b) / min(len(a), len(b)) if a and b else 0.0


def load_annot():
    d = pd.read_parquet(ANNOT)
    d["is_sym_rep"] = d.is_sym_rep.astype(bool)
    return d


def ipsae_tertile(x, edges=None):
    if edges is None:
        edges = np.quantile(x, [1 / 3, 2 / 3])
    return pd.cut(x, [-np.inf, edges[0], edges[1], np.inf], labels=["T1_low", "T2_mid", "T3_high"]), edges


def cliffs_delta(x, y):
    """P(X>Y) - P(X<Y) via ranks (O(n log n))."""
    from scipy.stats import rankdata
    x, y = np.asarray(x, float), np.asarray(y, float)
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan
    r = rankdata(np.concatenate([x, y]))
    u = r[:nx].sum() - nx * (nx + 1) / 2
    return 2 * u / (nx * ny) - 1


def cluster_boot(df, unit, stat, n_boot=N_BOOT, seed=0):
    """Cluster bootstrap: resample `unit` values with replacement, stat(df_resampled) -> float. Returns (est, lo, hi)."""
    rng = np.random.default_rng(seed)
    est = stat(df)
    groups = df.groupby(unit).indices
    keys = np.array(list(groups.keys()), dtype=object)
    idx_list = [groups[k] for k in keys]
    vals = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(keys), len(keys))
        idx = np.concatenate([idx_list[i] for i in pick])
        vals.append(stat(df.iloc[idx]))
    vals = np.array(vals, float)
    return est, np.nanpercentile(vals, 2.5), np.nanpercentile(vals, 97.5)
