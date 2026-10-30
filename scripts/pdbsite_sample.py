#!/usr/bin/env python
"""Stratified non-redundant sample of PDB ligand sites for pocket description (H2 secondary + EXPLORATORY).

Unit = (cluster30 x ccd_key), organic (non-metal) oligomer units of the Part 1 primary set
(data/part1/analysis/units_primary_v2.parquet, scripts/part1_analysis.py). One site per unit: a symmetric-collapsed
site in the unit's representative entry with the same ccd_key whose MBS label equals the unit label
(unit MBS = any site MBS), chosen at random (seed). Strata = MBS/SBS x ligand stratum (drug_like / cofactor / other);
equal allocation per stratum (all units if a stratum is smaller), remainder redistributed.
Output: data/interim/pdbsite/sample_n{N}_seed{S}.parquet (+ results/part1/pocket_descriptors/sample_strata_*.tsv)
Run with .venvs/part1/bin/python (needs scripts/part1_analysis.py loaders).
Status: first attempt, not yet verified.
"""
import argparse, sys
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import part1_analysis as P1
from pdbsite_common import UNITS, SITES, INTERIM, OUT, SEED, MAX_SITES, lig_stratum, write_tsv

SITE_COLS = ["ligand_uid", "pdb_id", "ccd", "ccd_key", "components", "label_asym", "auth_seq", "op", "centroid",
             "n_heavy_atoms_instance", "chain_counts", "contacts", "is_on_axis", "is_on_axis_geom", "is_mbs_min3",
             "dominant_chain", "assembly_id"]


def allocate(sizes, n):
    alloc = {k: 0 for k in sizes}
    left, open_ = n, set(sizes)
    while left > 0 and open_:
        share = max(left // len(open_), 1)
        for k in sorted(open_):
            add = min(share, sizes[k] - alloc[k], left)
            alloc[k] += add; left -= add
        open_ = {k for k in open_ if alloc[k] < sizes[k]}
    return alloc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=MAX_SITES)
    ap.add_argument("--seed", type=int, default=SEED)
    a = ap.parse_args()
    assert a.n <= MAX_SITES
    s = P1.load(P1.CENSUS_DIR)
    sites = P1.site_filter(s, dict(P1.PRIMARY))          # collapsed oligomer sites, primary variant
    u = pd.read_parquet(UNITS)
    u = u[~u.is_metal].copy()
    u["stratum_lig"] = [lig_stratum(c, m) for c, m in zip(u.cls, u.marsh_class)]
    u["stratum"] = np.where(u.y_mbs == 1, "MBS", "SBS") + "|" + u.stratum_lig
    # candidate sites per unit (representative entry, same ccd_key, same MBS label as the unit)
    c = sites.merge(u[["cluster30", "ccd_key", "pdb_id", "y_mbs"]].rename(columns={"y_mbs": "unit_mbs"}),
                    on=["cluster30", "ccd_key", "pdb_id"])
    c = c[c.y_mbs.astype(int) == c.unit_mbs]
    rng = np.random.default_rng(a.seed)
    c["r"] = rng.random(len(c))
    pick = c.sort_values("r").drop_duplicates(["cluster30", "ccd_key"])[["cluster30", "ccd_key", "ligand_uid", "symkey"]]
    u = u.merge(pick, on=["cluster30", "ccd_key"], how="left")
    n_nosite = int(u.ligand_uid.isna().sum())
    u = u.dropna(subset=["ligand_uid"])
    sizes = u.stratum.value_counts().to_dict()
    alloc = allocate(sizes, a.n)
    parts = [u[u.stratum == k].sample(n=alloc[k], random_state=a.seed) for k in sorted(alloc)]
    smp = pd.concat(parts, ignore_index=True)
    # site-level detail from the census
    d = pd.read_parquet(SITES, columns=SITE_COLS, filters=[("ligand_uid", "in", smp.ligand_uid.tolist())])
    d = d.drop_duplicates("ligand_uid")
    smp = smp.merge(d.drop(columns=["pdb_id", "ccd", "ccd_key"]), on="ligand_uid", how="left")
    smp = smp.sample(frac=1, random_state=a.seed).reset_index(drop=True)   # random processing order
    smp["site_id"] = [f"{p}_{i:05d}" for i, p in enumerate(smp.pdb_id)]
    out = INTERIM / f"sample_n{a.n}_seed{a.seed}.parquet"
    smp.to_parquet(out, index=False)
    tab = pd.DataFrame({"stratum": sorted(sizes), "units_available": [sizes[k] for k in sorted(sizes)],
                        "sampled": [alloc[k] for k in sorted(sizes)],
                        "clusters_sampled": [smp[smp.stratum == k].cluster30.nunique() for k in sorted(sizes)]})
    write_tsv(tab, OUT / f"sample_strata_n{a.n}_seed{a.seed}.tsv", note=f"units_without_matching_site={n_nosite}")
    print(tab.to_string(index=False)); print("units without a matching site:", n_nosite, "->", out)


if __name__ == "__main__":
    main()
