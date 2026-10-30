"""Independent re-derivation of the Part 1 unit table (verifier; does NOT import scripts/part1_analysis.py).

Definitions are taken from DESIGN.md v3 §0/§1 and docs/PART1_PIPELINE.md (v2 changes) only.
Status: verifier's first attempt; numbers are compared against the producer's ledger in the review.

Steps (each count written to the verifier ledger):
  instance filters -> oligomer -> cluster30 (re-derived from RCSB clusters-by-entity-30.txt using
  pdb_id + dominant_entity) -> symmetric-site collapse (pdb_id, site_group) -> unit (cluster30 x ccd_key)
  with representative entry = best resolution, ties -> earliest PDB ID.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SITES = ROOT / "data/part1/ligand_sites.parquet"
ENTRIES = ROOT / "data/part1/entries.parquet"
DRUGLIKE = ROOT / "data/processed/ccd_druglike.parquet"
CLUSTERS30 = ROOT / "data/raw/part1/clusters-by-entity-30.txt"
OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

OLIGOMER_TYPES = ("homomer", "heteromer")  # DESIGN §1 family A: oligomers only; same_uniprot_hetero excluded (PIPELINE v2 #6)

COLS = ["pdb_id", "ligand_uid", "method_class", "resolution", "release_date", "assembly_type",
        "n_protein_chains", "ccd", "ccd_key", "n_heavy_atoms_instance", "is_glycan", "is_dup_copy",
        "is_bound", "is_bound_4A", "is_mbs", "is_mbs_min3", "is_mbs_4A", "is_on_axis", "is_on_axis_geom",
        "on_axis_jaccard", "c2_axis_dist", "dominant_entity", "dominant_uniprot", "contacts_antibody",
        "biolip_covered", "is_additive", "plinder_artifact", "site_group", "ligand_class", "marsh_class",
        "cluster30", "in_primary_set", "is_photosynthetic"]


def load_cluster30():
    m = {}
    with open(CLUSTERS30) as fh:
        for i, line in enumerate(fh):
            for tok in line.split():
                m[tok.upper()] = i
    return m


def main():
    d = pd.read_parquet(SITES, columns=COLS)
    ledger = [("all ligand instances", len(d), None)]

    # --- independent cluster30 from the RCSB file
    c30 = load_cluster30()
    key = d.pdb_id.str.upper() + "_" + d.dominant_entity.astype(str)
    d["my_cluster30"] = key.map(c30)

    steps = [
        ("not duplicate copy", ~d.is_dup_copy),
        ("bound (>=2 res with >=2 atom contacts)", d.is_bound),
        ("not glycan", ~d.is_glycan),
        ("BioLiP2-covered entry", d.biolip_covered.fillna(False)),
        ("not additive", ~d.is_additive.fillna(True)),
        ("no antibody contact", ~d.contacts_antibody.fillna(False)),
        ("oligomer", d.assembly_type.isin(OLIGOMER_TYPES)),
        ("cluster30 assigned (verifier, RCSB file)", d.my_cluster30.notna()),
    ]
    mask = pd.Series(True, index=d.index)
    for name, m in steps:
        mask &= m.astype(bool)
        ledger.append((name, int(mask.sum()), int(d.loc[mask, "my_cluster30"].nunique())))
    p = d[mask].copy()
    p["my_cluster30"] = p.my_cluster30.astype(int)

    # agreement checks with census columns
    agree_c30 = (p.my_cluster30 == p.cluster30).mean()
    pre = d[mask.values | d.in_primary_set.fillna(False).values]
    print(f"cluster30 agreement (verifier vs census) on primary: {agree_c30:.5f}")
    base = (~d.is_dup_copy) & d.is_bound & ~d.is_glycan & d.biolip_covered.fillna(False) & ~d.is_additive.fillna(True)
    print("in_primary_set vs my base filter disagreements:", int((base != d.in_primary_set.fillna(False)).sum()))

    # --- drug-like flags (primary) via ccd_key
    dl = pd.read_parquet(DRUGLIKE, columns=["ccd", "name", "druglike_primary", "druglike_ro5qed"])
    p = p.merge(dl.rename(columns={"ccd": "ccd_key"}), on="ccd_key", how="left")
    p["druglike_primary"] = p.druglike_primary.fillna(False).astype(bool)
    p["druglike_ro5qed"] = p.druglike_ro5qed.fillna(False).astype(bool)
    # fine class: drug-like (ChEMBL primary) overrides; census Ro5-based 'drug_like' that fails primary -> other_organic
    fine = p.ligand_class.where(p.ligand_class != "drug_like", "other_organic")
    p["fine"] = np.where(p.druglike_primary, "drug_like", fine)

    # --- symmetric-site collapse: (pdb_id, site_group)
    g = p.groupby(["pdb_id", "site_group"], sort=False)
    sites = g.agg(my_cluster30=("my_cluster30", "first"), n_c30=("my_cluster30", "nunique"),
                  ccd_key=("ccd_key", "first"), fine=("fine", "first"), marsh_class=("marsh_class", "first"),
                  name=("name", "first"), druglike_ro5qed=("druglike_ro5qed", "first"),
                  assembly_type=("assembly_type", "first"), resolution=("resolution", "first"),
                  method_class=("method_class", "first"), release_date=("release_date", "first"),
                  n_protein_chains=("n_protein_chains", "first"),
                  heavy=("n_heavy_atoms_instance", "max"),
                  is_mbs=("is_mbs", "any"), is_mbs_min3=("is_mbs_min3", "any"),
                  on_axis=("is_on_axis", "any"), on_axis_geom=("is_on_axis_geom", "any"),
                  is_photosynthetic=("is_photosynthetic", "first"),
                  dominant_entity=("dominant_entity", "first"), dominant_uniprot=("dominant_uniprot", "first"),
                  n_inst=("ligand_uid", "size")).reset_index()
    print("site groups spanning >1 cluster30:", int((sites.n_c30 > 1).sum()))
    # on-axis only meaningful for MBS sites
    sites["on_axis_mbs"] = sites.is_mbs & sites.on_axis
    sites["on_axis_geom_mbs"] = sites.is_mbs & sites.on_axis_geom
    ledger.append(("sites (pdb_id x site_group)", len(sites), int(sites.my_cluster30.nunique())))

    # --- units: cluster30 x ccd_key; representative entry = best resolution, then earliest pdb_id
    ent = sites.groupby(["my_cluster30", "ccd_key", "pdb_id"], sort=False).agg(res=("resolution", "first")).reset_index()
    ent["res_sort"] = ent.res.fillna(np.inf)
    ent = ent.sort_values(["my_cluster30", "ccd_key", "res_sort", "pdb_id"])
    rep = ent.drop_duplicates(["my_cluster30", "ccd_key"])[["my_cluster30", "ccd_key", "pdb_id"]]
    # "MBS in any entry" sensitivity
    any_mbs = sites.groupby(["my_cluster30", "ccd_key"]).is_mbs.any().rename("mbs_any_entry").reset_index()
    rs = sites.merge(rep, on=["my_cluster30", "ccd_key", "pdb_id"])
    units = rs.groupby(["my_cluster30", "ccd_key"]).agg(
        pdb_id=("pdb_id", "first"), fine=("fine", "first"), marsh_class=("marsh_class", "first"),
        name=("name", "first"), druglike_ro5qed=("druglike_ro5qed", "first"),
        assembly_type=("assembly_type", "first"), resolution=("resolution", "first"),
        method_class=("method_class", "first"), release_date=("release_date", "first"),
        n_protein_chains=("n_protein_chains", "first"), heavy=("heavy", "max"),
        mbs=("is_mbs", "any"), mbs_strict=("is_mbs_min3", "any"),
        on_axis=("on_axis_mbs", "any"), on_axis_geom=("on_axis_geom_mbs", "any"),
        is_photosynthetic=("is_photosynthetic", "first"),
        dominant_entity=("dominant_entity", "first"), dominant_uniprot=("dominant_uniprot", "first"),
        n_sites=("site_group", "size")).reset_index()
    units = units.merge(any_mbs, on=["my_cluster30", "ccd_key"])
    e = pd.read_parquet(ENTRIES, columns=["pdb_id", "is_structural_genomics", "pisa_agrees"])
    units = units.merge(e, on="pdb_id", how="left")
    units["year"] = pd.to_datetime(units.release_date).dt.year
    ledger.append(("units (cluster30 x ccd_key, rep entry)", len(units), int(units.my_cluster30.nunique())))
    org = units[~units.fine.isin(["metal_ion", "metal_cluster", "inorganic"])]
    ledger.append(("organic units (excl metal_ion/metal_cluster/inorganic)", len(org), int(org.my_cluster30.nunique())))
    org2 = units[units.marsh_class != "metal"]
    ledger.append(("units with marsh_class != metal", len(org2), int(org2.my_cluster30.nunique())))
    org3 = units[units.fine != "metal_ion"]
    ledger.append(("units with fine != metal_ion", len(org3), int(org3.my_cluster30.nunique())))
    print(units.assembly_type.value_counts())

    led = pd.DataFrame(ledger, columns=["step", "n", "n_clusters"])
    led.to_csv(OUT / "verifier_ledger.tsv", sep="\t", index=False)
    print(led.to_string())
    units.to_parquet(OUT / "units_verifier.parquet")
    sites.to_parquet(OUT / "sites_verifier.parquet")


if __name__ == "__main__":
    sys.exit(main())
