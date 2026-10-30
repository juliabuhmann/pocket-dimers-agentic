#!/usr/bin/env python
"""H4 / Q2.1 step 1 feasibility gate (run BEFORE scripts/h4_recovery.py). First attempt, NOT verified.

From data/part2/h4/<census>/sites.parquet (h4_select.py) apply the mapping filters and count non-redundant sites
(unit = model x cluster30 x CCD x MBS label; symmetric copies collapsed by the census site group) per homo/hetero x
MBS/SBS x ligand set, and the number of bootstrap clusters. Hetero < HETERO_MIN_MBS organic MBS units -> hetero is
descriptive only (DESIGN §3). Writes results/part2/h4/<census>/feasibility.tsv and gate.json.
"""
import argparse
import json

import pandas as pd

from h4_common import OUT_DATA, OUT_RES, RESNAME_IDENTITY_MIN, HETERO_MIN_MBS, POWER_GATE_CLUSTERS, write_tsv, provenance


def mapping_filter(s, ledger):
    ledger.append(("site x model rows from selection", len(s)))
    s = s[s.status == "ok"]
    ledger.append(("chain pair resolved in assembly 1 (status ok)", len(s)))
    s = s[(s.chain1_sifts_acc == s.acc1) & (s.chain2_sifts_acc == s.acc2)]
    ledger.append(("both chains map to the model accessions in SIFTS XML", len(s)))
    s = s[s.resname_identity >= RESNAME_IDENTITY_MIN]
    ledger.append((f"binding-residue name identity >= {RESNAME_IDENTITY_MIN}", len(s)))
    return s


def collapse_and_units(s):
    """One row per (entry, model, census site group); then one representative entry per unit (best resolution, id)."""
    s = s.sort_values(["pdb_id", "model", "site_group", "ligand_uid"]).drop_duplicates(["pdb_id", "model", "site_group"])
    return s


def boot_clusters(s):
    """Connected components over (model, cluster30): sites sharing a model or a 30% cluster are one resampling unit."""
    parent = {}
    def f(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for m, c in zip(s.model, s.cluster30.astype(str)):
        parent[f("m:" + m)] = f("c:" + c)
    return s.model.map(lambda m: f("m:" + m))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True)
    a = ap.parse_args()
    s = pd.read_parquet(OUT_DATA / a.census / "sites.parquet")
    led = []
    s = mapping_filter(s, led)
    s = collapse_and_units(s)
    led.append(("after symmetric-copy collapse (entry x model x site group)", len(s)))
    s = s.copy()
    s["boot_cluster"] = boot_clusters(s)
    rows = []
    for (kind, ls, mbs), g in s.groupby(["kind", "lig_set", "is_mbs"]):
        rows.append(dict(kind=kind, lig_set=ls, site_type="MBS" if mbs else "SBS", n_rows=len(g), n_units=g.unit.nunique(),
                         n_models=g.model.nunique(), n_cluster30=g.cluster30.nunique(), n_boot_clusters=g.boot_cluster.nunique(),
                         n_entries=g.pdb_id.nunique()))
    ft = pd.DataFrame(rows)
    write_tsv(ft, OUT_RES / a.census / "feasibility.tsv", {"census": a.census, "note": "counted_before_recovery"})
    write_tsv(pd.DataFrame(led, columns=["step", "n_rows"]), OUT_RES / a.census / "mapping_ledger.tsv", {"census": a.census})
    org = ft[ft.lig_set == "organic"].set_index(["kind", "site_type"])
    get = lambda k, t, c: int(org.loc[(k, t), c]) if (k, t) in org.index else 0
    gate = dict(census=a.census, **provenance(),
                homo_mbs_units=get("homo", "MBS", "n_units"), homo_sbs_units=get("homo", "SBS", "n_units"),
                hetero_mbs_units=get("hetero", "MBS", "n_units"), hetero_sbs_units=get("hetero", "SBS", "n_units"),
                homo_mbs_clusters=get("homo", "MBS", "n_boot_clusters"), hetero_mbs_clusters=get("hetero", "MBS", "n_boot_clusters"))
    gate["hetero_inferential"] = gate["hetero_mbs_units"] >= HETERO_MIN_MBS and gate["hetero_mbs_clusters"] >= POWER_GATE_CLUSTERS
    gate["homo_inferential"] = gate["homo_mbs_clusters"] >= POWER_GATE_CLUSTERS
    json.dump(gate, open(OUT_RES / a.census / "gate.json", "w"), indent=1)
    print(pd.DataFrame(led).to_string()); print(ft.to_string()); print(json.dumps(gate, indent=1))


if __name__ == "__main__":
    main()
