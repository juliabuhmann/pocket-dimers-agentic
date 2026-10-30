#!/usr/bin/env python
"""Q2.3a (EXPLORATORY): human proteins whose only confident pocket of drug-like size appears in a dimer (dimer-only).

Qualifying dimer pocket (for each protein P of the model): interface (>=3 lining residues per chain) AND dimer-only
(max overlap coefficient with any isolated-chain pocket < 0.5) AND atlas filter (burrow_score >= 0.03, lining pLDDT
>= 70) AND volume_A3 in DRUGLIKE_VOL (300-1500 A^3; task-specified fallback because no burrow volumes of PDB
drug-like ligand pockets exist yet in this project).
P is listed if it has >= 1 qualifying dimer pocket in any model AND none of the following has a pocket that passes the
atlas filter with volume_A3 in DRUGLIKE_VOL:
  (i) P's monomer atlas pockets (human atlas snapshot, any fragment),
  (ii) isolated-chain pockets of P's chain(s) in any of P's dimer models,
  (iii) single-chain dimer pockets lying entirely on P's chain in any of P's models.
Primary list requires P to be present in the atlas structure table (status ok); others are flagged in_atlas=False.
Ranking: best qualifying pocket burrow_score (desc), then max_ipSAE.
Annotation: gene / protein name (UniProt REST, scripts/q22_fetch_uniprot.py), ChEMBL status (scripts/q22_chembl_targets.py).
Output: results/part2/q23a_dimer_only_ligandable.tsv (+ results/part2/q22/q23a_n_steps.tsv)
Status: first attempt, not yet verified.
"""
import glob, sys
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from q22_common import (WS, INTERIM, OUT, ATLAS, HUMAN_TAXID, DRUGLIKE_VOL, ATLAS_SCORE, ATLAS_PLDDT,
                        DIMER_ONLY_PRIMARY, write_tsv)

import argparse
ap = argparse.ArgumentParser()
ap.add_argument("--vol", nargs=2, type=float, default=DRUGLIKE_VOL, help="volume_A3 range (sensitivity runs)")
args = ap.parse_args()
lo, hi = args.vol
PRIMARY = tuple(args.vol) == tuple(DRUGLIKE_VOL)
SFX = "" if PRIMARY else f"_vol{int(lo)}-{int(hi)}"
d = pd.read_parquet(INTERIM / "dimer_pocket_features.parquet")
iso = pd.read_parquet(INTERIM / "iso_pocket_features.parquet")
am = pd.read_parquet(INTERIM / "atlas_human_meta.parquet")
genes = pd.read_csv(INTERIM / "uniprot_genes.tsv", sep="\t").set_index("Entry")
cf = INTERIM / "chembl_target_status.tsv"  # scripts/q22_chembl_targets.py (run after the candidate list exists)
chem = pd.read_csv(cf, sep="\t").set_index("acc") if cf.exists() else pd.DataFrame(
    columns=["chembl_status", "single_protein_targets", "n_act_pchembl6", "mech_max_phase", "chembl_release"])
st = pd.concat([pd.read_parquet(f, columns=["accession", "taxid", "status", "n_pockets_kept"])
                for f in glob.glob(str(ATLAS.parent / "structures" / "tier[12]" / "*.parquet"))])
st = st[(st.taxid == HUMAN_TAXID)]
in_atlas = set(st.loc[st.status == "ok", "accession"])
steps = []
sizeok = lambda x: x.volume_A3.between(lo, hi)
passes = lambda x: (x.burrow_score >= ATLAS_SCORE) & (x.plddt_mean >= ATLAS_PLDDT)

q = d[d.is_interface & (d.max_overlap_iso < DIMER_ONLY_PRIMARY) & passes(d) & sizeok(d)]
steps.append(("qualifying_dimer_only_druglike_pockets", len(q)))
long = pd.concat([q.assign(acc=q.uniprot_A, partner=q.uniprot_B, chain="A"),
                  q.assign(acc=q.uniprot_B, partner=q.uniprot_A, chain="B")])
long = long.drop_duplicates(["pocket_id", "acc"])  # homodimers: one row per pocket
cand = set(long.acc)
steps.append(("proteins_with_qualifying_pocket", len(cand)))

ex_atlas = set(am.loc[passes(am) & sizeok(am), "accession"])
# isolated-chain pockets: map chain -> accession via the model's accessions
macc = pd.concat([d[["modelEntityId", "uniprot_A"]].rename(columns={"uniprot_A": "acc"}).assign(chain="A"),
                  d[["modelEntityId", "uniprot_B"]].rename(columns={"uniprot_B": "acc"}).assign(chain="B")]).drop_duplicates()
iso = iso.merge(macc, on=["modelEntityId", "chain"], how="left")
ex_iso = set(iso.loc[passes(iso) & sizeok(iso), "acc"])
sc = d[passes(d) & sizeok(d) & ((d.n_res_A == 0) | (d.n_res_B == 0))]
ex_sc = set(sc.loc[sc.n_res_B == 0, "uniprot_A"]) | set(sc.loc[sc.n_res_A == 0, "uniprot_B"])
keep = cand - ex_atlas - ex_iso - ex_sc
steps += [("excluded_atlas_druglike_pocket", len(cand & ex_atlas)), ("excluded_isolated_chain_druglike_pocket", len(cand & ex_iso)),
          ("excluded_single_chain_dimer_druglike_pocket", len(cand & ex_sc)), ("listed_total", len(keep)),
          ("listed_in_atlas", len(keep & in_atlas)), ("listed_not_in_atlas", len(keep - in_atlas))]

if PRIMARY:
    (INTERIM / "q23a_candidates.txt").write_text("\n".join(sorted(keep)) + "\n")
L = long[long.acc.isin(keep)].copy()
L = L.sort_values(["acc", "burrow_score", "max_ipSAE"], ascending=[True, False, False])
agg = L.groupby("acc").agg(n_qualifying_pockets=("pocket_id", "size"), n_models=("modelEntityId", "nunique"),
                           partners=("partner", lambda x: ";".join(sorted(set(x)))),
                           kinds=("kind", lambda x: ";".join(sorted(set(x)))))
best = L.groupby("acc").head(1).set_index("acc")[["pocket_id", "modelEntityId", "kind", "partner", "volume_A3",
                                                   "burrow_score", "plddt_mean", "max_ipSAE", "max_pDockQ2",
                                                   "max_overlap_iso", "n_res_A", "n_res_B", "buriedness", "residues"]]
best.columns = ["best_" + c for c in best.columns]
out = agg.join(best)
out["gene"] = out.index.map(genes["Gene Names (primary)"]) if "Gene Names (primary)" in genes else np.nan
out["protein_name"] = out.index.map(genes["Protein names"]) if "Protein names" in genes else np.nan
out["partner_genes"] = out.partners.map(lambda s: ";".join(str(genes["Gene Names (primary)"].get(a, a)) for a in s.split(";")))
for c in ["chembl_status", "single_protein_targets", "n_act_pchembl6", "mech_max_phase", "chembl_release"]:
    out[c] = out.index.map(chem[c]) if len(chem) else pd.NA
out["in_atlas"] = out.index.isin(in_atlas)
out["n_atlas_pockets_any_size"] = out.index.map(am.accession.value_counts()).fillna(0).astype(int)
out = out.sort_values(["in_atlas", "best_burrow_score", "best_max_ipSAE"], ascending=[False, False, False])
out.insert(0, "rank", np.arange(1, len(out) + 1))
out = out.reset_index().rename(columns={"acc": "accession"})
write_tsv(out, (WS / "results/part2/q23a_dimer_only_ligandable.tsv") if PRIMARY else (OUT / f"q23a_sensitivity{SFX}.tsv"),
          f"Q2.3a; drug-like size = volume_A3 {lo}-{hi} (fallback); dimer-only OC<{DIMER_ONLY_PRIMARY}; rank by best burrow_score")
write_tsv(pd.DataFrame(steps, columns=["step", "n"]), OUT / f"q23a_n_steps{SFX}.tsv", "Q2.3a")
print(pd.DataFrame(steps).to_string())
print(out.chembl_status.value_counts())
print(out.head(15)[["rank", "accession", "gene", "partner_genes", "best_kind", "best_volume_A3", "best_burrow_score",
                    "best_max_ipSAE", "chembl_status", "in_atlas"]].to_string())
