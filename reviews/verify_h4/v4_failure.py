"""Verifier: failure-mode classification of misses, from the verifier's own recomputation (dcc_recompute_all.tsv).
Rules are the verifier's (named constants below), applied in order; first-attempt heuristics, UNVERIFIED.
Also writes the PI spot-check table (10 MBS dimer hits + 10 MBS dimer misses, seed-0 sample from v2_dcc.py)."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
OUT = Path(__file__).resolve().parent
DCC_CUT = 4.0
NEAR = 5.0            # near-miss band (4, 5] Å
GEOM_RMSD = 2.0       # local CA RMSD above which model geometry is considered locally wrong
REC_FOUND = 0.5       # best pocket covers >= half of the PDB binding residues -> pocket "found"
OVERSIZE = 2.0        # pocket size >= OVERSIZE x n binding residues -> oversized/merged
IFACE_PARTIAL = 0.5   # interface Jaccard below this flagged as partial interface match

def classify(r, lab):
    if r[f"v_local_{lab}_hit"]: return "hit"
    pk = json.loads(r[f"v_local_{lab}_pockets"]) if isinstance(r[f"v_local_{lab}_pockets"], str) else []
    nb = r[f"v_local_{lab}_nbind"]
    if r["v_local_rmsd"] > GEOM_RMSD: return "model geometry (local CA RMSD > 2 A)"
    if r[f"v_local_{lab}_dcc_any"] <= DCC_CUT: return "pocket exists but ranked outside top-(n+2)"
    if not pk: return "no pockets"
    best = max(pk, key=lambda p: (p[4], -p[5]))   # [rank, score, size, chains, recall, dcc]
    rec = sorted([p[4] for p in pk], reverse=True)
    if r[f"v_local_{lab}_dcc_topk"] <= NEAR: return "near miss (DCC 4-5 A)"
    if best[4] >= REC_FOUND:
        if best[2] >= OVERSIZE * nb: return "pocket found but oversized/merged (centre displaced)"
        return "pocket found, centre offset > 5 A"
    if len(rec) > 1 and rec[0] + rec[1] >= REC_FOUND: return "site split across pockets"
    return "no pocket at site (< 50% of binding residues in any pocket)"

def main():
    R = pd.read_csv(OUT / "dcc_recompute_all.tsv", sep="\t")
    R = R[R.v_status == "ok"].copy()
    for lab in ("dimer", "iso"):
        R[f"fm_{lab}"] = [classify(r, lab) for _, r in R.iterrows()]
    R["partial_iface"] = R.iface_jaccard < IFACE_PARTIAL
    M = R[R.is_mbs]
    tab = []
    for lab in ("dimer", "iso"):
        t = M[f"fm_{lab}"].value_counts(); tab.append(t.rename(f"MBS {lab}"))
        tab.append(R[~R.is_mbs][f"fm_{lab}"].value_counts().rename(f"SBS {lab}"))
    # discordant MBS: dimer miss & iso hit vs dimer hit & iso miss
    tab.append(M[~M.v_local_dimer_hit & M.v_local_iso_hit].fm_dimer.value_counts().rename("MBS dimer-miss/iso-hit (dimer reason)"))
    tab.append(M[M.v_local_dimer_hit & ~M.v_local_iso_hit].fm_iso.value_counts().rename("MBS dimer-hit/iso-miss (iso reason)"))
    T = pd.concat(tab, axis=1).fillna(0).astype(int)
    T.to_csv(OUT / "failure_modes.tsv", sep="\t")
    print(T.to_string())
    print("partial iface among MBS dimer misses:", int((M.partial_iface & ~M.v_local_dimer_hit).sum()), "/", int((~M.v_local_dimer_hit).sum()),
          "; hit rate partial vs full iface:", M.groupby("partial_iface").v_local_dimer_hit.mean().round(3).to_dict())
    # per-family top contributors to dimer-miss/iso-hit at MBS
    x = M[~M.v_local_dimer_hit & M.v_local_iso_hit]
    print(x.groupby(["model_acc_A", "fm_dimer"]).size().sort_values(ascending=False).head(12).to_string())
    # agreement with the producer over all sites
    for lab in ("dimer", "iso"):
        a = (R[f"{lab}_hit"] == R[f"v_local_{lab}_hit"]); g = (R[f"{lab}_hit"] == R[f"v_global_{lab}_hit"])
        dd = (R[f"{lab}_dcc"] - R[f"v_local_{lab}_dcc_topk"]).abs()
        print(lab, "agree local", int(a.sum()), "/", len(R), "global", int(g.sum()), "| |dDCC| local max", round(dd.max(), 3), "p99", round(dd.quantile(.99), 3))
    print("centroid shift > 1 A:", int((R.v_centroid_shift > 1).sum()), "seqid<0.9:", int((R.v_local_seqid < 0.9).sum()))
    # own headline from own recomputation
    for lab, c in [("MBS", True), ("SBS", False)]:
        s = R[R.is_mbs == c]; print(lab, "own recompute dimer", round(s.v_local_dimer_hit.mean(), 4), "iso", round(s.v_local_iso_hit.mean(), 4),
                                   "| global-frame dimer", round(s.v_global_dimer_hit.mean(), 4), "iso", round(s.v_global_iso_hit.mean(), 4))
    gm = lambda s, f: (s[f"v_{f}_dimer_hit"].astype(float) - s[f"v_{f}_iso_hit"].astype(float)).mean()
    for f in ("local", "global"):
        print(f, "DiD:", round(gm(R[R.is_mbs], f) - gm(R[~R.is_mbs], f), 4))
    R[["ligand_uid", "pdb_id", "ccd", "is_mbs", "kind", "acc1", "acc2", "model", "fm_dimer", "fm_iso", "v_local_dimer_dcc_topk", "v_local_iso_dcc_topk"]].to_csv(OUT / "failure_modes_sites.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
