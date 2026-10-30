#!/usr/bin/env python
"""Q2.3a curated list v2 (EXPLORATORY, post-hoc relative to the frozen design). Status: first attempt, not verified.

A-PRIORI FILTERS (declared before looking at which proteins pass; also in docs/Q22_RESULTS.md "Curated list v2"):
  F1 burrow_score >= 60th percentile of ALL atlas-filtered dimer pockets (atlas_pass rows of dimer_pocket_features)
  F2 lining pLDDT mean >= 80 on BOTH chains (plddt_lining_A and plddt_lining_B)
  F3 >= 5 lining residues on EACH chain (n_res_A, n_res_B)
  F4 helical fraction of lining residues < 0.7 (DSSP, pydssp numpy backend, run per chain on the dimer model;
     helix = H/G/I-like helix class of pydssp C3 output)
  F5 volume_A3 in 300-1500
  F6 max_ipSAE >= 0.7
  F7 dimer-only: max_overlap_iso < 0.5 (overlap coefficient vs isolated-chain pockets of both chains) AND residue overlap
     coefficient < 0.5 between each chain's lining residues (UniProt numbering) and every AFDB monomer-atlas pocket of that
     protein (if the atlas has the protein; otherwise the test is vacuous and flagged atlas_has_protein=False)
Dedup: symmetry representatives only (is_sym_rep); then pockets of the same unordered protein pair whose UniProt-position
residue sets overlap (OC >= 0.5) are merged (best burrow_score kept; n_models_merged recorded). One row per dimer pocket.
Output: results/part2/q23a_curated.tsv
"""
import sys, glob
from pathlib import Path
import importlib.util
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WS / "scripts"))
from q22_common import INTERIM, DIMER_ROOT, write_tsv, oc

# ---- declared thresholds (names, not magic numbers) ----
SCORE_PCTL = 60
MIN_PLDDT_LINING = 80.0
MIN_RES_PER_CHAIN = 5
MAX_HELIX_FRAC = 0.7
VOL_RANGE = (300.0, 1500.0)
MIN_IPSAE = 0.7
MAX_OC = 0.5

# pydssp's package __init__ needs torch (broken partial install in .venvs/q22): load only its numpy backend.
_sp = glob.glob(str(WS / ".venvs/q22/lib/python*/site-packages/pydssp/pydssp_numpy.py"))[0]
_spec = importlib.util.spec_from_file_location("pydssp_numpy", _sp)
pydssp_np = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(pydssp_np)

d = pd.read_parquet(INTERIM / "dimer_pocket_features.parquet")
steps = []
step = lambda name, m: steps.append((name, int(m.sum()) if hasattr(m, "sum") else int(m)))

thr_score = np.percentile(d.loc[d.atlas_pass, "burrow_score"], SCORE_PCTL)
print(f"burrow_score {SCORE_PCTL}th pct of {d.atlas_pass.sum()} atlas-pass dimer pockets = {thr_score:.4f}")
m = d.atlas_pass.copy(); step("atlas_pass", m)
for name, f in [("F1_score", d.burrow_score >= thr_score),
                ("F2_lining_plddt80_both", (d.plddt_lining_A >= MIN_PLDDT_LINING) & (d.plddt_lining_B >= MIN_PLDDT_LINING)),
                ("F3_ge5_res_each_chain", (d.n_res_A >= MIN_RES_PER_CHAIN) & (d.n_res_B >= MIN_RES_PER_CHAIN)),
                ("F5_volume", d.volume_A3.between(*VOL_RANGE)),
                ("F6_ipSAE", d.max_ipSAE >= MIN_IPSAE),
                ("F7a_iso_overlap", d.max_overlap_iso < MAX_OC),
                ("sym_rep", d.is_sym_rep.astype(bool))]:
    m &= f; step(name, m)
c = d[m].copy()
# isoform accessions (e.g. P00390-3): look up atlas / SIFTS / ChEMBL / genes by the canonical base accession
for s_ in "AB":
    c[f"base_{s_}"] = c[f"uniprot_{s_}"].str.split("-").str[0]
c["is_isoform_model"] = (c.base_A != c.uniprot_A) | (c.base_B != c.uniprot_B)

# ---- F7b: monomer atlas overlap ----
am = pd.read_parquet(INTERIM / "atlas_human_meta.parquet")
atlas_sets = {}
for acc, res in zip(am.accession, am.residues):
    atlas_sets.setdefault(acc, []).append({int(x.split(":")[1]) for x in res.split()})

def chain_atlas_oc(acc, pos):
    pos = {int(p) for p in pos}
    if acc not in atlas_sets:
        return np.nan, False
    if not pos:
        return np.nan, True
    return max(oc(pos, s) for s in atlas_sets[acc]), True

rows = []
for r in c.itertuples():
    oa, ha = chain_atlas_oc(r.base_A, r.uniprot_pos_A)
    ob, hb = chain_atlas_oc(r.base_B, r.uniprot_pos_B)
    rows.append((oa, ob, ha, hb))
c[["atlas_oc_A", "atlas_oc_B", "atlas_has_A", "atlas_has_B"]] = pd.DataFrame(rows, index=c.index)
c["atlas_has_protein"] = c.atlas_has_A & c.atlas_has_B
c["atlas_oc_max"] = c[["atlas_oc_A", "atlas_oc_B"]].max(axis=1)
c["atlas_testable"] = ~(c.atlas_has_A & c.atlas_oc_A.isna()) & ~(c.atlas_has_B & c.atlas_oc_B.isna())
okatlas = ~(c.atlas_oc_max >= MAX_OC)  # NaN (absent from atlas) passes; flagged
c = c[okatlas]; step("F7b_atlas_overlap", len(c))

# ---- F4: helical fraction via DSSP ----
def parse(pdb):
    ch = {}
    for ln in open(pdb):
        if ln.startswith("ATOM") and ln[12:16].strip() in ("N", "CA", "C", "O"):
            ch.setdefault(ln[21], {}).setdefault(int(ln[22:26]), {})[ln[12:16].strip()] = (
                float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return ch

helix_cache = {}
def helix_map(model, kind):
    if model in helix_cache: return helix_cache[model]
    out = {}
    for chain, res in parse(DIMER_ROOT / {"homo": "homodimers", "hetero": "heterodimers"}[kind] / f"{model}.pdb").items():
        nums = sorted(res)
        xyz = np.array([[res[n].get(a, (np.nan,) * 3) for a in ("N", "CA", "C", "O")] for n in nums])
        if len(nums) < 6 or np.isnan(xyz).any():
            out.update({(chain, n): np.nan for n in nums}); continue
        oh = pydssp_np.assign(xyz)  # [L, 3] = loop, helix, strand
        out.update({(chain, n): float(oh[i, 1]) for i, n in enumerate(nums)})
    helix_cache[model] = out
    return out

hf = []
for r in c.itertuples():
    hm = helix_map(r.modelEntityId, r.kind)
    v = [hm.get((x.split(":")[0], int(x.split(":")[1])), np.nan) for x in r.residues.split()]
    hf.append(np.nanmean(v) if not np.all(np.isnan(v)) else np.nan)
c["helix_frac_lining"] = hf
c = c[c.helix_frac_lining < MAX_HELIX_FRAC]; step("F4_helix_lt0.7", len(c))

# ---- dedup by pocket (pair-level merge of overlapping pockets) ----
c["pair"] = [";".join(sorted([a, b])) for a, b in zip(c.uniprot_A, c.uniprot_B)]
c["posset"] = [frozenset(zip("A" * len(a), a)) | frozenset(zip("B" * len(b), b)) if len(a) + len(b) else frozenset(x for x in r.split())
               for a, b, r in zip(c.uniprot_pos_A, c.uniprot_pos_B, c.residues)]
c = c.sort_values("burrow_score", ascending=False)
keep, merged = [], {}
for pair, g in c.groupby("pair", sort=False):
    kept = []
    for idx, r in g.iterrows():
        for k in kept:
            if oc(set(r.posset), set(c.loc[k, "posset"])) >= MAX_OC:
                merged[k] = merged.get(k, 1) + 1; break
        else:
            kept.append(idx); merged[idx] = 1
    keep += kept
c = c.loc[keep].copy(); c["n_models_merged"] = [merged[i] for i in c.index]
step("dedup_pocket", len(c))

# ---- annotation ----
genes = pd.read_csv(INTERIM / "uniprot_genes.tsv", sep="\t").set_index("Entry")
gn = genes["Gene Names (primary)"]
q = pd.read_csv(WS / "results/part2/q23a_dimer_only_ligandable.tsv", sep="\t", comment="#").set_index("accession")
for s in "AB":
    c[f"gene_{s}"] = c[f"base_{s}"].map(gn).fillna(c[f"gene_{'1' if s == 'A' else '2'}"])
    c[f"chembl_{s}"] = c[f"base_{s}"].map(q.chembl_status).fillna("not_in_q23a_table")
    c[f"in_prev_q23a_{s}"] = c[f"base_{s}"].isin(q.index)

sf = pd.read_parquet(WS / "data/processed/sifts_chain_uniprot.parquet", columns=["pdb", "chain", "sp_primary"])
ent = pd.read_parquet(WS / "data/part1/entries.parquet", columns=["pdb_id", "resolution", "skip_reason"])
res = ent.set_index("pdb_id").resolution
inpart1 = set(ent.pdb_id)
need = set(c.base_A) | set(c.base_B)
sf = sf[sf.sp_primary.isin(need)].drop_duplicates()
nch = sf.groupby(["sp_primary", "pdb"]).chain.nunique()
pdbs = sf.groupby("sp_primary").pdb.apply(set).to_dict()
def complex_info(a, b):
    if a == b:
        sh = {p for p in pdbs.get(a, set()) if nch.get((a, p), 0) >= 2}
    else:
        sh = pdbs.get(a, set()) & pdbs.get(b, set())
    p1 = sorted(p for p in sh if p in inpart1)
    best = min(sh, key=lambda p: res.get(p, 99)) if sh else ""
    return len(sh), len(p1), best
ci = [complex_info(a, b) for a, b in zip(c.base_A, c.base_B)]
c["n_pdb_complex"], c["n_pdb_complex_in_part1"], c["best_pdb_complex"] = zip(*ci) if ci else ([], [], [])
c["has_pdb_complex"] = c.n_pdb_complex > 0

cols = ["pocket_id", "modelEntityId", "kind", "uniprot_A", "gene_A", "uniprot_B", "gene_B", "chembl_A", "chembl_B",
        "is_isoform_model", "has_pdb_complex", "n_pdb_complex", "n_pdb_complex_in_part1", "best_pdb_complex",
        "burrow_score", "volume_A3", "max_ipSAE", "max_pDockQ2", "plddt_lining_A", "plddt_lining_B", "n_res_A", "n_res_B",
        "helix_frac_lining", "max_overlap_iso", "atlas_oc_A", "atlas_oc_B", "atlas_has_protein", "atlas_testable",
        "buriedness", "n_models_merged", "in_prev_q23a_A", "in_prev_q23a_B", "residues"]
out = c[cols].sort_values("burrow_score", ascending=False).reset_index(drop=True)
out.insert(0, "rank", np.arange(1, len(out) + 1))
out["n_unique_proteins_in_row"] = [len({a.split("-")[0], b.split("-")[0]}) for a, b in zip(out.uniprot_A, out.uniprot_B)]
write_tsv(out, WS / "results/part2/q23a_curated.tsv",
          f"Q2.3a curated v2 (post-hoc filters, see scripts/q23a_curate.py); score_thr={thr_score:.4f}")
sd = pd.DataFrame(steps, columns=["step", "n"])
sd.to_csv(WS / "results/part2/q22/q23a_curated_steps.tsv", sep="\t", index=False)
print(sd.to_string())
print("rows", len(out), "unique proteins", len({x.split("-")[0] for x in set(out.uniprot_A) | set(out.uniprot_B)}),
      "with PDB complex", int(out.has_pdb_complex.sum()))
print(out.head(10)[["rank", "pocket_id", "gene_A", "gene_B", "volume_A3", "burrow_score", "max_ipSAE", "helix_frac_lining",
                    "n_res_A", "n_res_B", "has_pdb_complex", "chembl_A", "chembl_B"]].to_string())
