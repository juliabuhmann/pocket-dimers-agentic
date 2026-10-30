#!/usr/bin/env python3
"""EXPLORATORY post-hoc extension of deep dive #2 (docs/COFOLD_H4.md): more heterodimer multi-chain (MBS) ligand sites.

Motivation: the original hetero stratum (n = 6) showed 4/6 dimer vs 2/6 single-chain successes; this script selects
additional heterodimer MBS cases with the same input protocol, to test that small-n hint. Post-hoc, not pre-registered.

Status: first attempt, UNVERIFIED.

Rerun:
    source scripts/env_threads.sh
    .venvs/part1/bin/python scripts/cofold_select_hetero_ext.py   # -> results/cofold/hetero_ext/cases.tsv,
                                                                   #    data/cofold/hetero_ext/inputs/*.json

Reuses (imports, unmodified) scripts/cofold_select_h4.py: uniprot_seqs, cached_msa, entry_info (chain pair, contacts,
ligand copies contacting both chains), write_inputs (AF3 JSONs, MSA remap, no templates, seed 0). Only the module's
output directory (WORK) is redirected. Identical per-case protocol: run A = PDB pair + 1 or 2 ligand copies (2 if
>= 2 copies contact both chains in assembly 1; reduced to 1 above SAFE_TOKENS); run B = chain with most contacts + 1.

Selection rules (all here):
- Sites: H4 v2 primary set, lig_set organic, is_mbs, kind == hetero; census: one CCD component, not covalent
  (as in the original).
- Ligand classes: TIER1 = the original classes (drug_like, cofactor). Because TIER1 alone has too few unused
  clusters for the target, TIER2 = all other ligand classes with >= MIN_HEAVY_TIER2 heavy atoms (excludes small
  additives/fragments such as EDO, MPD, PYR, BEN). The tier is recorded per case; results are reported both ways.
- Unit = unordered pair of RCSB 30% sequence clusters (clusters-by-entity-30, entities of the two H4 chains).
  One case per cluster pair; cluster pairs of the 6 original hetero cases are excluded. One case per PDB entry
  (case ids are <pdb>_<ccd>; an entry can hold two different pairs, e.g. 8umv).
- Within a pair: TIER1 first, post-cutoff release first, fewer estimated tokens, better resolution (up to
  MAX_TRIES_PER_PAIR tried). Across pairs: TIER1 first, then fewer tokens (wall-clock budget).
- Size: dimer tokens <= MAX_TOKENS (task limit, = largest AF3 bucket without unified memory, docs/COFOLDING.md).
- At most N_TARGET cases.
- Run order (file-name priority): ascending dimer tokens, dimer before monomer of the same case, so that an early
  stop leaves complete pairs.
"""
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import cofold_select_h4 as h4  # noqa: E402  (original selection helpers, reused unmodified)

H4_CASES = ROOT / "results/cofold/h4/cases.tsv"
CLUSTERS30 = ROOT / "data/raw/part1/clusters-by-entity-30.txt"   # fetched by scripts/fetch_part1_sources.py
TIER1 = h4.LIG_CLASSES                 # ("drug_like", "cofactor"), original deep dive
MIN_HEAVY_TIER2 = 10                   # tier-2 ligands below this are crystallisation additives / fragments
MAX_TOKENS = 1280                      # task limit
N_TARGET = 20                          # task target
MAX_TRIES_PER_PAIR = 5
OUT_DIR = ROOT / "results/cofold/hetero_ext"
WORK = ROOT / "data/cofold/hetero_ext"
h4.WORK = WORK                         # redirect the reused helpers' output dir (uniprot cache, inputs)


def log(*a):
    print(*a, flush=True)


def load_clusters():
    out = {}
    with open(CLUSTERS30) as f:
        for i, line in enumerate(f):
            for m in line.split():
                out[m] = i
    return out


_ENT = {}


def chain_entities(pdb_id):
    if pdb_id not in _ENT:
        st = gemmi.read_structure(str(h4.PDB_MIRROR / pdb_id[1:3] / f"{pdb_id}.cif.gz"))
        st.setup_entities()
        m = {}
        for ch in st[0]:
            pol = ch.get_polymer()
            if len(pol):
                m[ch.name] = pol[0].entity_id
        _ENT[pdb_id] = m
    return _ENT[pdb_id]


def cluster_pair(r, cl):
    ent = chain_entities(r.pdb_id)
    ks = []
    for c in (r.chain1_auth, r.chain2_auth):
        e = ent.get(c)
        ks.append(cl.get(f"{r.pdb_id.upper()}_{e}", -1) if e is not None else -1)
    return tuple(sorted(ks))


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    b = pd.read_parquet(h4.H4_TABLE)
    b = b[(b.lig_set == "organic") & b.is_mbs & (b.kind == "hetero")].copy()
    log("H4 primary hetero MBS organic sites", len(b))
    cen = pd.read_parquet(h4.CENSUS, columns=["ligand_uid", "n_components", "is_covalent"],
                          filters=[("ligand_uid", "in", sorted(b.ligand_uid.unique()))])
    b = b.merge(cen, on="ligand_uid", how="left")
    b = b[(b.n_components == 1) & ~b.is_covalent.fillna(False).astype(bool)]
    b["tier"] = np.where(b.ligand_class.isin(TIER1), 1, 2)
    b = b[(b.tier == 1) | (b.n_heavy >= MIN_HEAVY_TIER2)]
    log("after component/covalent/tier filters", len(b), b.tier.value_counts().to_dict())
    fm = pd.read_csv(h4.FAILURE_TABLE, sep="\t", usecols=["ligand_uid", "model", "fm_dimer", "fm_iso"])
    b = b.merge(fm, on=["ligand_uid", "model"], how="left")
    b["post_cutoff"] = b.release_date > h4.AF3_CUTOFF
    useq = h4.uniprot_seqs(set(b.acc1) | set(b.acc2))
    b["est_tokens"] = b.acc1.map(lambda a: len(useq.get(a, ""))) + b.acc2.map(lambda a: len(useq.get(a, ""))) + b.n_heavy
    cl = load_clusters()
    b["cl_pair"] = [cluster_pair(r, cl) for r in b.itertuples()]
    b = b[[-1 not in p for p in b.cl_pair]]

    orig = pd.read_csv(H4_CASES, sep="\t")
    orig = orig[orig.stratum == "hetero"]
    orig_pairs = set()
    for r in orig.itertuples():
        m = b[b.ligand_uid == r.ligand_uid]
        if len(m):
            orig_pairs.add(m.cl_pair.iloc[0])
        else:   # original site filtered out above (cannot happen for tier 1); compute from the H4 table
            row = pd.read_parquet(h4.H4_TABLE, filters=[("ligand_uid", "==", r.ligand_uid)]).iloc[0]
            orig_pairs.add(cluster_pair(row, cl))
    log("original hetero cluster pairs (excluded):", sorted(orig_pairs))
    b = b[~b.cl_pair.isin(orig_pairs)]
    log("eligible sites", len(b), "cluster pairs", b.cl_pair.nunique(),
        "tier1 pairs", b[b.tier == 1].cl_pair.nunique())

    b = b.sort_values(["tier", "post_cutoff", "est_tokens", "resolution", "pdb_id"],
                      ascending=[True, False, True, True, True])
    reps = b.groupby("cl_pair", sort=False).head(1).sort_values(["tier", "est_tokens"])
    msas, rows, rejects, used_pdb = {}, [], [], set()

    def get_msa(acc):
        if acc not in msas:
            p = h4.cached_msa(acc)
            msas[acc] = p if (p is not None and p["sequence"] == useq.get(acc)) else None
        return msas[acc]

    for pair in reps.cl_pair:
        if len(rows) >= N_TARGET:
            break
        for _, r in b[(b.cl_pair == pair) & ~b.pdb_id.isin(used_pdb)].head(MAX_TRIES_PER_PAIR).iterrows():
            if r.est_tokens > 2 * MAX_TOKENS:
                rejects.append((r.ligand_uid, str(pair), f"est_tokens={r.est_tokens}")); continue
            if get_msa(r.acc1) is None or get_msa(r.acc2) is None:
                rejects.append((r.ligand_uid, str(pair), "no_cached_msa_or_seq_ne_uniprot")); continue
            try:
                info = h4.entry_info(r, useq)
            except Exception as e:  # noqa: BLE001
                rejects.append((r.ligand_uid, str(pair), f"error:{e}")); continue
            if isinstance(info, str):
                rejects.append((r.ligand_uid, str(pair), info)); continue
            ch = info["chains"]
            nlig = 2 if info["n_copies_both"] >= 2 else 1
            tok = len(ch[0]["seq"]) + len(ch[1]["seq"]) + nlig * int(r.n_heavy)
            reduced = False
            if nlig == 2 and tok > h4.SAFE_TOKENS:
                nlig, tok, reduced = 1, tok - int(r.n_heavy), True
            if tok > MAX_TOKENS:
                rejects.append((r.ligand_uid, str(pair), f"tokens={tok}>{MAX_TOKENS}")); continue
            rows.append(dict(case_id=f"hext_{r.pdb_id}_{r.ccd}", priority=0, stratum="hetero_ext", tier=int(r.tier),
                             pdb_id=r.pdb_id, ligand_uid=r.ligand_uid, ccd=r.ccd, ligand_class=r.ligand_class,
                             n_heavy=int(r.n_heavy), kind=r.kind, acc1=r.acc1, acc2=r.acc2, model=r.model,
                             cluster30=int(r.cluster30), cl30_pair=f"{pair[0]}|{pair[1]}", resolution=r.resolution,
                             release_date=r.release_date, post_cutoff=bool(r.post_cutoff),
                             post_cutoff_pair=bool(r.post_cutoff_pair),
                             pair_complex_pre_cutoff=not bool(r.post_cutoff_pair),
                             iface_jaccard=r.iface_jaccard, apo_dimer_dcc=r.dimer_dcc, apo_dimer_hit=bool(r.dimer_hit),
                             apo_iso_dcc=r.iso_dcc, apo_iso_hit=bool(r.iso_hit), fm_dimer=r.fm_dimer, fm_iso=r.fm_iso,
                             centroid=h4.json.dumps([round(float(x), 3) for x in r.centroid]),
                             chainA_pdb=ch[0]["chain"], chainB_pdb=ch[1]["chain"], accA=ch[0]["acc"], accB=ch[1]["acc"],
                             loA=ch[0]["lo"], loB=ch[1]["lo"], lenA=len(ch[0]["seq"]), lenB=len(ch[1]["seq"]),
                             identA=round(ch[0]["identity"], 3), identB=round(ch[1]["identity"], 3),
                             n_contacts_A=ch[0]["n_contacts"], n_contacts_B=ch[1]["n_contacts"],
                             n_copies_both_pdb=info["n_copies_both"], n_lig_dimer_run=nlig,
                             lig_copies_reduced=reduced, tokens_dimer=tok,
                             tokens_monomer=len(ch[0]["seq"]) + int(r.n_heavy), _chains=ch))
            used_pdb.add(r.pdb_id)
            log(rows[-1]["case_id"], "tier", r.tier, "tokens", tok, r.release_date, pair)
            break

    order = sorted(rows, key=lambda c: (c["tokens_dimer"], c["case_id"]))
    for i, c in enumerate(order):
        c["priority"] = i
        c.update({f"input_{k}": v for k, v in h4.write_inputs(c, c["_chains"], msas, c["n_lig_dimer_run"]).items()})
    pd.DataFrame(rejects, columns=["ligand_uid", "cl30_pair", "reason"]).to_csv(OUT_DIR / "rejects.tsv", sep="\t",
                                                                               index=False)
    df = pd.DataFrame(order).drop(columns=["_chains"])
    df.to_csv(OUT_DIR / "cases.tsv", sep="\t", index=False)
    log(df.groupby(["tier", "post_cutoff", "pair_complex_pre_cutoff"]).size())
    log("n cases", len(df), "total dimer tokens", df.tokens_dimer.sum())


if __name__ == "__main__":
    main()
