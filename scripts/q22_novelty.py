#!/usr/bin/env python
"""Q2.2 analysis 4 (EXPLORATORY): novelty of dimer pockets vs the human monomer atlas (pocketcompare fingerprints).

For each dimer pocket: max cosine (inner product of 384-d unit fingerprints) to human atlas pockets (taxid 9606,
tiers 1+2), EXCLUDING atlas pockets of both proteins of the model and of their paralogs (mmseqs2, >=30% identity,
e-value <= 1e-3; scripts/q22_paralogs.sh). Also per-query z-score = (max_excl - mean_excl) / sd_excl over the
non-excluded atlas row, and max cosine to the same protein's own atlas pockets (reference only).
Comparison interface vs single-chain (atlas filter, sym reps, marginal pockets excluded) and dimer-only vs
pre-existing interface pockets: raw, z-score, and size-matched (within volume_A3 deciles of the compared pockets;
weighted mean of within-decile median differences, cluster bootstrap over paralog-cluster units).
Matrix products: torch on cuda:1 (never GPU 0); CPU fallback.
CAVEAT: pocketcompare recognises PPI-interface pockets poorly (docs/TOOLS.md §3) -> low max cosine for interface
pockets may reflect method weakness rather than biological novelty.
Positive control (PDB interface ligand sites with fingerprints): not available (Part 1 has no pocketcompare
fingerprints of PDB sites yet) -> skipped, logged.
Outputs: data/interim/q22/novelty.parquet, results/part2/q22/q22_novelty.tsv
Status: first attempt, not yet verified.
"""
import glob, json, sys
import numpy as np, pandas as pd, pyarrow.parquet as pq
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from q22_common import (INTERIM, OUT, POCKETS_DIMER, PARALOG_MIN_ID, PARALOG_MAX_EVALUE, SEEDS, N_BOOT,
                        DIMER_ONLY_PRIMARY, write_tsv)

DEVICE = "cuda:1"
N_DECILES = 10
CHUNK = 4096


def main():
    import torch
    dev = torch.device(DEVICE if torch.cuda.is_available() else "cpu")
    print("device", dev, flush=True)
    from q22_compare import paralog_clusters
    from q22_common import ipsae_tertile
    d = pd.read_parquet(INTERIM / "dimer_pocket_features.parquet")
    d["unit"] = d.uniprot_A.map(paralog_clusters(sorted(set(d.uniprot_A) | set(d.uniprot_B))))
    _, edges = ipsae_tertile(d.drop_duplicates("modelEntityId").max_ipSAE)
    d["tert"], _ = ipsae_tertile(d.max_ipSAE, edges)
    # fingerprints of dimer pockets
    fps, ids, methods = [], [], set()
    for kind in ("homo", "hetero"):
        for f in sorted(glob.glob(str(POCKETS_DIMER / kind / "chunk_[0-9][0-9][0-9][0-9].parquet"))):
            md = pq.read_schema(f).metadata
            methods.add(json.loads(md[b"pocketcompare"])["method"])
            t = pq.read_table(f, columns=["pocket_id", "fingerprint"])
            ids += t.column("pocket_id").to_pylist()
            fps.append(np.asarray(t.column("fingerprint").combine_chunks().flatten(), np.float32).reshape(-1, 384))
    Q = np.concatenate(fps); qid = pd.Series(np.arange(len(ids)), index=ids)
    assert qid.index.is_unique
    d = d[d.pocket_id.isin(qid.index)].reset_index(drop=True)
    Qd = Q[qid[d.pocket_id].to_numpy()]
    am = pd.read_parquet(INTERIM / "atlas_human_meta.parquet")
    Af = np.load(INTERIM / "atlas_human_fp16.npy").astype(np.float32)
    print("dimer methods", methods, "dimer pockets", len(d), "atlas", Af.shape, flush=True)
    # exclusion sets
    acc_codes = {a: i for i, a in enumerate(sorted(am.accession.unique()))}
    a_idx = am.accession.map(acc_codes).to_numpy()
    cols_of = pd.Series(np.arange(len(am))).groupby(a_idx).apply(np.array).to_dict()
    par = {}
    pf = INTERIM / "paralogs_mmseqs.tsv"
    paralog_status = "mmseqs2"
    if pf.exists():
        h = pd.read_csv(pf, sep="\t", header=None, names=["q", "t", "fident", "alnlen", "evalue", "qcov", "tcov"])
        h = h[(h.q != h.t) & (h.fident >= PARALOG_MIN_ID) & (h.evalue <= PARALOG_MAX_EVALUE)]
        for q, t in zip(h.q, h.t):
            par.setdefault(q, set()).add(t); par.setdefault(t, set()).add(q)
    else:
        paralog_status = "FALLBACK_self_only"
    print("paralog status", paralog_status, "proteins with paralogs", len(par), flush=True)

    def excl_cols(accs, with_par=True):
        s = set(accs)
        if with_par:
            for a in accs: s |= par.get(a, set())
        cols = [cols_of[acc_codes[a]] for a in s if a in acc_codes and acc_codes[a] in cols_of]
        return np.concatenate(cols) if cols else np.zeros(0, int)

    A = torch.from_numpy(Af).to(dev)
    n = len(d)
    out = {k: np.full(n, np.nan) for k in ["max_excl", "mean_excl", "sd_excl", "max_self", "max_all", "n_excluded"]}
    best = np.full(n, -1)
    model_rows = d.groupby("modelEntityId").indices
    models = list(model_rows)
    # process model batches up to CHUNK rows
    batch, rows_in = [], 0
    def run(batch):
        rows = np.concatenate([model_rows[m] for m in batch])
        S = torch.from_numpy(Qd[rows]).to(dev) @ A.T
        out["max_all"][rows] = S.max(1).values.cpu().numpy()
        mask = torch.zeros_like(S, dtype=torch.bool)
        selfm = torch.zeros_like(S, dtype=torch.bool)
        off = 0
        for m in batch:
            r = model_rows[m]; k = len(r)
            accs = {d.uniprot_A.iat[r[0]], d.uniprot_B.iat[r[0]]}
            ce = torch.from_numpy(excl_cols(accs)).to(dev); cs = torch.from_numpy(excl_cols(accs, False)).to(dev)
            if len(ce): mask[off:off + k, ce] = True
            if len(cs): selfm[off:off + k, cs] = True
            off += k
        Sm = S.masked_fill(mask, float("nan"))
        mx = torch.nan_to_num(Sm, nan=-9.0).max(1)
        out["max_excl"][rows] = mx.values.cpu().numpy(); best[rows] = mx.indices.cpu().numpy()
        out["mean_excl"][rows] = torch.nanmean(Sm, 1).cpu().numpy()
        cnt = (~mask).sum(1).float()
        mu = torch.nanmean(Sm, 1, keepdim=True)
        out["sd_excl"][rows] = torch.sqrt(torch.nansum((Sm - mu) ** 2, 1) / (cnt - 1)).cpu().numpy()
        out["max_self"][rows] = torch.where(selfm.any(1), S.masked_fill(~selfm, -9.0).max(1).values,
                                            torch.full_like(cnt, float("nan"))).cpu().numpy()
        out["n_excluded"][rows] = mask.sum(1).cpu().numpy()
    for m in models:
        batch.append(m); rows_in += len(model_rows[m])
        if rows_in >= CHUNK:
            run(batch); batch, rows_in = [], 0
    if batch: run(batch)
    for k, v in out.items(): d[k] = v
    d["z_excl"] = (d.max_excl - d.mean_excl) / d.sd_excl
    d["best_atlas_pocket"] = am.pocket_id.to_numpy()[best]; d["best_atlas_acc"] = am.accession.to_numpy()[best]
    d["paralog_status"] = paralog_status
    keep = ["pocket_id", "modelEntityId", "kind", "uniprot_A", "uniprot_B", "unit", "tert", "is_interface", "is_interface_ge1",
            "atlas_pass", "is_sym_rep", "max_overlap_iso", "volume_A3", "burrow_score"] + list(out) + ["z_excl", "best_atlas_pocket", "best_atlas_acc", "paralog_status"]
    d[keep].to_parquet(INTERIM / "novelty.parquet")
    print("done similarity", flush=True)
    stats(d[keep])


def size_matched(df, metric, seed, n_boot=N_BOOT):
    """weighted (by g1 count) mean of within-decile median differences g1-g0; deciles of volume_A3 over df."""
    df = df.copy()
    df["dec"] = pd.qcut(df.volume_A3.rank(method="first"), N_DECILES, labels=False)
    def est(x):
        diffs, w = [], []
        for _, g in x.groupby("dec"):
            a, b = g.loc[g.grp, metric].dropna(), g.loc[~g.grp, metric].dropna()
            if len(a) and len(b):
                diffs.append(a.median() - b.median()); w.append(len(a))
        return np.average(diffs, weights=w) if w else np.nan
    e = est(df)
    rng = np.random.default_rng(seed)
    gi = df.groupby("unit").indices; keys = list(gi)
    bs = [est(df.iloc[np.concatenate([gi[keys[i]] for i in rng.integers(0, len(keys), len(keys))])]) for _ in range(n_boot)]
    return e, np.nanpercentile(bs, 2.5), np.nanpercentile(bs, 97.5)


def stats(d):
    from q22_compare import boot_compare
    S = d[d.atlas_pass & d.is_sym_rep]
    marginal = S.is_interface_ge1 & ~S.is_interface
    rows = []
    contrasts = [("interface_vs_single_chain", S[~marginal].assign(grp=lambda x: x.is_interface)),
                 ("dimer_only_vs_preexisting_interface",
                  S[S.is_interface & S.max_overlap_iso.notna()].assign(grp=lambda x: x.max_overlap_iso < DIMER_ONLY_PRIMARY))]
    for cname, C in contrasts:
        C = C.reset_index(drop=True)
        for label, sub in [("all", C), ("homo", C[C.kind == "homo"]), ("hetero", C[C.kind == "hetero"])]:
            sub = sub.reset_index(drop=True)
            for r in boot_compare((sub, label, 0, ["max_excl", "z_excl", "max_self"])):
                r["contrast"] = cname; r["adjust"] = "none"; rows.append(r)
            for metric in ["max_excl", "z_excl"]:
                for s in (SEEDS if label == "all" else SEEDS[:1]):
                    e, lo, hi = size_matched(sub, metric, s, n_boot=N_BOOT if s == 0 else 200)
                    rows.append(dict(contrast=cname, stratum=label, seed=s, metric=metric, adjust="volume_decile_matched",
                                     n_g1=int(sub.grp.sum()), n_g0=int((~sub.grp).sum()), median_diff=e, md_lo=lo, md_hi=hi))
            print(cname, label, flush=True)
    write_tsv(pd.DataFrame(rows), OUT / "q22_novelty.tsv",
              "analysis4; g1 = interface / dimer-only; max_excl excludes self + >=30% id paralogs; PPI-weakness caveat; seeds 1-2 size-matched use 200 resamples")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--stats-only":
        stats(pd.read_parquet(INTERIM / "novelty.parquet"))
    else:
        main()
