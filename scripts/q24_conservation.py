#!/usr/bin/env python
"""Q2.4 (EXPLORATORY): conservation of human dimer pockets in ortholog dimers of other species.

Inputs: human annotated pockets (data/part2/human_dimer_pockets_annotated.parquet), bootstrap units from Q2.2
(data/interim/q22/dimer_pocket_features_units.parquet: paralog cluster of chain A's protein), other species' annotated pockets,
chain sequences and mmseqs2 hits (scripts/q24_prepare.py, scripts/q24_rbh.sh), divergence times (scripts/q24_divergence.py).

Definitions (declared here, first attempt, not yet verified):
- Human pockets: atlas_pass, symmetry representatives, class `interface` (>=3 lining residues per chain) or `single`
  (0 residues on one chain); marginal pockets (1-2 on one chain) excluded, as in Q2.2. Interface pockets are split into
  `dimer_only` / `preexisting` with the annotated dimer_only flag (max_overlap_iso < 0.5).
- Orthologs: reciprocal best hit by bitscore (mmseqs2, e <= 1e-3) between human and species dimer-model sequences.
- Ortholog dimer: species model whose unordered accession pair equals {RBH(A), RBH(B)} (homodimers: A = B).
- Mapping: each lining residue (chain, UniProt pos) is carried through the RBH pairwise alignment to the ortholog chain.
- Conservation score = max over ortholog models / chain orientations (both for homodimers) / target pockets T of
  |M & T| / min(|L|, |T|), where M = mapped set and L = human lining set (unaligned residues count as lost; this is the
  overlap coefficient with |L| instead of |M| so that poorly aligned pockets cannot score high on a few residues).
  Conserved = score >= 0.5. Primary target set: atlas_pass ortholog pockets; sensitivity: all ortholog pockets.
- Negative control: the same mapped set M (ortholog numbering) scored against a random species dimer of the same kind
  (homo/hetero) with no mmseqs2 hit to either human protein and chain lengths within LEN_TOL of the ortholog's (nearest total
  length if none); seed 0.
- Lining sequence conservation: fraction of aligned lining residues identical to the ortholog residue, and the same minus the
  chain-wide alignment identity (fident) of the chains the pocket lines.
- Bootstrap: 1000 resamples of human protein clusters (Q2.2 units), seed 0, percentile 95% CI.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
INTERIM = WS / "data/interim/q24"
OUT = WS / "results/part2/q24"
SPECIES = ["mouse", "rat", "zebrafish", "fly", "worm", "yeast", "pombe", "arabidopsis", "dicty", "ecoli", "mtb", "paeruginosa"]
CONSERVED_OC = 0.5           # same threshold as dimer-only / symmetry (DESIGN.md v3)
LEN_TOL = 0.2                # negative control: chain lengths within +-20% of the ortholog chains (our choice)
ID_BINS = [0, 0.3, 0.5, 0.7, 0.9, 1.0001]   # alignment identity bins (min over the two chains), our choice
ID_LABELS = ["<30", "30-50", "50-70", "70-90", ">=90"]
N_BOOT, SEED = 1000, 0       # DESIGN.md §5


def load_chains(sp):
    ch = pd.read_parquet(INTERIM / f"chains_{sp}.parquet")
    seqs = {}
    for ln in open(INTERIM / f"seqs_{sp}.fasta"):
        if ln.startswith(">"):
            k = ln[1:].strip()
        else:
            seqs[k] = ln.strip()
    ch["seq_ok"] = [seqs.get(a) == s and c for a, s, c in zip(ch.accession, ch.seq, ch.contiguous)]
    return ch, seqs


def best_hits(path):
    h = pd.read_csv(path, sep="\t", header=None, names="query,target,fident,alnlen,evalue,bits,qstart,qend,tstart,tend,qlen,tlen,qaln,taln".split(","))
    return h.sort_values(["query", "bits"], ascending=[True, False]).drop_duplicates("query"), h


def aln_map(r):
    """human seq index (0-based) -> species seq index (0-based) and identity flag."""
    i, j, m = r.qstart - 1, r.tstart - 1, {}
    for a, b in zip(r.qaln, r.taln):
        if a != "-" and b != "-":
            m[i] = (j, a == b)
        if a != "-": i += 1
        if b != "-": j += 1
    return m


def pocket_sets(df):
    out = {}
    for stem, rs, ap, itf in zip(df.structure_stem, df.residues, df.atlas_pass, df.is_interface):
        out.setdefault(stem, []).append((frozenset(rs.split()), bool(ap), bool(itf)))
    return out


def score(M, nL, targets):
    """best |M&T|/min(nL,|T|) over targets; returns (all, atlas_pass, atlas_pass & interface target)."""
    best = [0.0, 0.0, 0.0]
    for T, ap, itf in targets:
        s = len(M & T) / min(nL, len(T))
        if s > best[0]: best[0] = s
        if ap and s > best[1]: best[1] = s
        if ap and itf and s > best[2]: best[2] = s
    return best


def main(species=None):
    species = species or SPECIES
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    H = pd.read_parquet(WS / "data/part2/human_dimer_pockets_annotated.parquet")
    units = pd.read_parquet(WS / "data/interim/q22/dimer_pocket_features_units.parquet", columns=["pocket_id", "unit"])
    H = H.merge(units, on="pocket_id", how="left")
    H_models = H.drop_duplicates("structure_stem")[["structure_stem", "uniprot_A", "uniprot_B", "kind"]]  # all human models with pockets
    single = (H.n_res_A == 0) | (H.n_res_B == 0)
    H = H[H.atlas_pass & H.is_sym_rep & (H.is_interface | single)].copy()
    H["cls"] = np.where(H.is_interface, "interface", "single")
    H["sub"] = np.where(~H.is_interface, "single", np.where(H.dimer_only.astype("boolean").fillna(False), "dimer_only", "preexisting"))
    H["score_q"] = pd.qcut(H.burrow_score, 5, labels=[f"Q{i}" for i in range(1, 6)])
    H["unit"] = H.unit.fillna(H.uniprot_A)
    hch, hseq = load_chains("human")
    h_ok = {(m, c) for m, c, ok in zip(hch.modelEntityId, hch.chain, hch.seq_ok) if ok}
    h_first = {(m, c): f for m, c, f in zip(hch.modelEntityId, hch.chain, hch.first_num)}
    div = pd.read_csv(OUT / "q24_divergence.tsv", sep="\t", comment="#").set_index("species").divergence_mya

    rows, pair_rows, lin_rows = [], [], []
    for sp in species:
        if not (INTERIM / f"hits_{sp}_fwd.tsv").exists() or not (INTERIM / f"annot_{sp}.parquet").exists():
            print("missing inputs for", sp, file=sys.stderr); continue
        fwd, fwd_all = best_hits(INTERIM / f"hits_{sp}_fwd.tsv")
        rev, _ = best_hits(INTERIM / f"hits_{sp}_rev.tsv")
        rb = dict(zip(rev["query"], rev["target"]))
        fwd = fwd[[rb.get(t) == q for q, t in zip(fwd["query"], fwd["target"])]]
        rbh = {r.query: r for r in fwd.itertuples()}
        amap = {q: aln_map(r) for q, r in rbh.items()}
        homologs = fwd_all.groupby("query").target.agg(set).to_dict()
        sch, sseq = load_chains(sp)
        s_ok = {(m, c) for m, c, ok in zip(sch.modelEntityId, sch.chain, sch.seq_ok) if ok}
        s_first = {(m, c): f for m, c, f in zip(sch.modelEntityId, sch.chain, sch.first_num)}
        mod = sch.pivot_table(index=["modelEntityId", "kind"], columns="chain", values="accession", aggfunc="first").reset_index()
        mlen = sch.assign(L=sch.seq.str.len()).pivot_table(index="modelEntityId", columns="chain", values="L", aggfunc="first")
        mod = mod.join(mlen.add_suffix("_len"), on="modelEntityId")
        by_pair = {}
        for r in mod.itertuples():
            by_pair.setdefault(frozenset((r.A, r.B)), []).append(r)
        S = pd.read_parquet(INTERIM / f"annot_{sp}.parquet", columns=["structure_stem", "residues", "atlas_pass", "is_interface"])
        tgt = pocket_sets(S)
        # human models -> ortholog models
        hm = H_models
        n_pairs = {"homo": 0, "hetero": 0}; orth_of = {}
        for r in hm.itertuples():
            a, b = rbh.get(r.uniprot_A), rbh.get(r.uniprot_B)
            if a is None or b is None:
                continue
            ms = [m for m in by_pair.get(frozenset((a.target, b.target)), []) if (m.modelEntityId, "A") in s_ok and (m.modelEntityId, "B") in s_ok]
            if ms and (r.structure_stem, "A") in h_ok and (r.structure_stem, "B") in h_ok:
                orth_of[r.structure_stem] = ms; n_pairs[r.kind] += 1
        # negative-control partner per ortholog model
        hm_idx = hm.set_index("structure_stem")
        pool = {k: g for k, g in mod.groupby("kind")}
        null_of = {}
        for stem, ms in orth_of.items():
            hr = hm_idx.loc[stem]
            excl = homologs.get(hr.uniprot_A, set()) | homologs.get(hr.uniprot_B, set())
            lst = []
            for m in ms:
                g = pool[m.kind]
                g = g[~g.A.isin(excl) & ~g.B.isin(excl) & (g.modelEntityId != m.modelEntityId)]
                ok = g[(abs(g.A_len - m.A_len) <= LEN_TOL * m.A_len) & (abs(g.B_len - m.B_len) <= LEN_TOL * m.B_len)]
                if len(ok) == 0:
                    ok = g.iloc[[np.argmin(abs(g.A_len + g.B_len - m.A_len - m.B_len).values)]]
                lst.append(ok.iloc[rng.integers(len(ok))])
            null_of[stem] = lst
        pair_rows.append((sp, div.get(sp), len(rbh), n_pairs["homo"], n_pairs["hetero"],
                          len({frozenset((r.uniprot_A, r.uniprot_B)) for r in hm.itertuples() if r.structure_stem in orth_of})))
        # score pockets
        for p in H[H.structure_stem.isin(orth_of)].itertuples():
            L = p.residues.split(); nL = len(L)
            acc = {"A": p.uniprot_A, "B": p.uniprot_B}
            mapped = []  # (orig chain, species pos) or None
            ident, n_al = 0, 0
            for x in L:
                ch, pos = x[0], int(x[2:])
                hit = amap[acc[ch]].get(pos - h_first[(p.structure_stem, ch)])
                mapped.append(None if hit is None else (ch, hit[0]))
                if hit is not None:
                    n_al += 1; ident += hit[1]
            fid = [rbh[acc[c]].fident for c in ("A", "B") if any(x[0] == c for x in L)]
            idmin = min(rbh[p.uniprot_A].fident, rbh[p.uniprot_B].fident)

            def best(models):
                b = [0.0, 0.0, 0.0]
                for m in models:
                    orients = [{"A": "A", "B": "B"}, {"A": "B", "B": "A"}] if p.kind == "homo" else \
                              [{"A": "A", "B": "B"} if m.A == rbh[p.uniprot_A].target else {"A": "B", "B": "A"}]
                    if p.kind == "hetero" and m.A == rbh[p.uniprot_A].target and m.B != rbh[p.uniprot_B].target:
                        continue
                    for o in orients:
                        M = frozenset(f"{o[c]}:{j + s_first.get((m.modelEntityId, o[c]), 1)}" for c, j in (t for t in mapped if t))
                        s = score(M, nL, tgt.get(m.modelEntityId, []))
                        b = [max(u, v) for u, v in zip(b, s)]
                return b

            def best_null(models, orth_models):
                b = [0.0, 0.0, 0.0]
                for m, om in zip(models, orth_models):
                    # carry ortholog chain labels onto the random model by chain order; both orientations for homo
                    orients = [{"A": "A", "B": "B"}, {"A": "B", "B": "A"}] if p.kind == "homo" else \
                              [{"A": "A", "B": "B"} if om.A == rbh[p.uniprot_A].target else {"A": "B", "B": "A"}]
                    for o in orients:
                        M = frozenset(f"{o[c]}:{j + s_first.get((om.modelEntityId, o[c]), 1)}" for c, j in (t for t in mapped if t))
                        s = score(M, nL, tgt.get(m.modelEntityId, []))
                        b = [max(u, v) for u, v in zip(b, s)]
                return b

            real = best(orth_of[p.structure_stem])
            null = best_null(null_of[p.structure_stem], orth_of[p.structure_stem])
            rows.append((sp, p.pocket_id, p.structure_stem, p.kind, p.cls, p.sub, p.score_q, p.burrow_score, p.unit, nL, n_al / nL,
                         ident / n_al if n_al else np.nan, float(np.mean(fid)), idmin, *real, *null))
        print(sp, "RBH", len(rbh), "pairs", n_pairs, "pockets scored", sum(1 for r in rows if r[0] == sp), flush=True)

    R = pd.DataFrame(rows, columns=["species", "pocket_id", "structure_stem", "kind", "cls", "sub", "score_q", "burrow_score", "unit", "n_lining",
                                    "frac_aligned", "lining_identity", "chain_identity", "pair_identity_min",
                                    "oc_all", "oc_atlas", "oc_atlas_itf", "null_oc_all", "null_oc_atlas", "null_oc_atlas_itf"])
    R["id_bin"] = pd.cut(R.pair_identity_min, ID_BINS, labels=ID_LABELS, right=False)
    R["divergence_mya"] = R.species.map(div)
    INTERIM.mkdir(parents=True, exist_ok=True)
    R.to_parquet(INTERIM / "q24_pocket_scores.parquet", index=False)
    P = pd.DataFrame(pair_rows, columns=["species", "divergence_mya", "n_rbh_human_proteins", "n_human_homo_models_with_ortholog_dimer",
                                         "n_human_hetero_models_with_ortholog_dimer", "n_human_protein_pairs_with_ortholog_dimer"])
    P.to_csv(OUT / "q24_ortholog_pairs.tsv", sep="\t", index=False)
    print(P.to_string())


if __name__ == "__main__":
    main()
