#!/usr/bin/env python
"""Annotate dimer pockets (output of run_dimer_pockets.py) with chain / interface / symmetry / dimer-only info.

Per pocket: n_res_A/B, is_interface (PRIMARY: >=3 lining residues on each chain), is_interface_ge1 (sensitivity: >=1 each),
UniProt accession per chain (candidates.csv: A = uniprotAccession_1, B = _2), residue numbers = UniProt positions
(AFDB full-length; see --verify), ipTM / max_ipSAE / max_pDockQ2, mean pLDDT of lining residues per chain,
atlas_pass (burrow_score>=0.03 & plddt_mean>=70; ALL pockets are kept).
Homodimers: sym_group (A<->B swapped residue set matches another pocket with overlap coefficient >= SYM_OVERLAP), is_sym_rep.
Dimer-only: max_overlap_iso = max overlap coefficient of the pocket's residue set with any pocket of the isolated
chain runs (both chains, same coordinates); dimer_only = is_interface & max_overlap_iso < DIMER_ONLY_OVERLAP.
Overlap coefficient = |X & Y| / min(|X|, |Y|). Also writes the exclusion/failure ledger.
Status: first attempt, not yet verified.
"""
import os
import argparse, glob, json, sys
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
DIMER_ROOT = Path(os.environ.get("DIMERS_ROOT", "/cluster/scratch/gmagnusson/dimers"))
KINDS = {"homo": "homodimers", "hetero": "heterodimers"}
MIN_RES_PRIMARY, MIN_RES_SENS = 3, 1          # DESIGN.md v3
SYM_OVERLAP, DIMER_ONLY_OVERLAP = 0.5, 0.5     # DESIGN.md v3 (sensitivities 0.2 / 0.8 from max_overlap_iso)
ATLAS_SCORE, ATLAS_PLDDT = 0.03, 70.0


def oc(a, b):
    return len(a & b) / min(len(a), len(b)) if a and b else 0.0


def read_chunks(d):
    fs = sorted(glob.glob(str(d / "chunk_[0-9][0-9][0-9][0-9].parquet")))
    return (pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()), fs


def plddt_map(path):
    import gemmi
    st = gemmi.read_structure(str(path))
    return {(ch.name, f"{r.seqid.num}{r.seqid.icode.strip()}"): r[0].b_iso if len(r) else np.nan
            for ch in st[0] for r in ch}, {ch.name: len(ch) for ch in st[0]}


def sym_groups(df):
    """df: pockets of ONE homodimer. Returns sym_group (int) and rep flag."""
    sets = [frozenset(s.split()) for s in df.residues]
    swap = {"A": "B", "B": "A"}
    sw = [frozenset(swap.get(x[0], x[0]) + x[1:] for x in s) for s in sets]
    n = len(sets); parent = list(range(n))
    def f(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for i in range(n):
        for j in range(n):
            if oc(sw[i], sets[j]) >= SYM_OVERLAP:
                parent[f(i)] = f(j)
    root = [f(i) for i in range(n)]
    score = df.burrow_score.values
    grp = {}
    for i, r in enumerate(root):
        grp.setdefault(r, []).append(i)
    gid = np.zeros(n, int); rep = np.zeros(n, bool)
    for k, (r, idx) in enumerate(sorted(grp.items())):
        best = max(idx, key=lambda i: (score[i], -i)); rep[best] = True
        for i in idx: gid[i] = k
    return gid, rep


def annotate(species, out_root):
    frames, ledger = [], []
    for kind in ("homo", "hetero"):
        d = out_root / species / kind
        p, fs = read_chunks(d)
        c = pd.read_csv(DIMER_ROOT / f"afdb_{species}_dimers" / KINDS[kind] / "candidates.csv").drop_duplicates("modelEntityId")
        # failures / exclusions
        sfs = sorted(glob.glob(str(d / "chunk_*.structures.parquet")))
        if sfs:
            s = pd.concat([pd.read_parquet(f) for f in sfs])
            for r in s[s.status != "ok"].itertuples():
                ledger.append((species, kind, Path(r.structure).stem, "describe_" + str(r.status), str(r.error)[:200]))
            for r in s[(s.status == "ok") & (s.n_pockets == 0)].itertuples():
                ledger.append((species, kind, Path(r.structure).stem, "no_pockets", ""))
        for m in c.modelEntityId:
            if not (DIMER_ROOT / f"afdb_{species}_dimers" / KINDS[kind] / f"{m}.pdb").exists():
                ledger.append((species, kind, m, "missing_file", "listed in candidates.csv, no PDB (download_error)"))
        n_chunks = -(-len(c.modelEntityId.unique()) // 500)
        for i in range(n_chunks):
            if not (d / f"chunk_{i:04d}.parquet").exists():
                ledger.append((species, kind, f"chunk_{i:04d}", "chunk_not_done", "failed or not yet run"))
        if len(p) == 0:
            continue
        p = p.merge(c, left_on="structure_stem", right_on="modelEntityId", how="left")
        p["kind"] = kind; p["species"] = species
        frames.append(p)
    if not frames:
        return pd.DataFrame(), ledger
    df = pd.concat(frames, ignore_index=True)
    # per-chain counts
    toks = df.residues.str.split()
    df["res_A"] = [[int(x[2:]) for x in t if x[0] == "A" and x[2:].lstrip("-").isdigit()] for t in toks]
    df["res_B"] = [[int(x[2:]) for x in t if x[0] == "B" and x[2:].lstrip("-").isdigit()] for t in toks]
    df["n_res_A"] = df.res_A.map(len); df["n_res_B"] = df.res_B.map(len)
    df["is_interface"] = (df.n_res_A >= MIN_RES_PRIMARY) & (df.n_res_B >= MIN_RES_PRIMARY)
    df["is_interface_ge1"] = (df.n_res_A >= MIN_RES_SENS) & (df.n_res_B >= MIN_RES_SENS)
    df["atlas_pass"] = (df.burrow_score >= ATLAS_SCORE) & (df.plddt_mean >= ATLAS_PLDDT)
    df = df.rename(columns={"uniprotAccession_1": "uniprot_A", "uniprotAccession_2": "uniprot_B"})
    df["uniprot_pos_A"] = df.res_A; df["uniprot_pos_B"] = df.res_B  # = resnum (full-length models; verify with --verify)
    # pLDDT of lining residues per chain
    pa, pb = np.full(len(df), np.nan), np.full(len(df), np.nan)
    for (kind, stem), g in df.groupby(["kind", "structure_stem"]):
        mp, _ = plddt_map(DIMER_ROOT / f"afdb_{species}_dimers" / KINDS[kind] / f"{stem}.pdb")
        for i, rs in zip(g.index, g.residues):
            for ch, arr in (("A", pa), ("B", pb)):
                v = [mp.get((ch, x[2:]), np.nan) for x in rs.split() if x[0] == ch]
                if v: arr[df.index.get_loc(i)] = np.nanmean(v)
    df["plddt_lining_A"], df["plddt_lining_B"] = pa, pb
    # symmetry groups (homodimers only)
    df["sym_group"] = -1; df["is_sym_rep"] = True
    for stem, g in df[df.kind == "homo"].groupby("structure_stem"):
        gid, rep = sym_groups(g)
        df.loc[g.index, "sym_group"] = gid; df.loc[g.index, "is_sym_rep"] = rep
    # dimer-only vs isolated chains
    iso, _ = read_chunks_iso(out_root / f"{species}_isolated")
    df["max_overlap_iso"] = np.nan; df["dimer_only"] = pd.NA
    if iso is not None:
        iso_sets = {}
        for stem_ch, rs in zip(iso.structure_stem, iso.residues):
            iso_sets.setdefault(stem_ch.rsplit("_", 1)[0], []).append(frozenset(rs.split()))
        for k, (stem, rs) in enumerate(zip(df.structure_stem, df.residues)):
            if not df.is_interface.iat[k] and not df.is_interface_ge1.iat[k]:
                continue
            s = frozenset(rs.split())
            df.iat[k, df.columns.get_loc("max_overlap_iso")] = max([oc(s, t) for t in iso_sets.get(stem, [])] + [0.0])
        has = df.max_overlap_iso.notna() & df.structure_stem.isin(iso_sets)
        df.loc[has, "dimer_only"] = df.loc[has, "is_interface"] & (df.loc[has, "max_overlap_iso"] < DIMER_ONLY_OVERLAP)
    df = df.drop(columns=["res_A", "res_B", "fingerprint", "structure"], errors="ignore")
    return df, ledger


def read_chunks_iso(d):
    if not d.exists():
        return None, []
    fs = []
    for kind in ("homo", "hetero"):
        x, f = read_chunks(d / kind)
        if len(x): fs.append(x)
    return (pd.concat(fs, ignore_index=True) if fs else None), []


def verify(species="human", n=2):
    """Check chain numbering is 1..N contiguous and sequence equals the UniProt sequence (needs internet)."""
    import gemmi, urllib.request
    c = pd.concat([pd.read_csv(DIMER_ROOT / f"afdb_{species}_dimers/{KINDS[k]}/candidates.csv").assign(kind=k) for k in KINDS])
    c = c.drop_duplicates("modelEntityId").sample(n, random_state=0)
    for r in c.itertuples():
        st = gemmi.read_structure(str(DIMER_ROOT / f"afdb_{species}_dimers/{KINDS[r.kind]}/{r.modelEntityId}.pdb"))
        for ch, acc in zip(st[0], (r.uniprotAccession_1, r.uniprotAccession_2)):
            seq = "".join(gemmi.find_tabulated_residue(x.name).one_letter_code.upper() for x in ch)
            nums = [x.seqid.num for x in ch]
            try:
                fa = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=30).read().decode().split("\n", 1)[1].replace("\n", "")
            except Exception as e:
                print(r.modelEntityId, ch.name, acc, "UniProt fetch failed:", e, "len", len(ch), "nums_1..N", nums == list(range(1, len(ch) + 1))); continue
            print(r.modelEntityId, ch.name, acc, "len", len(ch), "uniprot_len", len(fa), "seq_equal", seq == fa,
                  "nums_1..N", nums == list(range(1, len(ch) + 1)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="human")
    ap.add_argument("--pockets", default=str(WS / "data/part2/pockets"))
    ap.add_argument("--out", default=str(WS / "data/part2/human_dimer_pockets_annotated.parquet"))
    ap.add_argument("--ledger", default=str(WS / "results/part2/pocket_run_ledger.tsv"))
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args()
    if a.verify:
        verify(a.species); return
    df, ledger = annotate(a.species, Path(a.pockets))
    df.to_parquet(a.out, index=False)
    L = pd.DataFrame(ledger, columns=["species", "kind", "structure", "reason", "detail"])
    Path(a.ledger).parent.mkdir(parents=True, exist_ok=True)
    L.to_csv(a.ledger, sep="\t", index=False)
    print(len(df), "pockets", len(L), "ledger rows")


if __name__ == "__main__":
    main()
