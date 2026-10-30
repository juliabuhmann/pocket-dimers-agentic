"""Verifier: humsavar P/LP counts from raw humsavar.txt and wt-vs-AFDB-dimer-model residue consistency (not verified itself)."""
from pathlib import Path
import re, numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
ROOT = str(Path(__file__).resolve().parents[2])
DIMERS = "/cluster/scratch/gmagnusson/dimers/afdb_human_dimers"
OUT = f"{ROOT}/reviews/verify_h5/out"
THREE = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I", LEU="L", LYS="K",
             MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
HUMSAVAR = f"{ROOT}/data/raw/humsavar/humsavar.txt"
VAR_RE = re.compile(r"^(\S+)\s+(\S+)\s+(VAR_\S+)\s+p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\s+(LP/P|LB/B|US)\s+(\S+)\s*(.*)$")

def ca_seq(args):
    mid, path, a1, a2 = args
    out = {}
    with open(path) as fh:
        for l in fh:
            if l.startswith("ATOM") and l[12:16] == " CA ":
                ch = l[21]; acc = a1 if ch == "A" else a2
                out[(acc, int(l[22:26]))] = THREE.get(l[17:20], "X")
    return mid, out

def main():
    rows = []
    for l in open(HUMSAVAR):
        m = VAR_RE.match(l.rstrip("\n"))
        if m:
            g, acc, ft, w, pos, mu, cc, dbsnp, dis = m.groups()
            rows.append((acc, g, ft, int(pos), THREE.get(w.upper()), THREE.get(mu.upper()), cc))
    H = pd.DataFrame(rows, columns=["accession", "gene", "ftid", "position", "wt", "mut", "cc"])
    H = H.dropna(subset=["wt", "mut"])
    print("raw missense rows", len(H), H.cc.value_counts().to_dict())
    c = []
    for d in ["homodimers", "heterodimers"]:
        x = pd.read_csv(f"{DIMERS}/{d}/candidates.csv"); x["path"] = [f"{DIMERS}/{d}/{m}.pdb" for m in x.modelEntityId]; c.append(x)
    C = pd.concat(c).drop_duplicates("modelEntityId")
    seq = {}; conflicts = 0
    with ProcessPoolExecutor(6) as ex:
        for mid, o in ex.map(ca_seq, [(r.modelEntityId, r.path, r.uniprotAccession_1, r.uniprotAccession_2) for r in C.itertuples()], chunksize=20):
            for k, v in o.items():
                if k in seq and seq[k] != v: conflicts += 1
                seq[k] = v
    accs = set(C.uniprotAccession_1) | set(C.uniprotAccession_2)
    P = H[H.cc == "LP/P"].copy()
    P["in_model_protein"] = P.accession.isin(accs)
    P["model_aa"] = [seq.get((a, p)) for a, p in zip(P.accession, P.position)]
    Q = P[P.in_model_protein]
    res = dict(raw_plp_missense=len(P), plp_in_model_proteins=len(Q),
               plp_pos_absent_in_model=int(Q.model_aa.isna().sum()),
               plp_wt_mismatch=int(((Q.model_aa != Q.wt) & Q.model_aa.notna()).sum()),
               plp_unique_residues=int(Q[["accession", "position"]].drop_duplicates().shape[0]),
               model_residue_conflicts_between_models=conflicts, n_model_proteins=len(accs))
    # all humsavar rows, for comparison
    H["model_aa"] = [seq.get((a, p)) for a, p in zip(H.accession, H.position)]
    Hm = H[H.model_aa.notna()]
    res["all_rows_in_models"] = len(Hm); res["all_rows_mismatch"] = int((Hm.model_aa != Hm.wt).sum())
    for k, v in res.items(): print(k, v)
    pd.Series(res).to_csv(f"{OUT}/seqcheck.tsv", sep="\t", header=False)
    Q.to_parquet(f"{OUT}/plp_in_models.parquet")
    pd.DataFrame([(a, p, v) for (a, p), v in seq.items()], columns=["accession", "position", "aa"]).to_parquet(f"{OUT}/model_seq.parquet")

if __name__ == "__main__":
    main()
