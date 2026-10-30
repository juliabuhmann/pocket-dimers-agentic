#!/usr/bin/env python
"""Q2.4 step 1: annotate non-human dimer pockets and extract chain sequences from every dimer model.

- Pocket annotation: scripts/annotate_pocket_chains.annotate() unchanged (interface = >=3 lining residues per chain, primary),
  one parquet per species -> data/interim/q24/annot_<species>.parquet (+ ledger_<species>.tsv). Non-human species have no
  isolated-chain runs, so max_overlap_iso / dimer_only are empty there (not needed: the ortholog side is only a pocket set).
- Sequences: one-letter sequence and residue numbers of chains A and B read from each model PDB (gemmi), mapped to UniProt
  accessions via candidates.csv (A = uniprotAccession_1, B = uniprotAccession_2; docs/DATA_RECON.md).
  -> data/interim/q24/chains_<species>.parquet (one row per model chain) and seqs_<species>.fasta (one sequence per accession;
  the longest if an accession occurs with different sequences, which is counted in the log).
Status: first attempt, not yet verified.
"""
import sys, argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

WS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WS / "scripts"))
from annotate_pocket_chains import annotate, DIMER_ROOT, KINDS  # noqa: E402

POCKETS = WS / "data/part2/pockets"
INTERIM = WS / "data/interim/q24"
SPECIES = ["arabidopsis", "human", "mouse", "rat", "zebrafish", "fly", "worm", "yeast", "pombe", "dicty", "ecoli", "mtb", "paeruginosa"]


def chains(species):
    import gemmi
    rows = []
    for kind, sub in KINDS.items():
        c = pd.read_csv(DIMER_ROOT / f"afdb_{species}_dimers" / sub / "candidates.csv").drop_duplicates("modelEntityId")
        for r in c.itertuples():
            p = DIMER_ROOT / f"afdb_{species}_dimers" / sub / f"{r.modelEntityId}.pdb"
            if not p.exists():
                continue
            st = gemmi.read_structure(str(p))
            for ch, acc in zip(st[0], (r.uniprotAccession_1, r.uniprotAccession_2)):
                seq = "".join(gemmi.find_tabulated_residue(x.name).one_letter_code.upper() if gemmi.find_tabulated_residue(x.name) else "X" for x in ch)
                nums = [x.seqid.num for x in ch]
                rows.append((species, kind, r.modelEntityId, ch.name, acc, seq, nums[0], nums[-1], nums == list(range(nums[0], nums[0] + len(nums)))))
    return pd.DataFrame(rows, columns=["species", "kind", "modelEntityId", "chain", "accession", "seq", "first_num", "last_num", "contiguous"])


def run(species):
    INTERIM.mkdir(parents=True, exist_ok=True)
    msg = []
    ch = chains(species)
    ch.to_parquet(INTERIM / f"chains_{species}.parquet", index=False)
    per = ch.groupby("accession").seq.agg(lambda s: max(set(s), key=len))
    n_multi = int((ch.groupby("accession").seq.nunique() > 1).sum())
    with open(INTERIM / f"seqs_{species}.fasta", "w") as f:
        for acc, s in per.items():
            f.write(f">{acc}\n{s}\n")
    msg.append(f"{species}: {len(ch)} chains, {len(per)} accessions, {n_multi} with >1 distinct sequence, "
               f"{int((~ch.contiguous).sum())} non-contiguous chains, {int((ch.first_num != 1).sum())} not starting at 1")
    if species != "human":
        out = INTERIM / f"annot_{species}.parquet"
        if not out.exists():
            df, ledger = annotate(species, POCKETS)
            df.to_parquet(out, index=False)
            pd.DataFrame(ledger, columns=["species", "kind", "structure", "reason", "detail"]).to_csv(INTERIM / f"ledger_{species}.tsv", sep="\t", index=False)
            msg.append(f"{species}: {len(df)} pockets annotated, {len(ledger)} ledger rows")
    return "\n".join(msg)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)   # docs/COMPUTE_POLICY.md; task cap 12
    a = ap.parse_args()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for m in ex.map(run, SPECIES):
            print(m, flush=True)
