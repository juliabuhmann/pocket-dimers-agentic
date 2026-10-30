#!/usr/bin/env python
"""Annotate burrow / pocketcompare pockets with the chain(s) lining them.

burrow writes `residues` as space-separated `chain:resnum[icode]` (author numbering); it has no explicit
chain-count column in its CSV output. This script derives, per pocket:
  chains        sorted chain ids lining the pocket, e.g. "A,B"
  n_chains      number of distinct chains
  n_res_<X>     not written; instead `res_per_chain` = "A:7,B:5"
  min_chain_frac  fraction of lining residues on the least-represented chain (0 for single-chain pockets)
and prints a summary (fraction of pockets with >= 2 lining chains).

Inputs: a burrow output directory (reads <stem>.csv, skips *.points.csv and timing.csv) or a
pocketcompare describe parquet (column `residues`).
Note: "multi-chain" here means >= 1 lining residue from each of >= 2 chains (burrow's residue
assignment); a stricter interface definition (e.g. >= 3 residues per chain) can be applied downstream
via MIN_RES_PER_CHAIN.
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import pandas as pd

MIN_RES_PER_CHAIN = 3  # stricter interface definition, mirroring burrow's dimer benchmark (>= 3 residues on >= 2 chains)


def chain_stats(residues: str) -> dict:
    toks = residues.split() if isinstance(residues, str) else []
    c = Counter(t.split(":", 1)[0] for t in toks)
    n = sum(c.values())
    return {
        "chains": ",".join(sorted(c)),
        "n_chains": len(c),
        "res_per_chain": ",".join(f"{k}:{v}" for k, v in sorted(c.items())),
        "min_chain_frac": (min(c.values()) / n) if len(c) > 1 else 0.0,
        f"n_chains_ge{MIN_RES_PER_CHAIN}": sum(v >= MIN_RES_PER_CHAIN for v in c.values()),
    }


def load(path: Path) -> pd.DataFrame:
    if path.is_dir():
        dfs = []
        for f in sorted(path.glob("*.csv")):
            if f.name.endswith(".points.csv") or f.name == "timing.csv":
                continue
            d = pd.read_csv(f)
            d.insert(0, "stem", f.stem)
            dfs.append(d)
        return pd.concat(dfs, ignore_index=True)
    return pd.read_parquet(path, columns=None).drop(columns=["fingerprint", "n1"], errors="ignore")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("inp", type=Path, help="burrow output dir or pocketcompare parquet")
    ap.add_argument("-o", "--out", type=Path, help="annotated table (.tsv)")
    ap.add_argument("--min-score", type=float, default=None, help="only pockets with score >= this")
    a = ap.parse_args()
    df = load(a.inp)
    score_col = "score" if "score" in df.columns else "burrow_score"
    if a.min_score is not None:
        df = df[df[score_col] >= a.min_score]
    df = pd.concat([df.reset_index(drop=True), pd.DataFrame([chain_stats(r) for r in df["residues"]])], axis=1)
    k = f"n_chains_ge{MIN_RES_PER_CHAIN}"
    n_struct = df.iloc[:, 0].nunique()
    print(f"{len(df)} pockets in {n_struct} structures; "
          f">=2 chains (any residue): {(df.n_chains >= 2).mean():.3f}; "
          f">=2 chains with >={MIN_RES_PER_CHAIN} residues each: {(df[k] >= 2).mean():.3f}")
    if a.out:
        df.to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
