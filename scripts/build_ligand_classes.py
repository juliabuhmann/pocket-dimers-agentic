#!/usr/bin/env python
"""Build data/part1/ligand_classes.parquet from the fetched CCD + lists (run scripts/fetch_part1_sources.py first)."""
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from part1.config import OUT_DIR, SUMMARY_DIR  # noqa: E402
from part1.ligands import build_table  # noqa: E402

if __name__ == "__main__":
    t = time.time()
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    df = build_table(limit)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT_DIR / "ligand_classes.parquet", index=False)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    summ = df.groupby(["ligand_class", "marsh_class"]).size().rename("n_ccd").reset_index()
    summ.to_csv(SUMMARY_DIR / "ccd_class_counts.tsv", sep="\t", index=False)
    print(summ.to_string(index=False))
    print(f"{len(df)} CCD components in {time.time() - t:.0f}s")
