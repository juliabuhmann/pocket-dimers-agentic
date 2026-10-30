#!/usr/bin/env python
"""Part 1 driver: PDB ligand-site census (single- vs multi-chain binding sites).

  python scripts/run_part1_census.py --sample          # 1% random sample (seed in config) -> data/part1/sample/
  python scripts/run_part1_census.py --full            # whole mirror -> data/part1/
  python scripts/run_part1_census.py --full --merge-only
Sentinel known-answer tests (tests/test_part1.py) run first and abort the run on failure.
Resumable: shards whose output parquet exists are skipped. Source scripts/env_threads.sh before running.
"""
import argparse
import hashlib
import multiprocessing as mp
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from part1 import merge as MG  # noqa: E402
from part1.config import (CENSUS_VERSION, N_WORKERS, OUT_DIR, PDB_MIRROR, SAMPLE_FRACTION, SAMPLE_SEED, SHARD_SIZE,  # noqa: E402
                          SUMMARY_DIR, raw_path)
from part1.entry import process_entry  # noqa: E402

SIFTS = None  # filled before fork, shared copy-on-write


def load_sifts():
    df = pd.read_csv(raw_path("sifts_chain_uniprot"), sep="\t", comment="#", dtype=str)
    for c in ("RES_BEG", "RES_END", "SP_BEG"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["RES_BEG", "RES_END", "SP_BEG"])
    out = {}
    for p, ch, acc, rb, re_, sb in zip(df.PDB, df.CHAIN, df.SP_PRIMARY, df.RES_BEG, df.RES_END, df.SP_BEG):
        out.setdefault((p, ch), []).append((acc, int(rb), int(re_), int(sb)))
    return out


def list_entries():
    ids = []
    for d in sorted(os.scandir(PDB_MIRROR), key=lambda e: e.name):
        if d.is_dir():
            ids += [f.name[:-7] for f in os.scandir(d.path) if f.name.endswith(".cif.gz")]
    return sorted(ids)


def entry_path(pid):
    return PDB_MIRROR / pid[1:3] / f"{pid}.cif.gz"


def run_shard(args):
    i, ids, shard_dir = args
    eo = shard_dir / f"entries_{i:05d}.parquet"
    lo = shard_dir / f"ligands_{i:05d}.parquet"
    if eo.exists() and lo.exists():
        return i, 0, 0.0
    t = time.time()
    erows, lrows = [], []
    for pid in ids:
        r, l = process_entry(entry_path(pid), SIFTS, pid)
        erows.append(r)
        lrows += l
    pd.DataFrame(lrows).to_parquet(lo.with_suffix(".tmp"), index=False)
    pd.DataFrame(erows).to_parquet(eo.with_suffix(".tmp"), index=False)
    lo.with_suffix(".tmp").rename(lo)
    eo.with_suffix(".tmp").rename(eo)
    return i, len(ids), time.time() - t


def sentinel():
    log = ROOT / "logs" / "part1_sentinel.log"
    with open(log, "w") as f:
        r = subprocess.run([sys.executable, "-m", "pytest", "-q", str(ROOT / "tests" / "test_part1.py")],
                           stdout=f, stderr=subprocess.STDOUT, cwd=ROOT)
    print(f"[sentinel] pytest exit={r.returncode} (log {log})", flush=True)
    if r.returncode != 0:
        sys.exit("Sentinel known-answer tests FAILED - aborting run")


def do_merge(out_dir: Path, tag: str):
    t = time.time()
    shard_dir = out_dir / "shards"
    ent = MG.concat_shards(shard_dir, "entries").to_pandas()
    big = MG.concat_shards(shard_dir, "ligands")
    slim = big.drop([c for c in MG.SLIM_DROP if c in big.column_names]).to_pandas()
    classes = pd.read_parquet(OUT_DIR / "ligand_classes.parquet")
    slim = MG.annotate(slim, ent, classes)
    cutoff = slim.attrs.get("biolip_cutoff")
    new_cols = [c for c in slim.columns if c not in big.column_names]
    add = pa.Table.from_pandas(slim[new_cols].reset_index(drop=True), preserve_index=False)
    for c in new_cols:
        big = big.append_column(c, add[c])
    pq.write_table(big, out_dir / "ligand_sites.parquet")
    ent["census_version"] = CENSUS_VERSION
    ent["is_photosynthetic"] = ent.pdb_id.isin(slim.attrs.get("photosynthetic_entries", set()))
    ent.to_parquet(out_dir / "entries.parquet", index=False)
    summ_dir = SUMMARY_DIR if tag == "full" else SUMMARY_DIR / tag
    summ_dir.mkdir(parents=True, exist_ok=True)
    led = MG.ledger(ent, slim)
    led.to_csv(summ_dir / "exclusion_ledger.tsv", sep="\t", index=False)
    for name, df in MG.summaries(slim).items():
        df.to_csv(summ_dir / f"{name}.tsv", sep="\t", index=False, float_format="%.4f")
    ent.groupby(["skip_reason"], dropna=False).size().rename("n").reset_index().to_csv(
        summ_dir / "entry_skip_reasons.tsv", sep="\t", index=False)
    ok = ent[ent.skip_reason.isna()]
    ok.groupby(["method_class", "assembly_type"]).size().rename("n_entries").reset_index().to_csv(
        summ_dir / "entries_by_method_assembly.tsv", sep="\t", index=False)
    with open(summ_dir / "run_info.tsv", "w") as f:
        f.write(f"key\tvalue\ncensus_version\t{CENSUS_VERSION}\n"
                f"n_photosynthetic_entries\t{len(slim.attrs.get('photosynthetic_entries', []))}\nn_entries_listed\t{len(ent)}\nn_ligand_rows\t{len(slim)}\n"
                f"biolip_coverage_release_cutoff\t{cutoff}\n"
                f"sum_entry_runtime_s\t{ent.runtime_s.sum():.0f}\nmerge_s\t{time.time() - t:.0f}\n")
    print(led.to_string(index=False))
    print(f"[merge] done in {time.time() - t:.0f}s -> {out_dir}")


def main():
    global SIFTS
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--sample", action="store_true")
    g.add_argument("--full", action="store_true")
    g.add_argument("--ids", nargs="+", help="explicit PDB ids (debug), output data/part1/debug/")
    ap.add_argument("--workers", type=int, default=N_WORKERS)
    ap.add_argument("--merge-only", action="store_true")
    ap.add_argument("--no-sentinel", action="store_true", help="debug only; never for production runs")
    a = ap.parse_args()
    tag = "sample" if a.sample else ("full" if a.full else "debug")
    out_dir = OUT_DIR if tag == "full" else OUT_DIR / tag
    shard_dir = out_dir / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    if not a.merge_only:
        if not a.no_sentinel:
            sentinel()
        t = time.time()
        ids = list_entries() if not a.ids else sorted(a.ids)
        idlist_sha = hashlib.sha256("\n".join(ids).encode()).hexdigest()
        print(f"[list] {len(ids)} entries in mirror (sha256 of sorted id list {idlist_sha}) {time.time() - t:.0f}s", flush=True)
        if a.sample:
            rng = np.random.default_rng(SAMPLE_SEED)
            k = int(round(len(ids) * SAMPLE_FRACTION))
            ids = sorted(rng.choice(ids, size=k, replace=False).tolist())
        (out_dir / "entry_ids.txt").write_text("\n".join(ids) + "\n")
        t = time.time()
        SIFTS = load_sifts()
        print(f"[sifts] {len(SIFTS)} chains {time.time() - t:.0f}s", flush=True)
        size = SHARD_SIZE if a.full else max(1, min(SHARD_SIZE, len(ids) // (a.workers * 4)))  # small runs: spread work
        shards = [(i, ids[s:s + size], shard_dir) for i, s in enumerate(range(0, len(ids), size))]
        todo = [s for s in shards if not ((shard_dir / f"entries_{s[0]:05d}.parquet").exists()
                                          and (shard_dir / f"ligands_{s[0]:05d}.parquet").exists())]
        print(f"[run] {len(shards)} shards, {len(todo)} to do, {a.workers} workers", flush=True)
        t = time.time()
        done = 0
        ctx = mp.get_context("fork")
        with ctx.Pool(processes=a.workers, maxtasksperchild=10) as pool:
            for i, n, dt in pool.imap_unordered(run_shard, todo):
                done += 1
                el = time.time() - t
                print(f"[shard] {i} n={n} {dt:.0f}s | {done}/{len(todo)} elapsed {el / 60:.1f} min "
                      f"eta {el / done * (len(todo) - done) / 60:.1f} min", flush=True)
        print(f"[run] finished {len(todo)} shards in {(time.time() - t) / 60:.1f} min", flush=True)
    do_merge(out_dir, tag)


if __name__ == "__main__":
    main()
