#!/usr/bin/env python
"""Run `pocketcompare describe` (atlas settings) on predicted dimer models, chunked + resumable.

Atlas settings: --af-mode, unit merged, min-points 10, --backend fast. ALL pockets are kept; apply the atlas
filter downstream (burrow_score >= 0.03 and plddt_mean >= 70).

  python scripts/run_dimer_pockets.py --species human --set homo --workers 16
  sets: homo | hetero | isolated   (isolated = each chain of every homo+hetero model written alone, same coords)
Outputs: <out>/<species>/<homo|hetero>/chunk_NNNN.parquet (+ .structures.parquet); isolated ->
         <out>/<species>_isolated/<homo|hetero>/chunk_NNNN.parquet. Existing chunks are skipped.
Status: first attempt, not yet verified.
"""
import argparse, json, os, subprocess, sys, shutil, time
from pathlib import Path
import pandas as pd

WS = Path(__file__).resolve().parents[1]
DIMER_ROOT = Path(os.environ.get("DIMERS_ROOT", "/path/to/afdb_dimers"))
KINDS = {"homo": "homodimers", "hetero": "heterodimers"}


def load_models(species, kind):
    d = DIMER_ROOT / f"afdb_{species}_dimers" / KINDS[kind]
    c = pd.read_csv(d / "candidates.csv").drop_duplicates("modelEntityId").sort_values("modelEntityId")
    c["path"] = [str(d / f"{m}.pdb") for m in c.modelEntityId]
    c = c[[os.path.exists(p) for p in c.path]]
    return c.reset_index(drop=True)


def write_chains(paths, tmpdir):
    import gemmi
    out = []
    tmpdir.mkdir(parents=True, exist_ok=True)
    for p in paths:
        st = gemmi.read_structure(p)
        stem = Path(p).stem
        for ch in st[0]:
            s2 = gemmi.Structure()
            s2.spacegroup_hm = st.spacegroup_hm
            s2.cell = st.cell
            m = gemmi.Model("1")
            m.add_chain(ch)
            s2.add_model(m)
            f = tmpdir / f"{stem}_{ch.name}.pdb"
            s2.write_pdb(str(f))
            out.append(str(f))
    return out


def tool_commits():
    r = {}
    for t in ("burrow", "pocketcompare"):
        try:
            r[t] = subprocess.check_output(["git", "-C", str(WS / "tools" / t), "rev-parse", "HEAD"], text=True).strip()
        except Exception as e:
            r[t] = f"unknown ({e})"
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--set", required=True, choices=["homo", "hetero", "isolated"])
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--out", default=str(WS / "data/part2/pockets"))
    ap.add_argument("--chunk", type=int, default=500, help="models per chunk")
    ap.add_argument("--limit", type=int, default=None, help="debug: first N models per kind")
    a = ap.parse_args()
    if a.workers > 16:
        sys.exit("policy: max 16 workers")
    iso = a.set == "isolated"
    kinds = ["homo", "hetero"] if iso else [a.set]
    sp_out = Path(a.out) / (a.species + ("_isolated" if iso else ""))
    for kind in kinds:
        models = load_models(a.species, kind)
        if a.limit:
            models = models.iloc[:a.limit]
        od = sp_out / kind
        od.mkdir(parents=True, exist_ok=True)
        (od / "run_meta.json").write_text(json.dumps(dict(
            species=a.species, kind=kind, isolated=iso, n_models=len(models), chunk=a.chunk, commits=tool_commits(),
            cmd="pocketcompare describe <files> --af-mode --backend fast --unit merged --min-points 10",
            workers=a.workers, started=time.strftime("%FT%T")), indent=1))
        n_ch = (len(models) + a.chunk - 1) // a.chunk
        for i in range(n_ch):
            outp = od / f"chunk_{i:04d}.parquet"
            if outp.exists():
                continue
            t0 = time.time()
            paths = list(models.path.iloc[i * a.chunk:(i + 1) * a.chunk])
            tmp = None
            if iso:
                tmp = WS / "data/interim/isolated_chains" / a.species / kind / f"chunk_{i:04d}"
                paths = write_chains(paths, tmp)
            tmpout = od / f"chunk_{i:04d}.partial.parquet"
            cmd = ["pocketcompare", "describe", *paths, "--af-mode", "--unit", "merged", "--min-points", "10",
                   "--backend", "fast", "--workers", str(a.workers), "-o", str(tmpout)]
            r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if r.returncode != 0 or not tmpout.exists():
                print(f"[{a.species} {kind}] chunk {i} FAILED rc={r.returncode}: {r.stdout[-500:]}", flush=True)
                continue
            sp = Path(str(tmpout) + ".structures.parquet") if Path(str(tmpout) + ".structures.parquet").exists() else tmpout.with_suffix(".structures.parquet")
            if sp.exists():
                sp.rename(od / f"chunk_{i:04d}.structures.parquet")
            tmpout.rename(outp)  # atomic marker of completion
            if tmp:
                shutil.rmtree(tmp, ignore_errors=True)
            print(f"[{a.species} {kind}] chunk {i+1}/{n_ch} done in {time.time()-t0:.0f}s ({len(paths)} inputs)", flush=True)


if __name__ == "__main__":
    main()
