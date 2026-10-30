#!/usr/bin/env python
"""ChEMBL annotations for the CCD->ChEMBL mapping (UniChem src1src3).

1. build ccd_to_chembl.parquet from data/raw/unichem/src1src3.txt.gz (file header: From src:'1' (ChEMBL) To src:'3' (PDBe))
2. max_phase for all mapped ChEMBL ids via REST (molecule_chembl_id__in, batches of 200 (URL length limit), sequential, polite)
3. pChEMBL>=6 flag (function `fetch_pchembl`): per-molecule activity query with limit=1; restricted to CCDs listed in
   data/part1/ligand_classes.parquet (column 'ccd' or 'comp_id' if present) unless --all-pchembl.
Outputs: data/processed/ccd_to_chembl.parquet, data/raw/chembl/{status.json,max_phase.tsv,pchembl.tsv}
and rows in DATA.md.  Usage: python scripts/fetch_chembl.py [--maxphase] [--pchembl] [--all-pchembl]
"""
import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from extlib import PROCESSED, RAW, ROOT, UA, append_data_md, now, sha256  # noqa: E402

API = "https://www.ebi.ac.uk/chembl/api/data"
UNICHEM = RAW / "unichem" / "src1src3.txt.gz"
OUT = RAW / "chembl"
BATCH = 200  # GET request-line limit is 4094 bytes (1000 ids -> HTTP 400/500; POST override not accepted)
SLEEP = 0.25  # polite delay after each request (per worker; <=4 workers)
PCHEMBL_MIN = 6
LIG_CLASSES = ROOT / "data" / "part1" / "ligand_classes.parquet"


def get_json(url, params=None, retries=5):
    for i in range(retries):
        try:
            r = requests.get(url, params=params, headers=UA, timeout=120)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(5 * (i + 1))
                continue
            r.raise_for_status()
        except requests.RequestException:
            time.sleep(5 * (i + 1))
    raise RuntimeError(f"failed {url} {params and str(params)[:100]}")


def build_ccd_map():
    df = pd.read_csv(UNICHEM, sep="\t", skiprows=1, names=["chembl_id", "ccd"], dtype=str)
    df = df.drop_duplicates()
    PROCESSED.mkdir(parents=True, exist_ok=True)
    df.to_parquet(PROCESSED / "ccd_to_chembl.parquet", index=False)
    print("ccd_to_chembl", df.shape, "unique ccd", df.ccd.nunique(), "unique chembl", df.chembl_id.nunique())
    return df


def chembl_release():
    OUT.mkdir(parents=True, exist_ok=True)
    st = get_json(f"{API}/status.json")
    (OUT / "status.json").write_text(json.dumps(st, indent=1))
    return st


def fetch_max_phase(ids):
    rows = []
    ids = sorted(set(ids))
    for i in range(0, len(ids), BATCH):
        b = ids[i:i + BATCH]
        j = get_json(f"{API}/molecule.json", {"molecule_chembl_id__in": ",".join(b), "limit": 1000,
                                              "only": "molecule_chembl_id,max_phase,molecule_type"})
        rows += [(m["molecule_chembl_id"], m.get("max_phase"), m.get("molecule_type")) for m in j["molecules"]]
        print(f"max_phase {i + len(b)}/{len(ids)}", flush=True)
        time.sleep(SLEEP)
    return pd.DataFrame(rows, columns=["chembl_id", "max_phase", "molecule_type"])


def has_pchembl(cid):
    j = get_json(f"{API}/activity.json", {"molecule_chembl_id": cid, "pchembl_value__gte": PCHEMBL_MIN,
                                          "limit": 1, "only": "molecule_chembl_id"})
    time.sleep(SLEEP)
    return cid, j["page_meta"]["total_count"] > 0


def fetch_pchembl(ids, workers=3):
    """Rerunnable: resumes from data/raw/chembl/pchembl.tsv. Returns DataFrame(chembl_id, has_pchembl_ge6)."""
    f = OUT / "pchembl.tsv"
    done = {}
    if f.exists():
        d = pd.read_csv(f, sep="\t")
        done = dict(zip(d.chembl_id, d.has_pchembl_ge6))
    todo = [c for c in sorted(set(ids)) if c not in done]
    print("pchembl todo", len(todo), "done", len(done), flush=True)
    with ThreadPoolExecutor(min(workers, 4)) as ex, open(f, "a") as fh:
        if not done:
            fh.write("chembl_id\thas_pchembl_ge6\n")
        for n, (cid, ok) in enumerate(ex.map(has_pchembl, todo), 1):
            done[cid] = ok
            fh.write(f"{cid}\t{ok}\n")
            if n % 200 == 0:
                fh.flush()
                print("pchembl", n, "/", len(todo), flush=True)
    return pd.DataFrame({"chembl_id": list(done), "has_pchembl_ge6": list(done.values())})


def log(path, name, st, url):
    append_data_md(name, f"ChEMBL {st['chembl_db_version']} ({st['chembl_release_date']})", sha256(path),
                   "CC BY-SA 3.0 (ChEMBL)", url, str(Path(path).relative_to(ROOT)), Path(path).stat().st_size)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxphase", action="store_true")
    ap.add_argument("--pchembl", action="store_true")
    ap.add_argument("--all-pchembl", action="store_true")
    a = ap.parse_args()
    m = build_ccd_map()
    st = chembl_release()
    print(st)
    if a.maxphase or not a.pchembl:
        mp = fetch_max_phase(m.chembl_id)
        mp.to_csv(OUT / "max_phase.tsv", sep="\t", index=False)
        log(OUT / "max_phase.tsv", "chembl_max_phase_for_unichem_ids", st, f"{API}/molecule.json?molecule_chembl_id__in=...")
        ids = m.chembl_id
    if a.pchembl:
        ids = m.chembl_id
        if not a.all_pchembl:
            if not LIG_CLASSES.exists():
                sys.exit(f"{LIG_CLASSES} missing; rerun later or use --all-pchembl")
            lc = pd.read_parquet(LIG_CLASSES)
            col = next(c for c in ("ccd", "comp_id", "ligand_id", "het_id") if c in lc.columns)
            ids = m[m.ccd.isin(set(lc[col]))].chembl_id
        pc = fetch_pchembl(ids)
        log(OUT / "pchembl.tsv", "chembl_pchembl_ge6_flag", st, f"{API}/activity.json?pchembl_value__gte=6")
    # merged table
    mp = pd.read_csv(OUT / "max_phase.tsv", sep="\t")
    out = m.merge(mp, on="chembl_id", how="left")
    if (OUT / "pchembl.tsv").exists():
        out = out.merge(pd.read_csv(OUT / "pchembl.tsv", sep="\t"), on="chembl_id", how="left")
    out["chembl_release"] = st["chembl_db_version"]
    out.to_parquet(PROCESSED / "ccd_to_chembl.parquet", index=False)
    print(out.shape, out.max_phase.describe().to_dict())
