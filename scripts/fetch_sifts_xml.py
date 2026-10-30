#!/usr/bin/env python
"""Per-entry SIFTS XML fetcher + per-residue parser.

fetch_sifts_xml(pdb_ids, workers=4) -> downloads https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/<id>.xml.gz
into data/raw/sifts_xml/ (skips existing; <=4 parallel; polite retries). Failures are logged to
data/raw/sifts_xml/FAILED.txt and returned. A provenance row for the fetch (date, count) is
appended to data/raw/sifts_xml/FETCH_LOG.tsv; the SIFTS release date is the server Last-Modified per file
(stored in FETCH_LOG.tsv per entry batch row-wise: id, last_modified, sha256).

parse_sifts_xml(path) -> DataFrame with one row per (residue, UniProt xref):
  pdb_id, chain (author/PDB chain), auth_resnum, ins_code, label_resnum (PDBe entity seq numbering),
  pdb_resname (3-letter), uniprot_acc, uniprot_pos, uniprot_resname (1-letter), observed (bool)
Residues without a UniProt xref get uniprot_* = None.

CLI: python scripts/fetch_sifts_xml.py 1hvr 4hhb 2lyz [--parse out.parquet]  or --ids-file ids.txt
"""
import argparse
import gzip
import re
import sys
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from extlib import RAW, UA, now, sha256  # noqa: E402

URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{id}.xml.gz"
CACHE = RAW / "sifts_xml"
MAX_WORKERS = 4  # COMPUTE_POLICY
_log_lock = threading.Lock()
_RES = re.compile(r"^(-?\d+)([A-Za-z]?)$")


def _fetch_one(pid, retries=4):
    dest = CACHE / f"{pid}.xml.gz"
    if dest.exists() and dest.stat().st_size > 0:
        return pid, "cached", None
    for i in range(retries):
        try:
            r = requests.get(URL.format(id=pid), headers=UA, timeout=60)
            if r.status_code == 404:
                return pid, "404", None
            r.raise_for_status()
            tmp = dest.with_suffix(".part")
            tmp.write_bytes(r.content)
            tmp.replace(dest)
            return pid, "ok", r.headers.get("Last-Modified", "")
        except Exception as e:  # noqa: BLE001
            if i == retries - 1:
                return pid, f"error:{e!r}", None
            time.sleep(2 * (i + 1))


def fetch_sifts_xml(pdb_ids, workers=MAX_WORKERS):
    """Fetch (cached) SIFTS XMLs. Returns dict id -> status."""
    CACHE.mkdir(parents=True, exist_ok=True)
    ids = sorted({p.lower() for p in pdb_ids})
    workers = min(workers, MAX_WORKERS)
    status = {}
    with ThreadPoolExecutor(workers) as ex:
        for pid, st, lm in ex.map(_fetch_one, ids):
            status[pid] = st
            if st == "ok":
                with _log_lock, open(CACHE / "FETCH_LOG.tsv", "a") as f:
                    f.write(f"{pid}\t{now()}\t{lm}\t{sha256(CACHE / (pid + '.xml.gz'))}\n")
    bad = {k: v for k, v in status.items() if v not in ("ok", "cached")}
    if bad:
        with open(CACHE / "FAILED.txt", "a") as f:
            for k, v in bad.items():
                f.write(f"{k}\t{v}\n")
    return status


def parse_sifts_xml(path):
    path = Path(path)
    pid = path.name.split(".")[0]
    rows = []
    with gzip.open(path, "rb") as fh:
        for _, el in ET.iterparse(fh):
            tag = el.tag.split("}")[-1]
            if tag != "residue":
                continue
            label = int(el.get("dbResNum")) if (el.get("dbResNum") or "").lstrip("-").isdigit() else None
            pdb = None
            ups = []
            observed = True
            for c in el:
                ct = c.tag.split("}")[-1]
                if ct == "crossRefDb":
                    src = c.get("dbSource")
                    if src == "PDB" and pdb is None:
                        pdb = c
                    elif src == "UniProt":
                        ups.append(c)
                elif ct == "residueDetail" and (c.text or "").strip() == "Not_Observed":
                    observed = False
            chain = resnum = ins = resname = None
            if pdb is not None:
                chain = pdb.get("dbChainId")
                resname = pdb.get("dbResName")
                m = _RES.match(pdb.get("dbResNum") or "")
                if m:
                    resnum, ins = int(m.group(1)), m.group(2)
                else:
                    observed = False  # 'null' auth number
            base = dict(pdb_id=pid, chain=chain, auth_resnum=resnum, ins_code=ins,
                        label_resnum=label, pdb_resname=resname, observed=observed)
            if ups:
                for u in ups:
                    rows.append({**base, "uniprot_acc": u.get("dbAccessionId"),
                                 "uniprot_pos": int(u.get("dbResNum")), "uniprot_resname": u.get("dbResName")})
            else:
                rows.append({**base, "uniprot_acc": None, "uniprot_pos": None, "uniprot_resname": None})
            el.clear()
    df = pd.DataFrame(rows)
    df["auth_resnum"] = df["auth_resnum"].astype("Int64")
    df["label_resnum"] = df["label_resnum"].astype("Int64")
    df["uniprot_pos"] = df["uniprot_pos"].astype("Int64")
    return df


def parse_many(pdb_ids):
    return pd.concat([parse_sifts_xml(CACHE / f"{p.lower()}.xml.gz") for p in pdb_ids], ignore_index=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="*")
    ap.add_argument("--ids-file")
    ap.add_argument("--workers", type=int, default=MAX_WORKERS)
    ap.add_argument("--parse", help="write per-residue parquet here")
    a = ap.parse_args()
    ids = list(a.ids)
    if a.ids_file:
        ids += Path(a.ids_file).read_text().split()
    st = fetch_sifts_xml(ids, a.workers)
    print({k: v for k, v in st.items() if v not in ("ok", "cached")} or "all ok", len(st))
    if a.parse:
        df = parse_many([k for k, v in st.items() if v in ("ok", "cached")])
        df.to_parquet(a.parse, index=False)
        print(df.shape)
