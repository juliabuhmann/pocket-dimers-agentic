#!/usr/bin/env python
"""Fetch every external list used by the Part 1 ligand-site census into data/raw/part1/.

For each source the URL, UTC fetch date, byte size and sha256 are written to
data/raw/part1/MANIFEST.tsv and a row is appended to DATA.md (section "Part 1 sources").
Re-running re-downloads (use --skip-existing to keep files already present).

Sources (URLs are the single place they are defined; see src/part1/config.py for paths):
  biolip_ligand_list   BioLiP artifact/additive list (Zhang lab)
  plinder_badlist      PLINDER artifacts_badlist.csv (plinder-org/plinder, main branch)
  ccd_components       wwPDB Chemical Component Dictionary (components.cif.gz)
  pdbe_cofactors       PDBe compound/cofactors API (CCD codes per cofactor class, derived from
                       the EMBL-EBI CoFactor database, Thornton group)
  sifts_chain_uniprot  SIFTS pdb_chain_uniprot.tsv.gz (EBI)
  rcsb_clusters_30/95  RCSB entity sequence clusters (MMseqs2, 30% / 95% identity)

Usage: python scripts/fetch_part1_sources.py [--skip-existing] [--only name ...]
"""
import argparse
import datetime as dt
import hashlib
import sys
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from part1.config import RAW_DIR, SOURCES  # noqa: E402

DATA_MD = ROOT / "DATA.md"
MANIFEST = RAW_DIR / "MANIFEST.tsv"
DATA_MD_HEADER = (
    "\n## Part 1 sources (fetched by `scripts/fetch_part1_sources.py`)\n\n"
    "| Source | URL | Fetched (UTC) | Local file | Bytes | sha256 |\n|---|---|---|---|---|---|\n"
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_chembl_paginated(url: str, dest: Path) -> None:
    """Follow ChEMBL REST page_meta.next links; write molecule_chembl_id<TAB>max_phase."""
    base = "https://www.ebi.ac.uk"
    tmp = dest.with_suffix(dest.suffix + ".part")
    with open(tmp, "w") as f:
        f.write("molecule_chembl_id\tmax_phase\n")
        while url:
            r = requests.get(url, timeout=120)
            r.raise_for_status()
            js = r.json()
            for m in js["molecules"]:
                f.write(f"{m['molecule_chembl_id']}\t{m['max_phase']}\n")
            nxt = js["page_meta"]["next"]
            url = base + nxt if nxt else None
    tmp.rename(dest)


def fetch(name: str, url: str, dest: Path) -> None:
    if name == "chembl_max_phase":
        return fetch_chembl_paginated(url, dest)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    tmp.rename(dest)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-existing", action="store_true")
    ap.add_argument("--only", nargs="*")
    args = ap.parse_args()
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, (url, fname) in SOURCES.items():
        if args.only and name not in args.only:
            continue
        dest = RAW_DIR / fname
        if args.skip_existing and dest.exists():
            print(f"[skip] {name}: {dest}")
            continue
        print(f"[fetch] {name}: {url}", flush=True)
        fetch(name, url, dest)
        when = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        rows.append((name, url, when, str(dest.relative_to(ROOT)), dest.stat().st_size, sha256(dest)))
    if not rows:
        return
    new_manifest = not MANIFEST.exists()
    with open(MANIFEST, "a") as f:
        if new_manifest:
            f.write("source\turl\tfetched_utc\tpath\tbytes\tsha256\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    text = DATA_MD.read_text()
    with open(DATA_MD, "a") as f:
        if "## Part 1 sources" not in text:
            f.write(DATA_MD_HEADER)
        for r in rows:
            f.write(f"| {r[0]} | {r[1]} | {r[2]} | `{r[3]}` | {r[4]} | `{r[5]}` |\n")
    for r in rows:
        print("\t".join(map(str, r)))


if __name__ == "__main__":
    main()
