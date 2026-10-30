#!/usr/bin/env python
"""Fetch external datasets into data/raw/<source>/ and log sha256/version/license to DATA.md.

Sources (URLs are defined only here):
  humsavar        UniProt human variants
  alphamissense   AlphaMissense aa_substitutions (1.2 GB)
  unichem         UniChem PDBe(1)->ChEMBL(3)... (src1src3)
  sifts_uniprot / sifts_pfam   SIFTS chain-level flat files
  pfam_clans      Pfam-A.clans.tsv.gz
  ligysis         small LIGYSIS Zenodo files (record 14645504; md5 verified from Zenodo API)
  plinder         attempt anonymous GCS fetch (expected to fail; see docs/EXTERNAL_DATA.md)
Usage: python scripts/fetch_external.py [--only name ...] [--force]
"""
import argparse
import hashlib
import re
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).parent))
from extlib import RAW, UA, append_data_md, download, now, sha256  # noqa: E402

BASE = "https://ftp.ebi.ac.uk/pub/databases"
SIMPLE = {
    # name: (url, filename, license)
    "humsavar": ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt",
                 "humsavar.txt", "CC BY 4.0"),
    "alphamissense": ("https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz",
                      "AlphaMissense_aa_substitutions.tsv.gz", "CC BY-NC-SA 4.0 (per file header; Cheng et al. 2023, Zenodo 8208688)"),
    "unichem": (f"{BASE}/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src3.txt.gz",
                "src1src3.txt.gz", "CC0 (UniChem)"),
    "sifts_uniprot": (f"{BASE}/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz",
                      "pdb_chain_uniprot.tsv.gz", "CC0 (PDBe/SIFTS)"),
    "sifts_pfam": (f"{BASE}/msd/sifts/flatfiles/tsv/pdb_chain_pfam.tsv.gz",
                   "pdb_chain_pfam.tsv.gz", "CC0 (PDBe/SIFTS)"),
    "pfam_clans": (f"{BASE}/Pfam/current_release/Pfam-A.clans.tsv.gz",
                   "Pfam-A.clans.tsv.gz", "CC0 (Pfam)"),
    # ChEBI flat files: metabolite-role exclusion for the drug-like class (scripts/build_druglike_flags.py)
    "chebi_relation": (f"{BASE}/chebi/flat_files/relation.tsv.gz", "relation.tsv.gz", "CC BY 4.0 (ChEBI)"),
    "chebi_compounds": (f"{BASE}/chebi/flat_files/compounds.tsv.gz", "compounds.tsv.gz", "CC BY 4.0 (ChEBI)"),
    "chebi_dbacc": (f"{BASE}/chebi/flat_files/database_accession.tsv.gz", "database_accession.tsv.gz",
                    "CC BY 4.0 (ChEBI)"),
    "chebi_source": (f"{BASE}/chebi/flat_files/source.tsv.gz", "source.tsv.gz", "CC BY 4.0 (ChEBI)"),
    "chebi_relation_type": (f"{BASE}/chebi/flat_files/relation_type.tsv.gz", "relation_type.tsv.gz",
                            "CC BY 4.0 (ChEBI)"),
}
ZENODO_RECORD = "14645504"
LIGYSIS_FILES = [  # small files only (<500 MB total); selected from the Zenodo file list
    "LIGYSIS_human_sites_APRIL_2024.pkl",
    "LIGYSIS_human_sites_ion_prop_APRIL_2024.pkl",
    "LIGYSIS_human_chains_per_lig_MAY_2024.pkl",
    "LIGYSIS_sites_DEF_TRANS.pkl",
    "LIGYSIS_3448_chains.pkl",
    "LIGYSIS_ress_binary_labels.pkl",
]
PLINDER_URL = "https://storage.googleapis.com/plinder/2024-06/v2/index/annotation_table.parquet"
PLINDER_MAX_BYTES = 3 * 10**9


def version_of(name, path, headers):
    if name == "humsavar":
        m = re.search(r"Release:\s+(.+)", Path(path).read_text(errors="replace")[:2000])
        return "UniProt " + (m.group(1).strip() if m else "?")
    return "Last-Modified: " + headers.get("Last-Modified", headers.get("last-modified", "n/a"))


def fetch_simple(name, force):
    url, fn, lic = SIMPLE[name]
    dest = RAW / name / fn
    headers = {}
    if force or not dest.exists():
        headers = download(url, dest)
    else:
        r = requests.head(url, headers=UA, allow_redirects=True, timeout=60)
        headers = dict(r.headers)
    append_data_md(name, version_of(name, dest, headers), sha256(dest), lic, url,
                   str(dest.relative_to(RAW.parent.parent)), dest.stat().st_size)
    print(name, dest.stat().st_size)


def fetch_ligysis(force):
    meta = requests.get(f"https://zenodo.org/api/records/{ZENODO_RECORD}", headers=UA, timeout=60).json()
    files = {f["key"]: f for f in meta["files"]}
    version = f"Zenodo {ZENODO_RECORD} (published {meta.get('metadata', {}).get('publication_date')})"
    lic = meta.get("metadata", {}).get("license", {}).get("id", "see Zenodo record")
    for key in LIGYSIS_FILES:
        f = files[key]
        dest = RAW / "ligysis" / key
        url = f["links"]["self"]
        if force or not dest.exists():
            download(url, dest)
        md5 = hashlib.md5(dest.read_bytes()).hexdigest()
        assert f"md5:{md5}" == f["checksum"], f"md5 mismatch for {key}"
        append_data_md(f"ligysis/{key}", version, sha256(dest), lic, url,
                       str(dest.relative_to(RAW.parent.parent)), dest.stat().st_size)
        print("ligysis", key, dest.stat().st_size)


def fetch_plinder():
    """Anonymous access test (HTTPS JSON API honours the cluster proxy; gcsfs/aiohttp does not)."""
    r = requests.head(PLINDER_URL, headers=UA, timeout=60)
    print("plinder HEAD", r.status_code)
    if r.status_code != 200:
        print("PLINDER NOT ACCESSIBLE anonymously:", r.status_code)
        return False
    size = int(r.headers.get("Content-Length", 0))
    if size > PLINDER_MAX_BYTES:
        print("PLINDER too large", size)
        return False
    dest = RAW / "plinder" / "annotation_table.parquet"
    h = download(PLINDER_URL, dest)
    append_data_md("plinder/annotation_table", "2024-06/v2", sha256(dest), "CC BY 4.0 (PLINDER)",
                   PLINDER_URL, str(dest.relative_to(RAW.parent.parent)), dest.stat().st_size)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    names = list(SIMPLE) + ["ligysis", "plinder"]
    for n in (a.only or names):
        print("==", n, now(), flush=True)
        try:
            if n in SIMPLE:
                fetch_simple(n, a.force)
            elif n == "ligysis":
                fetch_ligysis(a.force)
            elif n == "plinder":
                fetch_plinder()
        except Exception as e:  # noqa: BLE001
            print("FAILED", n, repr(e))


if __name__ == "__main__":
    main()
