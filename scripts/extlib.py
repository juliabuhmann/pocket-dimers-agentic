"""Shared helpers for external-data fetch scripts (download, sha256, DATA.md rows)."""
import datetime as dt
import hashlib
import os
import re
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
DATA_MD = ROOT / "DATA.md"
SECTION = "## External data (scripts/fetch_external.py, fetch_chembl.py, fetch_sifts_xml.py)"
HEADER = ("| Source | Version | sha256 | License | URL | Retrieved (UTC) | Local file | Bytes |\n"
          "|---|---|---|---|---|---|---|---|\n")
# Optional contact for the HTTP User-Agent (set CONTACT_EMAIL in the environment; empty by default).
_CONTACT = os.environ.get("CONTACT_EMAIL", "")
UA = {"User-Agent": "pocket-dimers-study/0.1" + (f" ({_CONTACT})" if _CONTACT else "")}


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def download(url, dest, retries=4):
    """Stream url to dest (via .part, atomic rename). Returns response headers."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    for i in range(retries):
        try:
            with requests.get(url, stream=True, timeout=120, headers=UA) as r:
                r.raise_for_status()
                tmp = dest.with_suffix(dest.suffix + ".part")
                with open(tmp, "wb") as f:
                    for c in r.iter_content(1 << 20):
                        f.write(c)
                tmp.replace(dest)
                return dict(r.headers)
        except Exception as e:  # noqa: BLE001
            if i == retries - 1:
                raise
            import time
            time.sleep(5 * (i + 1))


def append_data_md(source, version, sha, license_, url, local, nbytes):
    """Insert a row at the end of the External-data table in DATA.md (idempotent per sha256)."""
    txt = DATA_MD.read_text() if DATA_MD.exists() else ""
    if f"| {source} |" in txt and sha in txt:
        return
    row = (f"| {source} | {version} | `{sha}` | {license_} | {url} | {now()} | "
           f"`{local}` | {nbytes} |")
    if SECTION not in txt:
        txt = txt.rstrip("\n") + f"\n\n{SECTION}\n\n{HEADER}"
    lines = txt.split("\n")
    i = lines.index(SECTION) + 1
    while i < len(lines) and not lines[i].startswith("|"):
        i += 1
    while i < len(lines) and lines[i].startswith("|"):
        i += 1
    lines.insert(i, row)
    DATA_MD.write_text("\n".join(lines))
