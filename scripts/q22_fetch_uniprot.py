#!/usr/bin/env python
"""Fetch UniProt sequences + gene names for all human accessions in the dimer set and the human atlas.

Sources (UniProt REST, release recorded from the X-UniProt-Release header):
  - reference proteome stream  https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000005640&format=fasta
    (+ fields accession,gene_primary,protein_name,length in TSV)
  - any accession not in the proteome stream: /uniprotkb/accessions?accessions=...  (batches of 100)
Outputs: data/interim/q22/human_seqs.fasta, data/interim/q22/uniprot_genes.tsv, data/interim/q22/uniprot_fetch.json
Status: first attempt, not yet verified.
"""
import datetime, hashlib, io, json, sys
import pandas as pd, requests
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from q22_common import ANNOT, INTERIM, HUMAN_PROTEOME

BASE = "https://rest.uniprot.org/uniprotkb"
FIELDS = "accession,gene_primary,protein_name,length,reviewed"


def get(url, **params):
    r = requests.get(url, params=params, timeout=300)
    r.raise_for_status()
    return r


d = pd.read_parquet(ANNOT, columns=["uniprot_A", "uniprot_B"])
atlas = pd.read_parquet(INTERIM / "atlas_human_meta.parquet", columns=["accession"])
need = set(d.uniprot_A) | set(d.uniprot_B) | set(atlas.accession)

q = f"proteome:{HUMAN_PROTEOME}"
r_fa = get(f"{BASE}/stream", query=q, format="fasta")
r_tsv = get(f"{BASE}/stream", query=q, format="tsv", fields=FIELDS)
release = r_fa.headers.get("X-UniProt-Release", "NA")
fasta = r_fa.text
tsv = pd.read_csv(io.StringIO(r_tsv.text), sep="\t")
have = set(tsv.Entry)
missing = sorted(need - have)
extra_fa, extra_tsv = [], []
for i in range(0, len(missing), 100):
    chunk = ",".join(missing[i:i + 100])
    extra_fa.append(get(f"{BASE}/accessions", accessions=chunk, format="fasta").text)
    extra_tsv.append(pd.read_csv(io.StringIO(get(f"{BASE}/accessions", accessions=chunk, format="tsv", fields=FIELDS).text), sep="\t"))
fasta_all = fasta + "".join(extra_fa)
genes = pd.concat([tsv] + extra_tsv, ignore_index=True).drop_duplicates("Entry")
# FASTA headers -> plain accession ids for mmseqs
lines = []
for ln in fasta_all.splitlines():
    if ln.startswith(">"):
        ln = ">" + ln[1:].split("|")[1]
    lines.append(ln)
(INTERIM / "human_seqs.fasta").write_text("\n".join(lines) + "\n")
genes.to_csv(INTERIM / "uniprot_genes.tsv", sep="\t", index=False)
still = sorted(need - set(genes.Entry))
info = {"retrieved_utc": datetime.datetime.utcnow().isoformat(), "uniprot_release": release, "query": q,
        "n_needed": len(need), "n_proteome": len(have), "n_extra_fetched": len(missing), "n_still_missing": len(still),
        "still_missing_examples": still[:20], "fasta_sha256": hashlib.sha256("\n".join(lines).encode()).hexdigest()}
(INTERIM / "uniprot_fetch.json").write_text(json.dumps(info, indent=1))
print(info)
