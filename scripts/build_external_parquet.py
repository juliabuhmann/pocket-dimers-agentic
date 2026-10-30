#!/usr/bin/env python
"""Convert raw external files (data/raw/...) to tidy parquet in data/processed/.

Outputs: variants_humsavar, alphamissense_residue_mean, sifts_chain_uniprot, sifts_chain_pfam (with clan),
pfam_clans, ligysis_human_sites. (ccd_to_chembl.parquet is written by fetch_chembl.py.)
Usage: python scripts/build_external_parquet.py [--only name ...]
"""
import argparse
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pacsv

sys.path.insert(0, str(Path(__file__).parent))
from extlib import PROCESSED, RAW  # noqa: E402

AA3 = dict(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E", Gly="G", His="H", Ile="I", Leu="L",
           Lys="K", Met="M", Phe="F", Pro="P", Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V")
CATEGORY_MAP = {"LP/P": "P/LP", "LB/B": "B/LB", "US": "VUS"}  # humsavar codes -> readable (raw code kept too)


def humsavar():
    rows = []
    started = False
    for line in open(RAW / "humsavar" / "humsavar.txt", errors="replace"):
        if line.startswith("_____"):
            started = True
            continue
        if not started:
            continue
        if not line.strip() or line.startswith("---") or line.startswith("Copyrighted"):
            break
        p = line.rstrip("\n").split(None, 6)
        if len(p) < 6:
            continue
        gene, acc, ftid, change, cat, dbsnp = p[:6]
        disease = p[6].strip() if len(p) > 6 else "-"
        m = re.fullmatch(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})", change)
        if not m:  # non-missense (e.g. Ter, del) are skipped
            continue
        wt, pos, mut = m.groups()
        if wt not in AA3 or mut not in AA3:
            continue
        rows.append((acc, gene, ftid, int(pos), AA3[wt], AA3[mut], cat, CATEGORY_MAP.get(cat, cat),
                     None if disease == "-" else disease, None if dbsnp == "-" else dbsnp))
    df = pd.DataFrame(rows, columns=["accession", "gene", "ftid", "position", "wt", "mut", "category_code",
                                     "category", "disease", "dbSNP"])
    df.to_parquet(PROCESSED / "variants_humsavar.parquet", index=False)
    print("humsavar", df.shape, df.category_code.value_counts().to_dict())


def alphamissense():
    """Stream the 216M-row file in Arrow batches; aggregate sum/count/max/min per (uniprot, position)."""
    src = RAW / "alphamissense" / "AlphaMissense_aa_substitutions.tsv.gz"
    fh = gzip.open(src, "rb")
    # header: skip '#' comment lines
    pos0 = 0
    while True:
        line = fh.readline()
        if not line.startswith(b"#"):
            break
    header = line.decode().rstrip("\n").split("\t")
    reader = pacsv.open_csv(fh, read_options=pacsv.ReadOptions(column_names=header, block_size=1 << 27),
                            parse_options=pacsv.ParseOptions(delimiter="\t"),
                            convert_options=pacsv.ConvertOptions(include_columns=["uniprot_id", "protein_variant",
                                                                                  "am_pathogenicity"]))
    parts = []
    n = 0
    for batch in reader:
        d = batch.to_pandas()
        d["position"] = d["protein_variant"].str.slice(1, -1).astype(np.int32)
        g = d.groupby(["uniprot_id", "position"], sort=False)["am_pathogenicity"].agg(["sum", "count", "max", "min"])
        parts.append(g.reset_index())
        n += len(d)
        if len(parts) % 10 == 0:
            parts = [pd.concat(parts).groupby(["uniprot_id", "position"], sort=False)
                     .agg({"sum": "sum", "count": "sum", "max": "max", "min": "min"}).reset_index()]
            print("rows", n, flush=True)
    g = pd.concat(parts).groupby(["uniprot_id", "position"], sort=False).agg(
        {"sum": "sum", "count": "sum", "max": "max", "min": "min"}).reset_index()
    out = pd.DataFrame({"accession": g.uniprot_id, "position": g.position,
                        "am_mean": g["sum"] / g["count"], "am_max": g["max"], "am_min": g["min"],
                        "n_subst": g["count"].astype("int16")}).sort_values(["accession", "position"])
    out.to_parquet(PROCESSED / "alphamissense_residue_mean.parquet", index=False)
    print("alphamissense", out.shape, "rows read", n, "proteins", out.accession.nunique())


def sifts():
    u = pd.read_csv(RAW / "sifts_uniprot" / "pdb_chain_uniprot.tsv.gz", sep="\t", skiprows=1, dtype=str)
    for c in ["RES_BEG", "RES_END", "SP_BEG", "SP_END"]:
        u[c] = pd.to_numeric(u[c], errors="coerce").astype("Int64")
    u.columns = [c.lower() for c in u.columns]
    u.to_parquet(PROCESSED / "sifts_chain_uniprot.parquet", index=False)
    pf = pd.read_csv(RAW / "sifts_pfam" / "pdb_chain_pfam.tsv.gz", sep="\t", skiprows=1, dtype=str)
    pf.columns = [c.lower() for c in pf.columns]
    cl = pd.read_csv(RAW / "pfam_clans" / "Pfam-A.clans.tsv.gz", sep="\t", header=None,
                     names=["pfam_id", "clan_id", "clan_name", "pfam_name", "pfam_description"], dtype=str)
    cl.to_parquet(PROCESSED / "pfam_clans.parquet", index=False)
    pf = pf.merge(cl[["pfam_id", "clan_id", "clan_name", "pfam_name"]], on="pfam_id", how="left")
    pf.to_parquet(PROCESSED / "sifts_chain_pfam.parquet", index=False)
    print("sifts uniprot", u.shape, "pfam", pf.shape, "pfam with clan", pf.clan_id.notna().mean().round(3),
          "clans", cl.shape)


def ligysis():
    s = pd.read_pickle(RAW / "ligysis" / "LIGYSIS_human_sites_APRIL_2024.pkl")
    s = s.rename(columns={"ACC": "accession", "SEG": "segment", "ID": "site_id", "SEG_ID": "seg_id",
                          "SITE_NAME": "site_name", "n_up_aas": "n_residues"})
    ion = pd.read_pickle(RAW / "ligysis" / "LIGYSIS_human_sites_ion_prop_APRIL_2024.pkl")  # (acc, seg, site)->frac
    s["ion_proportion"] = [ion.get((a, int(g), int(i))) for a, g, i in zip(s.accession, s.segment, s.site_id)]
    ex = s.explode("up_aas").rename(columns={"up_aas": "position"})
    ex["position"] = ex["position"].astype("Int64")
    ex.drop(columns=["n_residues"]).to_parquet(PROCESSED / "ligysis_human_site_residues.parquet", index=False)
    s.to_parquet(PROCESSED / "ligysis_human_sites.parquet", index=False)
    print("ligysis sites", s.shape, "residue rows", ex.shape, "proteins", s.accession.nunique())


STEPS = dict(humsavar=humsavar, alphamissense=alphamissense, sifts=sifts, ligysis=ligysis)
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=list(STEPS))
    PROCESSED.mkdir(parents=True, exist_ok=True)
    for n in ap.parse_args().only:
        STEPS[n]()
