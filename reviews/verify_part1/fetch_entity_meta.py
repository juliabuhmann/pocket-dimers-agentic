"""Read organism (scientific name + NCBI taxid) and entity description of the dominant entity for the
representative entries of the H3 units, from the local mmCIF mirror (_entity_src_gen / _entity_src_nat /
_pdbx_entity_src_syn, _entity.pdbx_description). Verifier script; first attempt.
"""
import gzip
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gemmi
import pandas as pd

MIRROR = Path("/cluster/scratch/jjaenes/pdb")  # local mmCIF mirror named in docs/PART1_PIPELINE.md
OUT = Path(__file__).resolve().parent / "out"
N_WORKERS = 6  # task budget: <=6 workers

SRC = [("_entity_src_gen.", "pdbx_gene_src_scientific_name", "pdbx_gene_src_ncbi_taxonomy_id"),
       ("_entity_src_nat.", "pdbx_organism_scientific", "pdbx_ncbi_taxonomy_id"),
       ("_pdbx_entity_src_syn.", "organism_scientific", "ncbi_taxonomy_id")]


def one(pdb_id):
    f = MIRROR / pdb_id.lower()[1:3] / f"{pdb_id.lower()}.cif.gz"
    out = {}
    try:
        b = gemmi.cif.read(str(f)).sole_block()
        desc = {}
        t = b.find("_entity.", ["id", "pdbx_description"])
        for r in t:
            desc[r[0]] = gemmi.cif.as_string(r[1])
        for pre, nm, tx in SRC:
            t = b.find(pre, ["entity_id", nm, tx])
            for r in t:
                out.setdefault(r[0], (gemmi.cif.as_string(r[1]), gemmi.cif.as_string(r[2])))
        return [(pdb_id, e, desc.get(e), *out.get(e, (None, None))) for e in desc]
    except Exception as ex:  # noqa
        return [(pdb_id, None, f"ERROR {ex}", None, None)]


def main():
    u = pd.read_parquet(OUT / "units_verifier.parquet", columns=["pdb_id", "assembly_type", "mbs"])
    ids = sorted(u[(u.assembly_type == "homomer") & u.mbs].pdb_id.unique())
    with ProcessPoolExecutor(N_WORKERS) as ex:
        res = [r for rr in ex.map(one, ids, chunksize=20) for r in rr]
    df = pd.DataFrame(res, columns=["pdb_id", "entity", "description", "organism", "taxid"])
    df.to_parquet(OUT / "entity_meta.parquet")
    print(len(ids), len(df), df.organism.isna().sum())


if __name__ == "__main__":
    sys.exit(main())
