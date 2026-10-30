#!/usr/bin/env python
"""Smoke/timing test of pocketcompare on the tool-test inputs (scripts/fetch_tool_test_inputs.py first).

1. describe (all merged pockets, --af-mode) the AFDB monomer of P31946 and its predicted homodimer and a
   heterodimer; time it.
2. Reproducibility vs the AFDB pocket atlas: the atlas holds fingerprints of AF-<ATLAS_CHECK_ACC>-F1-model_v6;
   compare our freshly computed fingerprints with the stored float16 ones (same pocket_id).
3. compare monomer vs homodimer pockets (similarity matrix); time per pair.
4. all-vs-atlas search: build an exact FAISS index over the human tier-1 atlas pockets (taxid HUMAN_TAXID)
   and search the dimer pockets; time per query.
Outputs: results/tool_tests/pocketcompare/
"""
from __future__ import annotations
import os

import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds

from pocketcompare import fpindex, method  # noqa: F401  (method.similarity used below)

WS = Path(__file__).resolve().parents[1]
IN = WS / "data/interim/tool_tests"
OUT = WS / "results/tool_tests/pocketcompare"
ATLAS = Path(os.environ.get("ATLAS_ROOT", "/cluster/scratch/gmagnusson/pocketcompare_afdb"))
ATLAS_CHECK_ACC = "P31946"
HUMAN_TAXID = 9606
# atlas filters (from the atlas parquet schema metadata 'pocketcompare_afdb'.filters)
ATLAS_MIN_SCORE, ATLAS_MIN_PLDDT = 0.03, 70.0
K = 10

INPUTS = {
    "monomer": IN / f"afdb/AF-{ATLAS_CHECK_ACC}-F1-model_v6.cif",
    "homodimer": IN / "dimers/AF-0000000065760335.pdb",
    "heterodimer": IN / "dimers/AF-0000000203452075.pdb",
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rep: dict = {}
    V, info = {}, {}
    for name, path in INPUTS.items():
        # CLI, as documented in the atlas POCKETOME.md ("pocketcompare describe <model> --af-mode")
        out = OUT / f"{name}.parquet"
        t = time.perf_counter()
        subprocess.run(["pocketcompare", "describe", str(path), "--af-mode", "-o", str(out)], check=True,
                       stdout=subprocess.DEVNULL)
        rep[f"describe_cli_{name}_wall_s"] = time.perf_counter() - t      # includes start-up / imports
        meta, V[name], rec = fpindex.read_table(out)
        info[name] = meta
        st = pd.read_parquet(OUT / f"{name}.structures.parquet")
        rep[f"describe_{name}_structure_timing"] = st.drop(columns=[c for c in st.columns if st[c].dtype == object
                                                                     and c not in ("status",)]).iloc[0].to_dict()
        meta.to_csv(OUT / f"{name}_pockets.tsv", sep="\t", index=False)
        rep[f"n_pockets_{name}"] = int(len(meta))
    rep["method_record"] = rec
    rep["describe_columns"] = list(info["monomer"].columns)

    # 2. atlas reproducibility
    dset = ds.dataset(ATLAS / "pockets", format="parquet")
    at = dset.to_table(columns=["pocket_id", "burrow_score", "plddt_mean", "residues", "fingerprint"],
                       filter=pc.field("accession") == ATLAS_CHECK_ACC).to_pandas()
    A = np.stack(at.fingerprint.to_numpy()).astype(np.float32)
    mine = info["monomer"].assign(pocket_id=lambda d: d.pocket_id.astype(str))
    m = at.merge(mine.reset_index(names="row"), on="pocket_id", suffixes=("_atlas", ""))
    diffs = np.abs(A[m.index] - V["monomer"][m.row.to_numpy()]).max(axis=1) if len(m) else np.array([])
    rep["atlas_check"] = {"atlas_pockets": int(len(at)), "matched_by_pocket_id": int(len(m)),
                          "max_abs_diff": float(diffs.max()) if len(diffs) else None,
                          "residues_identical": bool((m.residues_atlas == m.residues).all()) if len(m) else None,
                          "atlas_ids": at.pocket_id.tolist(), "our_ids": mine.pocket_id.tolist()}

    # 3. compare monomer vs homodimer
    t = time.perf_counter()
    for _ in range(1000):
        S = method.similarity(V["monomer"], V["homodimer"])
    rep["compare_matrix_us_per_pair"] = (time.perf_counter() - t) / 1000 / S.size * 1e6
    pd.DataFrame(S, index=info["monomer"].pocket_id, columns=info["homodimer"].pocket_id).to_csv(OUT / "S_monomer_vs_homodimer.tsv", sep="\t")

    # 4. search against the human atlas (exact inner product)
    import faiss
    faiss.omp_set_num_threads(1)
    t = time.perf_counter()
    ht = dset.to_table(columns=["pocket_id", "fingerprint"], filter=pc.field("taxid") == HUMAN_TAXID)
    H = np.asarray(ht["fingerprint"].combine_chunks().flatten()).reshape(-1, 384).astype(np.float32)
    rep["atlas_human_load_s"] = time.perf_counter() - t
    rep["atlas_human_pockets"] = int(len(H))
    index = faiss.IndexFlatIP(384); index.add(H)
    Q = np.concatenate([V["homodimer"], V["heterodimer"]]).astype(np.float32)
    qid = list(info["homodimer"].pocket_id) + list(info["heterodimer"].pocket_id)
    t = time.perf_counter(); D, I = index.search(Q, K)
    rep["search_flat_ms_per_query_1thread"] = (time.perf_counter() - t) / len(Q) * 1e3
    ids = ht["pocket_id"].to_numpy(zero_copy_only=False)
    hits = pd.DataFrame([(qid[i], r + 1, ids[I[i, r]], D[i, r]) for i in range(len(Q)) for r in range(K)],
                        columns=["query", "rank", "hit", "score"])
    hits.to_csv(OUT / "dimer_vs_human_atlas_top10.tsv", sep="\t", index=False)
    (OUT / "report.json").write_text(json.dumps(rep, indent=1, default=str))
    print(json.dumps({k: v for k, v in rep.items() if k != "describe_columns"}, indent=1, default=str))


if __name__ == "__main__":
    main()
