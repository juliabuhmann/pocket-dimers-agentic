#!/usr/bin/env python
"""Extract human (taxid 9606) pockets from the monomer atlas snapshot (tiers 1+2) -> interim parquet + fp16 npy.

Status: first attempt, not yet verified.
"""
import glob, json, sys
import numpy as np, pandas as pd, pyarrow.parquet as pq, pyarrow.compute as pc
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from q22_common import ATLAS, HUMAN_TAXID, INTERIM

COLS = ["accession", "fragment", "tier", "pocket_id", "burrow_score", "volume_A3", "site_volume_A3",
        "residues", "plddt_mean", "fingerprint"]
metas, fps, methods = [], [], set()
files = sorted(glob.glob(str(ATLAS / "tier[12]" / "*.parquet")))
for f in files:
    pf = pq.ParquetFile(f)
    md = pf.schema_arrow.metadata or {}
    methods.add(json.loads(md[b"pocketcompare"])["method"] if b"pocketcompare" in md else "NA")
    t = pq.read_table(f, columns=COLS + ["taxid"])
    t = t.filter(pc.equal(t["taxid"], HUMAN_TAXID))
    if t.num_rows == 0:
        continue
    fp = np.asarray(t.column("fingerprint").combine_chunks().flatten(), dtype=np.float16).reshape(-1, 384)
    fps.append(fp)
    metas.append(t.drop(["fingerprint", "taxid"]).to_pandas())
meta = pd.concat(metas, ignore_index=True)
fp = np.concatenate(fps)
meta.to_parquet(INTERIM / "atlas_human_meta.parquet")
np.save(INTERIM / "atlas_human_fp16.npy", fp)
print(f"files={len(files)} human_pockets={len(meta)} accessions={meta.accession.nunique()} tiers={meta.tier.value_counts().to_dict()} methods={methods}")
