#!/usr/bin/env python
"""Summary TSV + sanity checks for the annotated human dimer pockets (first attempt, not verified)."""
import numpy as np, pandas as pd, gemmi
from pathlib import Path
WS = Path(__file__).resolve().parents[1]
df = pd.read_parquet(WS / "data/part2/human_dimer_pockets_annotated.parquet")
rows = []
for kind, g in df.groupby("kind"):
    n_struct = g.structure_stem.nunique()
    a = g[g.atlas_pass]
    rows.append(dict(set=kind, n_structures_with_pockets=n_struct, n_pockets=len(g), pockets_per_structure=len(g) / n_struct,
        n_atlas_pass=len(a), frac_interface_ge3_all=g.is_interface.mean(), frac_interface_ge1_all=g.is_interface_ge1.mean(),
        frac_interface_ge3_atlaspass=a.is_interface.mean(), frac_dimer_only_of_interface=g[g.is_interface].dimer_only.astype(float).mean(),
        n_sym_groups=(g.groupby("structure_stem").sym_group.nunique().sum() if kind == "homo" else np.nan)))
out = pd.DataFrame(rows); (WS / "results/part2").mkdir(parents=True, exist_ok=True)
out.to_csv(WS / "results/part2/pocket_run_summary.tsv", sep="\t", index=False)
print(out.T.to_string())
# geometric check: 5 random interface pockets, lining CA within 12 A of centre on both chains
iface = df[df.is_interface].sample(5, random_state=1)
for r in iface.itertuples():
    st = gemmi.read_structure(f"/cluster/scratch/gmagnusson/dimers/afdb_human_dimers/{'homodimers' if r.kind=='homo' else 'heterodimers'}/{r.structure_stem}.pdb")
    c = gemmi.Position(r.center_x, r.center_y, r.center_z); d = {"A": [], "B": []}
    for x in r.residues.split():
        res = st[0][x[0]][x[2:]][0]
        d[x[0]].append(res["CA"][0].pos.dist(c))
    print(r.pocket_id, {k: (len(v), round(min(v), 1), round(max(v), 1)) for k, v in d.items()})
