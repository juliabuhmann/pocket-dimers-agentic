#!/usr/bin/env python
"""H4 / Q2.1 step 1-2 (part 1): benchmark selection of PDB ligand sites on human UniProt pairs with an AFDB dimer model.

Intended behaviour (first attempt, NOT verified):
  A. Part 1 census primary-set sites (additive-filtered, bound, non-duplicate, non-glycan) in homomer/heteromer
     assemblies. MBS: the two chains with most binding residues must map to a dimer model's UniProt pair
     (homo: same accession; hetero: unordered pair). SBS: the binding chain's protein must be in a dimer model;
     the partner chain (same accession for homo, the model partner for hetero) is resolved from the assembly.
  B. Non-redundant unit = (model, cluster30, ccd_key, MBS label). Up to MAX_ENTRIES_PER_UNIT best-resolution entries
     per unit are processed (fallbacks); per entry: SIFTS XML per-residue mapping (fetched only for these entries),
     assembly-1 chains rebuilt with the census chain naming, partner choice (max interface), PDB interface residues
     (heavy atoms <= 5 Å), residue-name identity of the binding residues vs UniProt, CA coordinates in UniProt numbering.
Outputs: data/part2/h4/<census>/sites.parquet, chain_ca.parquet, iface.parquet; results/part2/h4/<census>/selection_ledger.tsv.
Usage: python scripts/h4_select.py --census v1|v2 [--workers 6] [--limit-entries N]
"""
import argparse
import json
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from h4_common import (WS, PDB_MIRROR, OUT_DATA, OUT_RES, IFACE_DIST, MAX_ENTRIES_PER_UNIT, ORGANIC_EXCLUDE,
                       METAL_CLASSES, census_paths, min_atoms_primary, site_group_col, load_candidates, write_tsv)

COLS = ["pdb_id", "ligand_uid", "assembly_type", "is_mbs", "contacts", "ccd", "ccd_key", "ligand_class", "marsh_class",
        "cluster30", "resolution", "release_date", "centroid", "n_heavy", "is_on_axis", "is_approved_drug",
        "chembl_max_phase", "method_class"]


def stage_a(version, log):
    lig_p, ent_p = census_paths(version)
    cand = load_candidates()
    homo = set(cand.loc[cand.kind == "homo", "uniprotAccession_1"])
    het = cand[cand.kind == "hetero"]
    hetpair = {frozenset((a, b)): m for a, b, m in zip(het.uniprotAccession_1, het.uniprotAccession_2, het.modelEntityId)}
    homo_model = dict(zip(cand.loc[cand.kind == "homo", "uniprotAccession_1"], cand.loc[cand.kind == "homo", "modelEntityId"]))
    het_partners = defaultdict(set)
    for p in hetpair:
        a, b = tuple(p)
        het_partners[a].add(b); het_partners[b].add(a)
    model_prot = homo | set(het_partners)

    ent = pq.read_table(ent_p, columns=["pdb_id", "uniprot_accs", "assembly_type"]).to_pandas()
    ent = ent[ent.uniprot_accs.notna()]
    ent["accs"] = ent.uniprot_accs.map(lambda x: set(x))
    ent = ent[ent.accs.map(lambda s: bool(s & model_prot))]
    entry_accs = dict(zip(ent.pdb_id, ent.accs))
    cols = COLS + ([site_group_col(version)] if site_group_col(version) not in COLS else [])
    schema = pq.read_schema(lig_p).names
    if "census_version" in schema:
        cols.append("census_version")
    t = time.time()
    df = pq.read_table(lig_p, columns=cols, filters=[("in_primary_set", "==", True),
                                                     ("assembly_type", "in", ["homomer", "heteromer"]),
                                                     ("pdb_id", "in", sorted(entry_accs))]).to_pandas()
    log.append(("primary-set homomer/heteromer sites in entries containing a model protein", len(df), df.pdb_id.nunique()))
    print(f"[A] loaded {len(df)} sites in {time.time()-t:.0f}s", flush=True)
    if "census_version" in df:
        assert set(df.census_version.unique()) == {version}, df.census_version.unique()
    df = df.rename(columns={site_group_col(version): "site_group"})
    mina = min_atoms_primary(version)

    rows = []
    n_drop = defaultdict(int)
    for r in df.itertuples(index=False):
        per = defaultdict(list)
        acc_of = {}
        for c in r.contacts:
            if c["kind"] != "protein" or c["n_vdw"] < mina:
                continue
            per[c["chain"]].append((int(c["label_seq"]) if c["label_seq"] is not None else None, c["comp"]))
            acc_of[c["chain"]] = c["unp_acc"]
        if not per:
            n_drop["no_binding_protein_residue"] += 1
            continue
        chains = sorted(per, key=lambda k: (-len(per[k]), k))
        base = dict(pdb_id=r.pdb_id, ligand_uid=r.ligand_uid, is_mbs=bool(r.is_mbs), ccd=r.ccd, ccd_key=r.ccd_key,
                    ligand_class=r.ligand_class, marsh_class=r.marsh_class, cluster30=r.cluster30,
                    resolution=r.resolution, release_date=r.release_date, centroid=list(map(float, r.centroid)),
                    n_heavy=r.n_heavy, is_on_axis=r.is_on_axis, site_group=r.site_group, method_class=r.method_class,
                    is_approved_drug=r.is_approved_drug, chembl_max_phase=r.chembl_max_phase,
                    bind=json.dumps({k: v for k, v in per.items()}), n_bind_chains=len(chains))
        base["lig_set"] = ("organic" if r.ligand_class not in ORGANIC_EXCLUDE else
                           "metal" if r.ligand_class in METAL_CLASSES else "inorganic")
        c1 = chains[0]
        a1 = acc_of[c1]
        if r.is_mbs:
            if len(chains) < 2:
                n_drop["mbs_lt2_chains_primary_rule"] += 1
                continue
            c2 = chains[1]
            a2 = acc_of[c2]
            if a1 is None or a2 is None:
                n_drop["mbs_chain_without_uniprot"] += 1
                continue
            if a1 == a2 and a1 in homo:
                rows.append(dict(base, kind="homo", model=homo_model[a1], chain1=c1, chain2=c2, acc1=a1, acc2=a2))
            elif a1 != a2 and frozenset((a1, a2)) in hetpair:
                rows.append(dict(base, kind="hetero", model=hetpair[frozenset((a1, a2))], chain1=c1, chain2=c2, acc1=a1, acc2=a2))
            else:
                n_drop["mbs_pair_not_modelled"] += 1
        else:
            if a1 is None:
                n_drop["sbs_chain_without_uniprot"] += 1
                continue
            added = 0
            accs = entry_accs.get(r.pdb_id, set())
            if a1 in homo:
                rows.append(dict(base, kind="homo", model=homo_model[a1], chain1=c1, chain2=None, acc1=a1, acc2=a1)); added += 1
            for b in sorted(het_partners.get(a1, set()) & accs):
                if b == a1:
                    continue
                rows.append(dict(base, kind="hetero", model=hetpair[frozenset((a1, b))], chain1=c1, chain2=None, acc1=a1, acc2=b)); added += 1
            if not added:
                n_drop["sbs_no_model_pair_in_entry"] += 1
    for k, v in n_drop.items():
        log.append((f"dropped: {k}", v, None))
    sm = pd.DataFrame(rows)
    log.append(("site x model candidate rows", len(sm), sm.pdb_id.nunique()))
    return sm


def choose_entries(sm, log):
    """Up to MAX_ENTRIES_PER_UNIT best-resolution entries per non-redundant unit."""
    sm["unit"] = sm.model + "|" + sm.cluster30.astype(str) + "|" + sm.ccd_key + "|" + np.where(sm.is_mbs, "MBS", "SBS")
    e = sm[["unit", "pdb_id", "resolution"]].drop_duplicates(["unit", "pdb_id"]).sort_values(["unit", "resolution", "pdb_id"])
    e = e.groupby("unit").head(MAX_ENTRIES_PER_UNIT)
    keep = set(e.pdb_id)
    log.append(("non-redundant units (model x cluster30 x CCD x MBS label), before mapping", sm.unit.nunique(), None))
    log.append((f"entries selected (<= {MAX_ENTRIES_PER_UNIT} best-resolution entries per unit)", len(keep), None))
    return sm[sm.pdb_id.isin(keep)].copy(), sorted(keep)


# ------------------------------------------------------------------ per-entry worker
def _one_letter(comp):
    import gemmi
    t = gemmi.find_tabulated_residue(comp)
    return t.one_letter_code.upper() if t is not None and t.one_letter_code.strip() else "X"


def process_entry(args):
    pdb_id, sites = args
    try:
        return _process(pdb_id, sites)
    except Exception as ex:  # logged, never silent
        return {"pdb_id": pdb_id, "error": f"{type(ex).__name__}: {ex}", "sites": [], "ca": [], "iface": []}


def _process(pdb_id, sites):
    from scipy.spatial import cKDTree
    from src.part1 import structure as S
    from fetch_sifts_xml import parse_sifts_xml, CACHE
    sx = CACHE / f"{pdb_id}.xml.gz"
    if not sx.exists():
        return {"pdb_id": pdb_id, "error": "no_sifts_xml", "sites": [], "ca": [], "iface": []}
    sf = parse_sifts_xml(sx)
    sf = sf[sf.uniprot_acc.notna() & sf.label_resnum.notna()]
    smap = {}
    for ch, lr, acc, pos, rn in zip(sf.chain, sf.label_resnum, sf.uniprot_acc, sf.uniprot_pos, sf.uniprot_resname):
        smap.setdefault((ch, int(lr)), {})[acc] = (int(pos), rn)
    block = S.read_block(PDB_MIRROR / pdb_id[1:3] / f"{pdb_id}.cif.gz")
    ents = S.entities(block)
    a2e = S.asym_to_entity(block)
    asm = S.assembly(block, None)
    gens = asm.gens if (asm is not None and asm.gens) else [(list(a2e), [("1", np.eye(4))])]
    at = S.atoms(block)
    keep = at["first_alt"]
    by_asym = defaultdict(list)
    for i, a in enumerate(at["label_asym_id"]):
        if keep[i]:
            by_asym[a].append(i)
    poly = {e for e, E in ents.items() if E.type == "polymer"}
    chains = {}
    seen = set()
    for asyms, ops in gens:
        for opname, M in ops:
            ident = np.allclose(M, np.eye(4))
            for a in asyms:
                e = a2e.get(a)
                if e not in poly or a not in by_asym:
                    continue
                idx = np.array(by_asym[a])
                auth = at["auth_asym_id"][idx[0]]
                cid = auth if ident else f"{auth}-{opname}"
                while cid in seen:
                    cid += "'"
                seen.add(cid)
                if not ents[e].is_protein:
                    continue
                xyz = at["xyz"][idx] @ M[:3, :3].T + M[:3, 3]
                ls = np.array([int(x) if x not in (".", "?") else -10**6 for x in at["label_seq_id"][idx]])
                chains[cid] = dict(auth=auth, xyz=xyz, ls=ls, atom=at["label_atom_id"][idx], comp=at["label_comp_id"][idx])
    # majority accession per chain from SIFTS
    for cid, c in chains.items():
        cnt = defaultdict(int)
        for l in np.unique(c["ls"]):
            for acc in smap.get((c["auth"], int(l)), {}):
                cnt[acc] += 1
        c["acc"] = max(cnt, key=cnt.get) if cnt else None
        c["tree"] = cKDTree(c["xyz"])

    def iface(c1, c2):
        A, B = chains[c1], chains[c2]
        hits = B["tree"].query_ball_point(A["xyz"], IFACE_DIST)
        ra = {int(A["ls"][i]) for i, h in enumerate(hits) if h}
        rb = {int(B["ls"][j]) for h in hits for j in h}
        return ra, rb

    def to_unp(cid, labels, acc):
        auth = chains[cid]["auth"]
        out = []
        for l in labels:
            m = smap.get((auth, l), {}).get(acc)
            if m:
                out.append(m[0])
        return sorted(set(out))

    out_sites, iface_cache, ca_needed = [], {}, set()
    for s in sites:
        c1, acc1, acc2 = s["chain1"], s["acc1"], s["acc2"]
        if c1 not in chains:
            out_sites.append(dict(s, status="chain1_not_found")); continue
        c2 = s["chain2"]
        if c2 is None:
            cands = [k for k, c in chains.items() if k != c1 and c["acc"] == acc2]
            if not cands:
                out_sites.append(dict(s, status="no_partner_chain")); continue
            best, bn = None, -1
            for k in cands:
                key = (c1, k)
                if key not in iface_cache:
                    iface_cache[key] = iface(c1, k)
                n = len(iface_cache[key][0])
                if n > bn:
                    best, bn = k, n
            c2 = best
        if c2 not in chains:
            out_sites.append(dict(s, status="chain2_not_found")); continue
        key = (c1, c2)
        if key not in iface_cache:
            iface_cache[key] = iface(c1, c2)
        ia, ib = iface_cache[key]
        bind = json.loads(s["bind"])
        n_tot = n_ok = 0
        bind_unp = {"1": [], "2": []}
        for role, cid, acc in (("1", c1, acc1), ("2", c2, acc2)):
            for l, comp in bind.get(cid, []):
                n_tot += 1
                m = smap.get((chains[cid]["auth"], l), {}).get(acc) if l is not None else None
                if m:
                    bind_unp[role].append(m[0])
                    if m[1] == _one_letter(comp):
                        n_ok += 1
        # SBS binding residues on other chains (not in the pair) are ignored by construction
        out_sites.append(dict(s, status="ok", chain2=c2, chain1_auth=chains[c1]["auth"], chain2_auth=chains[c2]["auth"],
                              chain1_sifts_acc=chains[c1]["acc"], chain2_sifts_acc=chains[c2]["acc"],
                              n_bind_pair=n_tot, n_bind_identical=n_ok, resname_identity=(n_ok / n_tot if n_tot else np.nan),
                              bind1=sorted(set(bind_unp["1"])), bind2=sorted(set(bind_unp["2"])),
                              pdb_iface1=to_unp(c1, ia, acc1), pdb_iface2=to_unp(c2, ib, acc2), n_iface_res_pdb=len(ia) + len(ib)))
        ca_needed.add((c1, acc1)); ca_needed.add((c2, acc2))
    ca_rows = []
    for cid, acc in ca_needed:
        c = chains[cid]
        m = c["atom"] == "CA"
        pos, xyz = [], []
        for l, x in zip(c["ls"][m], c["xyz"][m]):
            mm = smap.get((c["auth"], int(l)), {}).get(acc)
            if mm:
                pos.append(mm[0]); xyz.append(x)
        ca_rows.append(dict(pdb_id=pdb_id, chain=cid, acc=acc, pos=pos, xyz=np.asarray(xyz, float).ravel().tolist()))
    return {"pdb_id": pdb_id, "error": None, "sites": out_sites, "ca": ca_rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True, choices=["v1", "v2"])
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit-entries", type=int, default=None)
    ap.add_argument("--no-fetch", action="store_true", help="do not download; process entries whose XML is cached")
    ap.add_argument("--stage-a-only", action="store_true", help="write the entry list and stop")
    ap.add_argument("--redo", action="store_true", help="recompute stage A even if cached")
    ap.add_argument("--entries", nargs="*", help="restrict to these PDB ids (testing)")
    a = ap.parse_args()
    od = OUT_DATA / a.census
    od.mkdir(parents=True, exist_ok=True)
    log = []
    cache = od / "stageA.parquet"
    if cache.exists() and not a.redo:
        sm = pd.read_parquet(cache)
        log += [tuple(x) for x in pd.read_csv(od / "stageA_log.tsv", sep="\t").itertuples(index=False)]
        log = [(s, n, (None if pd.isna(e) else e)) for s, n, e in log]
    else:
        sm = stage_a(a.census, log)
        sm.to_parquet(cache)
        pd.DataFrame(log, columns=["step", "n", "n_entries"]).to_csv(od / "stageA_log.tsv", sep="\t", index=False)
    sm, entries = choose_entries(sm, log)
    if a.stage_a_only:
        open(od / "sifts_entries.txt", "w").write("\n".join(entries) + "\n")
        return
    if a.entries:
        entries = [e for e in entries if e in set(a.entries)]
    if a.limit_entries:
        entries = entries[: a.limit_entries]
    sm = sm[sm.pdb_id.isin(entries)]
    sm.to_parquet(od / "site_candidates.parquet")
    # SIFTS XML: fetch only needed entries (<=4 parallel; cached). --no-fetch: use what the background fetch has cached.
    from fetch_sifts_xml import fetch_sifts_xml, CACHE
    t = time.time()
    if not a.no_fetch:
        st = fetch_sifts_xml(entries)
        bad = {k: v for k, v in st.items() if v not in ("ok", "cached")}
    else:
        bad = {e: "not_cached" for e in entries if not (CACHE / f"{e}.xml.gz").exists()}
    log.append(("SIFTS XML available", len(entries) - len(bad), None))
    log.append(("SIFTS XML missing/failed", len(bad), None))
    print(f"[B] SIFTS {len(entries)} entries in {time.time()-t:.0f}s, missing {len(bad)}", flush=True)
    recs = {k: g.drop(columns=["unit"]).to_dict("records") for k, g in sm.groupby("pdb_id")}
    for g in recs.values():
        for s in g:
            s["chain2"] = s["chain2"] if isinstance(s["chain2"], str) else None
    # per-entry result cache (resumable; an entry is processed once its SIFTS XML exists)
    ec = od / "entry_cache"
    ec.mkdir(exist_ok=True)
    import pickle
    results = {}
    todo = []
    for k in entries:
        if k not in recs or k in bad:
            continue
        f = ec / f"{k}.pkl"
        if f.exists():
            results[k] = pickle.load(open(f, "rb"))
        else:
            todo.append((k, recs[k]))
    print(f"[B] cached {len(results)}, to process {len(todo)}", flush=True)
    t = time.time()
    with ProcessPoolExecutor(max_workers=min(a.workers, 6)) as ex:
        for i, res in enumerate(ex.map(process_entry, todo, chunksize=4)):
            results[res["pdb_id"]] = res
            if res["error"] != "no_sifts_xml":
                pickle.dump(res, open(ec / f"{res['pdb_id']}.pkl", "wb"))
            if (i + 1) % 500 == 0:
                print(f"[B] {i+1}/{len(todo)} entries {time.time()-t:.0f}s", flush=True)
    out_sites, out_ca, errs = [], [], []
    for k in entries:
        res = results.get(k)
        if res is None:
            continue
        if res["error"]:
            errs.append((k, res["error"]))
        out_sites += res["sites"]; out_ca += res["ca"]
    sites = pd.DataFrame(out_sites)
    sites["unit"] = sites.model + "|" + sites.cluster30.astype(str) + "|" + sites.ccd_key + "|" + np.where(sites.is_mbs, "MBS", "SBS")
    sites.to_parquet(od / "sites.parquet")
    pd.DataFrame(out_ca).to_parquet(od / "chain_ca.parquet")
    pd.DataFrame(errs, columns=["pdb_id", "error"]).to_csv(od / "entry_errors.tsv", sep="\t", index=False)
    log.append(("entries with errors", len(errs), None))
    for k, v in sites.status.value_counts().items():
        log.append((f"site x model rows status={k}", int(v), None))
    led = pd.DataFrame(log, columns=["step", "n", "n_entries"])
    write_tsv(led, OUT_RES / a.census / "selection_ledger.tsv", {"census": a.census})
    print(led.to_string(), flush=True)


if __name__ == "__main__":
    main()
