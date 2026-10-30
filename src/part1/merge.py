"""Merge shards, attach annotations (ligand classes, BioLiP2 relevance, clusters, H1 representatives),
write ligand_sites.parquet / entries.parquet, the exclusion ledger and summary tables.
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import re

from .config import (BIOLIP_MAX_COPIES, CENSUS_VERSION, OUT_DIR, PHOTOSYNTH_CCD_NAME_RE, PHOTOSYNTH_KEYWORD_RE,
                     raw_path)

PRIMARY_ASSEMBLY_TYPES = ("monomer", "homomer", "heteromer")  # v2: same_uniprot_hetero reported separately

SLIM_DROP = ["contacts", "chain_counts", "chain_counts_1atom", "chain_counts_4A", "components",
             "covalent_partner_comps", "centroid"]


def load_biolip_index():
    """Set of (pdb, ccd, ligand auth chain, ligand auth_seq) and set of pdb ids present in BioLiP.txt."""
    cols = [0, 4, 5, 19]
    df = pd.read_csv(raw_path("biolip_annotation"), sep="\t", header=None, usecols=cols, dtype=str,
                     quoting=3, na_filter=False)
    df.columns = ["pdb", "ccd", "lig_chain", "lig_seq"]
    keys = set(zip(df.pdb.str.lower(), df.ccd, df.lig_chain, df.lig_seq.str.strip()))
    return keys, set(df.pdb.str.lower())


def load_clusters(name):
    """'PDBID_entity' -> cluster index (line number) for RCSB clusters-by-entity file."""
    out = {}
    with open(raw_path(name)) as f:
        for i, line in enumerate(f):
            for m in line.split():
                out[m] = i
    return out


def concat_shards(shard_dir: Path, prefix: str) -> pa.Table:
    files = sorted(shard_dir.glob(f"{prefix}_*.parquet"))
    tabs = [pq.read_table(f) for f in files]
    tabs = [t for t in tabs if t.num_rows]
    if not tabs:
        return pa.table({})
    return pa.concat_tables(unify(tabs), promote_options="permissive")


def unify(tabs):
    """Shards are written independently, so all-null columns may be typed null/double/bool differently.
    For each conflicting field take the type of the shard with the most non-null values and cast the others."""
    best = {}
    for t in tabs:
        for f in t.schema:
            nn = len(t[f.name]) - t[f.name].null_count
            if f.name not in best or nn > best[f.name][1]:
                best[f.name] = (f.type, nn)
    out = []
    for t in tabs:
        for i, f in enumerate(t.schema):
            tgt = best[f.name][0]
            if f.type != tgt and not pa.types.is_null(tgt):
                t = t.set_column(i, pa.field(f.name, tgt), t[f.name].cast(tgt))
        out.append(t)
    return out


def annotate(lig: pd.DataFrame, entries: pd.DataFrame, classes: pd.DataFrame):
    """Add label columns to the slim ligand frame (in place) and return it."""
    keys, biolip_pdbs = load_biolip_index()
    # BioLiP coverage: entries released after the newest BioLiP-covered release date are 'not_covered'
    cov_dates = entries.loc[entries.pdb_id.isin(biolip_pdbs), "release_date"].dropna()
    cutoff = cov_dates.max() if len(cov_dates) else None
    in_b = np.array([(p, c, ch, s) in keys for p, c, ch, s in
                     zip(lig.pdb_id, lig.ccd, lig.main_auth_chain, lig.main_auth_seq.astype(str))])
    entry_in = lig.pdb_id.isin(biolip_pdbs).to_numpy()
    covered = entry_in | (lig.release_date.fillna("9999") <= (cutoff or "0000")).to_numpy()
    lig["biolip_status"] = np.where(in_b, "relevant", np.where(covered, "additive", "not_covered"))
    lig.attrs["biolip_cutoff"] = cutoff

    c = classes.set_index("ccd")
    keepcols = ["ligand_class", "marsh_class", "in_biolip_artifact_list", "plinder_badlist", "plinder_generic_filter",
                "plinder_artifact", "chembl_max_phase", "is_approved_drug", "mw", "qed", "ro5_pass", "n_heavy",
                "cofactor_class"]
    for k in keepcols:
        lig[k] = lig.ccd.map(c[k]) if k in c else None
    lig["ligand_class"] = lig["ligand_class"].fillna("unknown_ccd")
    lig["marsh_class"] = lig["marsh_class"].fillna("small_molecule")
    lig.loc[lig.is_glycan, "ligand_class"] = "sugar_glycan"
    # BioLiP2 geometric rule (sensitivity / fallback): listed AND (>15 copies OR 2-consecutive-residue site)
    listed = lig.in_biolip_artifact_list.fillna(False).astype(bool)
    lig["biolip_list_geom_artifact"] = listed & ((lig.n_copies_in_au > BIOLIP_MAX_COPIES) | lig.biolip2_two_consecutive_site)
    lig["biolip_list_strict_artifact"] = listed
    fb = lig.biolip_list_geom_artifact | lig.plinder_artifact.fillna(False).astype(bool)
    # v2: relevant = BioLiP2-relevant OR wwPDB SUBJECT OF INVESTIGATION (DESIGN s.0); coverage is its own column
    lig["biolip_covered"] = lig.biolip_status != "not_covered"
    loi = lig.is_ligand_of_interest.fillna(False).astype(bool)
    lig["is_relevant"] = (lig.biolip_status == "relevant") | loi
    lig["relevance_source"] = np.where(lig.biolip_status == "relevant", np.where(loi, "biolip2+soi", "biolip2"),
                                       np.where(loi, "soi", None))
    lig["is_additive"] = np.where(lig.biolip_covered, ~lig.is_relevant, fb & ~loi)
    lig["additive_label_source"] = np.where(lig.biolip_covered, "biolip2_or_soi", "FALLBACK_geom_plinder_or_soi")
    # v2: photosynthetic complexes (pigment CCDs by CCD name, or entry keywords)
    pig = set(classes.loc[classes.name.fillna("").str.contains(PHOTOSYNTH_CCD_NAME_RE, flags=re.I), "ccd"])
    kw = entries.set_index("pdb_id").keywords.fillna("")
    ps_entries = set(kw[kw.str.contains(PHOTOSYNTH_KEYWORD_RE, flags=re.I)].index) | set(lig.pdb_id[lig.ccd.isin(pig)])
    lig["is_photosynthetic"] = lig.pdb_id.isin(ps_entries)
    lig.attrs["photosynthetic_entries"] = sorted(ps_entries)
    lig["census_version"] = CENSUS_VERSION

    # clusters for the dominant contacted chain's entity
    for lvl in (30, 95):
        cl = load_clusters(f"rcsb_clusters_{lvl}")
        k = lig.pdb_id.str.upper() + "_" + lig.dominant_entity.fillna("").astype(str)
        lig[f"cluster{lvl}"] = k.map(cl).astype("Int64")

    # primary set cascade
    base = (~lig.is_dup_copy) & lig.is_bound & (~lig.is_glycan) & (~lig.is_additive.astype(bool)) \
        & lig.assembly_type.isin(PRIMARY_ASSEMBLY_TYPES)
    lig["in_primary_set"] = base & lig.biolip_covered            # v2: primary = BioLiP2-covered entries only
    lig["in_postsnapshot_sensitivity_set"] = base & ~lig.biolip_covered
    # H1 representative per (cluster30 x CCD): best resolution, ties by earliest pdb_id
    p = lig[lig.in_primary_set & lig.cluster30.notna()][["pdb_id", "ccd", "cluster30", "resolution", "ligand_uid"]]
    best = p.sort_values(["resolution", "pdb_id", "ligand_uid"]).drop_duplicates(["cluster30", "ccd"])
    rep_entries = set(zip(best.cluster30, best.ccd, best.pdb_id))
    lig["is_h1_rep_entry"] = np.array([(cl, cc, pid) in rep_entries for cl, cc, pid in
                                       zip(lig.cluster30, lig.ccd, lig.pdb_id)], dtype=bool) & lig.in_primary_set
    lig["is_h1_rep_instance"] = lig.ligand_uid.isin(set(best.ligand_uid))
    return lig


def ledger(entries: pd.DataFrame, lig: pd.DataFrame) -> pd.DataFrame:
    rows = []
    n = len(entries)
    order = ["method", "no_resolution", "resolution", "no_protein_chain_ge30", "huge_assembly", "no_protein_in_assembly",
             "error"]
    for r in order:
        k = int((entries.skip_reason == r).sum())
        rows.append(dict(level="entry", step=f"drop:{r}", n_in=n, n_out=n - k, n_dropped=k))
        n -= k
    has_lig = entries[entries.skip_reason.isna()].pdb_id.isin(set(lig.pdb_id)).sum()
    rows.append(dict(level="entry", step="info:entries_with_ligand_instances", n_in=n, n_out=int(has_lig), n_dropped=n - int(has_lig)))
    m = len(lig)
    steps = [("drop:assembly_type same_uniprot_hetero (reported separately)", ~lig.assembly_type.isin(PRIMARY_ASSEMBLY_TYPES)),
             ("drop:symmetry/altloc/overlap duplicate copy", lig.is_dup_copy),
             ("drop:not bound (<2 protein residues with >=2 atom contacts)", ~lig.is_bound),
             ("drop:glycan (branched or covalent sugar)", lig.is_glycan),
             ("drop:entry released after BioLiP2 snapshot (sensitivity set)", ~lig.biolip_covered),
             ("drop:additive (not BioLiP2-relevant and not SUBJECT OF INVESTIGATION)", lig.is_additive.astype(bool))]
    mask = pd.Series(True, index=lig.index)
    for name, bad in steps:
        k = int((mask & bad).sum())
        rows.append(dict(level="ligand_instance", step=name, n_in=m, n_out=m - k, n_dropped=k))
        mask &= ~bad
        m -= k
    return pd.DataFrame(rows)


def summaries(lig: pd.DataFrame) -> dict:
    out = {}
    base = lig[(~lig.is_dup_copy) & lig.is_bound & (~lig.is_glycan)]
    base = base[base.biolip_covered]
    frames = {"before_additive_filter": base, "primary_set": lig[lig.in_primary_set],
              "primary_excl_photosynthetic": lig[lig.in_primary_set & ~lig.is_photosynthetic],
              "postsnapshot_sensitivity": lig[lig.in_postsnapshot_sensitivity_set],
              "same_uniprot_hetero": lig[(~lig.is_dup_copy) & lig.is_bound & (~lig.is_glycan) & lig.biolip_covered
                                         & (~lig.is_additive.astype(bool)) & (lig.assembly_type == "same_uniprot_hetero")]}
    rows = []
    for setname, d in frames.items():
        for cls_col in ("ligand_class", "marsh_class"):
            for (at, cl), g in d.groupby(["assembly_type", cls_col]):
                rows.append(dict(set=setname, scheme=cls_col, assembly_type=at, ligand_class=cl, n=len(g),
                                 frac_mbs=g.is_mbs.mean(), frac_mbs_1atom=g.is_mbs_1atom.mean(),
                                 frac_mbs_4A=g.is_mbs_4A.mean(), frac_mbs_min2=g.is_mbs_min2.mean(),
                                 frac_mbs_min3=g.is_mbs_min3.mean(), frac_mbs_excl_antibody=g.is_mbs_excl_antibody.mean()))
            for at, g in d.groupby("assembly_type"):
                rows.append(dict(set=setname, scheme=cls_col, assembly_type=at, ligand_class="ALL", n=len(g),
                                 frac_mbs=g.is_mbs.mean(), frac_mbs_1atom=g.is_mbs_1atom.mean(),
                                 frac_mbs_4A=g.is_mbs_4A.mean(), frac_mbs_min2=g.is_mbs_min2.mean(),
                                 frac_mbs_min3=g.is_mbs_min3.mean(), frac_mbs_excl_antibody=g.is_mbs_excl_antibody.mean()))
    out["mbs_by_assembly_class"] = pd.DataFrame(rows).drop_duplicates()
    b = base.copy()
    b["loi"] = b.is_ligand_of_interest
    agree = b.groupby(["biolip_status", "biolip_list_strict_artifact", "biolip_list_geom_artifact", "plinder_artifact",
                       "loi"], dropna=False).size().rename("n").reset_index()
    out["additive_label_agreement"] = agree
    # v2: site-level collapse on 2-copy homodimers (pdb x CCD groups with exactly 2 primary instances)
    p = lig[lig.in_primary_set & (lig.assembly_type == "homomer") & (lig.n_protein_chains == 2)]
    g = p.groupby(["pdb_id", "ccd"])
    two = p[g.ligand_uid.transform("size") == 2]
    gg = two.groupby(["pdb_id", "ccd"])
    out["site_collapse_homodimers"] = pd.DataFrame([dict(
        n_pairs=gg.ngroups, frac_collapsed_exact_hash=(gg.sym_group.nunique() == 1).mean(),
        frac_collapsed_jaccard_site_group=(gg.site_group.nunique() == 1).mean(),
        n_primary_homomer_instances=int((lig.in_primary_set & (lig.assembly_type == "homomer")).sum()),
        n_primary_homomer_site_groups=int(lig[lig.in_primary_set & (lig.assembly_type == "homomer")].site_group.nunique()))])
    hm = lig[lig.in_primary_set & (lig.assembly_type == "homomer") & lig.is_mbs]
    out["homomer_mbs_geometry"] = pd.DataFrame([dict(
        n=len(hm), frac_symmetry_copy=hm.mbs_contrib_symmetry_copy.mean(), frac_on_axis=hm.is_on_axis.mean(),
        median_c2_axis_dist=hm.c2_axis_dist.median(), frac_with_c2_axis=hm.c2_axis_dist.notna().mean(),
        frac_on_axis_geom=hm.is_on_axis_geom.mean(),
        agree_jaccard_vs_geom=(hm.is_on_axis.fillna(False).astype(bool) == hm.is_on_axis_geom.fillna(False).astype(bool)).mean(),
        n_jaccard_only=int((hm.is_on_axis.fillna(False).astype(bool) & ~hm.is_on_axis_geom.fillna(False).astype(bool)).sum()),
        n_geom_only=int((~hm.is_on_axis.fillna(False).astype(bool) & hm.is_on_axis_geom.fillna(False).astype(bool)).sum()))])
    b2 = base.copy()
    out["relevance_source_counts"] = b2.groupby(["biolip_status", "is_ligand_of_interest", "is_additive"],
                                                 dropna=False).size().rename("n").reset_index()
    return out
