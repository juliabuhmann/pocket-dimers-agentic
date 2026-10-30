#!/usr/bin/env python
"""Per-CCD drug-likeness flags from the local ChEMBL SQLite (DESIGN.md v3 §0, "Drug-like (primary)").

Status: first attempt, not independently verified.

Inputs (all produced by scripts in this repo):
  data/processed/ccd_to_chembl.parquet   UniChem PDBe(src3) -> ChEMBL(src1) mapping (scripts/fetch_chembl.py)
  data/part1/ligand_classes.parquet      CCD classes, cofactor_class (PDBe CoFactor list), RDKit Ro5/QED
  data/raw/chembl_sqlite/chembl_<REL>/chembl_<REL>_sqlite/chembl_<REL>.db   (scripts/fetch_chembl_sqlite.sh)

Definition implemented (DESIGN §0):
  drug-like (primary) = CCD maps to ChEMBL (UniChem) AND (max_phase >= 1 OR any pChEMBL >= 6 on any target)
                        AND not a cofactor/metabolite.
  - max_phase and pChEMBL are read for the mapped molecule AND its parent (molecule_hierarchy), max taken
    (salt forms / parents share activities; this choice is ours, flagged in column `used_parent`).
  - "not a cofactor/metabolite": the PDBe CoFactor list (cofactor_class) plus the ligand-class precedence of
    src/part1/config.py CLASS_ORDER (heme, cofactor, nucleotide-like, sugar, amino acid, lipid classes win
    before drug_like), BioLiP/PLINDER solvent lists, and ChEBI metabolite/buffer roles (chebi_exclusions).
  - heavy-atom floor DRUGLIKE_MIN_HEAVY_ATOMS (src/part1/config.py, = census choice) is kept.
  Secondary: Ro5 (<= DRUGLIKE_MAX_RO5_VIOLATIONS) + QED >= DRUGLIKE_QED_MIN (same floor).
FALLBACK: if the SQLite is absent/unreadable, primary := Ro5+QED and `druglike_fallback` = True (DESIGN §0).
CCDs without a ChEMBL mapping are never drug-like under the primary definition.
Exclusions (cofactor/metabolite/buffer): census class precedence + PDBe CoFactor list; BioLiP artifact list and PLINDER
badlist (buffers/solvents); ChEBI roles (metabolite / buffer subtrees, see chebi_exclusions).

Output: data/processed/ccd_druglike.parquet (+ .json sidecar with ChEMBL release, sqlite sha256-of-tarball, counts).
"""
import argparse
import json
import sqlite3
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from part1.config import (CLASS_ORDER, DRUGLIKE_MAX_RO5_VIOLATIONS, DRUGLIKE_MIN_CHEMBL_PHASE,  # noqa: E402
                          DRUGLIKE_MIN_HEAVY_ATOMS, DRUGLIKE_QED_MIN)

CHEMBL_RELEASE = "chembl_37"          # must match scripts/fetch_chembl_sqlite.sh REL
PCHEMBL_MIN = 6.0                     # DESIGN §0: pChEMBL >= 6 (1 uM) on any target
SQLITE_DIR = ROOT / "data/raw/chembl_sqlite"
MAP = ROOT / "data/processed/ccd_to_chembl.parquet"
CLASSES = ROOT / "data/part1/ligand_classes.parquet"
OUT = ROOT / "data/processed/ccd_druglike.parquet"
# classes that precede drug_like in CLASS_ORDER = "cofactor or metabolite-like or non-organic" -> never drug-like
# Census classes the drug-like label may override. "amino_acid" is the census class for CCD type "peptide-linking",
# which also captures peptidomimetic drugs (e.g. ritonavir RIT); free amino acids are excluded via ChEBI metabolite
# roles instead. Decision 2026-10-29 (PI sentinel RIT must be drug-like). All other earlier classes win.
DRUGLIKE_OVERRIDABLE_CLASSES = ["amino_acid", "drug_like", "other_organic"]
NON_DRUG_CLASSES = set(CLASS_ORDER[:CLASS_ORDER.index("drug_like")]) - set(DRUGLIKE_OVERRIDABLE_CLASSES)
CHEBI_DIR = ROOT / "data/raw"   # chebi_* flat files, fetched by scripts/fetch_external.py --only chebi_*
# ChEBI roles (ascii names, resolved to ids at run time) defining an endogenous metabolite (our choice: model-organism
# metabolism + "fundamental metabolite"); xenobiotic/drug-metabolite roles are NOT counted as metabolites.
METABOLITE_ROLES = ["mammalian metabolite", "Escherichia coli metabolite", "Saccharomyces cerevisiae metabolite",
                    "fundamental metabolite"]
XENO_ROLES = ["xenobiotic metabolite", "drug metabolite"]
BUFFER_ROLES = ["buffer"]
# known-answer sentinels requested by the PI (2026-10-29): imatinib, ritonavir, methotrexate, tafamidis must be
# drug-like; MES, HEPES, citrate, cholesterol must not be.
KNOWN_DRUGS = ["STI", "RIT", "MTX", "3MI"]
KNOWN_NON_DRUGS = ["MES", "EPE", "CIT", "CLR"]


def find_db(override=None):
    if override:   # e.g. a byte-identical copy on node-local disk (network FS random reads are very slow)
        return Path(override)
    hits = sorted(SQLITE_DIR.rglob(f"{CHEMBL_RELEASE}.db"))
    return hits[0] if hits else None


def query_chembl(db, chembl_ids):
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    ver = con.execute("SELECT name, creation_date FROM version").fetchall()
    con.execute("CREATE TEMP TABLE q(chembl_id TEXT PRIMARY KEY)")
    con.executemany("INSERT OR IGNORE INTO q VALUES (?)", [(c,) for c in chembl_ids])
    # mapped molecule + its parent
    mol = pd.read_sql_query("""
        SELECT q.chembl_id, md.molregno, md.max_phase, mh.parent_molregno, pm.max_phase AS parent_max_phase
        FROM q JOIN molecule_dictionary md ON md.chembl_id = q.chembl_id
        LEFT JOIN molecule_hierarchy mh ON mh.molregno = md.molregno
        LEFT JOIN molecule_dictionary pm ON pm.molregno = mh.parent_molregno""", con)
    regs = set(mol.molregno.dropna().astype(int)) | set(mol.parent_molregno.dropna().astype(int))
    con.execute("CREATE TEMP TABLE r(molregno INTEGER PRIMARY KEY)")
    con.executemany("INSERT OR IGNORE INTO r VALUES (?)", [(int(x),) for x in regs])
    act = pd.read_sql_query("""
        SELECT a.molregno, MAX(a.pchembl_value) AS max_pchembl, COUNT(*) AS n_pchembl,
               SUM(a.pchembl_value >= ?) AS n_pchembl_ge,
               COUNT(DISTINCT CASE WHEN a.pchembl_value >= ? THEN ass.tid END) AS n_targets_ge
        FROM r JOIN activities a ON a.molregno = r.molregno
        JOIN assays ass ON ass.assay_id = a.assay_id
        WHERE a.pchembl_value IS NOT NULL GROUP BY a.molregno""", con, params=(PCHEMBL_MIN, PCHEMBL_MIN))
    con.close()
    a = act.set_index("molregno")
    rows = []
    for r in mol.itertuples():
        own = a.loc[r.molregno] if r.molregno in a.index else None
        par = a.loc[r.parent_molregno] if (pd.notna(r.parent_molregno) and r.parent_molregno in a.index
                                           and r.parent_molregno != r.molregno) else None
        pch = [x.max_pchembl for x in (own, par) if x is not None]
        tg = [x.n_targets_ge for x in (own, par) if x is not None]
        ph = [p for p in (r.max_phase, r.parent_max_phase) if pd.notna(p)]
        rows.append(dict(chembl_id=r.chembl_id, molregno=r.molregno,
                         max_phase_sqlite=max(ph) if ph else None,
                         max_pchembl=max(pch) if pch else None,
                         n_targets_pchembl_ge6=int(max(tg)) if tg else 0,
                         used_parent=bool(par is not None)))
    return pd.DataFrame(rows), ver


def chebi_exclusions(ccds, ccd_chembl):
    """CCD -> (is_endogenous_metabolite, is_buffer, chebi_ids) from ChEBI roles.

    CCD->ChEBI: ChEBI database_accession rows from PDBeChem (CCD code) + ChEBI's ChEMBL xrefs of the UniChem-mapped
    ChEMBL id (ChEMBL 37 SQLite has no chebi_par_id column).
    Entity expansion: conjugate acid/base, tautomer links (symmetric, closure) then is_a ancestors (roles are
    inherited along is_a). Role sets: is_a-descendants of METABOLITE_ROLES minus descendants of XENO_ROLES; BUFFER_ROLES.
    """
    comp = pd.read_csv(CHEBI_DIR / "chebi_compounds/compounds.tsv.gz", sep="\t", usecols=["id", "ascii_name", "parent_id"],
                       dtype={"parent_id": "Int64"})
    rel = pd.read_csv(CHEBI_DIR / "chebi_relation/relation.tsv.gz", sep="\t", usecols=["relation_type_id", "init_id", "final_id"])
    rtyp = pd.read_csv(CHEBI_DIR / "chebi_relation_type/relation_type.tsv.gz", sep="\t").set_index("code").id
    src = pd.read_csv(CHEBI_DIR / "chebi_source/source.tsv.gz", sep="\t")
    acc = pd.read_csv(CHEBI_DIR / "chebi_dbacc/database_accession.tsv.gz", sep="\t",
                      usecols=["compound_id", "accession_number", "source_id"])
    primary_of = comp.set_index("id").parent_id.dropna().astype(int).to_dict()   # secondary -> primary id
    canon = lambda i: primary_of.get(int(i), int(i))  # noqa: E731
    name2id = comp.drop_duplicates("ascii_name").set_index("ascii_name").id
    pdbechem = src.loc[src.name == "PDBeChem", "id"].iloc[0]
    a = acc[acc.source_id == pdbechem]
    ccd2chebi = {}
    for c, i in zip(a.accession_number, a.compound_id):
        ccd2chebi.setdefault(str(c), set()).add(canon(i))
    chembl_src = src.loc[src.name == "ChEMBL", "id"].iloc[0]
    ac = acc[acc.source_id == chembl_src]
    chembl2chebi = {}
    for c, i in zip(ac.accession_number, ac.compound_id):
        chembl2chebi.setdefault(str(c), set()).add(canon(i))
    for c, ch in ccd_chembl:
        if pd.notna(ch):
            for i in chembl2chebi.get(ch, ()):
                ccd2chebi.setdefault(c, set()).add(i)
    from collections import defaultdict
    isa_up, isa_down, sym, roles = defaultdict(set), defaultdict(set), defaultdict(set), defaultdict(set)
    for t, x, y in zip(rel.relation_type_id, rel.init_id, rel.final_id):
        x, y = canon(x), canon(y)
        if t == rtyp["is_a"]:
            isa_up[x].add(y); isa_down[y].add(x)
        elif t in (rtyp["is_conjugate_acid_of"], rtyp["is_conjugate_base_of"], rtyp["is_tautomer_of"]):
            sym[x].add(y); sym[y].add(x)
        elif t == rtyp["has_role"]:
            roles[x].add(y)

    def closure(seeds, nbr):
        seen, st = set(seeds), list(seeds)
        while st:
            for n in nbr.get(st.pop(), ()):
                if n not in seen:
                    seen.add(n); st.append(n)
        return seen

    rid = lambda names: {int(name2id[n]) for n in names}  # noqa: E731
    met = closure(rid(METABOLITE_ROLES), isa_down) - closure(rid(XENO_ROLES), isa_down)
    buf = closure(rid(BUFFER_ROLES), isa_down)
    out = {}
    for c in ccds:
        ids = ccd2chebi.get(c, set())
        if not ids:
            out[c] = (False, False, "")
            continue
        ents = closure(closure(ids, sym), isa_up)
        rs = set().union(*(roles.get(e, set()) for e in ents)) if ents else set()
        out[c] = (bool(rs & met), bool(rs & buf), ";".join(f"CHEBI:{i}" for i in sorted(ids)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", default=str(CLASSES), help="ligand_classes.parquet of the census version used")
    ap.add_argument("--db", default=None, help="path to a copy of chembl_<REL>.db (default: under data/raw/chembl_sqlite)")
    args = ap.parse_args()
    cls = pd.read_parquet(args.classes, columns=["ccd", "name", "n_heavy", "cofactor_class", "ligand_class",
                                                 "ro5_violations", "qed", "chembl_max_phase",
                                                 "in_biolip_artifact_list", "plinder_badlist"])
    m = pd.read_parquet(MAP)
    df = cls.merge(m[["ccd", "chembl_id"]], on="ccd", how="left")
    db = find_db(args.db)
    meta = {"chembl_release_expected": CHEMBL_RELEASE, "pchembl_min": PCHEMBL_MIN, "classes_file": args.classes,
            "min_phase": DRUGLIKE_MIN_CHEMBL_PHASE, "min_heavy_atoms": DRUGLIKE_MIN_HEAVY_ATOMS,
            "qed_min": DRUGLIKE_QED_MIN, "max_ro5_violations": DRUGLIKE_MAX_RO5_VIOLATIONS,
            "non_drug_classes": sorted(NON_DRUG_CLASSES),
            "druglike_overridable_classes": DRUGLIKE_OVERRIDABLE_CLASSES, "chebi_metabolite_roles": METABOLITE_ROLES,
            "chebi_xeno_roles_not_excluded": XENO_ROLES, "chebi_buffer_roles": BUFFER_ROLES}
    fallback = db is None
    if not fallback:
        try:
            q, ver = query_chembl(db, m.chembl_id.dropna().unique().tolist())
            meta["sqlite_path"] = str(db)
            meta["sqlite_version_table"] = [list(map(str, v)) for v in ver]
            tgz_sha = SQLITE_DIR / f"{CHEMBL_RELEASE}_sqlite.tar.gz.sha256"
            meta["tarball_sha256"] = tgz_sha.read_text().split()[0] if tgz_sha.exists() else None
            df = df.merge(q, on="chembl_id", how="left")
        except Exception as e:  # noqa: BLE001
            print("SQLite query failed -> FALLBACK:", e)
            fallback = True
            meta["sqlite_error"] = repr(e)
    if fallback:
        for c in ("molregno", "max_phase_sqlite", "max_pchembl", "n_targets_pchembl_ge6", "used_parent"):
            df[c] = None
    ex = chebi_exclusions(df.ccd.tolist(), list(zip(df.ccd, df.chembl_id)))
    df["chebi_ids"] = [ex[c][2] for c in df.ccd]
    df["chebi_endogenous_metabolite"] = [ex[c][0] for c in df.ccd]
    df["chebi_buffer"] = [ex[c][1] for c in df.ccd]
    big = df.n_heavy >= DRUGLIKE_MIN_HEAVY_ATOMS
    df["excl_cofactor_class"] = df.cofactor_class.notna() | df.ligand_class.isin(NON_DRUG_CLASSES)
    df["excl_buffer_solvent_list"] = df.in_biolip_artifact_list.fillna(False) | df.plinder_badlist.fillna(False) | df.chebi_buffer
    df["excl_metabolite"] = df.chebi_endogenous_metabolite
    excl = df.excl_cofactor_class | df.excl_buffer_solvent_list | df.excl_metabolite
    df["has_phase_ge1"] = df.max_phase_sqlite.astype(float).ge(DRUGLIKE_MIN_CHEMBL_PHASE)
    df["has_pchembl_ge6"] = df.max_pchembl.astype(float).ge(PCHEMBL_MIN)
    # secondary definition (DESIGN §0): Ro5 + QED, same exclusions; never used to fill unmapped CCDs in the primary
    df["druglike_ro5qed"] = (df.ro5_violations.le(DRUGLIKE_MAX_RO5_VIOLATIONS) & df.qed.ge(DRUGLIKE_QED_MIN) & big & ~excl)
    if fallback:   # DESIGN §0 declared FALLBACK: whole primary label replaced, flagged
        df["druglike_primary"] = df.druglike_ro5qed
    else:          # unmapped CCDs are NOT drug-like under the primary definition
        df["druglike_primary"] = df.chembl_id.notna() & (df.has_phase_ge1 | df.has_pchembl_ge6) & big & ~excl
    df["druglike_fallback"] = fallback
    df["druglike_reason"] = (df.has_phase_ge1.map({True: "phase", False: ""}) + df.has_pchembl_ge6.map(
        {True: "+pchembl", False: ""})).where(df.druglike_primary, None)
    df["chembl_release"] = None if fallback else CHEMBL_RELEASE
    pre = df.chembl_id.notna() & (df.has_phase_ge1 | df.has_pchembl_ge6) & big
    meta.update(fallback=fallback, n_ccd=int(len(df)), n_mapped=int(df.chembl_id.notna().sum()),
                n_found_in_sqlite=int(df.molregno.notna().sum()),
                n_phase_ge1=int(df.has_phase_ge1.sum()), n_pchembl_ge6=int(df.has_pchembl_ge6.sum()),
                n_chembl_rule_before_exclusions=int(pre.sum()),
                n_removed_cofactor_class=int((pre & df.excl_cofactor_class).sum()),
                n_removed_buffer_solvent=int((pre & ~df.excl_cofactor_class & df.excl_buffer_solvent_list).sum()),
                n_removed_metabolite=int((pre & ~df.excl_cofactor_class & ~df.excl_buffer_solvent_list & df.excl_metabolite).sum()),
                n_ccd_with_chebi=int((df.chebi_ids != "").sum()),
                n_druglike_primary=int(df.druglike_primary.sum()), n_druglike_ro5qed=int(df.druglike_ro5qed.sum()))
    both = df.chembl_max_phase.notna() | df.max_phase_sqlite.notna()
    meta["phase_rest_vs_sqlite_agree_frac"] = float(
        (df.chembl_max_phase.fillna(-9).astype(float) == df.max_phase_sqlite.astype(float).fillna(-9))[both].mean()) if both.any() else None
    # known-answer sentinels (PI list, 2026-10-29)
    d = df.set_index("ccd")
    sent = {c: bool(d.druglike_primary.get(c, False)) for c in KNOWN_DRUGS + KNOWN_NON_DRUGS}
    meta["sentinels"] = sent
    meta["sentinels_pass"] = all(sent[c] for c in KNOWN_DRUGS) and not any(sent[c] for c in KNOWN_NON_DRUGS)
    out = df.drop(columns=["ro5_violations", "qed", "chembl_max_phase"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT, index=False)
    OUT.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))
    if not meta["sentinels_pass"]:
        print("SENTINEL FAILURE:", {c: v for c, v in sent.items()})
        sys.exit(2)


if __name__ == "__main__":
    main()
