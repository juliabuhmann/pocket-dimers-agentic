"""Ligand (CCD component) classification -> data/part1/ligand_classes.parquet.

One row per CCD component. Kept separate from the site table so it can be revised without rerunning contacts.
First-attempt rules (documented in docs/PART1_PIPELINE.md); class labels have not been manually audited.
"""
import gzip
import json
import re
from collections import defaultdict

import gemmi
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import QED, Crippen, Descriptors, Lipinski

from .config import (CLASS_ORDER, DRUGLIKE_MAX_RO5_VIOLATIONS, DRUGLIKE_MIN_CHEMBL_PHASE, DRUGLIKE_MIN_HEAVY_ATOMS, DRUGLIKE_QED_MIN,
                     LIPID_MIN_CHAIN, NONMETALS, PLINDER_MAX_ABS_CHARGE, PLINDER_MAX_UNBRANCHED_CHAIN,
                     PLINDER_MIN_CARBONS, PLINDER_MIN_HEAVY_ATOMS, raw_path)

RDLogger.DisableLog("rdApp.*")

NUCLEOSIDE = Chem.MolFromSmarts("[#7;R]!@[#6;R]1[#8;R][#6;R][#6;R][#6;R]1")  # N-glycosidic base-furanose link
LIPID_CHAIN = Chem.MolFromSmarts("-".join(["[CX4;!R]"] * LIPID_MIN_CHAIN))
UNBRANCHED_HC = Chem.MolFromSmarts("-".join(["[CX4;!R;D1,D2]"] * (PLINDER_MAX_UNBRANCHED_CHAIN + 1)))
HEME_NAME = re.compile(r"PORPHYRIN|HEME|\bHAEM|CHLOROPHYLL|CHLORIN|BACTERIOPHEOPHYTIN|PHEOPHYTIN|CORRIN|COBALAMIN", re.I)
SMILES_PREF = [("SMILES_CANONICAL", "CACTVS"), ("SMILES_CANONICAL", "OpenEye"), ("SMILES", "CACTVS"),
               ("SMILES", "OpenEye"), ("SMILES", "ACDLabs")]


def load_biolip_list():
    codes = []
    with open(raw_path("biolip_ligand_list")) as f:
        for line in f:
            if line.strip():
                codes.append(line.split("\t")[0].strip())
    return set(codes)


def load_plinder_badlist():
    with open(raw_path("plinder_badlist")) as f:
        return {l.strip() for l in f if l.strip() and not l.startswith("#")}


def load_cofactors():
    """CCD code -> PDBe/CoFactor class name."""
    d = json.load(open(raw_path("pdbe_cofactors")))
    out = {}
    for cls, entries in d.items():
        for e in entries:
            for code in e.get("cofactors", []):
                out.setdefault(code, cls)
    return out


def load_chembl():
    """CCD -> (sorted ChEMBL ids, max_phase over mapped molecules or None)."""
    m = defaultdict(set)
    with gzip.open(raw_path("unichem_chembl_pdbe"), "rt") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                m[p[1]].add(p[0])
    ph = pd.read_csv(raw_path("chembl_max_phase"), sep="\t")
    phase = dict(zip(ph.molecule_chembl_id, ph.max_phase.astype(float)))
    out = {}
    for ccd, ids in m.items():
        ps = [phase[i] for i in ids if i in phase]
        out[ccd] = (sorted(ids), max(ps) if ps else None)
    return out


def _smiles(block):
    t = block.find("_pdbx_chem_comp_descriptor.", ["type", "program", "descriptor"])
    rows = [(gemmi.cif.as_string(r[0]), gemmi.cif.as_string(r[1]), gemmi.cif.as_string(r[2])) for r in t]
    for typ, prog in SMILES_PREF:
        for a, b, c in rows:
            if a == typ and b.startswith(prog):
                return c
    return None


def rdkit_props(smi):
    if not smi:
        return {}
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            return {"rdkit_ok": False}
        return {"rdkit_ok": False}
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    viol = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
    try:
        q = QED.qed(mol)
    except Exception:
        q = None
    return dict(rdkit_ok=True, mw=mw, clogp=logp, hbd=hbd, hba=hba, ro5_violations=viol,
                ro5_pass=viol <= DRUGLIKE_MAX_RO5_VIOLATIONS, qed=q,
                has_nucleoside=mol.HasSubstructMatch(NUCLEOSIDE), has_lipid_chain=mol.HasSubstructMatch(LIPID_CHAIN),
                has_long_unbranched_hc=mol.HasSubstructMatch(UNBRANCHED_HC),
                n_rings=Chem.rdMolDescriptors.CalcNumRings(mol))


def classify(r):
    """Return the first matching class in CLASS_ORDER (see config)."""
    els = r["elements"]
    tests = {
        "metal_ion": r["n_heavy"] == 1 and els and els[0] not in NONMETALS,
        "metal_cluster": r["n_carbon"] == 0 and r["n_heavy"] > 1 and any(e not in NONMETALS for e in els),
        "inorganic": r["n_carbon"] == 0,
        "heme_porphyrin": bool(HEME_NAME.search(r["name"] or "")) or (r["cofactor_class"] or "") == "Heme",
        "cofactor": r["cofactor_class"] is not None,
        "nucleotide_like": bool(r.get("has_nucleoside")) or "NUCLEOTIDE" in r["type"] or "DNA" in r["type"]
                           or "RNA" in r["type"],
        "sugar_glycan": "SACCHARIDE" in r["type"],
        "amino_acid": "PEPTIDE" in r["type"] or "AMINO ACID" in r["type"],
        "lipid_fatty_acid": bool(r.get("has_lipid_chain")),
        "drug_like": r["is_druglike_chembl"] or r["is_druglike_ro5qed"],
        "other_organic": True,
    }
    for c in CLASS_ORDER:
        if tests[c]:
            return c
    return "other_organic"


MARSH3 = {"metal_ion": "metal", "metal_cluster": "cofactor", "heme_porphyrin": "cofactor", "cofactor": "cofactor"}


def build_table(limit=None):
    biolip = load_biolip_list()
    plinder = load_plinder_badlist()
    cof = load_cofactors()
    chembl = load_chembl()
    doc = gemmi.cif.read(str(raw_path("ccd_components")))
    rows = []
    for n, b in enumerate(doc):
        if limit and n >= limit:
            break
        cc = b.get_mmcif_category("_chem_comp.")
        g = lambda k: (cc.get(k) or [None])[0] or None
        code = b.name
        els = [e.upper().capitalize() for e in (b.find_values("_chem_comp_atom.type_symbol") or [])]
        heavy = [e for e in els if e not in ("H", "D")]
        try:
            charge = int(g("pdbx_formal_charge"))
        except (TypeError, ValueError):
            charge = None
        r = dict(ccd=code, name=g("name"), type=(g("type") or "").upper(), formula=g("formula"),
                 release_status=g("pdbx_release_status"), formal_charge=charge,
                 n_heavy=len(heavy), n_carbon=heavy.count("C"), elements=sorted(set(heavy)),
                 has_unknown_atom=any(e in ("X", "") for e in heavy) or len(els) == 0,
                 cofactor_class=cof.get(code))
        ids, ph = chembl.get(code, ([], None))
        r["chembl_ids"] = ids
        r["chembl_max_phase"] = ph
        r["is_approved_drug"] = ph is not None and ph >= 4
        smi = _smiles(b)
        r["smiles"] = smi
        r.update(rdkit_props(smi))
        big = r["n_heavy"] >= DRUGLIKE_MIN_HEAVY_ATOMS
        r["is_druglike_chembl"] = bool(ph is not None and ph >= DRUGLIKE_MIN_CHEMBL_PHASE and big)
        r["is_druglike_ro5qed"] = bool(r.get("ro5_pass") and (r.get("qed") or 0) >= DRUGLIKE_QED_MIN and big)
        r["druglike_source"] = ("chembl+ro5qed" if r["is_druglike_chembl"] and r["is_druglike_ro5qed"] else
                                "chembl" if r["is_druglike_chembl"] else
                                "ro5qed_FALLBACK" if r["is_druglike_ro5qed"] else None)
        r["ligand_class"] = classify(r)
        r["marsh_class"] = MARSH3.get(r["ligand_class"], "small_molecule")
        r["in_biolip_artifact_list"] = code in biolip
        is_ion = len(heavy) == 1 and heavy[0] not in ("C", "H", "N", "O", "P", "S")
        generic = (not is_ion) and (r["n_carbon"] < PLINDER_MIN_CARBONS or r["n_heavy"] < PLINDER_MIN_HEAVY_ATOMS
                                    or (charge is not None and abs(charge) > PLINDER_MAX_ABS_CHARGE)
                                    or bool(r.get("has_long_unbranched_hc")) or r["has_unknown_atom"])
        r["plinder_badlist"] = code in plinder
        r["plinder_generic_filter"] = generic
        r["plinder_artifact"] = r["plinder_badlist"] or generic
        rows.append(r)
    return pd.DataFrame(rows)
