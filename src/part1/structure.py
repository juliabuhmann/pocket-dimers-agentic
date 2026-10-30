"""mmCIF parsing for the Part 1 census: entry metadata, atoms (numpy), entities, assembly operators.

Everything is read directly from the mmCIF block with gemmi.cif (no gemmi.Structure) so that the
provenance of each assembly copy (label_asym_id, operator) is kept explicitly.
"""
import re
from dataclasses import dataclass, field
from itertools import product

import gemmi
import numpy as np

from .config import MODEL_INDEX


def cat(block, name):
    """mmCIF category as dict of lists; '?' -> None and '.' -> None. Missing category -> {}."""
    d = block.get_mmcif_category(name)
    if not d:
        return {}
    return {k: [None if (v is None or v is False) else v for v in vals] for k, vals in d.items()}


def first(d, key):
    vals = d.get(key) or []
    for v in vals:
        if v is not None:
            return v
    return None


def to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------ metadata
def entry_metadata(block) -> dict:
    exptl = cat(block, "_exptl.")
    methods = sorted({m.upper() for m in exptl.get("method", []) if m})
    refine = cat(block, "_refine.")
    reflns = cat(block, "_reflns.")
    em = cat(block, "_em_3d_reconstruction.")
    res_xray = to_float(first(refine, "ls_d_res_high")) or to_float(first(reflns, "d_resolution_high"))
    res_em = None
    for v in em.get("resolution", []) or []:
        f = to_float(v)
        if f is not None:
            res_em = f if res_em is None else min(res_em, f)
    rev = cat(block, "_pdbx_audit_revision_history.")
    dates = sorted(d for d in rev.get("revision_date", []) if d)
    status = cat(block, "_pdbx_database_status.")
    kw = cat(block, "_struct_keywords.")
    kw_text = " ".join(x for x in (first(kw, "pdbx_keywords"), first(kw, "text")) if x)
    sg = cat(block, "_pdbx_SG_project.")
    details = cat(block, "_pdbx_entry_details.")
    return dict(
        methods=";".join(methods),
        resolution_xray=res_xray,
        resolution_em=res_em,
        release_date=dates[0] if dates else None,
        deposition_date=first(status, "recvd_initial_deposition_date"),
        keywords=kw_text or None,
        is_structural_genomics=bool(sg) or ("STRUCTURAL GENOMICS" in kw_text.upper()),
        has_ligand_of_interest=first(details, "has_ligand_of_interest"),
    )


# ------------------------------------------------------------------ entities
@dataclass
class Entity:
    id: str
    type: str                 # polymer / non-polymer / branched / water / macrolide
    poly_type: str = None     # e.g. polypeptide(L), polydeoxyribonucleotide
    seq: str = ""
    description: str = ""

    @property
    def is_protein(self):
        return self.type == "polymer" and (self.poly_type or "").startswith("polypeptide")

    @property
    def is_nucleic(self):
        return self.type == "polymer" and "nucleotide" in (self.poly_type or "")


def entities(block) -> dict:
    ent = cat(block, "_entity.")
    out = {}
    for i, eid in enumerate(ent.get("id", [])):
        out[eid] = Entity(eid, (ent["type"][i] or "").lower(), description=(ent.get("pdbx_description") or [None] * (i + 1))[i] or "")
    poly = cat(block, "_entity_poly.")
    for i, eid in enumerate(poly.get("entity_id", [])):
        if eid in out:
            out[eid].poly_type = poly["type"][i]
            seq = (poly.get("pdbx_seq_one_letter_code_can") or [None] * (i + 1))[i] or ""
            out[eid].seq = re.sub(r"\s+", "", seq)
    # fall back to entity_poly_seq length when canonical sequence is missing
    eps = cat(block, "_entity_poly_seq.")
    if eps:
        counts = {}
        for eid in eps.get("entity_id", []):
            counts[eid] = counts.get(eid, 0) + 1
        for eid, n in counts.items():
            if eid in out and not out[eid].seq:
                out[eid].seq = "X" * n
    return out


def asym_to_entity(block) -> dict:
    sa = cat(block, "_struct_asym.")
    return dict(zip(sa.get("id", []), sa.get("entity_id", [])))


def chem_comp_types(block) -> dict:
    cc = cat(block, "_chem_comp.")
    return {k: (t or "").upper() for k, t in zip(cc.get("id", []), cc.get("type", []))}


def struct_ref_uniprot(block) -> dict:
    """Fallback UniProt mapping from the entry: auth chain -> list of (acc, seq_align_beg, seq_align_end, db_align_beg)."""
    sr = cat(block, "_struct_ref.")
    unp = {rid: acc for rid, db, acc in zip(sr.get("id", []), sr.get("db_name", []), sr.get("pdbx_db_accession", []))
           if db == "UNP" and acc}
    srs = cat(block, "_struct_ref_seq.")
    out = {}
    for i, rid in enumerate(srs.get("ref_id", [])):
        if rid not in unp:
            continue
        ch = srs["pdbx_strand_id"][i]
        try:
            seg = (unp[rid], int(srs["seq_align_beg"][i]), int(srs["seq_align_end"][i]), int(srs["db_align_beg"][i]))
        except (TypeError, ValueError, KeyError):
            continue
        out.setdefault(ch, []).append(seg)
    return out


def ligand_of_interest(block) -> list:
    """(comp_id, label_asym_id or None, auth_seq_num or None) for 'SUBJECT OF INVESTIGATION' features."""
    f = cat(block, "_pdbx_entity_instance_feature.")
    out = []
    for i, ft in enumerate(f.get("feature_type", [])):
        if ft and ft.upper() == "SUBJECT OF INVESTIGATION":
            out.append(((f.get("comp_id") or [None] * (i + 1))[i],
                        (f.get("asym_id") or [None] * (i + 1))[i],
                        (f.get("auth_seq_num") or [None] * (i + 1))[i]))
    return out


def covalent_links(block) -> list:
    """struct_conn rows of type covale*: ((asym, auth_seq, comp), (asym, auth_seq, comp))."""
    sc = cat(block, "_struct_conn.")
    out = []
    for i, t in enumerate(sc.get("conn_type_id", [])):
        if not t or not t.lower().startswith("covale"):
            continue
        a = (sc["ptnr1_label_asym_id"][i], sc["ptnr1_auth_seq_id"][i], sc["ptnr1_label_comp_id"][i])
        b = (sc["ptnr2_label_asym_id"][i], sc["ptnr2_auth_seq_id"][i], sc["ptnr2_label_comp_id"][i])
        out.append((a, b))
    return out


# ------------------------------------------------------------------ assemblies
@dataclass
class Assembly:
    id: str
    details: str
    oligomeric_count: str
    gens: list = field(default_factory=list)  # list of (list_of_asym_ids, list_of_(op_name, 4x4 matrix))
    all_assemblies: list = field(default_factory=list)  # [(id, details, oligomeric_count)]


def _oper_matrices(block) -> dict:
    ol = cat(block, "_pdbx_struct_oper_list.")
    mats = {}
    for i, oid in enumerate(ol.get("id", [])):
        m = np.eye(4)
        for r in range(3):
            for c in range(3):
                m[r, c] = float(ol[f"matrix[{r+1}][{c+1}]"][i])
            m[r, 3] = float(ol[f"vector[{r+1}]"][i])
        mats[oid] = m
    return mats


def _expand_group(g: str) -> list:
    ids = []
    for part in g.split(","):
        part = part.strip()
        m = re.fullmatch(r"(\d+)-(\d+)", part)
        if m:
            ids += [str(k) for k in range(int(m.group(1)), int(m.group(2)) + 1)]
        elif part:
            ids.append(part)
    return ids


def parse_oper_expression(expr: str) -> list:
    """'1', '1,2', '1-60', '(1-60)(61-88)', '(X0)(1-5)' -> list of operator-id tuples (left-to-right product)."""
    expr = expr.replace(" ", "").replace("'", "")
    groups = re.findall(r"\(([^)]*)\)", expr) if "(" in expr else [expr]
    return list(product(*[_expand_group(g) for g in groups]))


def assembly(block, assembly_id=None):
    asm = cat(block, "_pdbx_struct_assembly.")
    if not asm:
        return None
    ids = asm["id"]
    allas = list(zip(ids, asm.get("details", [None] * len(ids)), asm.get("oligomeric_count", [None] * len(ids))))
    aid = assembly_id or ids[0]
    k = ids.index(aid)
    a = Assembly(aid, (asm.get("details") or [None] * len(ids))[k], (asm.get("oligomeric_count") or [None] * len(ids))[k],
                 all_assemblies=allas)
    mats = _oper_matrices(block)
    gen = cat(block, "_pdbx_struct_assembly_gen.")
    for i, gid in enumerate(gen.get("assembly_id", [])):
        if gid != aid:
            continue
        asyms = [x.strip() for x in gen["asym_id_list"][i].split(",") if x.strip()]
        ops = []
        for combo in parse_oper_expression(gen["oper_expression"][i]):
            m = np.eye(4)
            for oid in combo:
                m = m @ mats[oid]
            ops.append(("x".join(combo), m))
        a.gens.append((asyms, ops))
    return a


# ------------------------------------------------------------------ atoms
ATOM_TAGS = ["group_PDB", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id", "label_asym_id",
             "label_entity_id", "label_seq_id", "auth_seq_id", "auth_asym_id", "pdbx_PDB_ins_code",
             "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num"]


def atoms(block) -> dict:
    """Heavy atoms of the first model, all altlocs (see first_alt). Returns dict of numpy arrays."""
    tab = block.find("_atom_site.", ATOM_TAGS)
    n = len(tab)
    cols = {t: np.array(list(tab.column(j)), dtype=object) for j, t in enumerate(ATOM_TAGS)}
    cols["label_atom_id"] = np.array([gemmi.cif.as_string(x) for x in cols["label_atom_id"]], dtype=object)
    model = cols["pdbx_PDB_model_num"]
    first_model = model[0] if n else None
    alt = cols["label_alt_id"]
    alts = sorted(set(alt) - {".", "?"})
    keep_alt = alts[0] if alts else None
    el = np.array([e.upper() for e in cols["type_symbol"]], dtype=object)
    # all altlocs are kept (ligand alt conformers become duplicate copies, v2); first_alt marks '.' + lowest alt id
    mask = (model == first_model) & (el != "H") & (el != "D")
    out = {t: cols[t][mask] for t in ATOM_TAGS if t not in ("Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num", "label_alt_id")}
    out["xyz"] = np.stack([cols[c][mask].astype(float) for c in ("Cartn_x", "Cartn_y", "Cartn_z")], axis=1) if mask.any() \
        else np.zeros((0, 3))
    out["element"] = np.array([e.capitalize() for e in el[mask]], dtype=object)
    out["n_altloc_ids"] = len(alts)
    out["alt"] = alt[mask]
    out["first_alt"] = (out["alt"] == ".") | (out["alt"] == "?") | (out["alt"] == keep_alt)
    return out


def read_block(path):
    return gemmi.cif.read(str(path)).sole_block()
