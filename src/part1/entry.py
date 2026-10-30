"""Per-entry processing: selection -> assembly -> chain identity -> ligand instances -> contacts.

process_entry(path, sifts) returns (entry_row: dict, ligand_rows: list[dict]). Never raises: errors are
recorded in entry_row['skip_reason'] = 'error' and entry_row['error'].
First-attempt implementation; validated only by tests/test_part1.py known-answer cases.
"""
import hashlib
import re
import time
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree

from . import structure as S
from .config import (ASSEMBLY_ID, C2_ANGLE_TOL, DEDUP_CENTROID_MAX, EM_MAX_RESOLUTION, EM_METHODS, INCLUDE_BRANCHED, MAX_ASSEMBLY_CHAINS,
                     MBS_STRICT_LEVELS, MIN_ATOM_CONTACTS_SENS, MIN_ATOM_CONTACTS_PRIMARY, MIN_BOUND_RESIDUES,
                     MIN_CA_FOR_AXIS, ON_AXIS_JACCARD_MIN, DEDUP_ATOM_OVERLAP_DIST, DEDUP_ATOM_OVERLAP_FRAC,
                     SITE_JACCARD_MIN, ON_AXIS_GEOM_MAX_DIST, HOMOMER_MIN_IDENTITY, HOMOMER_MIN_COVERAGE,
                     MIN_PROTEIN_CHAIN_LEN, MIN_RES_PER_CHAIN_MBS, XRAY_MAX_RESOLUTION, XRAY_METHODS)
from .contacts import atom_pairs, radii, residue_contacts

ANTIBODY_RE = re.compile(r"immunoglobulin|\bfab\b|antibody|nanobody|\bvhh\b|\bscfv\b|\bfv\b", re.I)


def _uf_find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def select_entry(meta):
    """Return (method_class, resolution, skip_reason or None)."""
    methods = meta["methods"].split(";") if meta["methods"] else []
    if any(m in XRAY_METHODS for m in methods):
        res = meta["resolution_xray"]
        if res is None:
            return "xray", None, "no_resolution"
        return "xray", res, None if res <= XRAY_MAX_RESOLUTION else "resolution"
    if any(m in EM_METHODS for m in methods):
        res = meta["resolution_em"]
        if res is None:
            return "em", None, "no_resolution"
        return "em", res, None if res <= EM_MAX_RESOLUTION else "resolution"
    return "other", None, "method"


def pisa_agreement(asm):
    """True if assembly 1 is author_and_software_defined, or an author-only assembly 1 has a software-defined
    assembly with the same oligomeric count; None if no software-defined assembly exists."""
    if asm is None:
        return None
    if asm.details and "author_and_software" in asm.details:
        return True
    soft = [c for (_, d, c) in asm.all_assemblies if d and "software" in d]
    if not soft:
        return None
    if asm.details and asm.details.startswith("software"):
        return None  # assembly 1 itself is software-only: no author assembly to compare
    return asm.oligomeric_count in soft


def unp_segments(pdb, auth_chain, sifts, fallback):
    segs = sifts.get((pdb, auth_chain)) if sifts else None
    if segs:
        return segs, "sifts"
    if fallback.get(auth_chain):
        return [(a, b, e, db) for (a, b, e, db) in fallback[auth_chain]], "struct_ref"
    return [], None


def unp_pos(segs, label_seq):
    if label_seq is None:
        return None, None
    for acc, rb, re_, sb in segs:
        if rb <= label_seq <= re_:
            return acc, sb + (label_seq - rb)
    return None, None


def process_entry(path, sifts=None, pdb_id=None):
    t0 = time.time()
    pdb = pdb_id or path.name.split(".")[0].lower()
    row = dict(pdb_id=pdb, skip_reason=None, error=None)
    try:
        ligs = _process(path, pdb, sifts, row)
    except Exception as e:  # recorded, never fatal
        row["skip_reason"] = "error"
        row["error"] = f"{type(e).__name__}: {e}"[:300]
        ligs = []
    row["runtime_s"] = round(time.time() - t0, 3)
    row["n_ligand_rows"] = len(ligs)
    return row, ligs


def _process(path, pdb, sifts, row):
    block = S.read_block(path)
    meta = S.entry_metadata(block)
    row.update(meta)
    mclass, res, skip = select_entry(meta)
    row.update(method_class=mclass, resolution=res)
    if skip:
        row["skip_reason"] = skip
        return []
    ents = S.entities(block)
    a2e = S.asym_to_entity(block)
    prot_ents = {e for e, E in ents.items() if E.is_protein and len(E.seq) >= MIN_PROTEIN_CHAIN_LEN}
    row["n_protein_chains_au"] = sum(1 for a, e in a2e.items() if e in prot_ents)
    if not prot_ents:
        row["skip_reason"] = "no_protein_chain_ge30"
        return []
    asm = S.assembly(block, ASSEMBLY_ID)
    if asm is None or not asm.gens:
        asm_gens = [(list(a2e), [("1", np.eye(4))])]
        row.update(assembly_id="AU", assembly_details=None, assembly_oligomeric_count=None, n_assemblies=0)
    else:
        asm_gens = asm.gens
        row.update(assembly_id=asm.id, assembly_details=asm.details, assembly_oligomeric_count=asm.oligomeric_count,
                   n_assemblies=len(asm.all_assemblies))
    row["pisa_agrees"] = pisa_agreement(asm)
    poly_ents = {e for e, E in ents.items() if E.type == "polymer"}
    n_poly = sum(len(ops) for asyms, ops in asm_gens for a in asyms if a2e.get(a) in poly_ents)
    row["n_polymer_chains_asm"] = n_poly
    if n_poly > MAX_ASSEMBLY_CHAINS:
        row["skip_reason"] = "huge_assembly"
        return []

    at = S.atoms(block)
    row["n_altloc_ids"] = at["n_altloc_ids"]
    by_asym, by_asym_all = defaultdict(list), defaultdict(list)
    for i, (a, fa) in enumerate(zip(at["label_asym_id"], at["first_alt"])):
        by_asym_all[a].append(i)
        if fa:
            by_asym[a].append(i)
    by_asym = {a: np.array(v) for a, v in by_asym.items()}
    auth_of_asym = {a: at["auth_asym_id"][v[0]] for a, v in by_asym.items()}
    fallback = S.struct_ref_uniprot(block)
    cc_types = S.chem_comp_types(block)

    # ---------------- polymer chain instances in the assembly
    chains = []  # dicts
    res_rows = []  # residue table
    poly_xyz, poly_el, poly_res = [], [], []
    seen_chain_ids = set()
    for asyms, ops in asm_gens:
        for opname, M in ops:
            ident = np.allclose(M, np.eye(4))
            for a in asyms:
                e = a2e.get(a)
                if e not in poly_ents or a not in by_asym:
                    continue
                E = ents[e]
                auth = auth_of_asym[a]
                cid = auth if ident else f"{auth}-{opname}"
                while cid in seen_chain_ids:
                    cid += "'"
                seen_chain_ids.add(cid)
                kind = ("protein" if e in prot_ents else "peptide") if E.is_protein else ("nucleic" if E.is_nucleic else "other")
                segs, src = unp_segments(pdb, auth, sifts, fallback) if E.is_protein else ([], None)
                accs = {}
                for acc, rb, re_, sb in segs:
                    accs[acc] = accs.get(acc, 0) + (re_ - rb + 1)
                primary = max(accs, key=accs.get) if accs else None
                ci = len(chains)
                chains.append(dict(chain=cid, auth_chain=auth, label_asym=a, op=opname, entity=e, kind=kind,
                                   uniprot=primary, uniprot_src=src, chimeric=len(accs) > 1,
                                   key=primary or ("seq:" + hashlib.md5(E.seq.encode()).hexdigest()[:12]),
                                   segs=segs,
                                   antibody=bool(ANTIBODY_RE.search(E.description or ""))))
                idx = by_asym[a]
                xyz = at["xyz"][idx] @ M[:3, :3].T + M[:3, 3]
                # residue index within this chain
                lseq = at["label_seq_id"][idx]
                aseq = at["auth_seq_id"][idx]
                icode = at["pdbx_PDB_ins_code"][idx]
                comp = at["label_comp_id"][idx]
                keys = list(zip(lseq, aseq, icode))
                rmap = {}
                ridx = np.empty(len(idx), dtype=int)
                for k, key in enumerate(keys):
                    if key not in rmap:
                        rmap[key] = len(res_rows)
                        ls = int(key[0]) if key[0] not in (".", "?") else None
                        acc, up = unp_pos(segs, ls)
                        res_rows.append(dict(ci=ci, label_seq=ls,
                                             auth_seq=key[1] + ("" if key[2] in (".", "?") else key[2]),
                                             comp=comp[k], unp_acc=acc, unp_pos=up))
                    ridx[k] = rmap[key]
                if kind == "protein":
                    cam = at["label_atom_id"][idx] == "CA"
                    chains[ci]["ca"] = {int(l): x for l, x in zip(lseq[cam], xyz[cam]) if l not in (".", "?")}
                poly_xyz.append(xyz)
                poly_el.append(at["element"][idx])
                poly_res.append(ridx)
    prot_chains = [c for c in chains if c["kind"] == "protein"]
    row["n_protein_chains_asm"] = len(prot_chains)
    row["n_peptide_chains_asm"] = sum(c["kind"] == "peptide" for c in chains)
    row["n_nucleic_chains_asm"] = sum(c["kind"] == "nucleic" for c in chains)
    row["has_nucleic"] = row["n_nucleic_chains_asm"] > 0
    keys = sorted({c["key"] for c in prot_chains})
    row["n_distinct_proteins"] = len(keys)
    row["uniprot_accs"] = sorted({c["uniprot"] for c in prot_chains if c["uniprot"]})
    row["n_chains_without_uniprot"] = sum(c["uniprot"] is None for c in prot_chains)
    row["has_antibody"] = any(c["antibody"] for c in chains)
    row["any_chimeric"] = any(c["chimeric"] for c in prot_chains)
    if not prot_chains:
        row["skip_reason"] = "no_protein_in_assembly"
        return []
    atype, multi_ent = assembly_type(prot_chains, ents)
    row["assembly_type"] = atype
    row["has_same_uniprot_multi_entity"] = multi_ent

    if poly_xyz:
        P_xyz = np.concatenate(poly_xyz)
        P_el = np.concatenate(poly_el)
        P_res = np.concatenate(poly_res)
    else:
        P_xyz, P_el, P_res = np.zeros((0, 3)), np.zeros(0, dtype=object), np.zeros(0, dtype=int)
    P_r = radii(P_el)
    tree = cKDTree(P_xyz)
    max_pr = float(P_r.max()) if len(P_r) else 0.0

    # ---------------- ligand groups (AU level): residues of non-polymer (+branched) entities, merged by covalent links
    lig_ent_types = {"non-polymer"} | ({"branched"} if INCLUDE_BRANCHED else set())
    lig_res = {}  # (asym, auth_seq) -> comp
    for a, idx in by_asym.items():
        E = ents.get(a2e.get(a))
        if E is None or E.type not in lig_ent_types:
            continue
        for k in idx:
            lig_res.setdefault((a, at["auth_seq_id"][k]), at["label_comp_id"][k])
    parent = {k: k for k in lig_res}
    for a, idx in by_asym.items():  # branched entity: whole subchain is one ligand
        E = ents.get(a2e.get(a))
        if E is not None and E.type == "branched":
            ks = [k for k in lig_res if k[0] == a]
            for k in ks[1:]:
                parent[_uf_find(parent, k)] = _uf_find(parent, ks[0])
    cov_partner = defaultdict(set)  # lig residue -> protein residue comps covalently linked
    for (a1, s1, c1), (a2, s2, c2) in S.covalent_links(block):
        k1, k2 = (a1, s1), (a2, s2)
        if k1 in lig_res and k2 in lig_res:
            parent[_uf_find(parent, k1)] = _uf_find(parent, k2)
        elif k1 in lig_res and a2e.get(a2) in poly_ents:
            cov_partner[k1].add(c2)
        elif k2 in lig_res and a2e.get(a1) in poly_ents:
            cov_partner[k2].add(c1)
    groups = defaultdict(list)
    for k in lig_res:
        groups[_uf_find(parent, k)].append(k)
    au_count = defaultdict(int)
    for k, c in lig_res.items():
        au_count[c] += 1
    loi = S.ligand_of_interest(block)
    row["n_loi_features"] = len(loi)

    # atom indices per ligand residue, all altlocs: (atom index, alt id)
    lig_atoms = defaultdict(list)
    for a, idx in by_asym_all.items():
        E = ents.get(a2e.get(a))
        if E is None or E.type not in lig_ent_types:
            continue
        for k in idx:
            lig_atoms[(a, at["auth_seq_id"][k])].append(k)

    asym_ops = defaultdict(list)
    for asyms, ops in asm_gens:
        for a in asyms:
            for opname, M in ops:
                asym_ops[a].append((opname, M))

    chain_by_id = {c["chain"]: c for c in chains}

    def get_contacts(xyz, idx):
        lr = radii(at["element"][idx])
        li, pj, d, isv, iss = atom_pairs(xyz, lr, tree, P_r, max_pr)
        rc = residue_contacts(pj, d, isv, iss, P_res)
        out = []
        for ri, (nv, ns, md) in sorted(rc.items()):
            rr = res_rows[ri]
            c = chains[rr["ci"]]
            out.append(dict(chain=c["chain"], auth_chain=c["auth_chain"], label_asym=c["label_asym"], op=c["op"],
                            entity=c["entity"], kind=c["kind"], label_seq=rr["label_seq"], auth_seq=rr["auth_seq"],
                            comp=rr["comp"], unp_acc=rr["unp_acc"], unp_pos=rr["unp_pos"],
                            n_vdw=nv, n_simple=ns, min_dist=round(md, 2)))
        return out

    # ---- pass 1: all instances (group x operator x alt conformer)
    inst = []
    for root, members in groups.items():
        members = sorted(members, key=lambda k: (k[0], k[1]))
        comps = [lig_res[m] for m in members]
        main = max(members, key=lambda m: len(lig_atoms[m]))
        all_idx = np.array([k for m in members for k in lig_atoms[m]])
        if len(all_idx) == 0:
            continue
        alts_here = sorted({x for x in at["alt"][all_idx] if x not in (".", "?")})
        variants = alts_here if alts_here else [None]
        info = dict(members=members, comps=comps, main=main, ccd_main=lig_res[main], ccd_key="+".join(sorted(set(comps))),
                    is_branched=any(ents[a2e[m[0]]].type == "branched" for m in members),
                    is_sugar=any("SACCHARIDE" in cc_types.get(c, "") for c in comps),
                    cov=sorted(set().union(*[cov_partner[m] for m in members])))
        for opname, M in asym_ops.get(main[0], []):
            first_vi = None
            for v in variants:
                alt = at["alt"][all_idx]
                idx = all_idx[(alt == ".") | (alt == "?") | (alt == v)] if v else all_idx
                xyz = at["xyz"][idx] @ M[:3, :3].T + M[:3, 3]
                inst.append(dict(info=info, op=opname, ident=np.allclose(M, np.eye(4)), alt=v, idx=idx, xyz=xyz,
                                 cen=xyz.mean(axis=0), contacts=get_contacts(xyz, idx),
                                 alt_parent=first_vi))
                if first_vi is None:
                    first_vi = len(inst) - 1

    # ---- dedup: alt conformers -> first conformer; same-CCD copies overlapping (centroid or atom overlap) -> first copy
    reps = []
    for i, x in enumerate(inst):
        x["dup_of"], x["dup_reason"] = None, None
        if x["alt_parent"] is not None:
            x["dup_of"], x["dup_reason"] = x["alt_parent"], "altloc"
            continue
        for r in reps:
            y = inst[r]
            if y["info"]["ccd_key"] != x["info"]["ccd_key"]:
                continue
            if np.linalg.norm(y["cen"] - x["cen"]) < DEDUP_CENTROID_MAX:
                x["dup_of"], x["dup_reason"] = r, "centroid"
                break
            dm = np.linalg.norm(x["xyz"][:, None, :] - y["xyz"][None, :, :], axis=2)
            if (dm.min(axis=1) <= DEDUP_ATOM_OVERLAP_DIST).mean() > DEDUP_ATOM_OVERLAP_FRAC:
                x["dup_of"], x["dup_reason"] = r, "atom_overlap"
                break
        if x["dup_of"] is None:
            reps.append(i)
    # resolve chains of duplicates to their representative and merge contacts into it
    for x in inst:
        r = x["dup_of"]
        while r is not None and inst[r]["dup_of"] is not None:
            r = inst[r]["dup_of"]
        x["rep"] = r
    merged = {i: {} for i in reps}
    for i, x in enumerate(inst):
        tgt = i if x["rep"] is None else x["rep"]
        for c in x["contacts"]:
            k = (c["chain"], c["label_seq"], c["auth_seq"])
            m = merged[tgt].get(k)
            if m is None:
                merged[tgt][k] = dict(c)
            else:
                m["n_vdw"] = max(m["n_vdw"], c["n_vdw"])
                m["n_simple"] = max(m["n_simple"], c["n_simple"])
                m["min_dist"] = min(m["min_dist"], c["min_dist"])
    n_copies = defaultdict(int)
    for x in inst:
        if x["rep"] is not None:
            n_copies[x["rep"]] += 1

    # ---- pass 2: rows (representatives carry the union of their copies' contacts)
    ligs = []
    for i, x in enumerate(inst):
        f = x["info"]
        contacts = sorted(merged[i].values(), key=lambda c: (c["chain"], c["label_seq"] or -1, c["auth_seq"])) \
            if x["rep"] is None else x["contacts"]
        is_glycan = f["is_branched"] or (f["is_sugar"] and bool(f["cov"]))
        lrow = _ligand_row(pdb, row, f["members"], f["comps"], f["ccd_main"], f["ccd_key"], x["op"], x["ident"],
                           len(x["idx"]), contacts, chain_by_id, f["is_branched"], f["is_sugar"], is_glycan,
                           bool(f["cov"]), f["cov"], au_count, loi, x["cen"])
        if x["alt"]:
            lrow["ligand_uid"] += f":alt{x['alt']}"
        lrow["altloc"] = x["alt"]
        lrow["main_auth_chain"] = auth_of_asym[f["main"][0]]
        lrow["main_auth_seq"] = f["main"][1]
        lrow["is_dup_copy"] = x["rep"] is not None
        lrow["dup_reason"] = x["dup_reason"]
        lrow["n_merged_copies"] = n_copies[i] if x["rep"] is None else 0
        ligs.append(lrow)
    for i, x in enumerate(inst):
        if x["rep"] is not None:
            ligs[i]["dup_of"] = ligs[x["rep"]]["ligand_uid"]
        else:
            ligs[i]["dup_of"] = None
    # ---- symmetric-site groups: union-find over representatives with same CCD and residue-set Jaccard >= threshold
    par = list(range(len(ligs)))
    rep_idx = [i for i in range(len(ligs)) if not ligs[i]["is_dup_copy"]]
    for a_i, i in enumerate(rep_idx):
        si = ligs[i]["_symset"]
        if not si:
            continue
        for j in rep_idx[a_i + 1:]:
            if ligs[j]["ccd_key"] != ligs[i]["ccd_key"]:
                continue
            sj = ligs[j]["_symset"]
            if sj and len(si & sj) / len(si | sj) >= SITE_JACCARD_MIN:
                par[_uf_find(par, j)] = _uf_find(par, i)
    for i, l in enumerate(ligs):
        root = _uf_find(par, i if not l["is_dup_copy"] else inst[i]["rep"])
        l["site_group"] = ligs[root]["ligand_uid"]
        del l["_symset"]
    row["n_ligand_instances_asm"] = len(ligs)
    return ligs


def _unp_aa(chain, seq):
    """UniProt position -> one-letter residue for a chain, via its SIFTS / struct_ref segments."""
    out = {}
    for acc, rb, re_, sb in chain["segs"]:
        if acc != chain["uniprot"]:
            continue
        for ls in range(rb, min(re_, len(seq)) + 1):
            out[sb + ls - rb] = seq[ls - 1]
    return out


def entities_equivalent(c1, c2, ents):
    """Same entity, or same UniProt with identity >= HOMOMER_MIN_IDENTITY over >= HOMOMER_MIN_COVERAGE of the shorter chain."""
    if c1["entity"] == c2["entity"]:
        return True
    if c1["key"] != c2["key"] or c1["uniprot"] is None:
        return False
    s1, s2 = ents[c1["entity"]].seq, ents[c2["entity"]].seq
    a1, a2 = _unp_aa(c1, s1), _unp_aa(c2, s2)
    common = set(a1) & set(a2)
    if not common:
        return False
    # coverage is required over BOTH chains (i.e. over the longer one): with shorter-chain coverage only, a processed
    # prefix such as HIV RT p51 (subset of p66) would be called a homomer, contrary to the intended typing (PI v2 item 6)
    cov = len(common) / max(1, max(len(s1), len(s2)))
    ident = sum(a1[p] == a2[p] for p in common) / len(common)
    return cov >= HOMOMER_MIN_COVERAGE and ident >= HOMOMER_MIN_IDENTITY


def assembly_type(prot_chains, ents):
    """monomer / homomer / heteromer / same_uniprot_hetero (v2). Returns (type, has_same_uniprot_multi_entity)."""
    by_ent = {}
    for c in prot_chains:
        by_ent.setdefault(c["entity"], c)
    reps = list(by_ent.values())
    par = list(range(len(reps)))
    for i in range(len(reps)):
        for j in range(i + 1, len(reps)):
            if entities_equivalent(reps[i], reps[j], ents):
                par[_uf_find(par, j)] = _uf_find(par, i)
    comps_by_key = defaultdict(set)
    for i, c in enumerate(reps):
        comps_by_key[c["key"]].add(_uf_find(par, i))
    multi = any(len(v) > 1 for v in comps_by_key.values())
    if len(prot_chains) == 1:
        return "monomer", False
    if len(comps_by_key) > 1:
        return "heteromer", multi
    return ("same_uniprot_hetero" if multi else "homomer"), multi


def kabsch(X, Y):
    """R, t minimising |R X + t - Y|."""
    cx, cy = X.mean(0), Y.mean(0)
    H = (X - cx).T @ (Y - cy)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, cy - R @ cx


def c2_axis_distance(ca1, ca2, cen):
    """Superpose chain 1 onto chain 2 (common label_seq CAs). Returns (rotation angle deg, distance of cen to the
    rotation axis if the angle is within C2_ANGLE_TOL of 180 deg else None)."""
    common = sorted(set(ca1) & set(ca2))
    if len(common) < MIN_CA_FOR_AXIS:
        return None, None
    X = np.array([ca1[k] for k in common])
    Y = np.array([ca2[k] for k in common])
    R, t = kabsch(X, Y)
    ang = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
    if abs(ang - 180.0) > C2_ANGLE_TOL:
        return round(ang, 1), None
    w, v = np.linalg.eig(R)
    u = np.real(v[:, np.argmin(np.abs(w - 1))])
    u /= np.linalg.norm(u)
    s = float(u @ t)
    p = np.linalg.lstsq(np.eye(3) - R, t - s * u, rcond=None)[0]
    dvec = (np.asarray(cen) - p)
    dist = float(np.linalg.norm(dvec - (dvec @ u) * u))
    return round(ang, 1), round(dist, 2)


def _per_chain(contacts, kinds, min_atoms, simple=False, exclude_antibody=False, chain_by_id=None):
    out = defaultdict(int)
    for c in contacts:
        if c["kind"] not in kinds:
            continue
        if exclude_antibody and chain_by_id[c["chain"]]["antibody"]:
            continue
        ok = (c["n_simple"] >= 1) if simple else (c["n_vdw"] >= min_atoms)
        if ok:
            out[c["chain"]] += 1
    return dict(out)


def _ligand_row(pdb, erow, members, comps, ccd_main, ccd_key, opname, ident, n_heavy, contacts, chain_by_id,
                is_branched, is_sugar, is_glycan, is_cov, cov, au_count, loi, cen):
    prim = _per_chain(contacts, ("protein",), MIN_ATOM_CONTACTS_PRIMARY)
    b2 = _per_chain(contacts, ("protein",), MIN_ATOM_CONTACTS_SENS)
    simp = _per_chain(contacts, ("protein",), 1, simple=True)
    noab = _per_chain(contacts, ("protein",), MIN_ATOM_CONTACTS_PRIMARY, exclude_antibody=True, chain_by_id=chain_by_id)
    pep = _per_chain(contacts, ("peptide",), MIN_ATOM_CONTACTS_PRIMARY)
    nuc = _per_chain(contacts, ("nucleic",), MIN_ATOM_CONTACTS_PRIMARY)

    def mbs(d, k=MIN_RES_PER_CHAIN_MBS):
        return sum(v >= k for v in d.values()) >= 2

    n_prim = sum(prim.values())
    contrib = [chain_by_id[c] for c, v in prim.items() if v >= MIN_RES_PER_CHAIN_MBS]
    # BioLiP2 geometric artifact rule components (sensitivity)
    pres = [c for c in contacts if c["kind"] == "protein" and c["n_vdw"] >= MIN_ATOM_CONTACTS_PRIMARY]
    two_consec = (len(pres) == 2 and pres[0]["chain"] == pres[1]["chain"] and pres[0]["label_seq"] is not None
                  and pres[1]["label_seq"] is not None and abs(pres[0]["label_seq"] - pres[1]["label_seq"]) == 1)
    # symmetry-equivalence key: CCD + multiset of (protein identity, residue position) ignoring chain labels
    sym_items = sorted((chain_by_id[c["chain"]]["key"], c["unp_pos"] if c["unp_pos"] is not None else c["label_seq"] or -1)
                       for c in pres)
    sym_key = ccd_key + "|" + ";".join(f"{k}:{p}" for k, p in sym_items)
    symset = set(sym_items)
    # on-axis homomer flag: two top contributing chains of the same protein contact the same residues
    on_axis_j, axis_ang, axis_dist = None, None, None
    top = sorted(contrib, key=lambda c: -prim[c["chain"]])[:2]
    if len(top) == 2 and top[0]["key"] == top[1]["key"]:
        sets = [{(c["unp_pos"] if c["unp_pos"] is not None else c["label_seq"]) for c in pres if c["chain"] == t["chain"]}
                for t in top]
        on_axis_j = round(len(sets[0] & sets[1]) / max(1, len(sets[0] | sets[1])), 3)
        if "ca" in top[0] and "ca" in top[1]:
            axis_ang, axis_dist = c2_axis_distance(top[0]["ca"], top[1]["ca"], cen)
    asyms = {m[0] for m in members}
    seqs = {m[1] for m in members}
    is_loi = any(c in comps and (a is None or a in asyms) and (s is None or s in seqs) for c, a, s in loi)
    lab = members[0]
    return dict(
        pdb_id=pdb,
        ligand_uid=f"{pdb}:{lab[0]}:{lab[1]}:{ccd_main}:{opname}",
        method_class=erow["method_class"], resolution=erow["resolution"], release_date=erow["release_date"],
        assembly_id=erow["assembly_id"], assembly_type=erow["assembly_type"],
        n_protein_chains=erow["n_protein_chains_asm"], has_nucleic=erow["has_nucleic"], has_antibody=erow["has_antibody"],
        ccd=ccd_main, ccd_key=ccd_key, n_components=len(members), components=comps,
        label_asym=lab[0], auth_seq=lab[1], op=opname, is_identity_op=ident, n_heavy_atoms_instance=n_heavy,
        centroid=[round(float(x), 2) for x in cen],
        is_branched=is_branched, is_sugar_comp=is_sugar, is_glycan=is_glycan,
        is_covalent=is_cov, covalent_partner_comps=cov,
        n_copies_in_au=max(au_count[c] for c in comps), biolip2_two_consecutive_site=two_consec,
        is_ligand_of_interest=is_loi,
        n_res_total=n_prim, n_res_total_1atom=sum(b2.values()), n_res_total_4A=sum(simp.values()),
        is_bound=n_prim >= MIN_BOUND_RESIDUES,
        is_bound_1atom=sum(b2.values()) >= MIN_BOUND_RESIDUES,
        is_bound_4A=sum(simp.values()) >= MIN_BOUND_RESIDUES,
        chain_counts=[dict(chain=c, n_res=v) for c, v in sorted(prim.items())],
        chain_counts_1atom=[dict(chain=c, n_res=v) for c, v in sorted(b2.items())],
        chain_counts_4A=[dict(chain=c, n_res=v) for c, v in sorted(simp.items())],
        n_chains_contacted=len(prim), n_chains_contacted_1atom=len(b2), n_chains_contacted_4A=len(simp),
        is_mbs=mbs(prim), is_mbs_1atom=mbs(b2), is_mbs_4A=mbs(simp), is_mbs_excl_antibody=mbs(noab),
        **{f"is_mbs_min{k}": mbs(prim, k) for k in MBS_STRICT_LEVELS},
        mbs_contrib_same_entity=(len({c["entity"] for c in contrib}) == 1) if len(contrib) >= 2 else None,
        mbs_contrib_same_protein=(len({c["key"] for c in contrib}) == 1) if len(contrib) >= 2 else None,
        mbs_contrib_symmetry_copy=(len({c["label_asym"] for c in contrib}) < len(contrib)) if len(contrib) >= 2 else None,
        on_axis_jaccard=on_axis_j, is_on_axis=(on_axis_j >= ON_AXIS_JACCARD_MIN) if on_axis_j is not None else None,
        chainpair_rot_angle=axis_ang, c2_axis_dist=axis_dist,
        is_on_axis_geom=(axis_dist is not None and axis_dist <= ON_AXIS_GEOM_MAX_DIST) if len(top) == 2 else None,
        _symset=symset,
        dominant_chain=max(prim, key=prim.get) if prim else None,
        dominant_entity=chain_by_id[max(prim, key=prim.get)]["entity"] if prim else None,
        dominant_uniprot=chain_by_id[max(prim, key=prim.get)]["uniprot"] if prim else None,
        dominant_protein_key=chain_by_id[max(prim, key=prim.get)]["key"] if prim else None,
        contacts_antibody=any(chain_by_id[c]["antibody"] for c in prim),
        n_res_peptide=sum(pep.values()), n_res_nucleic=sum(nuc.values()),
        sym_group=pdb + ":" + hashlib.sha1(sym_key.encode()).hexdigest()[:12],
        contacts=contacts,
    )
