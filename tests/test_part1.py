"""Known-answer (sentinel) tests for the Part 1 census. Run automatically at the start of every census run.

Real-entry expectations were checked by inspecting the pipeline output against the mmCIF files and textbook
knowledge (see docs/PART1_PIPELINE.md "Validation"):
  1HVR HIV-1 protease + XK2 inhibitor: homodimer, one site on the C2 axis formed by both chains (classic MBS)
  5CFL anemone STING + c-di-GMP: Marsh 2018 Fig. 1 MBS example
  1STP streptavidin + biotin: biotin pocket in one subunit plus Trp120 from the neighbouring subunit
  3PTB trypsin + benzamidine: monomer, single-chain site; BioLiP2 drops BEN (listed additive), keeps CA
  4HHB deoxy-haemoglobin: heteromer (alpha2 beta2), each heme within one subunit; PO4 F has no protein contact
  3Q92 ketohexokinase (Marsh 2018 Fig. 1 SBS example): SO4 not in BioLiP2 -> additive
  1NEK succinate dehydrogenase: heteromer MBS heme b between SdhC and SdhD, ubiquinone site across B/C/D
  1HEW lysozyme + tri-NAG (branched entity) -> glycan flag
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from part1 import merge as MG  # noqa: E402
from part1.config import OUT_DIR, PDB_MIRROR, SIMPLE_CUTOFF, VDW_TOLERANCE  # noqa: E402
from part1.contacts import atom_pairs, radii, vdw_radius  # noqa: E402
from part1.entry import process_entry  # noqa: E402
from part1.structure import parse_oper_expression  # noqa: E402


def run(pid):
    return process_entry(PDB_MIRROR / pid[1:3] / f"{pid}.cif.gz", None, pid)


def by(ligs, **kw):
    return [l for l in ligs if all(l[k] == v for k, v in kw.items())]


# ---------------------------------------------------------------- toy / unit tests
def test_contact_thresholds():
    from scipy.spatial import cKDTree
    rc = vdw_radius("C")
    vdw_lim = 2 * rc + VDW_TOLERANCE
    prot = np.array([[0.0, 0, 0]])
    tree = cKDTree(prot)
    pr = radii(["C"])
    for d, exp_vdw, exp_simple in [(vdw_lim - 0.01, True, vdw_lim - 0.01 <= SIMPLE_CUTOFF),
                                   (vdw_lim + 0.01, False, vdw_lim + 0.01 <= SIMPLE_CUTOFF),
                                   (SIMPLE_CUTOFF + 0.05, (SIMPLE_CUTOFF + 0.05) <= vdw_lim, False)]:
        li, pj, dist, isv, iss = atom_pairs(np.array([[d, 0, 0]]), radii(["C"]), tree, pr, float(pr.max()))
        got_v = bool(isv[0]) if len(isv) else False
        got_s = bool(iss[0]) if len(iss) else False
        assert (got_v, got_s) == (exp_vdw, exp_simple), d


def test_oper_expression():
    assert parse_oper_expression("1") == [("1",)]
    assert parse_oper_expression("1,2") == [("1",), ("2",)]
    assert len(parse_oper_expression("(1-60)")) == 60
    assert parse_oper_expression("(1-2)(X0)") == [("1", "X0"), ("2", "X0")]


TOY_HEADER = """data_TOY1
_entry.id TOY1
_exptl.method 'X-RAY DIFFRACTION'
_refine.ls_d_res_high 2.0
_pdbx_audit_revision_history.revision_date 2000-01-01
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'toy protein'
2 non-polymer 'toy ligand'
3 non-polymer 'toy ligand 2'
_entity_poly.entity_id 1
_entity_poly.type 'polypeptide(L)'
_entity_poly.pdbx_seq_one_letter_code_can {seq}
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
C 2
D 2
E 2
F 3
loop_
_chem_comp.id
_chem_comp.type
GLY 'PEPTIDE LINKING'
LIG NON-POLYMER
LG2 NON-POLYMER
_pdbx_struct_assembly.id 1
_pdbx_struct_assembly.details author_defined_assembly
_pdbx_struct_assembly.oligomeric_count 2
_pdbx_struct_assembly_gen.assembly_id 1
_pdbx_struct_assembly_gen.oper_expression '1,2'
_pdbx_struct_assembly_gen.asym_id_list A,C,D,E,F
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 0
2 'twofold about x' 1 0 0 0 0 -1 0 0 0 0 -1 0
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
"""


def make_toy(tmp_path):
    """Chain A: 30 Gly CA atoms at (10*(i-1), 3, i%2) (zig-zag so the chain superposition is defined); assembly op 2 (2-fold about x) makes chain B at y=-3.
    Ligands (single C atoms unless noted):
      C at origin       -> contacts A1 and B1 (3.0 A each): bound, MBS; its op-2 copy coincides -> duplicate
      D at (500,500,500)-> no contacts (negative control)
      E at (100,6,0)    -> only A11: 1 residue -> not bound
      F (LG2, 2 atoms) at (200,6,0),(210,6,0) -> A21 + A22: bound, single-chain site"""
    lines, k = [], 1
    for i in range(1, 31):
        lines.append(f"ATOM {k} C CA . GLY A 1 {i} ? {10.0 * (i - 1):.3f} 3.000 {float(i % 2):.3f} {i} A 1"); k += 1
    # every ligand position carries 2 atoms (z and z+0.5) so contacted residues have 2 atom pairs (v2 primary rule)
    for asym, xyz in [("C", (0, 0, 0)), ("D", (500, 500, 500)), ("E", (100, 6, 0))]:
        for j, dz in enumerate((0.0, 0.5)):
            lines.append(f"HETATM {k} C C{j + 1} . LIG {asym} 2 . ? {xyz[0]:.3f} {xyz[1]:.3f} {xyz[2] + dz:.3f} 101 {asym} 1"); k += 1
    for j, (x, dz) in enumerate([(200, 0.0), (200, 0.5), (210, 0.0), (210, 0.5)]):
        lines.append(f"HETATM {k} C C{j + 1} . LG2 F 3 . ? {x:.3f} 6.000 {dz:.3f} 101 F 1"); k += 1
    txt = TOY_HEADER.format(seq="G" * 30) + "\n".join(lines) + "\n"
    p = tmp_path / "toy1.cif"
    p.write_text(txt)
    return p


def test_toy_assembly(tmp_path):
    row, ligs = process_entry(make_toy(tmp_path), None, "toy1")
    assert row["skip_reason"] is None, row
    assert row["assembly_type"] == "homomer" and row["n_protein_chains_asm"] == 2
    c1 = by(ligs, label_asym="C", op="1")[0]
    c2 = by(ligs, label_asym="C", op="2")[0]
    assert c1["is_bound"] and c1["is_mbs"] and c1["n_chains_contacted"] == 2 and not c1["is_dup_copy"]
    assert c1["mbs_contrib_symmetry_copy"] is True  # A and its op-2 copy
    assert c2["is_dup_copy"] and c2["dup_of"] == c1["ligand_uid"]
    assert c1["is_on_axis"] and c1["is_on_axis_geom"] and c1["c2_axis_dist"] < 0.5
    for l in by(ligs, label_asym="D"):
        assert l["n_res_total"] == 0 and not l["is_bound"]
    for l in by(ligs, label_asym="E"):
        assert l["n_res_total"] == 1 and not l["is_bound"]
    f = by(ligs, label_asym="F", op="1")[0]
    assert f["is_bound"] and not f["is_mbs"] and f["n_res_total"] == 2
    assert f["biolip2_two_consecutive_site"]  # residues 21,22 on one chain
    assert by(ligs, label_asym="F", op="2")[0]["sym_group"] == f["sym_group"]
    assert by(ligs, label_asym="F", op="2")[0]["site_group"] == f["site_group"]


# ---------------------------------------------------------------- real entries
def test_selection_filters():
    assert run("1d3z")[0]["skip_reason"] == "method"        # NMR ubiquitin
    assert run("5a40")[0]["skip_reason"] == "resolution"    # X-ray 3.6 A


def test_hiv_protease_1hvr():
    row, ligs = run("1hvr")
    assert row["assembly_type"] == "homomer" and row["n_protein_chains_asm"] == 2
    (x,) = by(ligs, ccd="XK2")
    assert x["is_mbs"] and x["is_mbs_min3"] and x["is_mbs_4A"] and x["is_mbs_1atom"]
    assert x["mbs_contrib_same_entity"] and x["mbs_contrib_symmetry_copy"] is False  # NCS copies A/B
    assert x["is_on_axis"] and x["c2_axis_dist"] < 1.0


def test_ttr_2rox_on_axis_after_dedup():
    """T4 in transthyretin: each channel ligand is disordered across the crystallographic 2-fold. v1 called channel 2
    off-axis (J=0.40); after merging the overlapping copy both channels are on-axis by both rules (audit M-3)."""
    row, ligs = run("2rox")
    reps = [l for l in ligs if l["ccd"] == "T44" and not l["is_dup_copy"]]
    dups = [l for l in ligs if l["ccd"] == "T44" and l["is_dup_copy"]]
    assert len(reps) == 2 and len(dups) == 2
    for l in reps:
        assert l["is_mbs"] and l["is_on_axis"] and l["is_on_axis_geom"] and l["c2_axis_dist"] < 1.0
    assert len({l["site_group"] for l in reps}) == 1


def test_ttr_1bm7_altloc_dedup():
    row, ligs = run("1bm7")
    reps = [l for l in ligs if l["ccd"] == "FLF" and not l["is_dup_copy"]]
    assert len(reps) == 2 and all(l["n_merged_copies"] == 3 for l in reps)
    assert all(l["is_on_axis_geom"] for l in reps)
    assert {l["dup_reason"] for l in ligs if l["ccd"] == "FLF" and l["is_dup_copy"]} >= {"altloc"}


def test_hiv_rt_3hvt_same_uniprot_hetero():
    row, ligs = run("3hvt")  # p66/p51, both P03366
    assert row["assembly_type"] == "same_uniprot_hetero" and row["has_same_uniprot_multi_entity"]


def test_glutathione_reductase_3grs_off_axis():
    row, ligs = run("3grs")
    fad = [l for l in ligs if l["ccd"] == "FAD" and not l["is_dup_copy"]]
    assert len(fad) == 2 and all(l["is_mbs"] and not l["is_on_axis"] and not l["is_on_axis_geom"] for l in fad)
    assert len({l["site_group"] for l in fad}) == 1


def test_sting_5cfl():
    row, ligs = run("5cfl")
    (c,) = by(ligs, ccd="C2E")
    assert c["is_mbs"] and c["is_on_axis"]
    assert all(not l["is_mbs"] for l in by(ligs, ccd="FLC"))


def test_streptavidin_1stp():
    row, ligs = run("1stp")
    assert row["assembly_type"] == "homomer" and row["n_protein_chains_asm"] == 4
    b = by(ligs, ccd="BTN")
    assert len(b) == 4
    for l in b:
        assert l["is_mbs"] and not l["is_mbs_min2"]          # 1 residue from the neighbour
        minor = [c for c in l["contacts"] if c["kind"] == "protein" and c["n_vdw"] >= 1
                 and c["chain"] != l["dominant_chain"]]
        assert [(c["comp"], c["auth_seq"]) for c in minor] == [("TRP", "120")]
        assert l["mbs_contrib_symmetry_copy"] is True


def test_hemoglobin_4hhb():
    row, ligs = run("4hhb")
    assert row["assembly_type"] == "heteromer" and set(row["uniprot_accs"]) == {"P68871", "P69905"}
    hem = by(ligs, ccd="HEM")
    assert len(hem) == 4 and all(l["is_bound"] and not l["is_mbs"] for l in hem)
    (po4f,) = by(ligs, ccd="PO4", label_asym="F")
    assert po4f["n_res_total"] == 0 and not po4f["is_bound"]  # real negative control


def test_heteromer_mbs_1nek():
    row, ligs = run("1nek")
    assert row["assembly_type"] == "heteromer"
    (h,) = by(ligs, ccd="HEM")
    chains = {d["chain"] for d in h["chain_counts"]}
    assert {"C", "D"} <= chains and h["is_mbs_min2"]
    assert by(ligs, ccd="FAD")[0]["is_mbs"] is False


def test_glycan_1hew():
    row, ligs = run("1hew")
    assert row["assembly_type"] == "monomer"
    assert all(l["is_glycan"] and l["is_branched"] for l in ligs if l["ccd"] == "NAG")


@pytest.fixture(scope="module")
def annotated():
    classes = pd.read_parquet(OUT_DIR / "ligand_classes.parquet")
    erows, lrows = [], []
    for pid in ["3ptb", "3q92", "2vb1", "1hvr", "1rb2"]:
        r, l = run(pid)
        erows.append(r)
        lrows += l
    ent = pd.DataFrame(erows)
    lig = pd.DataFrame(lrows).drop(columns=MG.SLIM_DROP)
    return MG.annotate(lig, ent, classes), classes.set_index("ccd")


def test_biolip_join_and_additives(annotated):
    lig, classes = annotated
    s = lambda pid, ccd: set(lig[(lig.pdb_id == pid) & (lig.ccd == ccd)].biolip_status)
    assert s("3ptb", "CA") == {"relevant"}
    assert s("3ptb", "BEN") == {"additive"}      # BioLiP2 dropped benzamidine in 3PTB
    assert s("3q92", "SO4") == {"additive"}
    assert s("3q92", "XNB") == {"relevant"}
    assert s("1hvr", "XK2") == {"relevant"}
    so4 = lig[(lig.pdb_id == "3q92") & (lig.ccd == "SO4")]
    assert so4.is_additive.all() and not so4.in_primary_set.any()
    # 1RB2 folate: absent from BioLiP2; relevant only if the entry flags it SUBJECT OF INVESTIGATION (it does not)
    fol = lig[(lig.pdb_id == "1rb2") & (lig.ccd == "FOL")]
    assert len(fol) and (fol.is_relevant == fol.is_ligand_of_interest).all()
    assert (fol.is_additive == ~fol.is_ligand_of_interest).all()
    assert lig.biolip_covered.all()
    for code in ["GOL", "SO4", "PEG", "EDO"]:
        assert classes.loc[code, "in_biolip_artifact_list"] and classes.loc[code, "plinder_artifact"], code


def test_ligand_classes(annotated):
    _, c = annotated
    exp = {"ZN": "metal_ion", "HEM": "heme_porphyrin", "FAD": "cofactor", "ATP": "nucleotide_like",
           "C2E": "nucleotide_like", "NAG": "sugar_glycan", "SO4": "inorganic", "SF4": "metal_cluster",
           "PLM": "lipid_fatty_acid", "NVP": "drug_like"}
    for code, cl in exp.items():
        assert c.loc[code, "ligand_class"] == cl, (code, c.loc[code, "ligand_class"])
    assert c.loc["ZN", "marsh_class"] == "metal" and c.loc["FAD", "marsh_class"] == "cofactor"
    assert c.loc["MK1", "is_druglike_chembl"]  # indinavir, approved
