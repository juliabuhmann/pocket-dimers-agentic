"""Known-answer checks for the SIFTS XML parser (4HHB, 1HVR, 2LYZ)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from fetch_sifts_xml import fetch_sifts_xml, parse_many

ids = ["1hvr", "4hhb", "2lyz"]
print(fetch_sifts_xml(ids))
df = parse_many(ids)
print(df.shape)
a = df[(df.pdb_id == "4hhb") & (df.chain == "A")].sort_values("label_resnum")
first = a.iloc[0]
print(first.to_dict())
assert first.uniprot_acc == "P69905" and first.uniprot_pos == 2 and first.pdb_resname == "VAL"
assert set(df[df.pdb_id == "4hhb"].uniprot_acc.dropna()) == {"P69905", "P68871"}
assert df[(df.pdb_id == "4hhb") & (df.chain == "B") & (df.label_resnum == 1)].uniprot_pos.iloc[0] == 2  # HBB Met1 cleaved
# 2LYZ: hen lysozyme P00698; mature chain starts Lys1 = UniProt 19
l = df[(df.pdb_id == "2lyz")].sort_values("label_resnum").iloc[0]
print(l.to_dict()); assert l.uniprot_acc == "P00698" and l.uniprot_pos == 19 and l.pdb_resname == "LYS"
# 1HVR: HIV-1 protease (P03366 or P04585); Pro1 -> UniProt pos 69 in Gag-Pol P03366 (check printed)
h = df[(df.pdb_id == "1hvr")].dropna(subset=["uniprot_acc"]).sort_values(["chain", "label_resnum"]).iloc[0]
print(h.to_dict())
# self-consistency: UniProt aa matches 1-letter of PDB residue for >=95% mapped standard residues
AA3 = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I",
           LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
m = df.dropna(subset=["uniprot_acc"]).copy()
m = m[m.pdb_resname.isin(AA3)]
agree = (m.pdb_resname.map(AA3) == m.uniprot_resname).mean()
print("PDB vs UniProt residue agreement", round(agree, 4), len(m)); assert agree > 0.95
print("ALL CHECKS PASSED")
