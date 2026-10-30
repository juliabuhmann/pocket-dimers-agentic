"""Part 1 (PDB ligand-site census) configuration: every path, URL and threshold in one place.

Each constant carries its rationale/source. Literature references point to LITERATURE.md
sections (a) (Abrusan & Marsh 2018 framework) and (f) (artifact lists).
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

CENSUS_VERSION = "v2"  # v2 = audit fixes (reviews/audit_part1_census.md), see docs/PART1_PIPELINE.md "v2 changes"

# ---------------------------------------------------------------- paths
PDB_MIRROR = Path(os.environ.get("PDB_MIRROR", "/cluster/scratch/jjaenes/pdb"))  # local mmCIF mirror, layout <id[1:3]>/<id>.cif.gz (docs/DATA_RECON.md s.3)
RAW_DIR = ROOT / "data" / "raw" / "part1"
OUT_DIR = ROOT / "data" / "part1"
SHARD_DIR = OUT_DIR / "shards"
SUMMARY_DIR = ROOT / "results" / "part1"  # small committed summary tables

# ---------------------------------------------------------------- external sources
# name -> (URL, local filename in RAW_DIR). Fetched by scripts/fetch_part1_sources.py.
SOURCES = {
    "biolip_ligand_list": ("https://zhanggroup.org/BioLiP/ligand_list", "biolip_ligand_list.tsv"),
    "plinder_badlist": (
        "https://raw.githubusercontent.com/plinder-org/plinder/main/src/plinder/data/annotations/static_files/artifacts_badlist.csv",
        "plinder_artifacts_badlist.csv",
    ),
    "ccd_components": ("https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz", "components.cif.gz"),
    # PDBe cofactor classes (CCD codes per class), derived from the EMBL-EBI CoFactor DB (Fischer/Thornton).
    "pdbe_cofactors": ("https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors", "pdbe_cofactors.json"),
    "sifts_chain_uniprot": (
        "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz",
        "pdb_chain_uniprot.tsv.gz",
    ),
    # BioLiP2 precomputed annotation: one row per biologically relevant ligand-protein site (artifacts removed by
    # BioLiP2's curation incl. the PubMed check). Primary additive/relevance label (PI design review round 1).
    "biolip_annotation": ("https://zhanggroup.org/BioLiP/download/BioLiP.txt.gz", "BioLiP.txt.gz"),
    # UniChem whole-source mapping ChEMBL (src 1) -> PDBe ligand / CCD (src 3)
    "unichem_chembl_pdbe": (
        "https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src3.txt.gz",
        "unichem_src1src3.txt.gz",
    ),
    # ChEMBL max clinical phase for molecules with max_phase > 0 (paginated REST API, concatenated to TSV)
    "chembl_max_phase": (
        "https://www.ebi.ac.uk/chembl/api/data/molecule.json?max_phase__gt=0&only=molecule_chembl_id,max_phase&limit=1000",
        "chembl_max_phase.tsv",
    ),
    "rcsb_clusters_30": ("https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-30.txt", "clusters-by-entity-30.txt"),
    "rcsb_clusters_95": ("https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-95.txt", "clusters-by-entity-95.txt"),
}


def raw_path(name: str) -> Path:
    return RAW_DIR / SOURCES[name][1]


# ---------------------------------------------------------------- entry selection
XRAY_METHODS = ("X-RAY DIFFRACTION",)
EM_METHODS = ("ELECTRON MICROSCOPY",)
XRAY_MAX_RESOLUTION = 3.0  # Angstrom; task spec (typical "medium-resolution" cut used in BioLiP-based surveys)
EM_MAX_RESOLUTION = 3.5    # Angstrom; task spec (side chains/ligands generally interpretable at <=3.5 A EM)
MIN_PROTEIN_CHAIN_LEN = 30  # SEQRES residues; polypeptide entities shorter than this are "peptides", not protein chains
                            # (task spec; also used to define protein chains for assembly typing and MBS calls)
ASSEMBLY_ID = None          # None = first assembly listed in _pdbx_struct_assembly (normally "1")
MAX_ASSEMBLY_CHAINS = 60    # polymer chains in the generated assembly; larger ones are skipped as "huge_assembly"
                            # (stand-in for Marsh's exclusion of viral capsids / fibrils / helical assemblies, LIT (a))
MODEL_INDEX = 0             # first model only (EM/X-ray entries normally have one)

# ---------------------------------------------------------------- contacts
VDW_TOLERANCE = 0.5          # A; BioLiP: d <= vdW_i + vdW_j + 0.5 (Yang et al. 2013; LIT (a),(f))
VDW_SOURCE = "gemmi.Element.vdw_r"  # van der Waals radii table used (gemmi built-in)
SIMPLE_CUTOFF = 4.0          # A; flat heavy-atom cutoff, sensitivity analysis (LIT (a) recommended defs)
MIN_ATOM_CONTACTS_PRIMARY = 2  # v2: atom pairs per residue for the primary residue call = BioLiP2 rule (DESIGN s.0, LIT (f))
MIN_ATOM_CONTACTS_SENS = 1     # v1 primary (BioLiP 2013 / Marsh); stored as *_1atom sensitivity columns
MIN_BOUND_RESIDUES = 2       # a ligand instance counts as bound if it contacts >=2 protein residues in total (BioLiP)
MIN_RES_PER_CHAIN_MBS = 1    # primary MBS: >=2 protein chains each with >=1 contacting residue (Marsh-like)
MBS_STRICT_LEVELS = (2, 3)   # also stored: MBS if >=2 chains each contribute >= k residues

DEDUP_CENTROID_MAX = 2.0     # A; same-CCD ligand copies (after assembly operators) with centroids closer than this are one
                             # physical ligand (special positions, overlapping altloc copies). DESIGN.md v3 (frozen).
DEDUP_ATOM_OVERLAP_DIST = 1.5  # A; v2: copies also merged if > DEDUP_ATOM_OVERLAP_FRAC of heavy atoms lie within this
DEDUP_ATOM_OVERLAP_FRAC = 0.5  #    distance of the other copy's atoms (ligands disordered across a 2-fold; audit M-3)
SITE_JACCARD_MIN = 0.5       # v2: same-CCD instances in one entry with chain-agnostic (protein, UniProt pos) residue-set
                             # Jaccard >= this are one symmetric site (union-find; audit H-2)
ON_AXIS_GEOM_MAX_DIST = 3.0  # A; v2 geometric on-axis flag: ligand centroid within this of the chain-pair C2 axis (audit M-3)
HOMOMER_MIN_IDENTITY = 0.90  # v2: chains of different entities count as the same protomer type only if identity >= 0.90
HOMOMER_MIN_COVERAGE = 0.80  #     over >= 80% of the chains (UniProt-position alignment; applied to the LONGER chain so that
                             #     p66/p51-type processed chains are not homomers); else same_uniprot_hetero
PHOTOSYNTH_CCD_NAME_RE = r"CHLOROPHYLL|PHEOPHYTIN"  # v2: CCD names marking photosynthetic pigments (from components.cif)
PHOTOSYNTH_KEYWORD_RE = r"PHOTOSYSTEM|LIGHT[- ]HARVESTING|PHOTOSYNTHE|REACTION CENT"  # entry keywords (struct_keywords)
ON_AXIS_JACCARD_MIN = 0.5    # homomer MBS "on-axis": UniProt residue sets contacted on the two top chains have Jaccard >= this
                             # (DESIGN.md v3)
C2_ANGLE_TOL = 15.0          # deg; chain-pair superposition rotation within 180 +/- tol is treated as a 2-fold axis
MIN_CA_FOR_AXIS = 10         # common CA atoms needed to superpose two chains for the axis calculation

# ---------------------------------------------------------------- artifacts
BIOLIP_MAX_COPIES = 15       # BioLiP2: listed ligand discarded if it appears >15 times in the entry (LIT (f))
# PLINDER generic artifact filters (LIT (f)); applied to organic (carbon-containing) ligands only
PLINDER_MIN_CARBONS = 2
PLINDER_MIN_HEAVY_ATOMS = 5
PLINDER_MAX_ABS_CHARGE = 2
PLINDER_MAX_UNBRANCHED_CHAIN = 12  # unbranched hydrocarbon chain longer than this -> artifact

# ---------------------------------------------------------------- ligand classes
INCLUDE_BRANCHED = True      # include branched (glycan) entities as ligand instances, flagged is_branched
NONMETALS = frozenset(
    "H He B C N O F Ne Si P S Cl Ar Ge As Se Br Kr Sb Te I Xe At Rn".split()
)  # elements not counted as metals (metalloids B/Si/Ge/As/Sb/Te treated as non-metals)
DRUGLIKE_QED_MIN = 0.5       # QED >= 0.5 ("attractive" compounds; Bickerton et al. 2012 Nat Chem, mean QED of oral drugs ~0.49-0.5)
DRUGLIKE_MAX_RO5_VIOLATIONS = 1  # Lipinski 1997: poor absorption more likely with >1 violation
DRUGLIKE_MIN_CHEMBL_PHASE = 1.0  # primary drug-like label: mapped (UniChem) ChEMBL molecule reached clinical phase >= 1
                                 # (DESIGN.md v3); Ro5+QED rule kept as FALLBACK flag
DRUGLIKE_MIN_HEAVY_ATOMS = 12    # our choice: excludes buffer/solvent-sized molecules from the drug-like class
LIPID_MIN_CHAIN = 8          # our choice: >=8 consecutive acyclic sp3 carbons -> lipid/fatty acid/detergent tail
# Class precedence (first match wins), see src/part1/ligands.py
CLASS_ORDER = ("metal_ion", "metal_cluster", "inorganic", "heme_porphyrin", "cofactor", "nucleotide_like",
               "sugar_glycan", "amino_acid", "lipid_fatty_acid", "drug_like", "other_organic")

# ---------------------------------------------------------------- run
N_WORKERS = 24               # docs/COMPUTE_POLICY.md: this task may use at most 24 CPU workers
SHARD_SIZE = 500             # entries per output shard (resumable unit)
SAMPLE_FRACTION = 0.01       # pilot run fraction
SAMPLE_SEED = 0
