# Manifest

Top-level layout of this repository. It was exported from the study workspace (branch `pocket-dimers-study`) with fresh history. Raw data, tool checkouts, environments and structure files are not included; see `project_record/DATA.md` for sources and `config/paths.env.example` for the data roots the scripts expect.

| Path | Contents |
|---|---|
| `briefing/` | Briefing material about the agentic-AI exercise (added separately) |
| `config/` | `paths.env.example`: data roots (PDB mirror, AFDB dimers, pocket atlas, AF3 MSAs and weights) read by the scripts |
| `docs/` | Per-step methods notes and results, each with the commands that produced its outputs (index below) |
| `envs/` | Environment lock files (pip freeze, conda lists, R lock) and pinned tool commits |
| `project_record/` | The AI PI's record: original prompt, frozen design, literature notes, research tree, lessons, data provenance |
| `report/` | `REPORT.md`, the main write-up; figure links point into `results/` |
| `results/` | Summary tables and figures per analysis (part1, part2, cofold, benchmarks), plus checksum listings of large inputs |
| `reviews/` | Design reviews, audits, report reviews and retrospectives (`*.md`), with the independent verifiers' code and small text outputs in `verify_*/` and `audit_*/` |
| `scripts/` | Pipeline scripts (Python, shell, R); each docs/*.md lists the commands for its step |
| `src/` | `part1` package for the PDB ligand-site census (structure parsing, contacts, ligand classes, site merging) |
| `tests/` | Known-answer tests for the Part 1 census (`pytest tests/test_part1.py`) |

## docs/ index

| File | Covers |
|---|---|
| `COMPUTE_POLICY.md` | Thread/process limits and job rules all agents followed |
| `DATA_RECON.md` | Reconnaissance of the mounted datasets (counts, layouts, commands) |
| `EXTERNAL_DATA.md` | Downloads of external tables (humsavar, AlphaMissense, SIFTS, Pfam, LIGYSIS, ChEBI, ChEMBL) |
| `TOOLS.md` | Install and smoke tests of the in-house tools burrow, pocketcompare and adock |
| `SEEDS.md` | Inventory of random seeds and RNG use across scripts |
| `PART1_PIPELINE.md` | Part 1 PDB ligand-site census pipeline: single- vs multi-chain binding sites |
| `PART1_RESULTS.md` | Part 1 results for pre-registered H1–H3 |
| `PDB_SITE_DESCRIPTORS.md` | Pocket descriptors at real PDB ligand sites, multi- vs single-chain |
| `FPOCKET_DRUGGABILITY.md` | fpocket druggability at PDB sites and predicted dimer-only pockets (exploratory) |
| `LIGYSIS_CALIBRATION.md` | Calibration of the burrow pocket finder against LIGYSIS (exploratory) |
| `H4_RESULTS.md` | H4: recovery of PDB ligand sites by predicted human dimers vs isolated chains |
| `Q22_RESULTS.md` | Q2.2/Q2.3a: dimer pockets vs monomer pockets in human, dimer-only ligandable pockets |
| `IFACE_LIGAND_MATCHED.md` | Ligand-matched vs unmatched predicted interface pockets (exploratory) |
| `H5_RESULTS.md` | H5: disease variants at dimer-only pockets |
| `Q24_RESULTS.md` | Q2.4: conservation of human dimer pockets in ortholog dimers (exploratory) |
| `COFOLDING.md` | AF3/Boltz co-folding feasibility, sizes and timings on the cluster node |
| `COFOLD_ONAXIS.md` | Deep dive 1: AF3 co-folding of on-axis drug sites in homodimers |
| `COFOLD_H4.md` | Deep dive 2: AF3 co-folding at H4 multi-chain ligand sites |
| `COFOLD_DIMERONLY.md` | Deep dive 3: does AF3 need the partner chain in curated dimer-only pockets |
| `COFOLD_HETERO_EXT.md` | Post-hoc extension with more heterodimer multi-chain sites (exploratory) |
| `COFOLD_SYNTHESIS.md` | Synthesis across co-folding runs: partner benefit vs how evenly the site is split |
