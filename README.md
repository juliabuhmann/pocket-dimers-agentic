# Pockets between protein chains — an exercise in agentic AI for bioinformatics

This repository documents an **experiment in agentic science**. An autonomous team of AI agents ran a computational structural-biology project end to end, in about 18 hours on one cluster node. The roles were:

- **an AI "principal investigator" (PI) agent** that planned, delegated, reviewed and kept the record;
- **AI worker agents** that ran the analyses;
- **AI reviewers and verifiers** that critiqued the design and independently re-derived results.

A human wrote the brief and made the data available, but did not intervene in the science.

**Scientific question:** how do ligand-binding pockets formed by *two* protein chains (multi-chain sites, MBS) differ from ordinary *single-chain* pockets (SBS)? The study looked at:

- experimental structures (the PDB);
- AlphaFold models of predicted human homo- and heterodimers, with co-folding re-evaluation in AlphaFold3.

> **Status.** These are research outputs produced by AI agents. Independent AI verifier agents re-derived most headline numbers, and every claim in the report is tagged VERIFIED or UNVERIFIED. No human has yet checked the results, and the pipeline has not been rerun from scratch by anyone else. Treat it as a well-documented first pass, not a peer-reviewed result.

## Start here

| If you want… | Read |
|---|---|
| A readable overview: the biology, the agentic process, figures | [`briefing/briefing.pdf`](briefing/briefing.pdf) (source: `briefing/briefing.md`) |
| The full technical report with all numbers, tags and deviations | [`report/REPORT.md`](report/REPORT.md) |
| How the AI team worked and decided | [`project_record/RESEARCH_TREE.md`](project_record/RESEARCH_TREE.md), [`project_record/LESSONS.md`](project_record/LESSONS.md) |
| What was planned before any data were analysed | [`project_record/DESIGN.md`](project_record/DESIGN.md) (frozen, pre-registered H1–H5) |
| The critiques and independent checks | [`reviews/`](reviews/) |

## Main findings (short version)

1. **About a quarter to a third of ligand sites in protein complexes are multi-chain.** By fpocket druggability and the in-house pocket score, they are **neither more nor less druggable** than single-chain sites.
2. **Symmetry matters more than "being multi-chain".** Ligands bound on the two-fold axis of homodimers are overwhelmingly drug-like (64 drug-like vs 8 cofactor protein families). This result is descriptive and between proteins.
3. **The major cofactors (NAD, FAD, SAM, PLP, CoA) are usually anchored in one chain.** A pre-registered "cofactor vs drug-like" comparison was therefore null.
4. **Predicted dimers recover real multi-chain sites somewhat more often than isolated chains do.** Exploratory: 68 families better vs 32 worse. The pre-registered test was **not** significant.
5. **AlphaFold3 co-folding is more accurate with the partner chain for symmetric drug sites.** For cofactors it rarely needs the partner. These are small case series with possible memorisation.
6. **Agentic process.** Independent verification caught several errors that would have changed conclusions, including a contaminated "drug-like" label, an invalid permutation test and a cropping bias. Pre-registration reduced, but did not remove, post-hoc choices.

## Repository layout

| Path | Contents |
|---|---|
| `briefing/` | Research briefing (Markdown, PDF, figures, figure script, PDF build script) |
| `report/` | Full technical report (`REPORT.md`) |
| `project_record/` | Original prompt, literature review, frozen design, research tree, process lessons, dataset provenance (`DATA.md`) |
| `reviews/` | AI design reviews (2 rounds), audit, independent verification reports with their own code, report reviews, retrospectives |
| `docs/` | Per-step methods and results notes, each with the commands that produced its outputs |
| `src/`, `scripts/`, `tests/` | Pipeline code (PDB census package, analysis and co-folding scripts, known-answer tests) |
| `results/` | Summary tables and figures per analysis; checksums of large inputs |
| `envs/` | Environment lock files and pinned tool commits |
| `config/` | `paths.env.example`: the data roots the scripts expect |

[`MANIFEST.md`](MANIFEST.md) has a per-file index of `docs/`.

## Reproducing

**Data and tools are not included.** This covers PDB mirror, AlphaFold DB models, the dimer model set, AlphaFold3 weights and alignments, and the in-house tools. Sources, versions, checksums and licences are listed in `project_record/DATA.md`. External tables (UniProt, ChEMBL, SIFTS, LIGYSIS, …) are downloaded by `scripts/fetch_*`.

1. **Set data roots:**
   ```bash
   cp config/paths.env.example config/paths.env
   # edit config/paths.env
   set -a; source config/paths.env; set +a
   ```
2. **Create the environments** from `envs/` and install the tools with `scripts/setup_tools.sh` and `scripts/setup_cofold.sh`. AlphaFold3 weights must be obtained from Google DeepMind under their terms.
3. **Limit threads** before any parallel job: `source scripts/env_threads.sh` (see `docs/COMPUTE_POLICY.md`).
4. **Run the steps in order:**
   1. Part 1 census and analysis;
   2. PDB site descriptors and fpocket;
   3. dimer pocket runs;
   4. H4 recovery;
   5. Q2.2/Q2.3a;
   6. H5;
   7. Q2.4;
   8. co-folding.

   Each `docs/*.md` gives the exact commands for its step, and `report/REPORT.md` §9 lists the order.
5. **Seeds:** [0, 1, 2], seed 0 primary; see `docs/SEEDS.md`.

The full pipeline needs a multi-core machine. The PDB census took about 16 min on 24 cores; pocket runs on ~10,600 dimer models took about 1 h on 16 cores. AlphaFold3 needs a GPU with ≥ 24 GB.

## What is deliberately not in this repository

- Raw and derived bulk data (`data/`), structure files, model outputs, model weights and environments, all excluded via `.gitignore`.
- **Per-residue tables derived from AlphaMissense** (licence CC BY-NC-SA). Only aggregate results are included.
- **Personal information.** Personal cluster paths were replaced by `/path/to/...` placeholders. The scripts read data locations from environment variables (`config/paths.env.example`); an optional contact e-mail for polite API use is read from `CONTACT_EMAIL`.

## About the agentic setup

- **Models:** the PI and most analysis/verifier agents ran on Claude Opus-class models; some data-engineering and review tasks ran on Sonnet-class models. All ran inside a Singularity container with no internet push rights.
- **Rules the agents followed:**
  - pre-register the hypotheses before analysis;
  - never change metrics to improve results;
  - have every headline number re-derived by a different agent with its own code;
  - log every deviation;
  - commit in small functional units.
- **What we learned about running research this way** is in `project_record/LESSONS.md` and in the second half of the briefing.

## Citation and licence

No licence has been chosen yet; please contact the repository owners before reuse. Third-party data keep their own licences (`project_record/DATA.md`).
