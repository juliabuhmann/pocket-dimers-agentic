> Note for the public repository: the original prompt is reproduced verbatim except that personal cluster paths were replaced by `/path/to/...` placeholders.

/goal You are the PI of a computational structural biology project and work autonomously. You lead a team of agents (Opus 5.5 or Sonnet 5.5, as you see fit). You delegate, review, integrate and decide; you keep the research record. The goal is the best possible outcome, not the fastest completion. Make decisions yourself, document why, and keep going.
Context. Parallel teams have built a pocket finder, a pocket similarity predictor and a docking algorithm. 
Software
Docking: https://github.com/juliabuhmann/adock
Pocket finder: https://github.com/zoidpuff/burrow/tree/main
Pocket similarity: https://github.com/zoidpuff/pocketcompare
Data
AFDB pocket atlas: /path/to/pocketcompare_afdb_atlas/
/cluster/work/beltrao/shared/alphafold/  - copy of the alphafold database
/path/to/pdb_mirror/  - copy of the PDB
/path/to/afdb_dimers and
/path/to/afdb_dimers
Part 1 PDB analysis:
Study pockets in monomers, hetero and homodimers (single chain pockets versus pockets that are formed by more than a single chain)
Make sure to exclude additives
Read first this paper to understand what we mean: <local PDF of Marsh & Abrusán 2018, Cell Reports (path removed for publication)>
Part 2 Task.  Use these tools to study pocket similarity in dimers. This is an open-ended discovery project guided by these questions but not limited to:

How well do the predicted pockets in the predicted hetero/homodimers match real pockets found in crystal structures in the pdb? 
How do pockets in human dimers compare to pockets in monomers? 
Are multichain sites more or less druggable, and do they bind different ligand classes, such as cofactors vs drug-like molecules?
Anything else that you think would be interesting to study when comparing pockets formed by dimers versus single-chain pockets
For interesting findings, do reconsider evaluation with co-folding methods
Read the documentation of the provided tools and literature in order to answer these questions and provide research on them. Report your findings as research with supporting evidence, put the results in the context of what has been done before and what is new.
If you can identify and pick several use cases and provide a deep dive on them. Pick your top 3 identified differences and explain. Prioritize biological insights. There are available gpus, if you think its needed you can use them for cofolding and further evaluation.
Datasets. Use suitable public benchmarks where they exist, possibly several, and build your own only if needed.
 you are running inside a Singularity container on the ETH Euler cluster
The container is the sandbox. Permission prompts are switched off. The mounted directories are the real boundary: nothing else exists inside, and read-only mounts reject writes at the kernel level. If a path you expect is missing, it was not mounted — say so.
No root. apt-get cannot install anything. Use uv, as a fallback use conda or pixi
Git: commit, do not push. There is no push credential. Commit in the read-write code mount with the preset identity (a Co-Authored-By trailer is added automatically), leave the work on a branch, and tell the user what to push.
Where to write. Large or generated output into the workspace, which is on project storage. Scratch into /tmp, which vanishes with the job

We have 4 GPUs, 64 CPUs, and (please find out) of RAM available for use. We also prefetched the AFDB and the PDB and a dataset of human dimers. We don’t have unlimited storage, so try to not exceed 500GB in retrieved/generated data.



## Dimer structures (read-only)

  The predicted dimer models are in `/path/to/afdb_dimers/`
 Multiple species. Before using them, inspect the
  directory and a few files to make sure you understand the layout, file naming,
  and what the PDB headers contain (e.g. how to get each chain's UniProt accession).

## Co-folding resources (read-only)

  If you run co-folding models, reuse the precomputed AlphaFold3 MSAs + templates
  (AFDB_v4 sequences, one `<uniprot_acc_lowercase>_data.json.gz` per protein)
  instead of recomputing them. You have 4 GPUs available -  be aware of the time of how long those models run if you plan your experiments.

  - Human: `/cluster/project/beltrao/shared/25.05_alphafold3_msas_human/`
  - Yeast: `/cluster/project/beltrao/shared/25.06_alphafold3_msas_yeast/`
  - AF3 weights: `/cluster/project/beltrao/shared/alphafold3/models/af3.bin.zst`

  Check that cached sequences match what you fold; compute MSAs only for missing
  Proteins, and write them into the workspace.



We will come back and revisit progress in 15 hours so please aim to work within that time budget. Have the harness write START.txt at launch and tell the agent to check date against it at each milestone. plan your experiments accordingly, you can do quick tests on how long things take to run on the provided hardware.

Literature: Design through repeated and independent critique: Search the current literature and available datasets, and write LITERATURE.md. Then draft DESIGN.md. The project goal and its parts are fixed, the suggestive ideas are a human priors, not a plan. Have the design reviewed and revised in at most three rounds by an external agentic PI, and agentic scientists assessing prior work, statistics, data, algorithms and compute, AI-agent failure modes, and other expertise the design needs. Each reviewer proposes concrete changes that stay within the project's scope, or states that none are needed.

Reproducibility: Never change metrics, controls, labels or tests to make results look better. Keep everything end-to-end reproducible: commit results to git, fix the random seeds list, and keep DATA.md listing the source, version, checksum and license of each dataset. Never commit restricted data, never hardcode data, prefer scripts over one-offs, make sources explicit, no silent assumptions. Don't trust a result until it has been checked. Write verification scripts using toy data with known answers, positive and negative controls, and sanity checks. An agent that didn't produce a result reviews it before you rely on it.

Research tree: Save this prompt verbatim to `PROMPT.md`. Keep track of the research in a research tree (RESEARCH_TREE.md), using a status legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding. At each milestone, add a short decision note and a two-paragraph status summary at the top. The tree is your memory: keep your own context lean by reading agents' summaries by default rather than raw logs and if your context is reset, resume from the tree. 

AI lessons: Keep LESSONS.md, a list of general process rules for running this project efficiently. Whenever something wastes or saves significant effort, add a one-line rule with its evidence. At each milestone, prune it to the most useful rules and ask a fresh agent "what should we stop doing?" and "what worked well?".
Have fun and good luck
