# R1 compute / methods review (DESIGN v1)

Basis: DESIGN, TOOLS.md, COMPUTE_POLICY, ledger. docs/COFOLDING.md does not exist yet. Nothing launched. Estimates are mine, not measured, unless stated.

## Runtime estimates (~19 h left)
- PDB census (261k gz mmCIF, gemmi parse + assembly 1 + contacts): ~0.3-1 s/entry on Lustre -> 16 workers ~3-5 h wall. Pre-filter by header (method, resolution, has non-polymer) to cut to ~100k first. Must be sharded + resumable.
- Pocket runs on ~50k dimers: 450 ms/structure/worker (measured) -> 16 workers ~25 min compute; plus isolated chains (2x, ~100k single-chain, ~150 ms) ~15 min. Cheap. Use describe (not predict) with --af-mode; run on CPU, not GPU.
- `describe --ligand` on PDB sites: cold start 15 s per call; must be one batch/process. Cost is dominated by per-site BURROW on the full assembly; expect ~0.5-2 s/site; limit to non-redundant (cluster x CCD) subset (~50-100k sites) -> 1-3 h at 16 workers. Not feasible for all sites.
- FAISS vs 3.23M: exact flat is ~0.25 s/query single-thread, batched matmul on GPU 1-3 is seconds for 1M queries. No bottleneck; build index in RAM (420 GB free). Copy atlas off scratch (15-day purge, ~5 GB) NOW.
- Co-folding: AF3 has NaN on Turing in all 5 try-logs, 3 of 6 size-sweep fails; Boltz-2 timing unfinished. Treat as the one real risk. Plan <=40-60 jobs, not 150, until a working model gives seconds/job. Budget 3 GPUs x ~4 h.

## Pitfalls
1. Numbering: dimer models use full-length UniProt numbering, but PDB auth numbering has insertion codes, tags and mutations. Use SIFTS per-residue (download XML/CSV, record URL/date); validate by residue-name identity after mapping, drop sites with <90% identity. Do not use struct_ref_seq offsets as primary.
2. Dimer PDB header chain->UniProt is unreliable; use DBREF-equivalent CSV columns (DATA_RECON). Chains >999 residues overflow PDB columns; parse robustly or use mmCIF.
3. Homodimer symmetry: a site in the PDB appears twice, and prediction gives two pockets. Collapse by entity-mapped residue sets before scoring (Jaccard >=0.5), and score recovery under both chain assignments (A<->B swap), taking max.
4. PDB assembly 1 may have >2 chains: define pair-level benchmark by chain pair that matches the predicted pair; sites at other interfaces must be excluded or flagged.
5. Pocket merge differences: match by residue overlap only; the merged unit is atlas-compatible, but burrow predict (score floor) is a different set (Gotcha 7). Use describe consistently for PDB and prediction.
6. Predicted-only comparison confound: PDB pockets include ligand-induced shape; ligand pocket by burrow locate may differ from the contacting set. Report recall against contact residues.

## Proposed changes
MUST
- Copy atlas fingerprints/metadata to workspace today; pin SIFTS/CCD download script.
- Decide co-folding engine by ~+4 h: if AF3 NaN unresolved, switch to Boltz-2 and cut scope to <=60 jobs.
- Sharded, resumable census jobs with ledger entries; all pocket batches via one process per worker list.
MUST (stats/numerics)
- Homodimer dedup and chain-swap rule in the Q2.1 metric.
SHOULD
- Subsample sites for `describe --ligand` (cluster x CCD, stratified MBS/SBS) instead of all.
- Run isolated chains with identical coordinates (already in design); cache dimer describe output once, reuse for Q2.2-2.4.
- Order: SIFTS+census -> dimer describe (parallel, 16 workers) -> ligand describe -> FAISS; cofold prep in parallel on GPUs.
COULD
- Use GPU matmul for atlas search; skip IVFPQ.
- Q2.4 ortholog conservation is the first thing to cut if time slips.
