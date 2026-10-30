# R1 — Agent failure modes (pocket_dimers_29oct)


## Project-specific risks
1. **Metric drift / garden of forking paths.** DESIGN.md has many tunable choices (Jaccard 0.3, recall 0.5, 1 vs 3 residues per chain, +0.5 Å vs 4 Å, two additive lists, two Marsh labels) and a "top 3 differences" step chosen by effect size after seeing results. That is an invitation to p-hack.
2. **Circular/leaky recovery benchmark (Q2.1).** AFDB dimers and AF-based co-folding may have seen the PDB sites; "recovery" may be memorisation. Selection bias: only confident (ipSAE ≥0.6) dimers.
3. **Silent fallbacks.** fpocket "if installable", ChEMBL/DrugBank mapping "where obtainable", OMA "if retrievable", AF3 to Boltz-2, atlas settings vs dimer settings: a swapped method may go unreported. Also SIFTS mapping failures and dropped structures (≤60 chains, parse errors) silently shrinking N.
4. **Unverified trust between agents.** Agent summaries (e.g. "3,378 files", DATA_RECON facts, tool README claims such as "pocketcompare weak on interfaces") get copied into DESIGN/REPORT without re-derivation.
5. **Runaway compute / context loss.** A 3.1 h agent overrun already occurred; 23.5 h budget, thread-limit crash history; long-session context compaction can drop frozen definitions.
6. **Over-claiming.** Druggability inferred from ligand class or burrow_score; cofactor-vs-drug-like causal language from observational, redundancy-confounded counts; co-folding on ≤150 jobs generalised.

## Guardrails
- **MUST** Freeze DESIGN.md at M1: tag in git, write SHA256 of DESIGN.md + config.yaml into every results file header; one primary definition per metric, all others labelled "sensitivity". Deviations logged in RESEARCH_TREE.md with timestamp before rerun.
- **MUST** Pre-register the deep-dive selection rule (rank by effect size with CI excluding 0 after BH, cluster-weighted, in primary definitions only), plus a "tried-and-not-reported" count of all comparisons run.
- **MUST** Sentinel known-answer cases in every pipeline stage (HIV-1 protease = MBS, lysozyme = SBS, a hand-checked heteromer; a planted fake interface pocket; shuffled labels give null). Pipeline aborts if sentinels fail.
- **MUST** Negative controls reported beside every headline (random pockets, isolated chains, label permutation), plus a PDB-release-date split (post AF training cutoff) for recovery.
- **MUST** Failure/exclusion ledger: every stage emits n_in, n_out, n_dropped by reason; report it in REPORT.md. Fallbacks must raise or write a `FALLBACK` flag into the output table, never degrade silently.
- **MUST** Headline numbers re-derived independently: a separate agent, with code written from the definitions only (not the original script), must reproduce each headline within tolerance before it enters REPORT.md.
- **SHOULD** Claim tags in REPORT.md: [VERIFIED: who/how], [UNVERIFIED], [LITERATURE-CLAIM]. No tag, no inclusion. Agent-sourced facts (tool behaviour, counts) need file+line or command evidence.
- **SHOULD** Spot audit: PI manually inspects ≥10 random sites per class (structure view, contact list) before trusting classification (additive filter, ligand class, MBS label).
- **SHOULD** Wall-clock and thread budget in every agent brief, with heartbeat in logs/compute_ledger.txt; auto-kill at 1.5x budget. Write STATE.md checkpoint (frozen definitions, hashes, open items) at each milestone for context recovery.
