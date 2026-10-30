#!/usr/bin/env python3
"""Deep dive #1 co-folding: select human homodimer PDB cases with an ON-AXIS drug-like ligand (primary set) plus
controls (off-axis MBS cofactor sites; SBS drug-like sites), and write AF3 inputs.

Status: first attempt, UNVERIFIED. Selection is fully rule-based (constants below); no hand-picked IDs.

Rerun:
    source scripts/env_threads.sh
    .venvs/part1/bin/python scripts/cofold_select_onaxis.py          # -> results/cofold/onaxis/cases.tsv, data/cofold/onaxis/inputs/

Rules (all in code):
- Sites = Part 1 primary collapsed sites (scripts/part1_analysis.py load + site_filter, variant "primary").
- Homodimer = assembly_type homomer and exactly N_CHAINS (=2) protein chains in assembly 1; all protein contacts map to
  one UniProt accession whose organism (UniProt REST) is Homo sapiens (taxon 9606).
- Strata: on_axis_druglike = MBS & census is_on_axis (Jaccard) & cls drug_like;
          offaxis_cofactor = MBS & not on-axis & cls cofactor;  sbs_druglike = not MBS & cls drug_like.
- Ligand: single CCD component, not covalent, given to AF3 by CCD code. Ligand copies for the dimer run = number of
  bound, non-duplicate instances of that CCD in the assembly (on-axis: forced to 1 as in the task; controls must have 2).
- AF3 protein sequence = PDB entity SEQRES trimmed to the first..last residue aligned to UniProt (drops tags);
  MSA columns are taken from the cached full-length UniProt AF3 MSA (cached query must equal the UniProt sequence),
  unmapped (linker/mutation-free insert) positions become gap columns. Templates are NOT used (avoid leaking holo
  structures; the cached templates go up to 2025).
- Tokens (dimer) = 2 * len(seq) + copies * heavy atoms <= MAX_TOKENS.
- One case per cluster30: post-cutoff release (> AF3_CUTOFF) preferred, then best resolution, then pdb_id.
- On-axis: up to N_ONAXIS clusters, ordered post-cutoff first, then fewer tokens (wall-clock budget).
- Controls: N_COF / N_SBS clusters (not used by another stratum), ordered post-cutoff first, then |tokens - median
  on-axis tokens| (rough size matching).
"""
import os
import argparse
import gzip
import io
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
import part1_analysis as P1  # noqa: E402
from part1.config import PDB_MIRROR  # noqa: E402

MSA_DIR = Path(os.environ.get("AF3_MSA_HUMAN", "/cluster/project/beltrao/shared/25.05_alphafold3_msas_human"))
AF3_CUTOFF = "2021-09-30"          # AF3 training data cutoff (Abramson et al. 2024)
N_CHAINS = 2
MAX_TOKENS = 1024                  # docs/COFOLDING.md: <=1280 fits; 1024 bucket ~12 min/seed -> cap for 4 h budget
N_ONAXIS, N_COF, N_SBS = 12, 6, 6  # task targets
HUMAN_TAXID = 9606
OUT_TSV = ROOT / "results/cofold/onaxis/cases.tsv"
WORK = ROOT / "data/cofold/onaxis"
CACHE_SITES = WORK / "candidate_sites.parquet"


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------------ candidates from census
def candidate_sites():
    if CACHE_SITES.exists():
        return pd.read_parquet(CACHE_SITES)
    s = P1.load(P1.CENSUS_DIR)
    v = dict(P1.PRIMARY)
    col = P1.site_filter(s, v)
    col = col[(col.assembly_type == "homomer") & (col.n_protein_chains == N_CHAINS)]
    col["onax"] = col.onax.fillna(False).astype(bool)
    col["stratum"] = None
    col.loc[col.y_mbs & col.onax & (col.c == "drug_like"), "stratum"] = "on_axis_druglike"
    col.loc[col.y_mbs & ~col.onax & (col.c == "cofactor"), "stratum"] = "offaxis_cofactor"
    col.loc[~col.y_mbs & (col.c == "drug_like"), "stratum"] = "sbs_druglike"
    col = col[col.stratum.notna()]
    keep = ["pdb_id", "ligand_uid", "ccd", "cluster30", "resolution", "release_date", "method_class", "stratum",
            "n_heavy_lig", "c2_axis_dist", "on_axis_jaccard", "c"]
    col = col[keep].copy()
    extra = pd.read_parquet(P1.CENSUS_DIR / "ligand_sites.parquet",
                            columns=["ligand_uid", "centroid", "is_covalent", "n_components", "dominant_uniprot"],
                            filters=[("pdb_id", "in", sorted(col.pdb_id.unique()))])
    col = col.merge(extra, on="ligand_uid", how="left")
    col["centroid"] = col.centroid.apply(lambda x: list(map(float, x)))
    WORK.mkdir(parents=True, exist_ok=True)
    col.to_parquet(CACHE_SITES)
    return col


def uniprot_meta(accs):
    """accession -> (taxid, sequence) from UniProt REST (batched, cached in WORK/uniprot_meta.tsv)."""
    cache = WORK / "uniprot_meta.tsv"
    have = pd.read_csv(cache, sep="\t") if cache.exists() else pd.DataFrame(columns=["Entry", "Organism (ID)", "Sequence"])
    accs = sorted(set(a for a in accs if isinstance(a, str) and a and "-" not in a) - set(have.Entry))
    new = [have]
    for i in range(0, len(accs), 100):
        url = ("https://rest.uniprot.org/uniprotkb/accessions?format=tsv&fields=accession,organism_id,sequence"
               "&accessions=" + ",".join(accs[i:i + 100]))
        with urllib.request.urlopen(url, timeout=120) as r:
            new.append(pd.read_csv(io.StringIO(r.read().decode()), sep="\t"))
    t = pd.concat(new, ignore_index=True).drop_duplicates("Entry")
    t.to_csv(cache, sep="\t", index=False)
    return {x["Entry"]: (int(x["Organism (ID)"]), x["Sequence"]) for _, x in t.iterrows()}


# ------------------------------------------------------------------ per-entry structure checks
def entry_info(row, unp_seq, acc):
    """Read assembly 1; return dict with protein seq for AF3, mapping, ligand copies; or reason string."""
    st = gemmi.read_structure(str(PDB_MIRROR / row.pdb_id[1:3] / f"{row.pdb_id}.cif.gz"))
    st.setup_entities()
    asm = st.assemblies[0] if st.assemblies else None
    if asm is None:
        return "no_assembly"
    a = gemmi.make_assembly(asm, st[0], gemmi.HowToNameCopiedChain.AddNumber)
    a.remove_hydrogens()
    pchains = [ch for ch in a if ch.get_polymer().check_polymer_type() in
               (gemmi.PolymerType.PeptideL,)]
    if len(pchains) != N_CHAINS:
        return f"n_protein_chains={len(pchains)}"
    ents = {ch.get_polymer()[0].entity_id for ch in pchains}
    if len(ents) != 1:
        return "two_entities"
    ent = st.get_entity(next(iter(ents)))
    seqres = gemmi.one_letter_code(ent.full_sequence).upper()
    if not seqres or seqres.count("X") > 2:
        return "nonstandard_seqres"
    al = gemmi.align_string_sequences(list(seqres), list(unp_seq), [])
    # map seqres index -> uniprot index (0-based), identical residues only count as mapped
    m, i, j = {}, 0, 0
    for op, n in _cigar(al.cigar_str()):
        if op == "M":
            for k in range(n):
                m[i + k] = j + k
            i += n; j += n
        elif op == "I":
            i += n
        else:
            j += n
    if not m:
        return "no_alignment"
    lo, hi = min(m), max(m)
    seq = seqres[lo:hi + 1]
    mapping = [m.get(k) for k in range(lo, hi + 1)]  # per AF3 residue -> uniprot 0-based or None
    ident = np.mean([seq[k] == unp_seq[u] for k, u in enumerate(mapping) if u is not None])
    if ident < 0.9:
        return f"low_identity={ident:.2f}"
    # ligand copies of this CCD bound in assembly: residues named ccd with >= 1 heavy atom within 4 A of protein
    ns = gemmi.NeighborSearch(a, st.cell, 5).populate()
    copies = []
    for ch in a:
        for res in ch:
            if res.name != row.ccd:
                continue
            if res.het_flag != "H":
                continue
            alt = {at.altloc for at in res}
            near = False
            for at in res:
                for mk in ns.find_atoms(at.pos, "\0", radius=4.0):
                    cra = mk.to_cra(a)
                    if cra.residue.het_flag != "H" and cra.chain.name in {c.name for c in pchains}:
                        near = True; break
                if near:
                    break
            if near:
                xyz = np.array([[at.pos.x, at.pos.y, at.pos.z] for at in res if at.altloc in ("\0", "A")])
                copies.append(xyz.mean(0))
    # collapse altloc/duplicate copies (centroids within 2 A)
    uniq = []
    for c in copies:
        if all(np.linalg.norm(c - u) > 2.0 for u in uniq):
            uniq.append(c)
    return dict(seq=seq, seqres_lo=lo, mapping=mapping, identity=float(ident), n_lig_copies=len(uniq),
                pchains=",".join(c.name for c in pchains), entity=ent.name)


def _cigar(s):
    import re
    return [(op, int(n)) for n, op in re.findall(r"(\d+)([MID])", s)]


def cached_msa(acc):
    p = MSA_DIR / f"{acc.lower()}_data.json.gz"
    if not p.exists():
        return None
    with gzip.open(p, "rt") as f:
        d = json.load(f)
    prots = [s["protein"] for s in d["sequences"] if "protein" in s]
    return prots[0] if len(prots) == 1 else None


def remap_a3m(a3m, mapping, new_query):
    """Keep match columns of a full-length a3m at `mapping` (uniprot 0-based or None -> gap)."""
    out, name = [], None
    lines = a3m.splitlines()
    first = True
    for ln in lines:
        if ln.startswith(">"):
            name = ln; continue
        if not ln:
            continue
        match = "".join(ch for ch in ln if not ch.islower())
        if first:
            new = new_query; first = False
        else:
            new = "".join(match[u] if u is not None and u < len(match) else "-" for u in mapping)
            if set(new) <= {"-"}:
                continue
        out += [name, new]
    return "\n".join(out) + "\n"


def write_inputs(case, prot, n_lig_dimer):
    """AF3 JSONs: (A) dimer + n_lig_dimer copies; (B) monomer + 1 copy. Seeds [0]. No templates."""
    unp = remap_a3m(prot["unpairedMsa"], case["mapping"], case["seq"])
    pair = remap_a3m(prot["pairedMsa"], case["mapping"], case["seq"]) if prot.get("pairedMsa") else ""
    paths = {}
    for cond, nprot, nlig in (("dimer", 2, n_lig_dimer), ("monomer", 1, 1)):
        ids = iter("ABCDEFGH")
        pid = [next(ids) for _ in range(nprot)]
        lid = [next(ids) for _ in range(nlig)]
        name = f"{case['case_id']}_{cond}"
        js = {"name": name, "modelSeeds": [0], "dialect": "alphafold3", "version": 2,
              "sequences": [{"protein": {"id": pid, "sequence": case["seq"], "modifications": [],
                                         "unpairedMsa": unp, "pairedMsa": pair, "templates": []}},
                            {"ligand": {"id": lid, "ccdCodes": [case["ccd"]]}}]}
        p = WORK / "inputs" / f"{name}.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(js))
        paths[cond] = str(p.relative_to(ROOT))
    return paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-tokens", type=int, default=MAX_TOKENS)
    a = ap.parse_args()
    c = candidate_sites()
    log("candidate collapsed sites", len(c), c.stratum.value_counts().to_dict())
    c = c[(c.n_components == 1) & ~c.is_covalent.fillna(False).astype(bool) & c.cluster30.notna()]
    c["post_cutoff"] = c.release_date > AF3_CUTOFF
    meta = uniprot_meta(c.dominant_uniprot.unique())
    c["taxid"] = c.dominant_uniprot.map(lambda x: meta.get(x, (None, None))[0])
    c = c[c.taxid == HUMAN_TAXID]
    c["unp_len"] = c.dominant_uniprot.map(lambda x: len(meta[x][1]))
    log("human single-CCD non-covalent", len(c), c.stratum.value_counts().to_dict(),
        "clusters", c.groupby("stratum").cluster30.nunique().to_dict())
    c = c.sort_values(["stratum", "cluster30", "post_cutoff", "resolution", "pdb_id"],
                      ascending=[True, True, False, True, True])
    rows, rejects, used_clusters = [], [], set()
    msa_cache = {}

    def try_stratum(stratum, n_target, order_key):
        sub = c[(c.stratum == stratum) & ~c.cluster30.isin(used_clusters)]
        # cluster-level order: best candidate per cluster determines order
        reps = sub.groupby("cluster30", sort=False).head(1).copy()
        reps["_k"] = reps.apply(order_key, axis=1)
        reps = reps.sort_values(["post_cutoff", "_k"], ascending=[False, True])
        got = 0
        for cl in reps.cluster30:
            if got >= n_target:
                break
            for _, r in sub[sub.cluster30 == cl].head(5).iterrows():  # up to 5 entries per cluster tried
                acc = r.dominant_uniprot
                if acc not in msa_cache:
                    msa_cache[acc] = cached_msa(acc)
                prot = msa_cache[acc]
                if prot is None:
                    rejects.append((r.ligand_uid, stratum, "no_cached_msa")); continue
                if prot["sequence"] != meta[acc][1]:
                    rejects.append((r.ligand_uid, stratum, "cached_msa_seq_ne_uniprot")); continue
                if 2 * r.unp_len > 6 * a.max_tokens:   # cheap prefilter, real check below on trimmed SEQRES
                    rejects.append((r.ligand_uid, stratum, "uniprot_too_long")); continue
                try:
                    info = entry_info(r, meta[acc][1], acc)
                except Exception as e:  # noqa: BLE001
                    rejects.append((r.ligand_uid, stratum, f"error:{e}")); continue
                if isinstance(info, str):
                    rejects.append((r.ligand_uid, stratum, info)); continue
                nlig = 1 if stratum == "on_axis_druglike" else info["n_lig_copies"]
                if stratum != "on_axis_druglike" and nlig != 2:
                    rejects.append((r.ligand_uid, stratum, f"lig_copies={nlig}")); continue
                tok = 2 * len(info["seq"]) + nlig * int(r.n_heavy_lig)
                if tok > a.max_tokens:
                    rejects.append((r.ligand_uid, stratum, f"tokens={tok}")); continue
                case = dict(case_id=f"{stratum.split('_')[0]}_{r.pdb_id}_{r.ccd}", stratum=stratum, pdb_id=r.pdb_id,
                            ligand_uid=r.ligand_uid, ccd=r.ccd, uniprot=acc, cluster30=int(r.cluster30),
                            resolution=r.resolution, release_date=r.release_date, post_cutoff=bool(r.post_cutoff),
                            method=r.method_class, n_heavy=int(r.n_heavy_lig), c2_axis_dist=r.c2_axis_dist,
                            on_axis_jaccard=r.on_axis_jaccard, centroid=json.dumps([round(x, 2) for x in r.centroid]),
                            seq_len=len(info["seq"]), seqres_lo=info["seqres_lo"], identity=round(info["identity"], 3),
                            n_lig_copies_pdb=info["n_lig_copies"], n_lig_dimer_run=nlig, tokens_dimer=tok,
                            tokens_monomer=len(info["seq"]) + int(r.n_heavy_lig), pchains=info["pchains"],
                            seq=info["seq"], mapping=info["mapping"],
                            unp_start=next(u for u in info["mapping"] if u is not None) + 1)
                case.update({f"input_{k}": v for k, v in write_inputs(case, prot, nlig).items()})
                rows.append(case)
                used_clusters.add(cl)
                got += 1
                log(stratum, case["case_id"], tok, r.release_date)
                break
        return got

    try_stratum("on_axis_druglike", N_ONAXIS, lambda r: 2 * r.unp_len)
    med = float(np.median([x["tokens_dimer"] for x in rows])) if rows else 600.0
    try_stratum("offaxis_cofactor", N_COF, lambda r: abs(2 * r.unp_len + 2 * r.n_heavy_lig - med))
    try_stratum("sbs_druglike", N_SBS, lambda r: abs(2 * r.unp_len + 2 * r.n_heavy_lig - med))

    pd.DataFrame(rejects, columns=["ligand_uid", "stratum", "reason"]).to_csv(WORK / "rejects.tsv", sep="\t", index=False)
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("no cases selected; see rejects.tsv")
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.drop(columns=["seq", "mapping"]).to_csv(OUT_TSV, sep="\t", index=False)
    df[["case_id", "seq", "mapping"]].assign(mapping=lambda d: d.mapping.map(json.dumps)).to_json(
        WORK / "case_seqmap.jsonl", orient="records", lines=True)
    log(df.groupby(["stratum", "post_cutoff"]).size())
    log("wrote", OUT_TSV)


if __name__ == "__main__":
    main()
