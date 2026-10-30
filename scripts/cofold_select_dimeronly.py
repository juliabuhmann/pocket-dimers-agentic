#!/usr/bin/env python3
"""Deep dive #3 co-folding: does AF3 need the partner chain to place the natural ligand / cofactor / drug in the
curated DIMER-ONLY pockets of human AF dimer models (results/part2/q23a_curated.tsv, docs/Q22_RESULTS.md section 6)?
Selects PDB ground-truth cases and writes AF3 inputs: (A) dimer + PDB ligand copies, (B) single chain + 1 ligand.

Status: first attempt, UNVERIFIED. Selection is rule-based (constants below); no hand-picked IDs.

Rerun:
    source scripts/env_threads.sh
    .venvs/part1/bin/python scripts/cofold_select_dimeronly.py   # -> results/cofold/dimeronly/cases.tsv, data/cofold/dimeronly/inputs/

Rules (all in code):
- Curated pocket -> canonical UniProt: isoform suffix stripped; an unreviewed (TrEMBL) accession is replaced by the
  unique reviewed human accession with the same primary gene name (UniProt REST). Rows without a gene are dropped.
- Pocket residues (dimer-model numbering) are mapped to canonical UniProt positions by aligning each model chain's
  sequence (from the AF dimer PDB file) to the canonical sequence; only identical aligned residues are mapped.
- PDB ligand sites = census v2 (data/part1/ligand_sites.parquet) rows with in_primary_set, n_components == 1, not
  covalent, ligand_class not in EXCLUDE_CLASSES, contacting both proteins of the pair (homodimer: >= 2 chains of
  the protein; heterodimer: both accessions), i.e. a multi-chain site between the two proteins of the curated pair.
- Overlap = overlap coefficient |L & P| / min(|L|, |P|) between the ligand's contact residues L and the curated pocket
  residues P, both as (accession, UniProt position); for a homodimer the two chains are pooled (chain-agnostic,
  because the symmetric copy of a site has the same positions). Keep overlap >= OVERLAP_MIN.
- Family = census cluster30 of the site; at most one case per family and per canonical protein pair.
  Within a family: post-AF3-cutoff release preferred, then higher overlap, then better resolution, then pdb_id.
  Families are taken in curated-list rank order (the curated score order), up to N_CASES.
- Reference chain pair (assembly 1): for the census instance (matched by CCD + centroid), the two protein chains of
  the pair (homo: two chains of the accession; hetero: one per accession) with most residues within CONTACT of it.
- Ligand stoichiometry (A) = bound copies of the same CCD (heavy atom within CONTACT of the chain pair, majority of
  contacting protein atoms on the pair; copies deduplicated by centroid within DEDUP) -> often 2 symmetric sites.
- AF3 protein sequence = PDB entity SEQRES trimmed to the first..last residue aligned to UniProt (drops tags), MSA
  columns remapped from the cached AF3 MSA (cached query must equal the UniProt sequence); no templates (avoid holo
  leakage), as in scripts/cofold_select_onaxis.py (helpers imported from there).
- Tokens (A) = sum of chain lengths + copies * heavy atoms <= MAX_TOKENS. (B) single chain = the pair chain with more
  contacts to the reference instance, + 1 ligand copy.
"""
import os
import gzip  # noqa: F401  (kept for parity with onaxis helpers)
import io
import json
import sys
import urllib.parse
import urllib.request
import urllib.error
import time
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
import cofold_select_onaxis as OA  # noqa: E402  reuse: cached_msa, remap_a3m, _cigar
from part1.config import PDB_MIRROR  # noqa: E402

CURATED = ROOT / "results/part2/q23a_curated.tsv"
CENSUS = ROOT / "data/part1/ligand_sites.parquet"
DIMER_ROOT = Path(os.environ.get("DIMERS_ROOT", "/path/to/afdb_dimers")) / "afdb_human_dimers"  # = scripts/annotate_pocket_chains.py
KIND_DIR = {"homo": "homodimers", "hetero": "heterodimers"}
AF3_CUTOFF = "2021-09-30"     # AF3 training data cutoff (Abramson et al. 2024)
OVERLAP_MIN = 0.3             # task definition
MAX_TOKENS = 1280             # docs/COFOLDING.md: largest bucket without unified memory
N_CASES = 12                  # task target
CONTACT = 4.0                 # A, heavy-atom contact (as census / onaxis)
DEDUP = 2.0                   # A, ligand copies with centroids closer than this are one copy (altloc/duplicates)
MIN_IDENT = 0.9               # PDB entity vs UniProt identity over aligned residues
EXCLUDE_CLASSES = {"metal_ion", "inorganic", "metal_cluster"}
HUMAN_TAXID = 9606
OUT = ROOT / "results/cofold/dimeronly"
WORK = ROOT / "data/cofold/dimeronly"


def log(*a):
    print(*a, flush=True)


def rest_tsv(query, fields="accession,reviewed,gene_primary,organism_id,sequence"):
    url = ("https://rest.uniprot.org/uniprotkb/search?format=tsv&size=500&fields=" + fields + "&query="
           + urllib.parse.quote(query))
    for attempt in range(4):   # UniProt REST returns sporadic HTTP 500s
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                txt = r.read().decode()
            return pd.read_csv(io.StringIO(txt), sep="\t") if txt.strip() else pd.DataFrame(columns=["Entry", "Reviewed"])
        except urllib.error.HTTPError:
            if attempt == 3:
                raise
            time.sleep(5 * (attempt + 1))


def canonical_map(cur):
    """curated accession -> (canonical reviewed accession, sequence); cached in WORK/canonical.tsv."""
    cache = WORK / "canonical.tsv"
    if cache.exists():
        t = pd.read_csv(cache, sep="\t")
    else:
        rows = []
        for acc, gene in sorted(set(zip(cur.uniprot_A, cur.gene_A)) | set(zip(cur.uniprot_B, cur.gene_B))):
            base = acc.split("-")[0]
            q = rest_tsv(f"accession:{base}")
            how = "base_reviewed"
            if q.empty or q.Reviewed.iloc[0] != "reviewed":
                if not isinstance(gene, str):
                    rows.append(dict(acc=acc, canon=None, how="no_gene")); continue
                q = rest_tsv(f"gene_exact:{gene} AND organism_id:{HUMAN_TAXID} AND reviewed:true")
                how = "gene_reviewed"
            if len(q) != 1:
                rows.append(dict(acc=acc, canon=None, how=f"n_hits={len(q)}")); continue
            rows.append(dict(acc=acc, canon=q.Entry.iloc[0], gene_canon=q["Gene Names (primary)"].iloc[0],
                             taxid=int(q["Organism (ID)"].iloc[0]), how=how, seq=q.Sequence.iloc[0]))
        t = pd.DataFrame(rows)
        WORK.mkdir(parents=True, exist_ok=True)
        t.to_csv(cache, sep="\t", index=False)
    return t.set_index("acc")


def model_chain_seqs(model, kind):
    st = gemmi.read_structure(str(DIMER_ROOT / KIND_DIR[kind] / f"{model}.pdb"))
    return {ch.name: ([r.seqid.num for r in ch], gemmi.one_letter_code([r.name for r in ch]).upper()) for ch in st[0]}


def align_map(qseq, tseq):
    """index in qseq -> index in tseq (0-based) for identical aligned residues."""
    al = gemmi.align_string_sequences(list(qseq), list(tseq), [])
    m, i, j = {}, 0, 0
    for op, n in OA._cigar(al.cigar_str()):
        if op == "M":
            for k in range(n):
                if qseq[i + k] == tseq[j + k]:
                    m[i + k] = j + k
            i += n; j += n
        elif op == "I":
            i += n
        else:
            j += n
    return m


def pocket_unp(row, canon):
    """curated pocket -> {'A': set(unp_pos 1-based), 'B': ...} on canonical sequences + model chain ids."""
    seqs = model_chain_seqs(row.modelEntityId, row.kind)
    chains = sorted(seqs)
    out = {}
    for role, ch, acc in (("A", chains[0], row.uniprot_A), ("B", chains[1], row.uniprot_B)):
        nums, s = seqs[ch]
        m = align_map(s, canon.loc[acc, "seq"])
        n2i = {n: i for i, n in enumerate(nums)}
        res = [int(x.split(":")[1]) for x in row.residues.split() if x.split(":")[0] == role]
        out[role] = {m[n2i[n]] + 1 for n in res if n in n2i and n2i[n] in m}
    return out


def census_sites(accs):
    cache = WORK / "census_sites.parquet"
    if cache.exists():
        t = pd.read_parquet(cache)
        t["contacts"] = t.contacts.map(json.loads)
        return t
    cols = ["pdb_id", "ligand_uid", "ccd", "dominant_uniprot", "is_mbs", "n_components", "is_covalent", "ligand_class",
            "cluster30", "resolution", "release_date", "centroid", "contacts", "n_heavy", "n_heavy_atoms_instance",
            "method_class"]
    f = (pc.field("in_primary_set") == True) & pc.field("dominant_uniprot").isin(sorted(accs))  # noqa: E712
    t = ds.dataset(CENSUS).to_table(columns=cols, filter=f).to_pandas()
    t["centroid"] = t.centroid.apply(lambda x: [float(v) for v in x])
    t["contacts"] = t.contacts.apply(lambda cs: [(c["chain"], c["unp_acc"], int(c["unp_pos"]) if c["unp_pos"] is not None
                                                  else None) for c in cs if c["kind"] == "protein"] if cs is not None else [])
    t.assign(contacts=t.contacts.map(json.dumps)).to_parquet(cache)
    return t


# ------------------------------------------------------------------ reference structure (assembly 1)
def entity_map(st, ent, seqs):
    """entity -> (accession, seqres_lo, trimmed seq, mapping) for the canonical sequence it matches best."""
    seqres = gemmi.one_letter_code(ent.full_sequence).upper()
    best = None
    for acc, us in seqs.items():
        m = align_map(seqres, us)
        if not m:
            continue
        lo, hi = min(m), max(m)
        mapping = [m.get(k) for k in range(lo, hi + 1)]
        ident = len(m) / (hi - lo + 1)
        if best is None or ident > best[4]:
            best = (acc, lo, seqres[lo:hi + 1], mapping, ident)
    return best


def entry_info(site, pair_accs, seqs):
    st = gemmi.read_structure(str(PDB_MIRROR / site.pdb_id[1:3] / f"{site.pdb_id}.cif.gz"))
    st.setup_entities()
    if not st.assemblies:
        return "no_assembly"
    a = gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.AddNumber)
    a.remove_hydrogens()
    pchains = {}
    for ch in a:
        pol = ch.get_polymer()
        if len(pol) and pol.check_polymer_type() == gemmi.PolymerType.PeptideL:
            ent = st.get_entity(pol[0].entity_id)
            em = entity_map(st, ent, {k: seqs[k] for k in set(pair_accs)})
            if em and em[4] >= MIN_IDENT and len(em[2]) >= 20:
                pchains[ch.name] = em
    # ligand instances of this CCD
    inst = []
    for ch in a:
        for r in ch:
            if r.name == site.ccd and r.het_flag == "H":
                xyz = np.array([at.pos.tolist() for at in r if at.altloc in ("\0", "A")])
                if len(xyz):
                    inst.append((ch.name, r, xyz))
    if not inst:
        return "no_ligand_instance"
    cen = np.array(site.centroid)
    ref = min(inst, key=lambda x: np.linalg.norm(x[2].mean(0) - cen))
    if np.linalg.norm(ref[2].mean(0) - cen) > DEDUP:
        return "centroid_mismatch"
    ns = gemmi.NeighborSearch(a, st.cell, 6).populate()

    def contacts(xyz):
        cnt = {}
        for p in xyz:
            for mk in ns.find_atoms(gemmi.Position(*p), "\0", radius=CONTACT):
                cra = mk.to_cra(a)
                if cra.chain.name in pchains and cra.residue.het_flag != "H":
                    cnt.setdefault(cra.chain.name, set()).add(cra.residue.seqid.num)
        return cnt

    cref = contacts(ref[2])
    if len(set(pair_accs)) == 1:
        cand = sorted(cref, key=lambda c: -len(cref[c]))[:2]
    else:
        cand = []
        for acc in pair_accs:
            cs = [c for c in cref if pchains[c][0] == acc]
            if cs:
                cand.append(max(cs, key=lambda c: len(cref[c])))
    if len(cand) != 2:
        return f"ref_instance_contacts_{len(cand)}_pair_chains"
    # order chains as (A acc, B acc); homo: dominant first
    if len(set(pair_accs)) == 2:
        cand = sorted(cand, key=lambda c: pair_accs.index(pchains[c][0]))
    copies = []
    for _, r, xyz in inst:
        c = contacts(xyz)
        on = sum(len(c.get(x, ())) for x in cand)
        tot = sum(len(v) for v in c.values())
        if on and on * 2 > tot and all(np.linalg.norm(xyz.mean(0) - u) > DEDUP for u in copies):
            copies.append(xyz.mean(0))
    mono = max(cand, key=lambda c: len(cref.get(c, ())))
    return dict(chains=cand, chain_info={c: pchains[c] for c in cand}, n_copies=len(copies), mono_chain=mono,
                n_ref_contacts=[len(cref.get(c, ())) for c in cand], n_heavy=int(len(ref[2])))


def write_inputs(case, chain_rows, prots):
    """(A) dimer + n copies; (B) the monomer chain + 1 copy. Seeds [0]. No templates."""
    msas = {}
    for cr in chain_rows:
        p = prots[cr["acc"]]
        msas[cr["chain"]] = (OA.remap_a3m(p["unpairedMsa"], cr["mapping"], cr["seq"]),
                             OA.remap_a3m(p["pairedMsa"], cr["mapping"], cr["seq"]) if p.get("pairedMsa") else "")
    paths = {}
    for cond in ("dimer", "monomer"):
        crs = chain_rows if cond == "dimer" else [c for c in chain_rows if c["chain"] == case["mono_chain"]]
        nlig = case["n_lig_dimer_run"] if cond == "dimer" else 1
        ids = iter("ABCDEFGHIJ")
        seqs = []
        # identical sequences -> one entity with several ids (as onaxis)
        groups = {}
        for cr in crs:
            groups.setdefault(cr["seq"], []).append(cr)
        for s, g in groups.items():
            u, pmsa = msas[g[0]["chain"]]
            seqs.append({"protein": {"id": [next(ids) for _ in g], "sequence": s, "modifications": [],
                                     "unpairedMsa": u, "pairedMsa": pmsa, "templates": []}})
        seqs.append({"ligand": {"id": [next(ids) for _ in range(nlig)], "ccdCodes": [case["ccd"]]}})
        name = f"{case['case_id']}_{cond}"
        js = {"name": name, "modelSeeds": [0], "dialect": "alphafold3", "version": 2, "sequences": seqs}
        p = WORK / "inputs" / f"{name}.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(js))
        paths[cond] = str(p.relative_to(ROOT))
    return paths


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    cur = pd.read_csv(CURATED, sep="\t", comment="#")
    canon = canonical_map(cur)
    ok = canon.canon.notna() & (canon.get("taxid") == HUMAN_TAXID)
    log("curated rows", len(cur), "accessions", len(canon), "canonical mapped", int(ok.sum()))
    cur = cur[cur.uniprot_A.map(ok) & cur.uniprot_B.map(ok)].copy()
    cur["canon_A"] = cur.uniprot_A.map(canon.canon)
    cur["canon_B"] = cur.uniprot_B.map(canon.canon)
    seqs = dict(zip(canon.canon.dropna(), canon.seq.dropna()))
    sites = census_sites(set(cur.canon_A) | set(cur.canon_B))
    sites = sites[(sites.n_components == 1) & ~sites.is_covalent.fillna(False).astype(bool)
                  & ~sites.ligand_class.isin(EXCLUDE_CLASSES) & sites.cluster30.notna()]
    log("census primary-set candidate sites on curated proteins", len(sites))
    rows = []
    for _, r in cur.iterrows():
        P = pocket_unp(r, canon)
        accA, accB = r.canon_A, r.canon_B
        Pset = {(accA, p) for p in P["A"]} | {(accB, p) for p in P["B"]}
        sub = sites[sites.dominant_uniprot.isin([accA, accB])]
        for _, s in sub.iterrows():
            L = {(a, p) for _, a, p in s.contacts if a in (accA, accB) and p is not None}
            chains_by_acc = {}
            for ch, a, _ in s.contacts:
                chains_by_acc.setdefault(a, set()).add(ch)
            pair_ok = (len(chains_by_acc.get(accA, ())) >= 2 if accA == accB
                       else bool(chains_by_acc.get(accA)) and bool(chains_by_acc.get(accB)))
            if not L or not Pset:
                continue
            ov = len(L & Pset) / min(len(L), len(Pset))
            rows.append(dict(rank=r["rank"], pocket_id=r.pocket_id, kind=r.kind, gene_A=r.gene_A, gene_B=r.gene_B,
                             canon_A=accA, canon_B=accB, pdb_id=s.pdb_id, ligand_uid=s.ligand_uid, ccd=s.ccd,
                             ligand_class=s.ligand_class, cluster30=int(s.cluster30), resolution=s.resolution,
                             release_date=s.release_date, method=s.method_class, overlap=ov, pair_contact=pair_ok,
                             n_lig_res=len(L), n_pocket_res=len(Pset), centroid=s.centroid,
                             post_cutoff=s.release_date > AF3_CUTOFF))
    cand = pd.DataFrame(rows)
    cand.drop(columns=["centroid"]).to_csv(WORK / "candidates_all.tsv", sep="\t", index=False, float_format="%.3f")
    log("site x pocket pairs", len(cand), "overlap>=min", int((cand.overlap >= OVERLAP_MIN).sum()),
        "and pair contact", int(((cand.overlap >= OVERLAP_MIN) & cand.pair_contact).sum()))
    cand = cand[(cand.overlap >= OVERLAP_MIN) & cand.pair_contact]
    cand = cand.sort_values(["rank", "post_cutoff", "overlap", "resolution", "pdb_id"],
                            ascending=[True, False, False, True, True])
    msa = {}
    for acc in set(cand.canon_A) | set(cand.canon_B):
        msa[acc] = OA.cached_msa(acc)
    cases, rejects, used_fam, used_pair = [], [], set(), set()
    for (rank, cl), g in cand.groupby(["rank", "cluster30"], sort=False):
        if len(cases) >= N_CASES:
            break
        pair = tuple(sorted((g.canon_A.iloc[0], g.canon_B.iloc[0])))
        if cl in used_fam or pair in used_pair:
            continue
        for _, s in g.head(8).iterrows():   # up to 8 entries per (pocket, family) tried
            accs = [s.canon_A, s.canon_B]
            bad = [a for a in set(accs) if msa[a] is None or msa[a]["sequence"] != seqs[a]]
            if bad:
                rejects.append((s.ligand_uid, s.pocket_id, f"msa_missing_or_ne_uniprot:{bad}")); continue
            try:
                info = entry_info(s, accs, seqs)
            except Exception as e:  # noqa: BLE001
                rejects.append((s.ligand_uid, s.pocket_id, f"error:{e}")); continue
            if isinstance(info, str):
                rejects.append((s.ligand_uid, s.pocket_id, info)); continue
            nlig = info["n_copies"]
            crs = []
            for c in info["chains"]:
                acc, lo, seq, mapping, ident = info["chain_info"][c]
                crs.append(dict(chain=c, acc=acc, seqres_lo=lo, seq=seq, mapping=mapping, identity=round(ident, 3)))
            tok = sum(len(c["seq"]) for c in crs) + nlig * info["n_heavy"]
            if tok > MAX_TOKENS:
                rejects.append((s.ligand_uid, s.pocket_id, f"tokens={tok}")); continue
            mono = next(c for c in crs if c["chain"] == info["mono_chain"])
            gene = s.gene_A if s.kind == "homo" else f"{s.gene_A}-{s.gene_B}"
            case = dict(case_id=f"{gene}_{s.pdb_id}_{s.ccd}", rank=int(rank), pocket_id=s.pocket_id, kind=s.kind,
                        gene_A=s.gene_A, gene_B=s.gene_B, uniprot_A=s.canon_A, uniprot_B=s.canon_B, pdb_id=s.pdb_id,
                        ligand_uid=s.ligand_uid, ccd=s.ccd, ligand_class=s.ligand_class, cluster30=int(cl),
                        resolution=s.resolution, release_date=s.release_date, post_cutoff=bool(s.post_cutoff),
                        method=s.method, overlap=round(s.overlap, 3), n_lig_res=s.n_lig_res,
                        n_pocket_res=s.n_pocket_res, n_heavy=info["n_heavy"], n_lig_copies_pdb=nlig,
                        n_lig_dimer_run=nlig, tokens_dimer=tok, tokens_monomer=len(mono["seq"]) + info["n_heavy"],
                        pchains=",".join(info["chains"]), mono_chain=info["mono_chain"],
                        ref_contacts=",".join(map(str, info["n_ref_contacts"])),
                        seq_len=",".join(str(len(c["seq"])) for c in crs),
                        centroid=json.dumps([round(x, 2) for x in s.centroid]))
            if nlig < 1:
                rejects.append((s.ligand_uid, s.pocket_id, "no_bound_copy")); continue
            case.update({f"input_{k}": v for k, v in write_inputs(case, crs, msa).items()})
            case["_chains"] = crs
            cases.append(case)
            used_fam.add(cl); used_pair.add(pair)
            log("case", case["case_id"], "tokens", tok, "copies", nlig, s.release_date, "overlap", case["overlap"])
            break
    pd.DataFrame(rejects, columns=["ligand_uid", "pocket_id", "reason"]).to_csv(WORK / "rejects.tsv", sep="\t", index=False)
    if not cases:
        sys.exit("no cases")
    df = pd.DataFrame(cases)
    OUT.mkdir(parents=True, exist_ok=True)
    df.drop(columns=["_chains"]).to_csv(OUT / "cases.tsv", sep="\t", index=False)
    with open(WORK / "case_chains.json", "w") as f:
        json.dump({c["case_id"]: c["_chains"] for c in cases}, f)
    # curated pocket residues per case on canonical positions (for the metrics' curated-pocket DCC)
    pk = {}
    for c in cases:
        r = cur[cur.pocket_id == c["pocket_id"]].iloc[0]
        P = pocket_unp(r, canon)
        pk[c["case_id"]] = {"A": [c["uniprot_A"], sorted(P["A"])], "B": [c["uniprot_B"], sorted(P["B"])]}
    (WORK / "case_pockets.json").write_text(json.dumps(pk))
    log(df[["case_id", "rank", "tokens_dimer", "tokens_monomer", "n_lig_dimer_run", "post_cutoff", "overlap"]].to_string())


if __name__ == "__main__":
    main()
