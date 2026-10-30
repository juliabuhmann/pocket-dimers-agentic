#!/usr/bin/env python3
"""Deep dive #2 co-folding: select H4 multi-chain (MBS) ligand sites that the apo dimer pocket finder MISSED or HIT,
plus heterodimer MBS sites, and write AF3 inputs (A: dimer pair + ligand; B: single chain + ligand).

Status: first attempt, UNVERIFIED. Selection is rule-based (constants below); no hand-picked IDs.

Rerun:
    source scripts/env_threads.sh
    .venvs/part1/bin/python scripts/cofold_select_h4.py   # -> results/cofold/h4/cases.tsv, data/cofold/h4/inputs/

Rules (all in code):
- Sites: H4 v2 primary set (data/part2/h4/v2/benchmark_filtered.parquet, lig_set == organic; interface-filtered),
  is_mbs, ligand_class in LIG_CLASSES; census (data/part1/ligand_sites.parquet): one CCD component, not covalent.
- Strata: homo_miss (homodimer, H4 dimer_hit False), homo_hit (dimer_hit True), hetero (heterodimer, either).
- Forced families: the FORCE_N_NEG MBS families with the most negative contribution to the H4 DiD
  (reviews/verify_h4/family_contributions.tsv; these are NAMPT and MAT2A) are tried first among the misses, as asked.
- One case per cluster30 (no cluster reused across strata). Within a cluster: post-cutoff release (> AF3_CUTOFF)
  first, then fewer estimated tokens, then resolution. Across clusters: post-cutoff first, then fewer tokens.
- Size: dimer tokens <= MAX_TOKENS (docs/COFOLDING.md limit). Non-forced cases additionally <= SOFT_TOKENS
  (wall-clock budget of this deep dive: 4 h incl. analysis, one GPU).
- AF3 sequence: PDB entity SEQRES of each chain trimmed to the first..last residue aligned to UniProt (drops tags);
  MSA columns remapped from the cached full-length AF3 MSA (cached query must equal the UniProt REST sequence);
  no templates (cached templates go up to 2025 and would include the holo entries).
- Chains: in assembly 1, the reference ligand instance is located by the census centroid; the two protein chains
  with the most residues within CONTACT_CUTOFF of it and mapping to (acc1, acc2) form the pair. Run B uses the
  chain with the most contacts.
- Ligand copies (run A): bound copies of the CCD in the assembly that contact BOTH pair chains (deduplicated by
  centroid); 2 if >= 2 such copies ("symmetric-paired"), else 1; reduced to 1 if 2 copies exceed SAFE_TOKENS.
  Run B: 1 copy.
"""
import os
import gzip
import io
import json
import re
import sys
import urllib.request
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from cofold_select_onaxis import remap_a3m  # noqa: E402  (deep dive #1 helper, reused unmodified)

H4_TABLE = ROOT / "data/part2/h4/v2/benchmark_filtered.parquet"
CENSUS = ROOT / "data/part1/ligand_sites.parquet"
FAMILY_TABLE = ROOT / "reviews/verify_h4/family_contributions.tsv"
FAILURE_TABLE = ROOT / "reviews/verify_h4/failure_modes_sites.tsv"
PDB_MIRROR = Path(os.environ.get("PDB_MIRROR", "/cluster/scratch/jjaenes/pdb"))  # = src/part1/config.py PDB_MIRROR
MSA_DIR = Path(os.environ.get("AF3_MSA_HUMAN", "/cluster/project/beltrao/shared/25.05_alphafold3_msas_human"))
AF3_CUTOFF = "2021-09-30"           # AF3 training cutoff (Abramson et al. 2024)
LIG_CLASSES = ("drug_like", "cofactor")
MAX_TOKENS = 1280                   # docs/COFOLDING.md: largest bucket without unified memory
SAFE_TOKENS = 1024                  # 1280 bucket OOMed on GPU 1 for 8y55 (NAMPT x2 + 2 ligands, 1034 tokens;
                                    # logs/cofold_h4/af3_gpu1.log, first attempt). Above this, use 1 ligand copy.
SOFT_TOKENS = 768                   # budget cap for non-forced cases (1024 bucket ~12 min per dimer run)
FORCE_N_NEG = 2                     # most negative MBS family contributions (verify_h4 §3: NAMPT, MAT2A)
N_TARGET = {"homo_miss": 10, "homo_hit": 8, "hetero": 6}
CONTACT_CUTOFF = 4.0
MIN_IDENTITY = 0.9
MAX_TRIES_PER_CLUSTER = 5
OUT_TSV = ROOT / "results/cofold/h4/cases.tsv"
WORK = ROOT / "data/cofold/h4"


def log(*a):
    print(*a, flush=True)


def uniprot_seqs(accs):
    cache = WORK / "uniprot_seq.tsv"
    have = pd.read_csv(cache, sep="\t") if cache.exists() else pd.DataFrame(columns=["Entry", "Sequence"])
    need = sorted(set(accs) - set(have.Entry))
    parts = [have]
    for i in range(0, len(need), 100):
        url = ("https://rest.uniprot.org/uniprotkb/accessions?format=tsv&fields=accession,sequence&accessions="
               + ",".join(need[i:i + 100]))
        with urllib.request.urlopen(url, timeout=120) as r:
            parts.append(pd.read_csv(io.StringIO(r.read().decode()), sep="\t"))
    t = pd.concat(parts, ignore_index=True).drop_duplicates("Entry")
    WORK.mkdir(parents=True, exist_ok=True)
    t.to_csv(cache, sep="\t", index=False)
    return dict(zip(t.Entry, t.Sequence))


def cached_msa(acc):
    p = MSA_DIR / f"{acc.lower()}_data.json.gz"
    if not p.exists():
        return None
    with gzip.open(p, "rt") as f:
        d = json.load(f)
    prots = [s["protein"] for s in d["sequences"] if "protein" in s]
    return prots[0] if len(prots) == 1 else None


def _cigar(s):
    return [(op, int(n)) for n, op in re.findall(r"(\d+)([MID])", s)]


def map_entity(seqres, unp):
    al = gemmi.align_string_sequences(list(seqres), list(unp), [])
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
        return None
    lo, hi = min(m), max(m)
    seq = seqres[lo:hi + 1]
    mapping = [m.get(k) for k in range(lo, hi + 1)]
    ident = float(np.mean([seq[k] == unp[u] for k, u in enumerate(mapping) if u is not None]))
    return dict(seq=seq, lo=lo, mapping=mapping, identity=ident)


def heavy_xyz(res):
    return np.array([a.pos.tolist() for a in res if a.element.name != "H" and a.altloc in ("\0", "A")])


def entry_info(r, useq):
    st = gemmi.read_structure(str(PDB_MIRROR / r.pdb_id[1:3] / f"{r.pdb_id}.cif.gz"))
    st.setup_entities()
    if not st.assemblies:
        return "no_assembly"
    a = gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.AddNumber)
    a.remove_hydrogens()
    pch = {ch.name: ch for ch in a if ch.get_polymer().check_polymer_type() == gemmi.PolymerType.PeptideL}
    cen = np.array(list(map(float, r.centroid)))
    ligs = []
    for ch in a:
        for res in ch:
            if res.name == r.ccd and res.het_flag == "H":
                x = heavy_xyz(res)
                if len(x):
                    ligs.append(x)
    if not ligs:
        return "ligand_not_found"
    d = [np.linalg.norm(x.mean(0) - cen) for x in ligs]
    if min(d) > 1.0:
        return f"ligand_centroid_off={min(d):.1f}"
    ref = ligs[int(np.argmin(d))]

    def contacts(chain, xyz):
        n = 0
        for res in chain.get_polymer():
            hx = heavy_xyz(res)
            if len(hx) and np.min(np.linalg.norm(hx[:, None] - xyz[None], axis=-1)) <= CONTACT_CUTOFF:
                n += 1
        return n

    # map chains to accessions by entity sequence
    ent_map = {}
    chain_acc = {}
    for cn, ch in pch.items():
        eid = ch.get_polymer()[0].entity_id
        if eid not in ent_map:
            seqres = gemmi.one_letter_code(st.get_entity(eid).full_sequence).upper()
            best = None
            for acc in {r.acc1, r.acc2}:
                mp = map_entity(seqres, useq[acc])
                if mp and mp["identity"] >= MIN_IDENTITY and (best is None or len(mp["seq"]) > len(best[1]["seq"])):
                    best = (acc, mp)
            ent_map[eid] = best
        if ent_map[eid] is not None:
            chain_acc[cn] = ent_map[eid]
    cnt = {cn: contacts(pch[cn], ref) for cn in chain_acc}
    cnt = {k: v for k, v in cnt.items() if v > 0}
    order = sorted(cnt, key=lambda k: -cnt[k])
    if r.acc1 == r.acc2:
        pair = order[:2]
    else:
        c1 = [c for c in order if chain_acc[c][0] == r.acc1][:1]
        c2 = [c for c in order if chain_acc[c][0] == r.acc2][:1]
        pair = sorted(c1 + c2, key=lambda k: -cnt[k])
    if len(pair) != 2:
        return "pair_not_contacting"
    # ligand copies contacting both chains (dedupe by centroid 2 A)
    both = []
    for x in ligs:
        if all(contacts(pch[c], x) > 0 for c in pair):
            c = x.mean(0)
            if all(np.linalg.norm(c - u) > 2.0 for u in both):
                both.append(c)
    chains = []
    for cn in pair:
        acc, mp = chain_acc[cn]
        chains.append(dict(chain=cn, acc=acc, seq=mp["seq"], lo=mp["lo"], mapping=mp["mapping"],
                           identity=mp["identity"], n_contacts=cnt[cn]))
    return dict(chains=chains, n_copies_both=len(both))


def write_inputs(case, chains, msas, n_lig):
    """AF3 JSONs: dimer (chains A,B) + n_lig copies; monomer (chain with most contacts = chains[0]) + 1 copy."""
    ents = []
    for c in chains:
        prot = msas[c["acc"]]
        ents.append(dict(seq=c["seq"], unp=remap_a3m(prot["unpairedMsa"], c["mapping"], c["seq"]),
                         pair=remap_a3m(prot["pairedMsa"], c["mapping"], c["seq"]) if prot.get("pairedMsa") else ""))
    paths = {}
    for cond, use in (("dimer", ents), ("monomer", ents[:1])):
        ids = iter("ABCDEFGH")
        seqs = []
        if cond == "dimer" and ents[0]["seq"] == ents[1]["seq"] and chains[0]["acc"] == chains[1]["acc"]:
            groups = [(ents[0], [next(ids), next(ids)])]
        else:
            groups = [(e, [next(ids)]) for e in use]
        for e, pid in groups:
            seqs.append({"protein": {"id": pid, "sequence": e["seq"], "modifications": [], "unpairedMsa": e["unp"],
                                     "pairedMsa": e["pair"], "templates": []}})
        nl = n_lig if cond == "dimer" else 1
        seqs.append({"ligand": {"id": [next(ids) for _ in range(nl)], "ccdCodes": [case["ccd"]]}})
        name = f"{case['case_id']}_{cond}"
        js = {"name": name, "modelSeeds": [0], "dialect": "alphafold3", "version": 2, "sequences": seqs}
        p = WORK / "inputs" / f"{case['priority']:02d}_{name}.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(js))
        paths[cond] = str(p.relative_to(ROOT))
    return paths


def main():
    b = pd.read_parquet(H4_TABLE)
    b = b[(b.lig_set == "organic") & b.is_mbs & b.ligand_class.isin(LIG_CLASSES)].copy()
    log("H4 primary MBS drug-like/cofactor sites", len(b), "clusters", b.cluster30.nunique())
    cen = pd.read_parquet(CENSUS, columns=["ligand_uid", "n_components", "is_covalent"],
                          filters=[("ligand_uid", "in", sorted(b.ligand_uid.unique()))])
    b = b.merge(cen, on="ligand_uid", how="left")
    b = b[(b.n_components == 1) & ~b.is_covalent.fillna(False).astype(bool)]
    fm = pd.read_csv(FAILURE_TABLE, sep="\t", usecols=["ligand_uid", "model", "fm_dimer", "fm_iso"])
    b = b.merge(fm, on=["ligand_uid", "model"], how="left")
    b["stratum"] = np.where(b.kind == "hetero", "hetero", np.where(b.dimer_hit, "homo_hit", "homo_miss"))
    b["post_cutoff"] = b.release_date > AF3_CUTOFF
    useq = uniprot_seqs(set(b.acc1) | set(b.acc2))
    b["est_tokens"] = b.acc1.map(lambda a: len(useq.get(a, ""))) + b.acc2.map(lambda a: len(useq.get(a, ""))) + b.n_heavy
    log("after single-component / non-covalent", len(b), b.groupby("stratum").cluster30.nunique().to_dict())

    fam = pd.read_csv(FAMILY_TABLE, sep="\t")
    fam = fam[fam.is_mbs & (fam.contrib < 0)].sort_values("contrib").head(FORCE_N_NEG)
    force_accs = set(",".join(fam.accs).split(","))
    log("forced families (most negative MBS DiD contribution):", fam[["boot_cluster", "accs", "gain"]].to_dict("records"))
    b["forced"] = b.acc1.isin(force_accs) & (b.stratum == "homo_miss")

    b = b.sort_values(["cluster30", "forced", "post_cutoff", "est_tokens", "resolution", "pdb_id"],
                      ascending=[True, False, False, True, True, True])
    msas, rows, rejects, used = {}, [], [], set()

    def get_msa(acc):
        if acc not in msas:
            p = cached_msa(acc)
            msas[acc] = p if (p is not None and p["sequence"] == useq.get(acc)) else None
        return msas[acc]

    for stratum in ("homo_miss", "homo_hit", "hetero"):
        sub = b[(b.stratum == stratum) & ~b.cluster30.isin(used)]
        reps = sub.groupby("cluster30", sort=False).head(1)
        reps = reps.sort_values(["forced", "post_cutoff", "est_tokens"], ascending=[False, False, True])
        got = 0
        for cl, forced in zip(reps.cluster30, reps.forced):
            if got >= N_TARGET[stratum]:
                break
            cap = MAX_TOKENS if forced else SOFT_TOKENS
            cand = sub[sub.cluster30 == cl]
            if forced:   # keep the forced protein (cluster may mix proteins), still post-cutoff first
                cand = cand[cand.forced]
            for _, r in cand.head(MAX_TRIES_PER_CLUSTER).iterrows():
                if r.est_tokens > 2 * cap:
                    rejects.append((r.ligand_uid, stratum, f"est_tokens={r.est_tokens}")); continue
                if get_msa(r.acc1) is None or get_msa(r.acc2) is None:
                    rejects.append((r.ligand_uid, stratum, "no_cached_msa_or_seq_ne_uniprot")); continue
                try:
                    info = entry_info(r, useq)
                except Exception as e:  # noqa: BLE001
                    rejects.append((r.ligand_uid, stratum, f"error:{e}")); continue
                if isinstance(info, str):
                    rejects.append((r.ligand_uid, stratum, info)); continue
                ch = info["chains"]
                nlig = 2 if info["n_copies_both"] >= 2 else 1
                tok = len(ch[0]["seq"]) + len(ch[1]["seq"]) + nlig * int(r.n_heavy)
                reduced = False
                if nlig == 2 and tok > SAFE_TOKENS:   # memory fallback (see SAFE_TOKENS)
                    nlig, tok, reduced = 1, tok - int(r.n_heavy), True
                if tok > cap:
                    rejects.append((r.ligand_uid, stratum, f"tokens={tok}>{cap}")); continue
                case = dict(case_id=f"{stratum}_{r.pdb_id}_{r.ccd}", priority=0, stratum=stratum, forced=bool(forced),
                            pdb_id=r.pdb_id, ligand_uid=r.ligand_uid, ccd=r.ccd, ligand_class=r.ligand_class,
                            n_heavy=int(r.n_heavy), kind=r.kind, acc1=r.acc1, acc2=r.acc2, model=r.model,
                            cluster30=int(r.cluster30), resolution=r.resolution, release_date=r.release_date,
                            post_cutoff=bool(r.post_cutoff), post_cutoff_pair=bool(r.post_cutoff_pair),
                            iface_jaccard=r.iface_jaccard, apo_dimer_dcc=r.dimer_dcc, apo_dimer_hit=bool(r.dimer_hit),
                            apo_iso_dcc=r.iso_dcc, apo_iso_hit=bool(r.iso_hit), apo_dimer_recall=r.dimer_recall,
                            apo_dimer_pksize=r.dimer_pksize, fm_dimer=r.fm_dimer, fm_iso=r.fm_iso,
                            centroid=json.dumps([round(float(x), 3) for x in r.centroid]),
                            chainA_pdb=ch[0]["chain"], chainB_pdb=ch[1]["chain"],
                            accA=ch[0]["acc"], accB=ch[1]["acc"],
                            loA=ch[0]["lo"], loB=ch[1]["lo"], lenA=len(ch[0]["seq"]), lenB=len(ch[1]["seq"]),
                            identA=round(ch[0]["identity"], 3), identB=round(ch[1]["identity"], 3),
                            n_contacts_A=ch[0]["n_contacts"], n_contacts_B=ch[1]["n_contacts"],
                            n_copies_both_pdb=info["n_copies_both"], n_lig_dimer_run=nlig, lig_copies_reduced=reduced, tokens_dimer=tok,
                            tokens_monomer=len(ch[0]["seq"]) + int(r.n_heavy), _chains=ch)
                rows.append(case)
                used.add(cl)
                got += 1
                log(stratum, case["case_id"], "tokens", tok, r.release_date, "forced" if forced else "")
                break

    # run priority: forced first, then round-robin over strata (balanced partial results if the budget runs out),
    # smaller cases first within stratum
    by = {s: sorted([c for c in rows if c["stratum"] == s and not c["forced"]], key=lambda c: c["tokens_dimer"])
          for s in N_TARGET}
    order = [c for c in rows if c["forced"]]
    while any(by.values()):
        for s in N_TARGET:
            if by[s]:
                order.append(by[s].pop(0))
    for i, c in enumerate(order):
        c["priority"] = i
        c.update({f"input_{k}": v for k, v in write_inputs(c, c["_chains"], msas, c["n_lig_dimer_run"]).items()})
    pd.DataFrame(rejects, columns=["ligand_uid", "stratum", "reason"]).to_csv(WORK / "rejects.tsv", sep="\t", index=False)
    df = pd.DataFrame(order).drop(columns=["_chains"])
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    log(df.groupby(["stratum", "post_cutoff"]).size())
    log("total dimer tokens", df.tokens_dimer.sum(), "wrote", OUT_TSV)


if __name__ == "__main__":
    main()
