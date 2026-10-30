#!/usr/bin/env python
"""ChEMBL target status per UniProt accession via the ChEMBL REST API (for the Q2.3a candidate proteins).

Why REST and not the local SQLite (data/raw/chembl_sqlite/.../chembl_37.db): the activities join over the 30 GB
SQLite on the shared filesystem ran >35 min IO-bound (2026-10-29) and was stopped; REST is the documented
alternative in DESIGN.md §0.
Endpoints (https://www.ebi.ac.uk/chembl/api/data):
  target.json?target_components__accession=<acc>           -> targets containing the accession (type, chembl id)
  activity.json?target_chembl_id=<id>&pchembl_value__gte=6  -> page_meta.total_count (activity records, SINGLE PROTEIN)
  mechanism.json?target_chembl_id__in=<ids>                 -> molecules with a mechanism on any containing target
  molecule.json?molecule_chembl_id__in=<ids>                -> max_phase of those molecules
status: drug_target (mechanism molecule max_phase >= 1) > potent_ligands (>=1 pChEMBL>=6 record on a single-protein
target) > chembl_target > not_in_chembl.  Release: from status.json (chembl_db_version).
Input: accessions on stdin-less mode = all proteins in data/interim/q22/q23a_candidates.txt
Output: data/interim/q22/chembl_target_status.tsv
Status: first attempt, not yet verified.
"""
import sys, time
from concurrent.futures import ThreadPoolExecutor
import pandas as pd, requests
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from q22_common import INTERIM

API = "https://www.ebi.ac.uk/chembl/api/data"
PCHEMBL_MIN = 6  # same potency threshold as the DESIGN.md §0 drug-like definition


def get(path, **params):
    for k in range(4):
        try:
            r = requests.get(f"{API}/{path}", params={**params, "format": "json"}, timeout=60)
            r.raise_for_status(); return r.json()
        except Exception:
            time.sleep(2 ** k)
    raise RuntimeError(f"failed {path} {params}")


def one(acc):
    t = get("target.json", target_components__accession=acc, limit=1000)["targets"]
    sp = [x["target_chembl_id"] for x in t if x["target_type"] == "SINGLE PROTEIN" and x.get("organism") == "Homo sapiens"]
    allt = [x["target_chembl_id"] for x in t]
    n_act = sum(get("activity.json", target_chembl_id=s, pchembl_value__gte=PCHEMBL_MIN, limit=1)["page_meta"]["total_count"] for s in sp)
    phase = None
    if allt:
        mech = get("mechanism.json", target_chembl_id__in=",".join(allt), limit=1000)["mechanisms"]
        mols = sorted({m["molecule_chembl_id"] for m in mech if m.get("molecule_chembl_id")})
        if mols:
            ph = [float(m["max_phase"]) for m in get("molecule.json", molecule_chembl_id__in=",".join(mols[:200]), limit=1000)["molecules"]
                  if m.get("max_phase") is not None]
            phase = max(ph) if ph else None
    st = ("drug_target" if phase is not None and phase >= 1 else "potent_ligands" if n_act > 0
          else "chembl_target" if allt else "not_in_chembl")
    return dict(acc=acc, single_protein_targets=";".join(sp), n_targets_any=len(allt), n_act_pchembl6=n_act,
                mech_max_phase=phase, chembl_status=st)


def main():
    accs = [a.strip() for a in open(INTERIM / "q23a_candidates.txt") if a.strip()]
    rel = get("status.json").get("chembl_db_version", "NA")
    with ThreadPoolExecutor(6) as ex:
        rows = list(ex.map(one, accs))
    out = pd.DataFrame(rows); out["chembl_release"] = rel
    out.to_csv(INTERIM / "chembl_target_status.tsv", sep="\t", index=False)
    print(rel, len(out)); print(out.chembl_status.value_counts())


if __name__ == "__main__":
    main()
