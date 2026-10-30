"""Ligand-protein atom contact rules (pure numpy/scipy; unit-tested on toy coordinates).

Primary rule (BioLiP): heavy-atom pair in contact if d <= vdW_i + vdW_j + VDW_TOLERANCE.
Sensitivity rule: d <= SIMPLE_CUTOFF.
"""
from functools import lru_cache

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .config import SIMPLE_CUTOFF, VDW_TOLERANCE


@lru_cache(maxsize=None)
def vdw_radius(element: str) -> float:
    """van der Waals radius (A) from gemmi's element table; unknown elements -> carbon radius."""
    el = gemmi.Element(element)
    r = el.vdw_r
    return float(r) if r and r > 0 else float(gemmi.Element("C").vdw_r)


def radii(elements) -> np.ndarray:
    return np.array([vdw_radius(e) for e in elements], dtype=float)


def atom_pairs(lig_xyz, lig_r, tree: cKDTree, poly_r, max_poly_r: float):
    """All ligand-atom x polymer-atom pairs that satisfy the vdW or the simple-cutoff rule.

    Returns arrays (lig_idx, poly_idx, dist, is_vdw, is_simple)."""
    if len(lig_xyz) == 0 or tree.n == 0:
        e = np.zeros(0, dtype=int)
        return e, e, np.zeros(0), np.zeros(0, bool), np.zeros(0, bool)
    rmax = max(SIMPLE_CUTOFF, float(lig_r.max()) + max_poly_r + VDW_TOLERANCE)
    hits = tree.query_ball_point(lig_xyz, r=rmax)
    li = np.fromiter((i for i, h in enumerate(hits) for _ in h), dtype=int)
    pj = np.fromiter((j for h in hits for j in h), dtype=int)
    if len(li) == 0:
        e = np.zeros(0, dtype=int)
        return e, e, np.zeros(0), np.zeros(0, bool), np.zeros(0, bool)
    d = np.linalg.norm(lig_xyz[li] - tree.data[pj], axis=1)
    is_vdw = d <= lig_r[li] + poly_r[pj] + VDW_TOLERANCE
    is_simple = d <= SIMPLE_CUTOFF
    keep = is_vdw | is_simple
    return li[keep], pj[keep], d[keep], is_vdw[keep], is_simple[keep]


def residue_contacts(pj, d, is_vdw, is_simple, poly_res_idx):
    """Aggregate atom pairs to residues: {res_idx: (n_vdw_pairs, n_simple_pairs, min_dist)}."""
    out = {}
    for j, dist, v, s in zip(pj, d, is_vdw, is_simple):
        r = int(poly_res_idx[j])
        nv, ns, md = out.get(r, (0, 0, 1e9))
        out[r] = (nv + int(v), ns + int(s), min(md, float(dist)))
    return out
