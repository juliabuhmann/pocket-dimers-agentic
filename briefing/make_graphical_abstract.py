#!/usr/bin/env python
"""Graphical abstract and README hero image for the pocket-dimers study.

First attempt; rendered and looked at by the producing agent only. Whether the
figure states the findings correctly has NOT been confirmed by the user.

Every number drawn is read from the same result / review tables that
briefing/make_figures.py uses (paths named at each read, see load_numbers()).
Phrases that are only in prose (caveats, cofactor names) are checked with
make_figures.require_text, so the build fails if the source wording changes.
Style constants and cartoon helpers (ligand, chain, notch) are imported from
make_figures.py so the visual language matches figures 1-8.

Run from anywhere:
    .venvs/briefing/bin/python briefing/make_graphical_abstract.py
Outputs: briefing/figures/graphical_abstract.png (16:9) and
         briefing/figures/readme_hero.png (~2.6:1), both PNG at dpi 200.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_figures as mf  # noqa: E402  (sets rcParams, palette, helpers)
from make_figures import (  # noqa: E402
    C_COF, C_DIMER, C_DRUG, C_SINGLE, CHAIN_A, CHAIN_B, GRID, INK, INK2, LIG, MUTED,
    SURFACE, chain, ligand, notch, require_text, tsv,
)

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Rectangle  # noqa: E402

OUT = mf.OUT
DPI = 200
REPO_NAME = mf.WS.name           # repository directory name, e.g. pocket-dimers-agentic
BRIEF = "briefing/briefing.md"
CARD_FC = "#f7f7f4"              # card background (neutral, lighter than GRID)
FONT = "DejaVu Sans"


# ---------------------------------------------------------------------------
# Numbers (all read from files; nothing typed in)
# ---------------------------------------------------------------------------
def load_numbers() -> dict:
    n = {}
    # Finding 1a: MBS fraction of symmetry-collapsed sites, after additive filter
    # (REPORT/briefing: "sites 26.7% -> 26.3%"), producer table, verified in verify_v3.
    d = tsv("results/part1/analysis/mbs_fraction_descriptives_seed0.tsv")
    s = d[(d.level == "site_symcollapsed") & (d["filter"] == "after_additive_filter")
          & (d.grouping == "all_oligomers")].iloc[0]
    n["site_frac"] = float(s.mbs_frac)
    n["one_in"] = round(1 / n["site_frac"])            # "about 1 in N"
    # Finding 1b: fpocket druggability (fig6 source), weighted, all sites
    fp = tsv("results/part1/fpocket/pdbsite_mbs_vs_sbs_cliffs_delta.tsv")
    f = fp[(fp.metric == "fp_drug") & (fp.stratum == "all_weighted")].iloc[0]
    n["fp_delta"], n["fp_lo"], n["fp_hi"] = f.cliffs_delta, f.cd_lo, f.cd_hi
    n["fp_mbs"], n["fp_sbs"] = f.frac_gt_thr_mbs, f.frac_gt_thr_sbs
    assert n["fp_lo"] < 0 < n["fp_hi"], "fpocket delta CI excludes 0; 'not more druggable' wording invalid"
    # Finding 1c: pre-registered burrow_score (H2-secondary), fig6 source
    b = tsv("results/part1/pocket_descriptors/mbs_vs_sbs_cliffs_delta.tsv")
    b = b[(b.metric == "burrow_score") & (b.stratum == "all_weighted") & (b.seed == 0)].iloc[0]
    n["bs_delta"], n["bs_lo"], n["bs_hi"] = b.cliffs_delta, b.cd_lo, b.cd_hi
    assert n["bs_lo"] < 0 < n["bs_hi"], "burrow_score CI excludes 0"

    # Finding 2: on/off-axis cluster counts (verifier, as fig4b)
    st = pd.read_csv(mf.WS / "reviews/verify_part1/out/h3_stress.tsv", sep="\t")
    prim = st[st.variant == "primary"].iloc[0]
    n["on_dl"], n["on_cof"] = int(prim.cl_on_dl), int(prim.cl_on_cof)
    n["off_dl"], n["off_cof"] = int(prim.cl_off_dl), int(prim.cl_off_cof)

    # Finding 3a: H1 strict-MBS sensitivity OR (fig4a source)
    ht = tsv("results/part1/analysis/hypothesis_tests.tsv")
    r = ht[(ht.test == "H1") & (ht.variant == "strict_mbs")].iloc[0]
    n["h1s_or"], n["h1s_lo"], n["h1s_hi"] = r.estimate, r.ci_lo, r.ci_hi
    assert n["h1s_hi"] < 1, "strict-MBS OR CI no longer below 1"
    # Finding 3b: off-axis cofactors placed without partner (verifier seed 0, as fig7)
    rs = pd.read_csv(mf.WS / "reviews/verify_cofold/rescored_seed0.tsv", sep="\t")
    rs["success"] = rs["success"].astype(str) == "True"
    off = rs[rs.case_id.str.startswith("offaxis_") & (rs.cond == "monomer")]
    n["off_mono_k"], n["off_mono_n"] = int(off.success.sum()), len(off)

    # Finding 4a: H4 hit rates, homodimer MBS, organic, filtered (fig5 source)
    hr = tsv("results/part2/h4/v2/hit_rates.tsv")
    hr = hr[(hr["filter"] == "filtered") & (hr.lig_set == "organic")
            & (hr.kind == "homo") & (hr.site_type == "MBS")].set_index("condition")
    n["h4_dimer"], n["h4_iso"] = float(hr.loc["dimer", "hit_rate"]), float(hr.loc["iso", "hit_rate"])
    # pre-registered DiD (verifier CI, as fig5b)
    sv = pd.read_csv(mf.WS / "reviews/verify_h4/stats_results.tsv", sep="\t")
    p = sv[sv.analysis.str.startswith("primary, unstratified cluster bootstrap (producer")].iloc[0]
    n["did"], n["did_lo"], n["did_hi"] = p.did, p.did_lo, p.did_hi
    n["h4_sig"] = not (p.did_lo < 0 < p.did_hi)
    # Finding 4b: AF3 on-axis drugs pooled over seeds (fig7 source)
    seeds = tsv("results/cofold/onaxis/seeds_summary.tsv")
    pooled = seeds[(seeds.table == "pooled") & (seeds.metric == "sym") & (seeds.group == "on_axis_druglike")]
    dd, mm = pooled[pooled.cond == "dimer"].iloc[0], pooled[pooled.cond == "monomer"].iloc[0]
    n["af_d_k"], n["af_d_n"] = int(dd.n_success), int(dd.n)
    n["af_m_k"], n["af_m_n"] = int(mm.n_success), int(mm.n)
    n["n_seeds"] = int(seeds[(seeds.table == "per_seed") & (seeds.group == "on_axis_druglike")].seed.nunique())

    # Prose anchors (wording used in the figure must exist in the briefing)
    require_text(BRIEF, "NAD, FAD, SAM, PLP and CoA")
    require_text(BRIEF, "touched at the rim by the partner")
    require_text(BRIEF, "pre-selected as confident predictions")
    require_text(BRIEF, "keep this result descriptive")
    require_text(BRIEF, "small case series")
    return n


def pct(x):
    return f"{x * 100:.0f}%"


def sdelta(x):
    return f"{x:+.2f}".replace("-", "−")


# ---------------------------------------------------------------------------
# Drawing helpers (data units of a single full-figure axes, aspect equal)
# ---------------------------------------------------------------------------
def canvas(w_in, h_in, scale=10):
    fig = plt.figure(figsize=(w_in, h_in))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w_in * scale)
    ax.set_ylim(0, h_in * scale)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def T(ax, x, y, s, size=10, color=INK, weight="normal", ha="left", va="baseline", family=FONT, **kw):
    return ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va,
                   family=family, zorder=10, **kw)


def card(ax, x, y, w, h):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                fc=CARD_FC, ec="none", zorder=0))


def badge(ax, x, y, k, r=1.3, size=10):
    ax.add_patch(Circle((x, y), r, fc=INK, ec="none", zorder=5))
    T(ax, x, y - 0.05, str(k), size=size, color="white", weight="bold", ha="center", va="center")


def track_bar(ax, x, y, w, h, frac, color):
    """Horizontal proportion bar: grey track + coloured fill (left-anchored)."""
    ax.add_patch(Rectangle((x, y), w, h, fc=GRID, ec="none", zorder=3))
    ax.add_patch(Rectangle((x, y), w * frac, h, fc=color, ec="none", zorder=4))


def split_bar(ax, x, y, w, h, a, b, ca, cb, size=9, labels=True):
    """100% stacked bar of two counts with a 2px-ish surface gap between segments."""
    fa = a / (a + b)
    gap = 0.25
    ax.add_patch(Rectangle((x, y), w * fa - gap / 2, h, fc=ca, ec="none", zorder=4))
    ax.add_patch(Rectangle((x + w * fa + gap / 2, y), w * (1 - fa) - gap / 2, h, fc=cb, ec="none", zorder=4))
    if labels:
        T(ax, x - 0.6, y + h / 2, f"{a}", size=size, color=INK, weight="bold", ha="right", va="center")
        T(ax, x + w + 0.6, y + h / 2, f"{b}", size=size, color=INK, weight="bold", ha="left", va="center")


def sbs_mbs_pair(ax, cx, cy, s=1.0, labels=True, lab_size=9, chain_size=9):
    """fig2a-style pair: SBS (ligand in chain A only) and MBS (ligand in the A/B cleft)."""
    for dx, lab in [(-3.0 * s, "single-chain\nsite"), (3.0 * s, "shared\n(multi-chain) site")]:
        x0 = cx + dx
        chain(ax, x0, cy + 0.93 * s, 3.0 * s, 1.9 * s, CHAIN_A)
        chain(ax, x0, cy - 0.93 * s, 3.0 * s, 1.9 * s, CHAIN_B)
        T(ax, x0 - 0.6 * s, cy + 0.98 * s, "A", size=chain_size, color=INK2, ha="center", va="center")
        T(ax, x0 - 0.6 * s, cy - 1.05 * s, "B", size=chain_size, color=INK2, ha="center", va="center")
        if labels:
            T(ax, x0, cy - 2.55 * s, lab, size=lab_size, ha="center", va="top", linespacing=1.15)
    notch(ax, cx - 3.0 * s + 1.15 * s, cy + 1.38 * s, 0.5 * s)
    ligand(ax, cx - 3.0 * s + 1.2 * s, cy + 1.38 * s, r=0.28 * s)
    notch(ax, cx + 3.0 * s + 1.0 * s, cy + 0.0 * s, 0.55 * s)
    ligand(ax, cx + 3.0 * s + 1.0 * s, cy + 0.0 * s, r=0.3 * s)


def axis_dimer(ax, cx, cy, s, on_axis: bool, lig_color):
    """fig2b-style homodimer seen down its twofold axis; on-axis = one ligand on the axis."""
    chain(ax, cx - 1.0 * s, cy, 2.0 * s, 3.8 * s, CHAIN_A)
    chain(ax, cx + 1.0 * s, cy, 2.0 * s, 3.8 * s, CHAIN_A)
    if on_axis:
        notch(ax, cx, cy, 0.55 * s)
        ligand(ax, cx, cy, r=0.32 * s, color=lig_color)
    else:
        for dx, dy in [(-0.3, 1.25), (0.3, -1.25)]:
            notch(ax, cx + dx * s, cy + dy * s, 0.48 * s)
            ligand(ax, cx + dx * s, cy + dy * s, r=0.29 * s, color=lig_color)
    ax.add_patch(Ellipse((cx, cy), 0.42 * s, 0.25 * s, fc=C_DIMER, ec="none", zorder=7))


def cofactor_one_chain(ax, cx, cy, s):
    """Cofactor held in a pocket of chain A; partner B only touches the rim."""
    chain(ax, cx, cy + 0.9 * s, 3.4 * s, 2.2 * s, CHAIN_A)          # bottom edge at cy - 0.2 s
    chain(ax, cx + 0.3 * s, cy - 1.05 * s, 3.4 * s, 1.9 * s, CHAIN_B)  # top edge at cy - 0.1 s
    notch(ax, cx + 0.4 * s, cy + 0.47 * s, 0.5 * s)                  # pocket inside A, rim reaches B
    ligand(ax, cx + 0.4 * s, cy + 0.47 * s, r=0.3 * s, color=C_COF)


def save(fig, name):
    path = OUT / name
    fig.savefig(path, dpi=DPI, facecolor=SURFACE)  # fixed canvas, no tight bbox
    plt.close(fig)
    print("wrote", path)


# ---------------------------------------------------------------------------
# Graphical abstract (16 x 9 in)
# ---------------------------------------------------------------------------
def graphical_abstract(n):
    fig, ax = canvas(16, 9)           # 160 x 90 units, 1 unit = 0.1 in
    # title
    T(ax, 4, 83, "Ligand pockets formed by one protein chain vs by two chains",
      size=24, weight="bold", va="center")
    T(ax, 4, 78.2, "PDB structures and AlphaFold dimer models", size=15, color=INK2, va="center")

    # left concept panel
    card(ax, 4, 12, 36, 60)
    T(ax, 6.5, 68, "Two kinds of ligand site", size=13, weight="bold", va="center")
    sbs_mbs_pair(ax, 22.3, 53, s=3.2, lab_size=11.5, chain_size=13)
    # legend
    ly = 32
    items = [("hex", LIG, "ligand"), ("ell", CHAIN_A, "protein chain A"), ("ell", CHAIN_B, "partner chain B"),
             ("hex", C_DRUG, "drug-like ligand"), ("hex", C_COF, "cofactor"),
             ("axis", C_DIMER, "homodimer two-fold axis")]
    for i, (kind, col, lab) in enumerate(items):
        yy = ly - i * 3.0
        if kind == "hex":
            ligand(ax, 9, yy, r=0.75, color=col)
        elif kind == "ell":
            ax.add_patch(Ellipse((9, yy), 2.2, 1.5, fc=col, ec=MUTED, lw=0.9, zorder=3))
        else:
            ax.add_patch(Ellipse((9, yy), 1.1, 0.65, fc=col, ec="none", zorder=3))
        T(ax, 11.5, yy, lab, size=11, color=INK2, va="center")

    # 2x2 finding cards
    X0, Y0, W, H, G = 43, 12, 55.5, 28.5, 3.0
    cells = [(X0, Y0 + H + G), (X0 + W + G, Y0 + H + G), (X0, Y0), (X0 + W + G, Y0)]

    def head(i, x, y, text):
        badge(ax, x + 3.2, y + H - 3.4, i + 1, r=1.5, size=12)
        T(ax, x + 5.8, y + H - 3.4, text, size=13, weight="bold", va="center")

    # --- 1: shared sites common, not more druggable
    x, y = cells[0]
    card(ax, x, y, W, H)
    head(0, x, y, "Shared sites: common, not more druggable")
    T(ax, x + 3, y + 15.5, pct(n["site_frac"]), size=34, weight="bold", va="baseline")
    T(ax, x + 3, y + 11.3, f"of PDB ligand sites in complexes\nare shared between chains\n(about 1 in {n['one_in']})",
      size=11, color=INK2, va="top", linespacing=1.25)
    bx, bw = x + 30, 21
    T(ax, bx, y + 19.8, "fpocket druggability > 0.5", size=10.5, color=INK, va="center")
    for j, (lab, v, col) in enumerate([("shared", n["fp_mbs"], C_DIMER), ("single-chain", n["fp_sbs"], C_SINGLE)]):
        yy = y + 15.2 - j * 4.2
        T(ax, bx, yy + 1.9, lab, size=10, color=INK2, va="center")
        track_bar(ax, bx, yy - 0.8, bw - 5, 1.3, v, col)
        T(ax, bx + bw - 4.3, yy - 0.15, pct(v), size=11, weight="bold", va="center")
    T(ax, bx, y + 3.4, f"fpocket Cliff's δ {sdelta(n['fp_delta'])}\npocket score δ {sdelta(n['bs_delta'])} (no difference)",
      size=9.5, color=INK2, va="center", linespacing=1.3)

    # --- 2: on-axis drugs
    x, y = cells[1]
    card(ax, x, y, W, H)
    head(1, x, y, "Drugs dominate the homodimer two-fold axis")
    for j, (on, a, b, lab) in enumerate([(True, n["on_dl"], n["on_cof"], "on-axis"),
                                         (False, n["off_dl"], n["off_cof"], "off-axis")]):
        cx = x + 14 + j * 27
        axis_dimer(ax, cx, y + 15.6, 2.3, on, C_DRUG if on else C_COF)
        T(ax, cx, y + 8.4, lab, size=11.5, weight="bold", ha="center", va="center")
        split_bar(ax, cx - 8, y + 4.6, 16, 1.7, a, b, C_DRUG, C_COF, size=11)
    T(ax, x + W / 2, y + 1.9, "protein families with drug-like vs cofactor ligands (descriptive, between proteins)",
      size=9.5, color=INK2, ha="center", va="center")

    # --- 3: cofactors anchored in one chain
    x, y = cells[2]
    card(ax, x, y, W, H)
    head(2, x, y, "Classic cofactors sit in one chain")
    cofactor_one_chain(ax, x + 10.5, y + 12.5, 2.5)
    T(ax, x + 11.2, y + 3.4, "partner touches\nonly the rim", size=10, color=INK2, ha="center", va="center",
      linespacing=1.2)
    tx = x + 22
    T(ax, tx, y + 20.2, "NAD, FAD, SAM, PLP, CoA", size=12, weight="bold", va="center")
    T(ax, tx, y + 16.6, f"strict multi-chain definition: cofactor\nOR {n['h1s_or']:.2f} "
      f"(95% CI {n['h1s_lo']:.2f}–{n['h1s_hi']:.2f}; sensitivity)", size=10.5, color=INK2, va="top",
      linespacing=1.25)
    T(ax, tx, y + 5.4, f"{n['off_mono_k']}/{n['off_mono_n']}", size=26, weight="bold", va="baseline")
    T(ax, tx + 10.5, y + 7.4, "off-axis cofactors placed\ncorrectly by AlphaFold3\nwithout the partner",
      size=10.5, color=INK2, va="center", linespacing=1.2)

    # --- 4: dimer models
    x, y = cells[3]
    card(ax, x, y, W, H)
    head(3, x, y, "Dimer models find shared sites more often")
    bx, bw = x + 4, 30
    rows = [("Pocket finding: real homodimer shared sites found",
             [("dimer model", n["h4_dimer"], pct(n["h4_dimer"]), C_DIMER),
              ("single chain", n["h4_iso"], pct(n["h4_iso"]), C_SINGLE)]),
            (f"AlphaFold3: on-axis drug poses correct ({n['n_seeds']} seeds)",
             [("with partner", n["af_d_k"] / n["af_d_n"], f"{n['af_d_k']}/{n['af_d_n']}", C_DIMER),
              ("single chain", n["af_m_k"] / n["af_m_n"], f"{n['af_m_k']}/{n['af_m_n']}", C_SINGLE)])]
    for k, (title, bars) in enumerate(rows):
        y0 = y + 21.4 - k * 8.4
        T(ax, bx, y0, title, size=10.5, color=INK, va="center")
        for j, (lab, v, vt, col) in enumerate(bars):
            yy = y0 - 3.0 - j * 2.7
            T(ax, bx, yy, lab, size=9.5, color=INK2, va="center")
            track_bar(ax, bx + 11, yy - 0.75, bw, 1.5, v, col)
            T(ax, bx + 11 + bw + 1.0, yy, vt, size=11, weight="bold", va="center")
    sig = "significant" if n["h4_sig"] else "not significant"
    T(ax, bx, y + 4.6,
      f"Pre-registered recovery test {sig} (DiD {sdelta(n['did'])},\n95% CI {sdelta(n['did_lo'])} to "
      f"{sdelta(n['did_hi'])}). Co-folding: small case series.", size=9.5, color=INK2, va="top", linespacing=1.3)

    # footer
    T(ax, 4, 5.5, "Descriptive / exploratory where marked; small co-folding sets; confident AlphaFold interfaces only.",
      size=10.5, color=INK2, va="center")
    T(ax, 156, 5.5, REPO_NAME, size=10.5, color=INK2, ha="right", va="center", family="DejaVu Sans Mono")
    save(fig, "graphical_abstract.png")


# ---------------------------------------------------------------------------
# README hero (10.4 x 4 in, ~2.6:1)
# ---------------------------------------------------------------------------
def readme_hero(n):
    fig, ax = canvas(10.4, 4.0)       # 104 x 40 units
    T(ax, 3, 35.6, "Ligand pockets shared\nbetween protein chains", size=16.5, weight="bold", va="top",
      linespacing=1.1)
    T(ax, 3, 28.4, "PDB structures and\nAlphaFold dimer models", size=11, color=INK2, va="top", linespacing=1.2)
    sbs_mbs_pair(ax, 16.0, 14.0, s=1.75, lab_size=9.5, chain_size=10)

    tiles = [
        ("mbs", pct(n["site_frac"]), "of PDB ligand sites\nare shared between\nchains, not more\ndruggable"),
        ("axis", f"{n['on_dl']} : {n['on_cof']}", "drug-like vs cofactor\nfamilies on the\nhomodimer two-fold\naxis"),
        ("cof", f"{n['off_mono_k']}/{n['off_mono_n']}", "off-axis cofactors\nco-folded correctly\nfrom one chain"),
        ("h4", f"{pct(n['h4_dimer'])} vs {pct(n['h4_iso'])}", "shared sites found in\ndimer models vs\nsingle chains"),
    ]
    X0, W, G, Y0, H = 35.0, 15.6, 1.3, 4.5, 32.5
    for i, (icon, big, line) in enumerate(tiles):
        x = X0 + i * (W + G)
        card(ax, x, Y0, W, H)
        cx, cy = x + W / 2, Y0 + 23.5
        if icon == "mbs":
            chain(ax, cx, cy + 1.3, 5.0, 2.8, CHAIN_A)
            chain(ax, cx, cy - 1.3, 5.0, 2.8, CHAIN_B)
            notch(ax, cx + 1.6, cy, 0.85)
            ligand(ax, cx + 1.6, cy, r=0.5)
        elif icon == "axis":
            axis_dimer(ax, cx, cy, 1.35, True, C_DRUG)
        elif icon == "cof":
            cofactor_one_chain(ax, cx - 0.6, cy + 0.4, 1.4)
        else:
            for j, (v, col) in enumerate([(n["h4_dimer"], C_DIMER), (n["h4_iso"], C_SINGLE)]):
                ax.add_patch(Rectangle((cx - 4.0 + j * 4.8, cy - 3.2), 3.2, 9.0 * v / 0.5, fc=col, ec="none",
                                       zorder=4))
            ax.plot([cx - 5.0, cx + 5.0], [cy - 3.2, cy - 3.2], color=MUTED, lw=0.9, zorder=5)
            for j, lab in enumerate(["dimer", "single"]):
                T(ax, cx - 2.4 + j * 4.8, cy - 4.0, lab, size=8.5, color=INK2, ha="center", va="top")
        T(ax, cx, Y0 + 15.0, big, size=21 if len(big) < 7 else 15.5, weight="bold", ha="center", va="center")
        T(ax, cx, Y0 + 11.3, line, size=10, color=INK2, ha="center", va="top", linespacing=1.2)
    sig = "not significant" if not n["h4_sig"] else "significant"
    T(ax, 101.5, 1.9, f"Descriptive / exploratory; pre-registered recovery test {sig}; small co-folding sets.",
      size=8.5, color=INK2, ha="right", va="center")
    save(fig, "readme_hero.png")


if __name__ == "__main__":
    nums = load_numbers()
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in nums.items()})
    graphical_abstract(nums)
    readme_hero(nums)
