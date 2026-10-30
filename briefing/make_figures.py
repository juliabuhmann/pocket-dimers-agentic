#!/usr/bin/env python
"""Briefing figures for the pocket-dimers study (fig1-fig8).

First attempt; the figures have been rendered and looked at once, but their
correctness against the report has NOT been confirmed by the user.

Every number drawn comes from a results table, docs/*.md or reviews/verify_*.md
in this workspace (paths are named next to each read). Where REPORT.md quotes a
verifier's number that differs from the producer table, the verifier's file is
read instead and the switch is noted in a comment.

Run from anywhere:
    .venvs/briefing/bin/python briefing_src/make_figures.py
Outputs: briefing_src/figures/fig{1..8}_*.png (PNG, dpi 200).
"""
from __future__ import annotations

import re
import subprocess
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import (  # noqa: E402
    Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle, RegularPolygon,
)

WS = Path(__file__).resolve().parents[1]  # repository root
OUT = WS / "briefing" / "figures"
# Files that live under project_record/ or report/ in the published repository
# (the script was first run in the study workspace, where they sat at the root).
_RELOCATED = {"REPORT.md": "report/REPORT.md", "RESEARCH_TREE.md": "project_record/RESEARCH_TREE.md",
              "LESSONS.md": "project_record/LESSONS.md"}


def rp(rel: str) -> Path:
    return WS / _RELOCATED.get(rel, rel)
OUT.mkdir(parents=True, exist_ok=True)
DPI = 200  # print resolution requested for the briefing

# ---------------------------------------------------------------------------
# Style. Palette = dataviz-skill reference palette (light mode); every set used
# together was run through scripts/validate_palette.js against #ffffff.
# Semantic roles are kept fixed across all figures.
# ---------------------------------------------------------------------------
INK = "#0b0b0b"
INK2 = "#52514e"      # secondary text, PDB/experimental series
MUTED = "#8a8984"
GRID = "#e6e5e0"
SURFACE = "#ffffff"
C_DIMER = "#2a78d6"   # with partner / dimer model / multi-chain
C_SINGLE = "#eb6834"  # isolated chain / single chain
C_RANDOM = "#a3a29c"  # random-patch null, controls
C_DRUG = "#eda100"    # drug-like ligands
C_COF = "#008300"     # cofactors
C_HOMO = "#4a3aa7"    # homomer / homodimer
C_HETERO = "#1baf7a"  # heteromer / heterodimer
C_PDB = INK2          # experimental PDB sites
C_PRED = C_DIMER      # predicted AF dimers

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 9.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 8.5,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "axes.linewidth": 0.8,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def tsv(rel: str) -> pd.DataFrame:
    """Read a workspace TSV, skipping '#' provenance header lines."""
    return pd.read_csv(rp(rel), sep="\t", comment="#")


def panel_label(ax, s, x=-0.02, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="bottom", ha="right", color=INK)


def ygrid(ax, axis="y"):
    ax.grid(axis=axis, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)


def save(fig, name):
    fig.savefig(OUT / name, dpi=DPI, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print("wrote", OUT / name)


def require_text(rel: str, phrase: str) -> None:
    """Fail loudly if an annotation's source phrase is not in the cited file."""
    txt = rp(rel).read_text(encoding="utf-8")
    if phrase not in txt:
        raise AssertionError(f"phrase not found in {rel}: {phrase!r}")


# ---------------------------------------------------------------------------
# Fig 1 - agentic workflow schematic (counts derived from files / git)
# ---------------------------------------------------------------------------
def project_counts() -> dict:
    # The published repository has a fresh git history and no START.txt, so the
    # counts computed in the study workspace are read from a frozen JSON file.
    frozen = WS / "briefing" / "figure_inputs" / "fig1_counts.json"
    if not (WS / "START.txt").exists() and frozen.exists():
        import json
        return json.loads(frozen.read_text())
    rv = WS / "reviews"
    c = {
        "r1": len(list(rv.glob("r1_*.md"))),
        "r2": len(list(rv.glob("r2_*.md"))),
        "audit": len(list(rv.glob("audit_*.md"))),
        "verify": len(list(rv.glob("verify_*.md"))),
        "report_review": len(list(rv.glob("report_review_*.md"))),
        "retro": len(list(rv.glob("retro_*.md"))),
    }
    c["design_reviews"] = c["r1"] + c["r2"]
    # Commit count and wall clock are frozen at the last commit that touched
    # RESEARCH_TREE.md (the project's final milestone), so that later briefing
    # commits do not change the figure.
    git = ["git", "-C", str(WS)]
    last = subprocess.check_output(git + ["log", "-1", "--format=%H", "--", "RESEARCH_TREE.md"],
                                   text=True).strip()
    c["commits"] = int(subprocess.check_output(git + ["rev-list", "--count", last], text=True))
    t_end = datetime.fromisoformat(subprocess.check_output(
        git + ["log", "-1", "--format=%aI", last], text=True).strip())
    start_txt = (WS / "START.txt").read_text().split()[0]  # first token = ISO start time
    t0 = datetime.fromisoformat(start_txt)
    c["hours"] = (t_end - t0).total_seconds() / 3600
    # Number of pre-registered hypotheses and human dimer models: parsed from RESEARCH_TREE.md
    tree = (WS / "RESEARCH_TREE.md").read_text()
    m = re.search(r"It has (\w+) pre-registered hypotheses", tree)
    c["n_hyp"] = {"five": 5}.get(m.group(1), m.group(1)) if m else "?"
    m = re.search(r"Human dimer pocket runs are done \(([\d,]+) models", tree)
    c["n_models"] = m.group(1) if m else "?"
    return c


def fig1():
    c = project_counts()
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 64)
    ax.axis("off")

    def box(x, y, w, h, title, sub="", fc="#f4f3ef", ec=MUTED, tc=INK, bold=True):
        p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.6",
                           fc=fc, ec=ec, lw=1.0)
        ax.add_patch(p)
        ax.text(x + w / 2, y + h / 2 + (2.0 if sub else 0), title, ha="center", va="center",
                fontsize=7.6, fontweight="bold" if bold else "normal", color=tc, wrap=True)
        if sub:
            ax.text(x + w / 2, y + h / 2 - 2.6, sub, ha="center", va="center", fontsize=6.5,
                    color=INK2, linespacing=1.15)
        return (x, y, w, h)

    def arrow(p, q, color=INK2, style="-|>", rad=0.0, ls="-"):
        a = FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=10, color=color, lw=1.1,
                            connectionstyle=f"arc3,rad={rad}", linestyle=ls)
        ax.add_patch(a)

    AGENT = "#e8f0fb"   # light tint of C_DIMER for AI-agent boxes
    CHECK = "#fdeee6"   # light tint of C_SINGLE for independent-check boxes
    HUMAN = "#f4f3ef"

    # Row 1 (top, left -> right)
    y1, h = 46, 12
    b1 = box(0.5, y1, 16, h, "Human prompt", "research question\n+ deadline", fc=HUMAN)
    b2 = box(20.5, y1, 16, h, 'AI "PI" agent', "plans, delegates,\nwrites report", fc=AGENT, ec=C_DIMER)
    b3 = box(40.5, y1, 17, h, "Literature &\nrecon", "4 parallel agents", fc=AGENT, ec=C_DIMER)
    b4 = box(61.5, y1, 17, h, "Design reviews",
             f"2 rounds, {c['design_reviews']} reviewer\nreports ({c['r1']} + {c['r2']})", fc=CHECK, ec=C_SINGLE)
    b5 = box(82.5, y1, 17, h, "Frozen design",
             f"{c['n_hyp']} pre-registered\nhypotheses, git tag", fc=HUMAN, ec=INK2)
    for a, b in [(b1, b2), (b2, b3), (b3, b4), (b4, b5)]:
        arrow((a[0] + a[2] + 0.6, y1 + h / 2), (b[0] - 0.6, y1 + h / 2))

    # Row 2: parallel producers
    prods = ["PDB census", f"Pocket runs\n({c['n_models']} dimers)", "Analyses\nH1-H5, Q2.x", "AF3 co-folding"]
    px = [2, 26, 50, 74]
    y2, h2 = 24, 11
    ax.text(50, y2 + h2 + 3.6, "parallel producer agents", ha="center", fontsize=7.6,
            color=C_DIMER, style="italic")
    for x, t in zip(px, prods):
        box(x, y2, 22, h2, t, fc=AGENT, ec=C_DIMER)
    ax.plot([2, 96], [y2 + h2 + 1.6] * 2, color=INK2, lw=1.0)
    arrow((91, y1 - 0.8), (91, y2 + h2 + 1.8))

    # Row 3 (bottom, right -> left)
    y3, h3 = 2, 12
    b6 = box(62, y3, 37.5, h3, "Audit + independent verifiers",
             f"{c['audit']} audit + {c['verify']} verification reports;\nnumbers re-derived without producer code",
             fc=CHECK, ec=C_SINGLE)
    b7 = box(36.5, y3, 21, h3, "Report reviews", f"{c['report_review']} reviewers\n(PI, statistics)",
             fc=CHECK, ec=C_SINGLE)
    b8 = box(9, y3, 23, h3, "Report", "claims tagged VERIFIED /\nPRE-REG / EXPLORATORY", fc=HUMAN, ec=INK2)
    arrow((81, y2 - 0.6), (81, y3 + h3 + 0.8))
    arrow((b6[0] - 0.6, y3 + h3 / 2), (b7[0] + b7[2] + 0.6, y3 + h3 / 2))
    arrow((b7[0] - 0.6, y3 + h3 / 2), (b8[0] + b8[2] + 0.6, y3 + h3 / 2))
    # feedback loop: verifiers -> producers (fixes / reruns)
    arrow((70, y3 + h3 + 0.8), (60, y2 - 0.8), color=C_SINGLE, rad=0.25, ls="--")
    ax.text(58.5, 19.0, "fixes &\nreruns", fontsize=7, color=C_SINGLE, ha="right", va="center")

    ax.text(1, 0.0 - 1.5,
            f"{c['commits']} git commits  ·  ~{c['hours']:.1f} h wall clock (start to final milestone)  ·  "
            f"{c['retro']} process retrospective; LESSONS.md",
            fontsize=7.6, color=INK2, va="top")
    lg = [Line2D([], [], marker="s", ls="", ms=8, mfc=AGENT, mec=C_DIMER, label="AI agent producing work"),
          Line2D([], [], marker="s", ls="", ms=8, mfc=CHECK, mec=C_SINGLE, label="independent AI check"),
          Line2D([], [], marker="s", ls="", ms=8, mfc=HUMAN, mec=INK2, label="input / frozen artefact")]
    ax.legend(handles=lg, loc="upper left", bbox_to_anchor=(0.0, 1.08), ncol=3, fontsize=7.5)
    save(fig, "fig1_agentic_workflow.png")
    return c


# ---------------------------------------------------------------------------
# Fig 2 - concept cartoon (no data)
# ---------------------------------------------------------------------------
CHAIN_A = "#cfe0f6"   # light tint of C_DIMER
CHAIN_B = "#fbd9c8"   # light tint of C_SINGLE
LIG = INK


def ligand(ax, x, y, r=0.32, color=LIG):
    ax.add_patch(RegularPolygon((x, y), 6, radius=r, fc=color, ec="white", lw=0.8, zorder=6))


def chain(ax, x, y, w, h, fc, label=None, angle=0, lx=None, ly=None):
    ax.add_patch(Ellipse((x, y), w, h, angle=angle, fc=fc, ec=MUTED, lw=0.9, zorder=2))
    if label:
        ax.text(lx if lx is not None else x, ly if ly is not None else y, label, ha="center",
                va="center", fontsize=8, color=INK2, zorder=7)


def notch(ax, x, y, r):
    """A pocket = a surface-coloured bite out of the chain outline."""
    ax.add_patch(Circle((x, y), r, fc=SURFACE, ec=MUTED, lw=0.9, zorder=3))


def fig2():
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), gridspec_kw={"wspace": 0.08})
    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 7.6)
        ax.set_aspect("equal")
        ax.axis("off")

    def title(ax, letter, text):
        ax.text(0.0, 7.55, letter, fontsize=11, fontweight="bold", color=INK, va="top")
        ax.text(0.8, 7.45, text, fontsize=8.3, fontweight="bold", color=INK, va="top")

    # (a) SBS vs MBS: the same dimer drawn twice, ligand position differs
    ax = axes[0]
    title(ax, "a", "Single- vs multi-chain site")
    for cx, lab in [(2.4, "single-chain\nsite (SBS)"), (7.4, "multi-chain\nsite (MBS)")]:
        chain(ax, cx, 5.1, 3.0, 1.9, CHAIN_A, "A", lx=cx - 0.6, ly=5.2)
        chain(ax, cx, 3.25, 3.0, 1.9, CHAIN_B, "B", lx=cx - 0.6, ly=3.1)
        ax.text(cx, 1.15, lab, ha="center", va="center", fontsize=7.2, color=INK)
    notch(ax, 3.55, 5.55, 0.5)          # SBS: pocket in chain A only, away from B
    ligand(ax, 3.6, 5.55, r=0.28)
    notch(ax, 8.4, 4.18, 0.55)          # MBS: pocket in the A/B cleft
    ligand(ax, 8.4, 4.18, r=0.3)
    ax.plot([4.9, 4.9], [0.6, 6.4], color=GRID, lw=1)

    # (b) twofold axis perpendicular to the page (lens symbol at the centre)
    ax = axes[1]
    title(ax, "b", "Homodimer twofold axis")
    for cx, lab in [(2.5, "on-axis:\none ligand"), (7.5, "off-axis:\ntwo copies")]:
        chain(ax, cx - 1.0, 4.1, 2.0, 3.8, CHAIN_A, "A", lx=cx - 1.35, ly=4.1)
        chain(ax, cx + 1.0, 4.1, 2.0, 3.8, CHAIN_A, "A'", lx=cx + 1.35, ly=4.1)
        ax.text(cx, 1.15, lab, ha="center", va="center", fontsize=7.2, color=INK)
    notch(ax, 2.5, 4.1, 0.55)
    ligand(ax, 2.5, 4.1, r=0.3)
    ax.add_patch(Ellipse((2.5, 4.1), 0.34, 0.2, fc=C_DIMER, ec="none", zorder=7))
    for dx, dy in [(-0.3, 1.25), (0.3, -1.25)]:
        notch(ax, 7.5 + dx, 4.1 + dy, 0.48)
        ligand(ax, 7.5 + dx, 4.1 + dy, r=0.27)
    ax.add_patch(Ellipse((7.5, 4.1), 0.5, 0.3, fc=C_DIMER, ec="none", zorder=7))
    ax.add_patch(Ellipse((4.0, 6.55), 0.5, 0.3, fc=C_DIMER, ec="none"))
    ax.text(4.4, 6.55, "twofold axis (⊥ page)", fontsize=6.8, color=INK2, va="center")

    # (c) pocket that only exists with the partner
    ax = axes[2]
    title(ax, "c", "Pocket needs the partner")
    chain(ax, 1.55, 4.2, 2.0, 3.6, CHAIN_A, "A", lx=1.15, ly=4.2)
    chain(ax, 3.45, 4.2, 2.0, 3.6, CHAIN_B, "B", lx=3.85, ly=4.2)
    notch(ax, 2.5, 4.2, 0.5)
    ax.add_patch(Circle((2.5, 4.2), 0.5, fc="none", ec=C_DIMER, lw=1.6, ls="--", zorder=6))
    ax.text(2.5, 1.15, "dimer:\ninterface pocket", ha="center", va="center", fontsize=7.4, color=INK)
    ax.add_patch(FancyArrowPatch((4.75, 4.2), (6.15, 4.2), arrowstyle="-|>", mutation_scale=10,
                                 color=INK2, lw=1.1))
    ax.text(5.45, 4.5, "remove\npartner", ha="center", va="bottom", fontsize=6.6, color=INK2)
    chain(ax, 7.6, 4.2, 2.0, 3.6, CHAIN_A, "A", lx=7.2, ly=4.2)
    ax.add_patch(Circle((8.6, 4.2), 0.5, fc="none", ec=C_SINGLE, lw=1.4, ls=":", zorder=6))
    ax.text(7.6, 1.15, "isolated chain:\nno pocket", ha="center", va="center",
            fontsize=7.4, color=INK)

    lg = [Line2D([], [], marker="h", ls="", ms=8, mfc=LIG, mec="white", label="ligand"),
          Line2D([], [], marker="o", ls="", ms=9, mfc=CHAIN_A, mec=MUTED, label="protein chain (A, A')"),
          Line2D([], [], marker="o", ls="", ms=9, mfc=CHAIN_B, mec=MUTED, label="partner chain B"),
          Line2D([], [], marker="o", ls="", ms=9, mfc="none", mec=C_DIMER, mew=1.4, label="predicted pocket")]
    fig.legend(handles=lg, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.02))
    save(fig, "fig2_concepts.png")


# ---------------------------------------------------------------------------
# Fig 3 - PDB census
# ---------------------------------------------------------------------------
CLASS_NAMES = {
    "lipid_fatty_acid": "Lipid / fatty acid", "heme_porphyrin": "Heme / porphyrin",
    "sugar_glycan": "Sugar / glycan", "cofactor": "Cofactor", "drug_like": "Drug-like",
    "nucleotide_like": "Nucleotide-like", "other_organic": "Other organic",
    "amino_acid": "Amino acid", "inorganic": "Inorganic", "metal_cluster": "Metal cluster",
    "metal_ion": "Metal ion",
}


def fig3():
    # Producer table; site/unit fractions were re-derived in reviews/verify_v3.md
    # (v0_part1_fractions.tsv) and match to <=0.3 pp, so producer values are used.
    d = tsv("results/part1/analysis/mbs_fraction_descriptives_seed0.tsv")
    a = d[(d.level == "unit_nonredundant") & (d["filter"] == "after_additive_filter")
          & (d.grouping == "assembly_type+cls")]
    homo = a[a.assembly_type == "homomer"].set_index("cls")
    het = a[a.assembly_type == "heteromer"].set_index("cls")
    order = homo.sort_values("mbs_frac").index.tolist()
    order = [c for c in order if c in CLASS_NAMES]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.9),
                                   gridspec_kw={"width_ratios": [1.35, 1], "wspace": 0.55})
    y = np.arange(len(order))
    for df, col, off, lab, mk in [(homo, C_HOMO, 0.17, "Homomer", "o"), (het, C_HETERO, -0.17, "Heteromer", "s")]:
        v = df.loc[order]
        ax1.errorbar(v.mbs_frac * 100, y + off,
                     xerr=[(v.mbs_frac - v.ci_lo) * 100, (v.ci_hi - v.mbs_frac) * 100],
                     fmt=mk, ms=5, color=col, mec="white", mew=0.6, elinewidth=1.2, capsize=0, label=lab,
                     zorder=3)
    allu = d[(d.level == "unit_nonredundant") & (d["filter"] == "after_additive_filter")
             & (d.grouping == "all_oligomers")].iloc[0]
    ax1.axvline(allu.mbs_frac * 100, color=MUTED, lw=0.9, ls="--", zorder=1)
    ax1.text(allu.mbs_frac * 100 + 1, len(order) - 0.35, f"all units {allu.mbs_frac*100:.1f}%",
             fontsize=7, color=INK2, va="center")
    ax1.set_yticks(y)
    ax1.set_yticklabels([CLASS_NAMES[c] for c in order])
    ax1.set_xlabel("Multi-chain sites (% of non-redundant units)")
    ax1.set_xlim(0, 85)
    ax1.set_ylim(-0.7, len(order) - 0.1)
    ygrid(ax1, "x")
    ax1.legend(loc="lower right")
    ax1.set_title("MBS fraction by ligand class", loc="left")
    panel_label(ax1, "a", x=-0.42)

    # (b) before vs after additive filter
    rows = [("site_symcollapsed", "all_oligomers", "", "Sites, all"),
            ("site_symcollapsed", "assembly_type", "homomer", "Sites, homomer"),
            ("site_symcollapsed", "assembly_type", "heteromer", "Sites, heteromer"),
            ("unit_nonredundant", "all_oligomers", "", "Units, all"),
            ("unit_nonredundant", "assembly_type", "homomer", "Units, homomer"),
            ("unit_nonredundant", "assembly_type", "heteromer", "Units, heteromer")]
    yy = np.arange(len(rows))[::-1]
    for yi, (lvl, grp, asm, lab) in zip(yy, rows):
        vals = {}
        for f in ["before_additive_filter", "after_additive_filter"]:
            s = d[(d.level == lvl) & (d["filter"] == f) & (d.grouping == grp)]
            if asm:
                s = s[s.assembly_type == asm]
            vals[f] = s.iloc[0]
        b, af = vals["before_additive_filter"], vals["after_additive_filter"]
        ax2.plot([b.mbs_frac * 100, af.mbs_frac * 100], [yi, yi], color=GRID, lw=3, zorder=1,
                 solid_capstyle="round")
        ax2.errorbar(b.mbs_frac * 100, yi + 0.12, xerr=[[(b.mbs_frac - b.ci_lo) * 100], [(b.ci_hi - b.mbs_frac) * 100]],
                     fmt="o", mfc="white", mec=INK2, color=INK2, ms=5, elinewidth=0.9, zorder=3)
        ax2.errorbar(af.mbs_frac * 100, yi - 0.12, xerr=[[(af.mbs_frac - af.ci_lo) * 100], [(af.ci_hi - af.mbs_frac) * 100]],
                     fmt="o", color=INK, ms=5, mec="white", elinewidth=0.9, zorder=3)
    ax2.set_yticks(yy)
    ax2.set_yticklabels([r[3] for r in rows])
    ax2.set_xlabel("Multi-chain (%)")
    ax2.set_xlim(15, 45)
    ax2.set_ylim(-0.5, len(rows) - 0.2)
    ygrid(ax2, "x")
    ax2.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=INK2, label="before additive filter"),
                        Line2D([], [], marker="o", ls="", color=INK, label="after (primary set)")],
               loc="upper center", bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=7.0, handletextpad=0.1, columnspacing=0.8)
    ax2.set_title("Effect of removing additives", loc="left")
    panel_label(ax2, "b", x=-0.5)
    save(fig, "fig3_pdb_census.png")


# ---------------------------------------------------------------------------
# Fig 4 - pre-registered H1-H3 forest + on/off-axis composition
# ---------------------------------------------------------------------------
def fig4():
    ht = tsv("results/part1/analysis/hypothesis_tests.tsv")

    def row(test, variant):
        r = ht[(ht.test == test) & (ht.variant == variant)].iloc[0]
        return r.estimate, r.ci_lo, r.ci_hi

    # H3 cluster level and the on/off-axis cluster counts are verifier numbers
    # (REPORT quotes them as VERIFIED; the producer table has unit-level only).
    st = pd.read_csv(WS / "reviews/verify_part1/out/h3_stress.tsv", sep="\t")
    pres = st[st.variant.str.startswith("presence-only")].iloc[0]
    prim = st[st.variant == "primary"].iloc[0]

    items = [
        ("H1  cofactor vs drug-like, MBS odds", *row("H1", "primary"), "primary"),
        ("H1 sensitivity: strict MBS", *row("H1", "strict_mbs"), "sens"),
        ("H2  drug-like occupancy, MBS vs SBS", *row("H2", "primary"), "primary"),
        ("H3  drug-like on- vs off-axis, units", *row("H3", "primary"), "primary"),
        ("H3 at cluster level (verifier)", pres.fisher_or, pres.robust_lo, pres.robust_hi, "sens"),
    ]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.1),
                                   gridspec_kw={"width_ratios": [1.6, 1], "wspace": 0.75})
    y = np.arange(len(items))[::-1]
    for yi, (lab, est, lo, hi, kind) in zip(y, items):
        filled = kind == "primary"
        ax1.plot([lo, hi], [yi, yi], color=INK, lw=1.3, zorder=2)
        ax1.plot(est, yi, "D" if filled else "o", ms=6, mfc=INK if filled else "white", mec=INK, zorder=3)
        ax1.text(1.02, yi, f"{est:.2f} ({lo:.2f}–{hi:.2f})" if est < 10 else f"{est:.1f} ({lo:.1f}–{hi:.1f})",
                 va="center", fontsize=7.4, color=INK2, transform=ax1.get_yaxis_transform())
    ax1.axvline(1, color=MUTED, lw=0.9, ls="--", zorder=1)
    ax1.set_xscale("log")
    ax1.set_xlim(0.15, 130)
    ax1.set_xticks([0.2, 0.5, 1, 2, 5, 10, 20, 50, 100])
    ax1.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10", "20", "50", "100"])
    ax1.set_yticks(y)
    ax1.set_yticklabels([i[0] for i in items])
    ax1.set_xlabel("Odds ratio (95% CI, log scale)")
    ax1.spines["left"].set_visible(False)
    ax1.tick_params(axis="y", length=0)
    ygrid(ax1, "x")
    ax1.set_title("Pre-registered tests, family A (PDB)", loc="left")
    ax1.legend(handles=[Line2D([], [], marker="D", ls="", color=INK, label="pre-registered primary"),
                        Line2D([], [], marker="o", ls="", mfc="white", mec=INK, label="sensitivity / descriptive")],
               loc="upper center", bbox_to_anchor=(0.35, -0.2), ncol=2)
    panel_label(ax1, "a", x=-0.9)

    # (b) composition, cluster counts
    cats = [("On-axis", prim.cl_on_dl, prim.cl_on_cof), ("Off-axis", prim.cl_off_dl, prim.cl_off_cof)]
    for i, (lab, dl, cf) in enumerate(cats):
        tot = dl + cf
        ax2.bar(i, dl / tot * 100, width=0.6, color=C_DRUG, edgecolor=SURFACE, lw=2)
        ax2.bar(i, cf / tot * 100, bottom=dl / tot * 100, width=0.6, color=C_COF, edgecolor=SURFACE, lw=2)
        ax2.text(i, dl / tot * 50, f"{int(dl)}", ha="center", va="center", fontsize=8, color=INK, fontweight="bold")
        ax2.text(i, dl / tot * 100 + cf / tot * 50, f"{int(cf)}", ha="center", va="center", fontsize=8,
                 color="white", fontweight="bold")
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels([c[0] for c in cats])
    ax2.set_ylabel("Homomer MBS clusters (%)")
    ax2.set_ylim(0, 100)
    ax2.set_xlim(-0.6, 1.6)
    ygrid(ax2)
    ax2.legend(handles=[Rectangle((0, 0), 1, 1, color=C_DRUG, label="drug-like"),
                        Rectangle((0, 0), 1, 1, color=C_COF, label="cofactor")],
               loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    ax2.set_title("Twofold-axis composition", loc="left")
    panel_label(ax2, "b", x=-0.3)
    save(fig, "fig4_ligand_classes.png")


# ---------------------------------------------------------------------------
# Fig 5 - H4 recovery
# ---------------------------------------------------------------------------
def fig5():
    hr = tsv("results/part2/h4/v2/hit_rates.tsv")
    hr = hr[(hr["filter"] == "filtered") & (hr.lig_set == "organic")]
    # random patches: dimer-side surface-patch null, seed 0 (primary seed, DESIGN §5);
    # REPORT §4.1 quotes these seed-0 values (0.005 / 0.000 / ...).
    conds = [("dimer", "Dimer model", C_DIMER), ("iso", "Isolated chains", C_SINGLE),
             ("dimer_null0", "Random patches", C_RANDOM)]
    groups = [("homo", "MBS"), ("homo", "SBS"), ("hetero", "MBS"), ("hetero", "SBS")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2),
                                   gridspec_kw={"width_ratios": [1.6, 1], "wspace": 0.45})
    w = 0.26
    for gi, (k, s) in enumerate(groups):
        for ci, (cond, lab, col) in enumerate(conds):
            r = hr[(hr.kind == k) & (hr.site_type == s) & (hr.condition == cond)].iloc[0]
            x = gi + (ci - 1) * w
            ax1.bar(x, r.hit_rate, width=w - 0.03, color=col, label=lab if gi == 0 else None, zorder=2)
            ax1.errorbar(x, r.hit_rate, yerr=[[r.hit_rate - r.ci_lo], [r.ci_hi - r.hit_rate]],
                         color=INK, lw=0.9, capsize=0, zorder=3)
    ns = [int(hr[(hr.kind == k) & (hr.site_type == s) & (hr.condition == "dimer")].n.iloc[0]) for k, s in groups]
    ax1.set_xticks(range(4))
    ax1.set_xticklabels([f"{'Homo' if k == 'homo' else 'Hetero'} {s}\nn={n:,}" for (k, s), n in zip(groups, ns)],
                        fontsize=7.4)
    ax1.set_ylabel("Sites recovered (DCC ≤4 Å, top n+2)")
    ax1.set_ylim(0, 0.65)
    ax1.axvline(1.5, color=GRID, lw=1)
    ygrid(ax1)
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
    ax1.set_title("Recovery of PDB sites", loc="left")
    panel_label(ax1, "a", x=-0.12)

    # (b) DiD: verifier values (reviews/verify_h4/stats_results.tsv). REPORT quotes
    # the verifier's unstratified CI (-0.006, 0.338) rather than the producer's
    # (-0.015, 0.332; results/part2/h4/v2/h4_tests.tsv).
    sv = pd.read_csv(WS / "reviews/verify_h4/stats_results.tsv", sep="\t")
    p = sv[sv.analysis.str.startswith("primary, unstratified cluster bootstrap (producer")].iloc[0]
    f = sv[sv.analysis.str.startswith("family-weighted (each cluster")].iloc[0]
    rows = [("Site-weighted\n(pre-registered)", p, "D", True), ("Family-weighted\n(exploratory)", f, "o", False)]
    for yi, (lab, r, mk, filled) in zip([1, 0], rows):
        ax2.plot([r.did_lo, r.did_hi], [yi, yi], color=INK, lw=1.4)
        ax2.plot(r.did, yi, mk, ms=7, mfc=INK if filled else "white", mec=INK, zorder=3)
        ax2.text(r.did_hi, yi + 0.22, f"{r.did:+.2f} ({r.did_lo:.2f} to {r.did_hi:.2f})", ha="right",
                 fontsize=7.2, color=INK2, bbox=dict(fc=SURFACE, ec="none", pad=0.6), zorder=4)
    ax2.axvline(0, color=C_SINGLE, lw=1.1, ls="--")
    ax2.set_yticks([1, 0])
    ax2.set_yticklabels([r[0] for r in rows])
    ax2.set_ylim(-0.6, 1.6)
    ax2.set_xlim(-0.08, 0.4)
    ax2.set_xlabel("Difference-in-differences\n(MBS gain − SBS gain)")
    ax2.spines["left"].set_visible(False)
    ax2.tick_params(axis="y", length=0)
    ygrid(ax2, "x")
    ax2.set_title("Dimer advantage at MBS", loc="left")
    panel_label(ax2, "b", x=-0.55)
    save(fig, "fig5_recovery.png")


# ---------------------------------------------------------------------------
# Fig 6 - pocket chemistry: predicted dimers vs PDB sites
# ---------------------------------------------------------------------------
def fig6():
    pat = tsv("results/part1/pocket_descriptors/pdb_vs_afdimer_pattern.tsv")
    pat = pat[pat.pdb_contrast == "PDB MBS vs SBS ligand sites"].set_index("metric")
    q = tsv("results/part2/q22/q22_props.tsv")
    q = q[(q.stratum == "all") & (q.seed == 0)].set_index("metric")
    desc = tsv("results/part1/pocket_descriptors/mbs_vs_sbs_cliffs_delta.tsv")
    bs = desc[(desc.metric == "burrow_score") & (desc.stratum == "all_weighted") & (desc.seed == 0)].iloc[0]
    fp = tsv("results/part1/fpocket/pdbsite_mbs_vs_sbs_cliffs_delta.tsv")
    fpr = fp[(fp.metric == "fp_drug") & (fp.stratum == "all_weighted")].iloc[0]

    metrics = [
        # label, AF (delta, lo, hi), PDB (delta, lo, hi)
        ("Pocket volume", (q.loc["volume_A3", "cliffs_delta"], q.loc["volume_A3", "cd_lo"], q.loc["volume_A3", "cd_hi"]),
         tuple(pat.loc["volume_A3", ["pdb_delta", "pdb_lo", "pdb_hi"]])),
        ("Charged fraction", (q.loc["frac_charged", "cliffs_delta"], q.loc["frac_charged", "cd_lo"], q.loc["frac_charged", "cd_hi"]),
         tuple(pat.loc["frac_charged", ["pdb_delta", "pdb_lo", "pdb_hi"]])),
        ("Hydrophobic fraction", (q.loc["frac_hydrophobic", "cliffs_delta"], q.loc["frac_hydrophobic", "cd_lo"], q.loc["frac_hydrophobic", "cd_hi"]),
         tuple(pat.loc["frac_hydrophobic", ["pdb_delta", "pdb_lo", "pdb_hi"]])),
        # PDB burrow_score = pre-registered H2-secondary (located pocket, uncropped), as in REPORT §3.3
        ("burrow_score", (q.loc["burrow_score", "cliffs_delta"], q.loc["burrow_score", "cd_lo"], q.loc["burrow_score", "cd_hi"]),
         (bs.cliffs_delta, bs.cd_lo, bs.cd_hi)),
        ("fpocket druggability", None, (fpr.cliffs_delta, fpr.cd_lo, fpr.cd_hi)),
    ]
    fig, ax = plt.subplots(figsize=(6.4, 3.3))
    y = np.arange(len(metrics))[::-1]
    off = 0.17
    for yi, (lab, af, pdb) in zip(y, metrics):
        if af is not None:
            ax.plot([af[1], af[2]], [yi + off] * 2, color=C_PRED, lw=1.6)
            ax.plot(af[0], yi + off, "o", ms=6, color=C_PRED, mec="white", mew=0.6, zorder=3)
        else:
            ax.text(0.03, yi + off, "(no matched AF-dimer contrast)", fontsize=6.8, color=MUTED, va="center")
        ax.plot([pdb[1], pdb[2]], [yi - off] * 2, color=C_PDB, lw=1.6)
        ax.plot(pdb[0], yi - off, "s", ms=6, color=C_PDB, mec="white", mew=0.6, zorder=3)
    ax.axvline(0, color=MUTED, lw=0.9, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([m[0] for m in metrics])
    ax.set_xlim(-0.3, 0.3)
    ax.set_xlabel("Cliff's δ (95% CI)   ←  lower in multi-chain / interface   |   higher  →")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ygrid(ax, "x")
    ax.legend(handles=[Line2D([], [], marker="o", color=C_PRED, label="Predicted AF dimers: interface vs single-chain pocket"),
                       Line2D([], [], marker="s", color=C_PDB, label="PDB ligand sites: MBS vs SBS (uncropped)")],
              loc="upper center", bbox_to_anchor=(0.4, -0.2), ncol=1)
    save(fig, "fig6_pocket_chemistry.png")


# ---------------------------------------------------------------------------
# Fig 7 - co-folding
# ---------------------------------------------------------------------------
def fig7():
    seeds = tsv("results/cofold/onaxis/seeds_summary.tsv")
    pooled = seeds[(seeds.table == "pooled") & (seeds.metric == "sym") & (seeds.group == "on_axis_druglike")]
    on_d = pooled[pooled.cond == "dimer"].iloc[0]
    on_m = pooled[pooled.cond == "monomer"].iloc[0]
    # Seed-0 success calls from the verifier's symmetry-aware RMSD (REPORT §5 uses
    # these; producer metric gives 5/8 for dimer-only single chain, verifier 6/8).
    rs = pd.read_csv(WS / "reviews/verify_cofold/rescored_seed0.tsv", sep="\t")
    rs["success"] = rs["success"].astype(str) == "True"

    def cnt(mask, cond):
        s = rs[mask & (rs.cond == cond)]
        return int(s.success.sum()), len(s)

    off = rs.case_id.str.startswith("offaxis_")
    dmo = rs.batch == "dimeronly"
    het = tsv("results/cofold/hetero_ext/summary.tsv")
    h = het[het.subset == "all"].iloc[0]
    groups = [
        ("On-axis\ndrug-like\n(12 × 3 seeds)", (int(on_d.n_success), int(on_d.n)), (int(on_m.n_success), int(on_m.n))),
        ("Off-axis\ncofactor\n(control)", cnt(off, "dimer"), cnt(off, "monomer")),
        ("Dimer-only\nactive sites", cnt(dmo, "dimer"), cnt(dmo, "monomer")),
        ("Hetero-\ndimer\nMBS", (int(h.dimer_success), int(h.n)), (int(h.single_success), int(h.n))),
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.5),
                                   gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.35})
    w = 0.38
    for gi, (lab, dd, mm) in enumerate(groups):
        for j, ((k, n), col, nm) in enumerate([(dd, C_DIMER, "with partner (dimer)"), (mm, C_SINGLE, "single chain")]):
            x = gi + (j - 0.5) * w
            ax1.bar(x, k / n * 100, width=w - 0.04, color=col, label=nm if gi == 0 else None, zorder=2)
            ax1.text(x, k / n * 100 + 1.5, f"{k}/{n}", ha="center", va="bottom", fontsize=7, color=INK)
    ax1.set_xticks(range(len(groups)))
    ax1.set_xticklabels([g[0] for g in groups], fontsize=7.4)
    ax1.set_ylabel("Correct poses (ligand RMSD ≤2 Å, %)")
    ax1.set_ylim(0, 110)
    ax1.set_yticks([0, 25, 50, 75, 100])
    ygrid(ax1)
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=1)
    ax1.set_title("AF3 pose success", loc="left")
    panel_label(ax1, "a", x=-0.12)

    # (b) paired RMSD, seed 0: verifier RMSD for on-axis/off-axis/SBS/dimer-only/H4 homo;
    # heterodimers (26 pairs) from hetero_ext/paired.tsv rmsd_ver_* (verifier code).
    piv = rs.pivot_table(index=["batch", "case_id"], columns="cond", values="lig_rmsd").reset_index()

    def grp(cid):
        if cid.startswith("on_"):
            return "on"
        if cid.startswith("offaxis_"):
            return "off"
        if cid.startswith("sbs_"):
            return "sbs"
        if cid.startswith("hetero_"):
            return None  # taken from paired.tsv
        return "homo"

    piv["g"] = piv.case_id.map(grp)
    piv = piv.dropna(subset=["g"])
    hp = tsv("results/cofold/hetero_ext/paired.tsv")
    hp = pd.DataFrame({"case_id": hp.case_id, "dimer": hp.rmsd_ver_dimer, "monomer": hp.rmsd_ver_monomer, "g": "het"})
    pts = pd.concat([piv[["case_id", "dimer", "monomer", "g"]], hp], ignore_index=True).dropna()
    style = {"on": ("On-axis drug-like", C_DRUG, "D"), "off": ("Off-axis cofactor", C_COF, "^"),
             "homo": ("Homodimer MBS (dimer-only, H4)", C_HOMO, "o"), "het": ("Heterodimer MBS", C_HETERO, "s"),
             "sbs": ("SBS drug-like (control)", C_RANDOM, "v")}
    lo, hi = 0.1, 60
    ax2.fill_between([lo, 2], [2, 2], [hi, hi], color="#eef4fc", zorder=0, lw=0)
    ax2.plot([lo, hi], [lo, hi], color=MUTED, lw=0.9, zorder=1)
    ax2.axvline(2, color=MUTED, lw=0.8, ls=":")
    ax2.axhline(2, color=MUTED, lw=0.8, ls=":")
    for g, (lab, col, mk) in style.items():
        s = pts[pts.g == g]
        ax2.scatter(s.monomer.clip(lo, hi), s.dimer.clip(lo, hi), s=26, c=col, marker=mk,
                    edgecolors="white", linewidths=0.6, label=f"{lab} (n={len(s)})", zorder=3)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(lo, hi)
    ax2.set_ylim(lo, hi)
    ax2.set_xticks([0.2, 0.5, 1, 2, 5, 10, 20])
    ax2.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10", "20"])
    ax2.set_yticks([0.2, 0.5, 1, 2, 5, 10, 20])
    ax2.set_yticklabels(["0.2", "0.5", "1", "2", "5", "10", "20"])
    ax2.set_xlabel("Ligand RMSD, single chain (Å)")
    ax2.set_ylabel("Ligand RMSD, dimer (Å)")
    ax2.text(0.13, 30, "only single\nchain correct", fontsize=6.6, color=INK2, va="top")
    ax2.text(50, 0.13, "only dimer\ncorrect", fontsize=6.6, color=INK2, ha="right")
    ax2.legend(loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=1, fontsize=7, handletextpad=0.2)
    ax2.set_title("Paired ligand RMSD (seed 0)", loc="left")
    panel_label(ax2, "b", x=-0.2)
    save(fig, "fig7_cofolding.png")
    return groups, pts.groupby("g").size().to_dict()


# ---------------------------------------------------------------------------
# Fig 8 - what independent checks caught (text table; each claim tied to a file)
# ---------------------------------------------------------------------------
CATCHES = [
    # (stage, source file, phrase that must occur in it, caught, changed)
    ("Audit, before tests", "reviews/audit_part1_census.md", "44,230 (75%) from `ro5qed_FALLBACK`",
     "75% of 'drug-like' labels came from a\nRo5+QED fallback that included buffers",
     "ChEMBL-based drug-like label;\ncensus v2 before any hypothesis test"),
    ("Audit, before tests", "reviews/audit_part1_census.md", "fails for about 43% of symmetric copies",
     "~43% of symmetric ligand copies in\nhomodimers were not collapsed",
     "Symmetry collapse fixed in census v2"),
    ("Verifier, H3", "reviews/verify_part1_results.md", "Cluster-level ORs are about 14",
     "Unit-level OR 35 driven by a few drugged\ntargets (HIV protease = 127 units)",
     "Cluster-level OR ≈14 reported;\nH3 labelled descriptive"),
    ("Verifier, H4", "reviews/verify_h4.md", "The permutation null is centred at **0.091, not 0**",
     "Permutation p = 1e-4 tested a null\ncentred at 0.091, not 0",
     "Permutation p dropped; cluster bootstrap:\nH4 'not established' (CI includes 0)"),
    ("Verifier, H5", "reviews/verify_h5.md", "1.09 [0.97, 1.23]",
     "OR 1.20 falls to 1.09 (0.97–1.23) after\nadjusting for amino-acid identity",
     "Reported as pocket-lining effect,\nnot a dimer-specific disease signal"),
    ("Verifier, descriptors", "reviews/verify_misc.md", "cropped fraction is higher for MBS (24.4%) than SBS (17.8%)",
     "30 Å crop removed chains at 24% of MBS\nvs 18% of SBS sites (differential bias)",
     "Uncropped rerun; volume δ +0.107 → +0.071"),
    ("Verifier, co-folding", "reviews/verify_cofold.md", "dimer-only single chain 6/8 instead of 5/8",
     "Symmetry-aware RMSD changes some calls\n(109/112 agree)",
     "Verifier counts used in report\n(dimer-only 6/8 vs 6/8)"),
    ("Report review", "reviews/report_review_pi.md", "H1 is therefore a new contrast that came out null",
     "H1 framed as a failed Marsh replication",
     "Reframed: cofactor > metal reproduced;\ncofactor vs drug-like is a new, null test"),
]
# additional anchors for numbers quoted in the CATCHES text
EXTRA_ANCHORS = [
    ("reviews/verify_part1_results.md", "gives HIV PR 127 units"),
    ("reviews/verify_cofold.md", "109/112 at 2 Å"),
    ("reviews/verify_v3.md", "Cropped merged volume +0.107"),
    ("REPORT.md", "+0.107 cropped → +0.071 uncropped"),
]


def fig8():
    for rel, phrase, *_ in [(c[1], c[2]) for c in CATCHES] + EXTRA_ANCHORS:
        require_text(rel, phrase)
    n = len(CATCHES)
    fig, ax = plt.subplots(figsize=(7.2, 0.66 * n + 0.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.1, n + 0.6)
    ax.axis("off")
    X_CAUGHT, X_CHANGED = 21, 61
    for x, t in [(1, "Check"), (X_CAUGHT, "What it caught"), (X_CHANGED, "What changed")]:
        ax.text(x, n + 0.25, t, fontsize=8.4, fontweight="bold", color=INK, va="center")
    ax.plot([0, 100], [n - 0.02, n - 0.02], color=INK2, lw=0.8)
    for i, (stage, src, _, caught, changed) in enumerate(CATCHES):
        yc = n - 1 - i + 0.5
        if i % 2 == 0:
            ax.add_patch(Rectangle((0, yc - 0.5), 100, 1.0, fc="#f6f5f1", ec="none", zorder=0))
        chip_col = C_HOMO if "review" in stage.lower() else C_SINGLE
        ax.add_patch(FancyBboxPatch((1, yc - 0.05), 17.5, 0.4, boxstyle="round,pad=0.04,rounding_size=0.12",
                                    fc="white", ec=chip_col, lw=1.0))
        ax.text(9.75, yc + 0.15, stage.replace("\n", " "), ha="center", va="center", fontsize=6.6, color=INK)
        ax.text(9.75, yc - 0.27, Path(src).name, ha="center", va="center", fontsize=5.6,
                color=INK2, family="DejaVu Sans Mono")
        ax.text(X_CAUGHT, yc, caught, va="center", fontsize=6.9, color=INK, linespacing=1.15)
        ax.annotate("", xy=(X_CHANGED - 0.8, yc), xytext=(X_CHANGED - 3.2, yc),
                    arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.9, mutation_scale=8))
        ax.text(X_CHANGED, yc, changed, va="center", fontsize=6.9, color=INK, linespacing=1.15)
    save(fig, "fig8_verification.png")


if __name__ == "__main__":
    import sys
    # optional: pass figure numbers to build a subset, e.g. `make_figures.py 1 2`
    wanted = set(sys.argv[1:]) or {str(i) for i in range(1, 9)}
    funcs = {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5, "6": fig6, "7": fig7, "8": fig8}
    for k in sorted(wanted):
        res = funcs[k]()
        if res is not None:
            print(f"fig{k} derived values:", res)
