#!/usr/bin/env bash
# Render a briefing markdown file to PDF with pandoc + tectonic (XeTeX engine).
# Usage:  briefing_src/build_pdf.sh [input.md] [output.pdf]
#   default input  = briefing_src/briefing.md
#   default output = same name with .pdf
# Toolchain: conda env conda_envs/pdf (created with
#   mamba create -p conda_envs/pdf -c conda-forge pandoc tectonic poppler).
# tectonic downloads LaTeX packages/fonts on first use (needs internet once);
# set TECTONIC_CACHE_DIR to relocate its cache.
# Fonts are named by file (TeX Gyre Pagella from the tectonic bundle) so the
# result does not depend on system fontconfig.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WS="$(dirname "$HERE")"
ENV="${PDF_ENV:-$WS/conda_envs/pdf}"
IN="$(realpath "${1:-$HERE/briefing.md}")"
OUT="${2:-${IN%.md}.pdf}"
OUT="$(realpath -m "$OUT")"

export PATH="$ENV/bin:$PATH"
command -v pandoc >/dev/null || { echo "pandoc not found in $ENV/bin" >&2; exit 1; }
command -v tectonic >/dev/null || { echo "tectonic not found in $ENV/bin" >&2; exit 1; }

# LaTeX header: every figure at full text width (capped at 70% of text height),
# captions in small type, figures kept near their reference.
HDR="$(mktemp --suffix=.tex)"
trap 'rm -f "$HDR"' EXIT
cat > "$HDR" <<'TEX'
\usepackage{float}
\floatplacement{figure}{htbp}
\makeatletter
\setkeys{Gin}{width=\linewidth,height=0.7\textheight,keepaspectratio}
\makeatother
\usepackage[font=small,labelfont=bf]{caption}
\usepackage{booktabs}
TEX

cd "$(dirname "$IN")"   # relative figure paths (figures/...) resolve from the md's folder
pandoc "$IN" \
  --from markdown+implicit_figures+pipe_tables+yaml_metadata_block \
  --pdf-engine=tectonic \
  --number-sections \
  --include-in-header="$HDR" \
  --resource-path=".:$HERE" \
  -V papersize=a4 \
  -V fontsize=11pt \
  -V geometry:margin=2.5cm \
  -V mainfont="texgyrepagella" \
  -V mainfontoptions="Extension=.otf" \
  -V mainfontoptions="UprightFont=*-regular" \
  -V mainfontoptions="BoldFont=*-bold" \
  -V mainfontoptions="ItalicFont=*-italic" \
  -V mainfontoptions="BoldItalicFont=*-bolditalic" \
  -V mathfont="texgyrepagella-math.otf" \
  -V linestretch=1.1 \
  -V colorlinks=true \
  -V linkcolor=black -V urlcolor=blue -V citecolor=black \
  -o "$OUT"
echo "wrote $OUT"
