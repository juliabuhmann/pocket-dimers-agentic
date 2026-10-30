#!/usr/bin/env bash
# Records the checked-out commit of each tool under tools/ (documented pins: docs/TOOLS.md, docs/COFOLDING.md).
cd "$(dirname "$0")/.."
for t in burrow pocketcompare adock alphafold3; do printf '%s\t%s\n' $t "$(git -C tools/$t rev-parse HEAD)"; done > envs/TOOL_COMMITS.txt
