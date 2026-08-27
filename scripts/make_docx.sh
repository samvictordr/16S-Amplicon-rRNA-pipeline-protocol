#!/bin/bash
# ==============================================================================
# scripts/make_docx.sh — Word version of the Frontiers manuscript
# ==============================================================================
# Produces paper/frontiers/manuscript.docx for people who need to read or check
# the text outside LaTeX: supervisors, similarity checking, tracked-change
# comments. It is NOT a submission artifact. Frontiers takes the LaTeX source,
# and this conversion drops formatting the class provides.
#
# Direct `pandoc manuscript.tex` does not give a usable result. Three problems:
#
#   1. The custom \bp unit is not an siunitx unit pandoc knows, and equations
#      carrying \label fail to parse as math. Both emit raw TeX into the output.
#   2. \ref{} has nothing to resolve against, so every cross-reference becomes
#      a placeholder. A document saying "see Section ??" 40 times is not
#      reviewable. LaTeX already resolved them; manuscript.aux holds the answers.
#   3. The Frontiers class defines the title and author block through macros
#      (\Authors, \corrAuthor) that pandoc does not evaluate, so the byline
#      vanishes.
#
# So the .tex is preprocessed first, then handed to pandoc.
#
# Usage:
#   bash scripts/make_docx.sh            # build the PDF first if stale, then convert
#   bash scripts/make_docx.sh --no-latex # skip the LaTeX build, reuse existing .aux
#
# Requires: pandoc. The .aux file must exist, so the manuscript must have been
# built at least once.
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DIR="$REPO_ROOT/paper/frontiers"
RUN_LATEX=1

[[ "${1:-}" == "--no-latex" ]] && RUN_LATEX=0

command -v pandoc >/dev/null || {
    echo "pandoc not found. Install it: brew install pandoc" >&2
    exit 1
}

cd "$DIR"

if [[ $RUN_LATEX -eq 1 ]]; then
    echo "[1/3] Building the PDF, so cross-references resolve..."
    command -v latexmk >/dev/null || {
        echo "  latexmk not found; pass --no-latex to reuse the existing .aux" >&2
        exit 1
    }
    latexmk -pdf -interaction=nonstopmode manuscript.tex >/dev/null 2>&1 || true
fi

[[ -f manuscript.aux ]] || {
    echo "manuscript.aux not found. Build the PDF once before converting." >&2
    exit 1
}

echo "[2/3] Preprocessing LaTeX for pandoc..."
python3 "$SCRIPT_DIR/latex_to_pandoc.py" manuscript.tex manuscript.aux \
        /tmp/manuscript-pandoc.tex

echo "[3/3] Converting..."
pandoc /tmp/manuscript-pandoc.tex \
    --from=latex \
    --to=docx \
    --citeproc \
    --bibliography=references.bib \
    --resource-path=".:figures" \
    --number-sections \
    --metadata reference-section-title="References" \
    --output=manuscript.docx

echo
echo "Wrote $DIR/manuscript.docx"
python3 - <<'PY'
import zipfile, re, os
p = "manuscript.docx"
z = zipfile.ZipFile(p)
xml = z.read("word/document.xml").decode("utf8", "replace")
text = re.sub(r"<[^>]+>", "", xml)
words = len(text.split())
imgs = len([n for n in z.namelist() if n.startswith("word/media/")])
print("  %.1f KB, ~%d words, %d embedded figures" % (os.path.getsize(p)/1024, words, imgs))
for bad, why in (("??", "unresolved cross-reference"),
                 ("\\SI{", "unexpanded siunitx"),
                 ("[@", "unresolved citation")):
    n = text.count(bad)
    if n:
        print("  WARNING: %d x %s (%s)" % (n, bad, why))
PY
