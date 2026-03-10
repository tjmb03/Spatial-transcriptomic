#!/usr/bin/env bash
# =============================================================================
# render.sh — Build and optionally publish the Quarto documentation site
#
# Usage:
#   bash render.sh              # preview locally
#   bash render.sh --publish    # render + push to GitHub Pages
#   bash render.sh --pdf        # render PDF version (requires LaTeX)
#
# Prerequisites:
#   brew install quarto          # macOS
#   # or: https://quarto.org/docs/get-started/
#
#   pip install jupyter           # for Python code blocks (if any)
#   R -e 'install.packages("rmarkdown")'
# =============================================================================
set -euo pipefail

PUBLISH=false
PDF=false
PREVIEW=false

for arg in "$@"; do
    case "$arg" in
        --publish) PUBLISH=true ;;
        --pdf)     PDF=true ;;
        --preview) PREVIEW=true ;;
    esac
done

SITE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "════════════════════════════════════════"
echo " Glaucoma ST Pipeline — Quarto Build"
echo "════════════════════════════════════════"
echo " Site root : $SITE_DIR"
echo " Publish   : $PUBLISH"
echo " PDF       : $PDF"
echo ""

# ── Check Quarto version ──────────────────────────────────────────────────────
if ! command -v quarto &> /dev/null; then
    echo "❌  Quarto not found."
    echo "    Install from: https://quarto.org/docs/get-started/"
    exit 1
fi

QUARTO_VERSION=$(quarto --version)
echo "✓  Quarto $QUARTO_VERSION"

# ── Render HTML site ─────────────────────────────────────────────────────────
echo ""
echo "Rendering HTML site..."
quarto render "$SITE_DIR" --to html

echo "✓  HTML site rendered → $SITE_DIR/_site/"

# ── Optional: render PDF ──────────────────────────────────────────────────────
if [[ "$PDF" == true ]]; then
    echo ""
    echo "Rendering PDF..."

    # Check for tinytex
    if ! command -v pdflatex &> /dev/null && ! quarto check 2>&1 | grep -q "tinytex"; then
        echo "  Installing TinyTeX..."
        quarto install tinytex
    fi

    quarto render "$SITE_DIR" --to pdf
    echo "✓  PDF rendered"
fi

# ── Preview locally ───────────────────────────────────────────────────────────
if [[ "$PREVIEW" == true ]]; then
    echo ""
    echo "Starting preview server at http://localhost:4321"
    quarto preview "$SITE_DIR" --port 4321 --no-browser
fi

# ── Publish to GitHub Pages ───────────────────────────────────────────────────
if [[ "$PUBLISH" == true ]]; then
    echo ""
    echo "Publishing to GitHub Pages..."

    # Requires: git remote origin pointing to GitHub repo
    # The gh-pages branch is created automatically by quarto publish
    quarto publish gh-pages "$SITE_DIR" --no-prompt

    echo "✓  Published → https://$(git -C "$SITE_DIR" remote get-url origin | \
        sed 's/.*github.com[:/]//' | sed 's/\.git$//' | \
        awk -F/ '{print $1".github.io/"$2}')"
fi

echo ""
echo "════════════════════════════════════════"
echo " Done."
echo "════════════════════════════════════════"
