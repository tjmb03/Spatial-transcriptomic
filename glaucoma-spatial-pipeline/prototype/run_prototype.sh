#!/usr/bin/env bash
# =============================================================================
# run_prototype.sh
# End-to-end R↔Python handoff prototype on Tran 2019 ONC data.
#
# Runs all four steps:
#   1. Download Tran 2019 data from GEO
#   2. Seurat v5 processing → h5seurat → h5ad  (R, glaucoma-r env)
#   3. Python validation + minimal scVI         (Python, glaucoma-scvi env)
#   4. Print summary
#
# Usage:
#   bash prototype/run_prototype.sh
#   bash prototype/run_prototype.sh --skip_download  # if data already present
# =============================================================================
set -euo pipefail

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
pass()  { echo -e "${GREEN}✓${NC}  $1"; }
fail()  { echo -e "${RED}✗${NC}  $1"; exit 1; }
header(){ echo -e "\n${YELLOW}════════════════════════════════════════${NC}"; \
          echo -e "${YELLOW} $1${NC}"; \
          echo -e "${YELLOW}════════════════════════════════════════${NC}"; }

SKIP_DOWNLOAD=false
for arg in "$@"; do
    [[ "$arg" == "--skip_download" ]] && SKIP_DOWNLOAD=true
done

DATA_DIR="prototype/data/tran2019"
SEURAT_OUT="prototype/output/seurat"
PYTHON_OUT="prototype/output/python"

mkdir -p "$DATA_DIR" "$SEURAT_OUT" "$PYTHON_OUT"

# ── Step 1: Download ──────────────────────────────────────────────────────────
if [[ "$SKIP_DOWNLOAD" == false ]]; then
    header "Step 1 — Download Tran 2019 ONC data"
    bash prototype/00_download_tran2019.sh "$DATA_DIR"
    pass "Download complete"
else
    header "Step 1 — Skipping download (--skip_download)"
    echo "  Data dir: $DATA_DIR"
    ls "$DATA_DIR/" 2>/dev/null || echo "  (empty)"
fi

# ── Step 2: R processing ──────────────────────────────────────────────────────
header "Step 2 — Seurat v5 processing (R)"

# Check we're in the right conda environment
if ! command -v Rscript &>/dev/null; then
    fail "Rscript not found. Run: conda activate glaucoma-r"
fi

R_ENV=$(conda info --active-env-prefix 2>/dev/null || echo "unknown")
echo "  R environment: $R_ENV"
echo "  Seurat version: $(Rscript -e 'cat(as.character(packageVersion("Seurat")))')"

Rscript prototype/01_seurat_processing.R \
    --datadir "$DATA_DIR" \
    --outdir  "$SEURAT_OUT"

H5AD_PATH=$(ls "$SEURAT_OUT"/*.h5ad 2>/dev/null | head -1)
MANIFEST_PATH="$SEURAT_OUT/handoff_manifest.json"

if [ -z "$H5AD_PATH" ]; then
    fail "No h5ad file produced by Seurat step"
fi

pass "h5ad produced: $H5AD_PATH"
echo "  Size: $(du -sh "$H5AD_PATH" | cut -f1)"

# ── Step 3: Python validation ─────────────────────────────────────────────────
header "Step 3 — Python handoff validation + scVI"

# Check we need to switch environments
CURRENT_ENV=$(conda info --json 2>/dev/null | python3 -c \
    "import sys,json; d=json.load(sys.stdin); print(d.get('active_prefix_name','unknown'))" \
    2>/dev/null || echo "unknown")

echo "  Current env: $CURRENT_ENV"

if [[ "$CURRENT_ENV" != "glaucoma-scvi" ]]; then
    echo ""
    echo "  ⚠  Switch to glaucoma-scvi environment and run:"
    echo ""
    echo "     conda activate glaucoma-scvi"
    echo "     python3 prototype/02_python_validation.py \\"
    echo "         --h5ad $H5AD_PATH \\"
    echo "         --manifest $MANIFEST_PATH \\"
    echo "         --outdir $PYTHON_OUT"
    echo ""
    echo "  Then check results at: $PYTHON_OUT/handoff_validation_results.json"
    echo ""
    echo "  (The R step completed successfully — only the Python step remains)"
    exit 0
fi

python3 prototype/02_python_validation.py \
    --h5ad     "$H5AD_PATH" \
    --manifest "$MANIFEST_PATH" \
    --outdir   "$PYTHON_OUT"

# ── Step 4: Summary ───────────────────────────────────────────────────────────
header "Prototype Complete"

RESULTS_FILE="$PYTHON_OUT/handoff_validation_results.json"
if [ -f "$RESULTS_FILE" ]; then
    ALL_PASS=$(python3 -c \
        "import json; d=json.load(open('$RESULTS_FILE')); print(d['all_pass'])")

    if [[ "$ALL_PASS" == "True" ]]; then
        pass "All handoff validation checks passed"
        echo ""
        echo "  The R↔Python bridge is working correctly."
        echo "  Proceed to Step 5 — Nextflow orchestration."
        echo ""
        echo "  Key outputs:"
        echo "    h5ad (R output):       $H5AD_PATH"
        echo "    Validated h5ad:        $PYTHON_OUT/tran2019_validated.h5ad"
        echo "    Validation results:    $RESULTS_FILE"
    else
        fail "Handoff validation failed — see $RESULTS_FILE for details"
    fi
else
    echo "  Run the Python step to complete validation."
fi
