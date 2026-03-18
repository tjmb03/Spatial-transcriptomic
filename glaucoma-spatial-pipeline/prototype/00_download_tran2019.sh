#!/usr/bin/env bash
# =============================================================================
# 00_download_tran2019.sh
# Download Tran et al. 2019 (GSE133382) RGC single-cell RNA-seq data.
# This is the ONC (optic nerve crush) time series used as handoff prototype.
#
# Dataset: 46 RGC subtypes, 6 timepoints post-ONC (0h, 12h, 1d, 3d, 7d, 14d)
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133382
# SCP: https://singlecell.broadinstitute.org/single_cell/study/SCP509
#
# Usage:
#   conda activate glaucoma-r
#   bash prototype/00_download_tran2019.sh
# =============================================================================
set -euo pipefail

OUTDIR="${1:-prototype/data/tran2019}"
mkdir -p "$OUTDIR"

echo "Downloading Tran 2019 ONC dataset to $OUTDIR..."
echo ""

# ── Option 1: GEO FTP (full count matrix) ─────────────────────────────────────
# The processed count matrix is available as a supplementary file on GEO.
# This is the fastest download for prototyping.

BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133382/suppl"

FILES=(
    "GSE133382_RAW.tar"
)

echo "Downloading from GEO FTP..."
for f in "${FILES[@]}"; do
    if [ ! -f "$OUTDIR/$f" ]; then
        echo "  Downloading: $f"
        curl -L "$BASE_URL/$f" -o "$OUTDIR/$f" --progress-bar
    else
        echo "  Already exists: $f"
    fi
done

# ── Extract ───────────────────────────────────────────────────────────────────
echo ""
echo "Extracting..."
cd "$OUTDIR"
tar -xf GSE133382_RAW.tar

echo ""
echo "Files downloaded:"
ls -lh "$OUTDIR/"

echo ""
echo "Done. Next step:"
echo "  conda activate glaucoma-r"
echo "  Rscript prototype/01_seurat_processing.R --datadir $OUTDIR"
