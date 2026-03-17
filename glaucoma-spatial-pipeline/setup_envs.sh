#!/usr/bin/env bash
# =============================================================================
# setup_envs.sh
# Build and verify all three conda environments for the glaucoma ST pipeline.
#
# Usage:
#   bash setup_envs.sh              # build all three environments
#   bash setup_envs.sh --verify     # verify existing environments only
#   bash setup_envs.sh --env scvi   # build one environment (r / scvi / velocity)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="$SCRIPT_DIR/envs"

# ── Colour output ─────────────────────────────────────────────────────────────
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

pass()  { echo -e "${GREEN}✓${NC}  $1"; }
warn()  { echo -e "${YELLOW}⚠${NC}  $1"; }
fail()  { echo -e "${RED}✗${NC}  $1"; }
header(){ echo -e "\n${YELLOW}════════════════════════════════════════${NC}"; \
          echo -e "${YELLOW} $1${NC}"; \
          echo -e "${YELLOW}════════════════════════════════════════${NC}"; }

# ── Parse arguments ───────────────────────────────────────────────────────────
VERIFY_ONLY=false
TARGET_ENV="all"

for arg in "$@"; do
    case "$arg" in
        --verify)     VERIFY_ONLY=true ;;
        --env)        shift; TARGET_ENV="$1" ;;
        --env=*)      TARGET_ENV="${arg#--env=}" ;;
    esac
done

# ── Check conda is available ──────────────────────────────────────────────────
if ! command -v conda &> /dev/null; then
    fail "conda not found. Install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

CONDA_VERSION=$(conda --version)
pass "Found $CONDA_VERSION"

# ── Build function ────────────────────────────────────────────────────────────
build_env() {
    local env_name="$1"
    local yml_file="$ENV_DIR/${env_name}.yml"

    header "Building: $env_name"

    if ! [ -f "$yml_file" ]; then
        fail "Environment file not found: $yml_file"
        return 1
    fi

    # Remove existing environment if present
    if conda env list | grep -q "^${env_name} "; then
        warn "Environment $env_name already exists — removing and rebuilding"
        conda env remove -n "$env_name" -y
    fi

    echo "Creating environment from $yml_file..."
    conda env create -f "$yml_file"

    pass "Built: $env_name"
}

# ── Verify function ───────────────────────────────────────────────────────────
verify_r_env() {
    header "Verifying: glaucoma-r"

    conda run -n glaucoma-r Rscript - << 'REOF'
pkgs <- c("Seurat", "SeuratDisk", "optparse", "jsonlite",
          "spdep", "clusterProfiler", "org.Mm.eg.db")
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
if (length(missing) > 0) {
    cat("MISSING:", paste(missing, collapse=", "), "\n")
    quit(status=1)
}

# Version checks
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("spdep version:", as.character(packageVersion("spdep")), "\n")

# Test Moran's I (confirms spdep works)
library(spdep)
set.seed(42)
coords <- cbind(runif(50), runif(50))
nb <- knn2nb(knearneigh(coords, k=5))
lw <- nb2listw(nb, style="W")
x  <- runif(50)
mi <- moran.test(x, lw)
cat("Moran's I test: PASS (I =", round(mi$estimate[1], 3), ")\n")
cat("glaucoma-r: ALL CHECKS PASSED\n")
REOF

    pass "glaucoma-r verified"
}

verify_scvi_env() {
    header "Verifying: glaucoma-scvi"

    conda run -n glaucoma-scvi python3 - << 'PYEOF'
import sys

checks = {}

# Core imports
try:
    import scvi; checks["scvi"] = scvi.__version__
except ImportError as e:
    print(f"FAIL: scvi import error: {e}"); sys.exit(1)

try:
    import scanpy as sc; checks["scanpy"] = sc.__version__
except ImportError as e:
    print(f"FAIL: scanpy import error: {e}"); sys.exit(1)

try:
    import anndata as ad; checks["anndata"] = ad.__version__
except ImportError as e:
    print(f"FAIL: anndata import error: {e}"); sys.exit(1)

try:
    import torch; checks["torch"] = torch.__version__
    checks["cuda"] = str(torch.cuda.is_available())
except ImportError as e:
    print(f"FAIL: torch import error: {e}"); sys.exit(1)

# Functional test: build and train a minimal scVI model
import numpy as np
import scipy.sparse as sp

n_cells, n_genes = 200, 100
counts = sp.csr_matrix(
    np.random.negative_binomial(2, 0.3, size=(n_cells, n_genes))
)
adata = ad.AnnData(X=counts)
adata.layers["counts"] = counts
adata.obs["batch"] = ["A"] * 100 + ["B"] * 100

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model = scvi.model.SCVI(adata, n_latent=5, n_layers=1)
model.train(max_epochs=3, progress_bar=False)
latent = model.get_latent_representation()
assert latent.shape == (n_cells, 5), f"Expected (200,5), got {latent.shape}"

print("Versions:", checks)
print("scVI functional test: PASS")
print("glaucoma-scvi: ALL CHECKS PASSED")
PYEOF

    pass "glaucoma-scvi verified"
}

verify_velocity_env() {
    header "Verifying: glaucoma-velocity"

    conda run -n glaucoma-velocity python3 - << 'PYEOF'
import sys

checks = {}

try:
    import scvelo as scv; checks["scvelo"] = scv.__version__
except ImportError as e:
    print(f"FAIL: scvelo import error: {e}"); sys.exit(1)

try:
    import cellrank as cr; checks["cellrank"] = cr.__version__
except ImportError as e:
    print(f"FAIL: cellrank import error: {e}"); sys.exit(1)

try:
    import libpysal; checks["libpysal"] = libpysal.__version__
except ImportError as e:
    print(f"FAIL: libpysal import error: {e}"); sys.exit(1)

try:
    import esda; checks["esda"] = esda.__version__
except ImportError as e:
    print(f"FAIL: esda import error: {e}"); sys.exit(1)

try:
    import gseapy; checks["gseapy"] = gseapy.__version__
except ImportError as e:
    print(f"FAIL: gseapy import error: {e}"); sys.exit(1)

# Functional test: spatial weights + Moran's I
import numpy as np
from libpysal.weights import DistanceBand
import esda.moran

np.random.seed(42)
coords = np.random.uniform(0, 1000, size=(100, 2))
w = DistanceBand(coords, threshold=200, binary=False)
w.transform = 'r'  # row-standardize
y = np.random.uniform(0, 1, 100)
mi = esda.moran.Moran(y, w, permutations=99)
assert 0 < mi.p_sim < 1, "Moran's I p-value out of range"

print("Versions:", checks)
print("Spatial Moran's I test: PASS (I =", round(mi.I, 3), ")")
print("glaucoma-velocity: ALL CHECKS PASSED")
PYEOF

    pass "glaucoma-velocity verified"
}

# ── Main ──────────────────────────────────────────────────────────────────────
echo ""
echo "  Glaucoma ST Pipeline — Environment Setup"
echo "  Conda prefix: $(conda info --base)"
echo ""

if [[ "$VERIFY_ONLY" == true ]]; then
    verify_r_env
    verify_scvi_env
    verify_velocity_env
    echo ""
    pass "All environments verified successfully."
    exit 0
fi

# Build selected environments
case "$TARGET_ENV" in
    "all")
        build_env "glaucoma-r"
        build_env "glaucoma-scvi"
        build_env "glaucoma-velocity"
        verify_r_env
        verify_scvi_env
        verify_velocity_env
        ;;
    "r")
        build_env "glaucoma-r"
        verify_r_env
        ;;
    "scvi")
        build_env "glaucoma-scvi"
        verify_scvi_env
        ;;
    "velocity")
        build_env "glaucoma-velocity"
        verify_velocity_env
        ;;
    *)
        fail "Unknown environment: $TARGET_ENV. Use: r / scvi / velocity / all"
        exit 1
        ;;
esac

echo ""
header "Setup Complete"
echo ""
echo "  Activate environments with:"
echo "    conda activate glaucoma-r         # for Seurat processing"
echo "    conda activate glaucoma-scvi      # for batch correction"
echo "    conda activate glaucoma-velocity  # for velocity + CellRank"
echo ""
echo "  Verify at any time with:"
echo "    bash setup_envs.sh --verify"
echo ""
