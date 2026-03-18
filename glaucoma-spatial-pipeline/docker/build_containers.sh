#!/usr/bin/env bash
# =============================================================================
# build_containers.sh
# Build, test, and optionally push all pipeline Docker images.
#
# Usage:
#   bash build_containers.sh              # build all, run smoke tests
#   bash build_containers.sh --push       # build + push to ghcr.io
#   bash build_containers.sh --image scvi # build one image
#   bash build_containers.sh --test       # smoke test existing images
# =============================================================================
set -euo pipefail

REGISTRY="ghcr.io/tjmb03"
VERSION="1.0.0"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOCKER_DIR="$REPO_ROOT/docker"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
pass()   { echo -e "${GREEN}✓${NC}  $1"; }
fail()   { echo -e "${RED}✗${NC}  $1"; exit 1; }
header() { echo -e "\n${YELLOW}════════════════════════════════════════${NC}"; \
           echo -e "${YELLOW} $1${NC}"; \
           echo -e "${YELLOW}════════════════════════════════════════${NC}"; }

PUSH=false
TARGET="all"
TEST_ONLY=false

for arg in "$@"; do
    case "$arg" in
        --push)        PUSH=true ;;
        --test)        TEST_ONLY=true ;;
        --image)       shift; TARGET="$1" ;;
        --image=*)     TARGET="${arg#--image=}" ;;
    esac
done

# ── Check Docker is running ───────────────────────────────────────────────────
if ! docker info &> /dev/null; then
    fail "Docker is not running. Start Docker Desktop and retry."
fi
pass "Docker is running"

# ── Build function ────────────────────────────────────────────────────────────
build_image() {
    local name="$1"
    local tag="${REGISTRY}/${name}:${VERSION}"

    header "Building: $name"

    docker build \
        --file "$DOCKER_DIR/${name}/Dockerfile" \
        --tag "${name}:${VERSION}" \
        --tag "${name}:latest" \
        --tag "$tag" \
        --progress=plain \
        "$REPO_ROOT"

    pass "Built: ${name}:${VERSION}"
}

# ── Smoke test function ───────────────────────────────────────────────────────
smoke_test() {
    local name="$1"
    header "Smoke test: $name"

    case "$name" in
        microenvi-star)
            docker run --rm "${name}:${VERSION}" \
                bash -c "STAR --version && samtools --version | head -1"
            ;;
        microenvi-seurat)
            docker run --rm "${name}:${VERSION}" \
                Rscript -e "
                    pkgs <- c('Seurat','SeuratDisk','spdep','DESeq2','limma','BiocParallel')
                    ok <- sapply(pkgs, requireNamespace, quietly=TRUE)
                    stopifnot(all(ok))
                    cat('Seurat:', as.character(packageVersion('Seurat')), '\n')
                    cat('All R packages OK\n')
                "
            ;;
        microenvi-scvi)
            docker run --rm "${name}:${VERSION}" \
                python -c "
import scvi, scanpy, anndata, torch, scipy.sparse as sp, numpy as np
adata = anndata.AnnData(X=sp.csr_matrix(
    np.random.negative_binomial(2, 0.3, size=(100, 50)).astype('float32')))
adata.layers['counts'] = adata.X.copy()
adata.obs['batch'] = ['A']*50 + ['B']*50
scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='batch')
m = scvi.model.SCVI(adata, n_latent=5)
m.train(max_epochs=2)
z = m.get_latent_representation()
assert z.shape == (100, 5)
print('scVI smoke test: PASS')
"
            ;;
        microenvi-scvelo)
            docker run --rm "${name}:${VERSION}" \
                python -c "
import scvelo as scv
print('scvelo:', scv.__version__)
print('scvelo smoke test: PASS')
"
            ;;
        microenvi-cellrank)
            docker run --rm "${name}:${VERSION}" \
                python -c "
import cellrank as cr, libpysal, esda, gseapy
from libpysal.weights import DistanceBand
import esda.moran, numpy as np
coords = np.random.uniform(0, 1000, size=(80, 2))
w = DistanceBand(coords, threshold=200, binary=False)
w.transform = 'r'
mi = esda.moran.Moran(np.random.uniform(0,1,80), w, permutations=49)
assert 0 < mi.p_sim < 1
print('CellRank + spatial Moran test: PASS')
"
            ;;
    esac

    pass "Smoke test passed: $name"
}

# ── Images to process ─────────────────────────────────────────────────────────
ALL_IMAGES=(
    microenvi-star
    microenvi-seurat
    microenvi-scvi
    microenvi-scvelo
    microenvi-cellrank
)

if [[ "$TARGET" == "all" ]]; then
    IMAGES=("${ALL_IMAGES[@]}")
else
    IMAGES=("glaucoma-${TARGET}")
fi

# ── Main ──────────────────────────────────────────────────────────────────────
header "Glaucoma ST Pipeline — Docker Build"
echo "  Registry : $REGISTRY"
echo "  Version  : $VERSION"
echo "  Images   : ${IMAGES[*]}"
echo "  Push     : $PUSH"

if [[ "$TEST_ONLY" == false ]]; then
    for img in "${IMAGES[@]}"; do
        build_image "$img"
    done
fi

for img in "${IMAGES[@]}"; do
    smoke_test "$img"
done

if [[ "$PUSH" == true ]]; then
    header "Pushing to $REGISTRY"

    if ! docker login ghcr.io &> /dev/null; then
        echo "Logging in to ghcr.io..."
        echo "Use a GitHub Personal Access Token with 'write:packages' scope"
        docker login ghcr.io -u tjmb03
    fi

    for img in "${IMAGES[@]}"; do
        echo "Pushing ${REGISTRY}/${img}:${VERSION}..."
        docker push "${REGISTRY}/${img}:${VERSION}"
        pass "Pushed: ${img}"
    done
fi

header "Done"
echo ""
echo "  Images built:"
for img in "${IMAGES[@]}"; do
    size=$(docker image inspect "${img}:${VERSION}" \
        --format='{{.Size}}' 2>/dev/null | \
        awk '{printf "%.1f GB", $1/1073741824}')
    echo "    ${img}:${VERSION}  (${size})"
done
echo ""
echo "  Run locally:"
echo "    docker compose -f docker/docker-compose.yml run scvi"
echo ""
echo "  Push to registry:"
echo "    bash docker/build_containers.sh --push"
