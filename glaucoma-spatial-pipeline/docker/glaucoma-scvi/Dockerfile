# ── Stage 1: builder ─────────────────────────────────────────────────────────
# Separate build stage keeps the final image lean
FROM python:3.10-slim AS builder

RUN apt-get update && apt-get install -y \
    build-essential gcc g++ \
    && rm -rf /var/lib/apt/lists/*

# Install all Python packages into a virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# PyTorch CPU (swap this line for GPU build — see bottom of file)
RUN pip install --no-cache-dir \
    torch==2.1.2 \
    --index-url https://download.pytorch.org/whl/cpu

RUN pip install --no-cache-dir \
    scvi-tools==1.2.0 \
    scanpy==1.9.8 \
    anndata==0.10.3 \
    mudata==0.2.3 \
    numpy==1.26.3 \
    scipy==1.11.4 \
    pandas==2.1.3 \
    scikit-learn==1.3.2 \
    matplotlib==3.8.2 \
    seaborn==0.13.0 \
    umap-learn==0.5.5 \
    leidenalg==0.10.2 \
    igraph==0.11.3 \
    h5py==3.10.0 \
    loompy==3.0.7 \
    zarr==2.16.1 \
    tqdm==4.66.1 \
    rich==13.7.0 \
    pyyaml==6.0.1 \
    requests==2.31.0

# ── Stage 2: runtime ──────────────────────────────────────────────────────────
FROM python:3.10-slim
LABEL maintainer="tjmb03"
LABEL version="1.0.0"
LABEL description="scVI batch correction + 7-test validation suite"

# Copy only the venv from builder — no build tools in final image
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# ── Runtime system deps ───────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# ── Pipeline scripts ──────────────────────────────────────────────────────────
RUN mkdir -p /pipeline/scripts /pipeline/config /pipeline/validation
COPY microenvi-spatial-pipeline/config/preregistration.py /pipeline/config/

# ── Smoke test ────────────────────────────────────────────────────────────────
RUN python -c "\
import scvi, scanpy, anndata, torch; \
print('scvi:', scvi.__version__); \
print('scanpy:', scanpy.__version__); \
print('torch:', torch.__version__); \
print('CUDA available:', torch.cuda.is_available())"

WORKDIR /pipeline
CMD ["python", "scripts/03_scvi_integration.py", "--help"]

# ── GPU variant ───────────────────────────────────────────────────────────────
# To build a GPU-enabled image for HPC:
#
#   1. Change the base image in Stage 1 to:
#      FROM pytorch/pytorch:2.1.0-cuda11.8-cudnn8-devel AS builder
#
#   2. Remove the CPU PyTorch pip install line
#      (PyTorch is already in the pytorch base image)
#
#   3. Change Stage 2 base to:
#      FROM pytorch/pytorch:2.1.0-cuda11.8-cudnn8-runtime
#
#   4. Build with:
#      docker build -t microenvi-scvi:gpu -f Dockerfile.gpu .
#
# The GPU image will be ~8GB vs ~3GB for CPU.
# Use the CPU image for local development and CI.
# Use the GPU image for HPC/SLURM execution.
