FROM python:3.10-slim AS builder
LABEL maintainer="tjmb03"

# GEOS/PROJ needed for libpysal — must install from source in slim image
RUN apt-get update && apt-get install -y \
    build-essential gcc g++ \
    libgeos-dev libproj-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install in dependency order to maximise Docker layer cache reuse
RUN pip install --no-cache-dir \
    numpy==1.26.3 \
    scipy==1.11.4 \
    pandas==2.1.3

RUN pip install --no-cache-dir \
    scanpy==1.9.8 \
    anndata==0.10.3 \
    scikit-learn==1.3.2

RUN pip install --no-cache-dir \
    scvelo==0.3.2 \
    cellrank==2.0.7

RUN pip install --no-cache-dir \
    libpysal==4.9.2 \
    esda==2.5.1 \
    splot==1.1.5.post1

RUN pip install --no-cache-dir \
    gseapy==1.1.1 \
    statsmodels==0.14.1 \
    matplotlib==3.8.2 \
    seaborn==0.13.0 \
    requests==2.31.0 \
    tqdm==4.66.1 \
    h5py==3.10.0 \
    loompy==3.0.7

# ── Runtime ───────────────────────────────────────────────────────────────────
FROM python:3.10-slim
LABEL version="1.0.0"
LABEL description="CellRank 2 trajectories + spatial statistics + lineage drivers"

# Runtime needs GEOS shared libraries
RUN apt-get update && apt-get install -y \
    libgeos-dev \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN mkdir -p /pipeline/scripts /pipeline/config
COPY microenvi-spatial-pipeline/config/preregistration.py /pipeline/config/

# ── Smoke test ────────────────────────────────────────────────────────────────
RUN python -c "\
import cellrank as cr; \
import libpysal; \
import esda; \
import gseapy; \
print('cellrank:', cr.__version__); \
print('libpysal:', libpysal.__version__); \
print('esda:', esda.__version__)"

WORKDIR /pipeline
CMD ["python", "scripts/05_cellrank_trajectories.py", "--help"]
