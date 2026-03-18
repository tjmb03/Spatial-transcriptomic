FROM python:3.10-slim AS builder
LABEL maintainer="tjmb03"

RUN apt-get update && apt-get install -y \
    build-essential gcc g++ \
    && rm -rf /var/lib/apt/lists/*

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --no-cache-dir \
    scvelo==0.3.2 \
    scanpy==1.9.8 \
    anndata==0.10.3 \
    numpy==1.26.3 \
    scipy==1.11.4 \
    pandas==2.1.3 \
    matplotlib==3.8.2 \
    loompy==3.0.7 \
    h5py==3.10.0 \
    tqdm==4.66.1

# ── Runtime ───────────────────────────────────────────────────────────────────
FROM python:3.10-slim
LABEL version="1.0.0"
LABEL description="scVelo RNA velocity — dynamical mode"

COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN apt-get update && apt-get install -y libgomp1 \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /pipeline/scripts /pipeline/config
COPY glaucoma-spatial-pipeline/config/preregistration.py /pipeline/config/

RUN python -c "\
import scvelo as scv; \
print('scvelo:', scv.__version__)"

WORKDIR /pipeline
CMD ["python", "scripts/04_scvelo_velocity.py", "--help"]
