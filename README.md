> **Disclaimer:** This repository contains a preregistered computational 
> pipeline design and associated documentation. All reference datasets 
> cited are publicly available (GEO accession numbers provided in 
> documentation). No proprietary, confidential, patient-derived, or 
> employer-affiliated data is included. Pipeline architecture and 
> analytical frameworks represent independent methodological development 
> conducted outside of any employment context and do not reflect the 
> proprietary methods, data, or intellectual property of any employer 
> or collaborator.
>
> This repository is released under the [MIT License](LICENSE).  
> © 2026 Bo Ma (tjmb03). Reuse with attribution.

# Spatial Transcriptomics — Microenvironment Pipeline

> A preregistered single-nucleus spatial transcriptomics pipeline for documenting pathological microenvironment development, combining Seurat v5 (R) and CellRank 2 (Python).

[![GitHub Pages](https://img.shields.io/badge/docs-live-f5c842?style=flat-square&logo=github)](https://tjmb03.github.io/Spatial-transcriptomic/)
[![License: MIT](https://img.shields.io/badge/license-MIT-3dd6c8?style=flat-square)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Slide--tag-9ba8bc?style=flat-square)]()
[![Preregistered](https://img.shields.io/badge/preregistered-OSF-f5c842?style=flat-square)](https://osf.io)

---

## 📖 Documentation

Full methods, validation suite, and therapeutic target scoring framework:

**[https://tjmb03.github.io/Spatial-transcriptomic/](https://tjmb03.github.io/Spatial-transcriptomic/)**

---

## Overview

This pipeline processes Slide-tag snRNA-seq data from the C57BL/6J microbead occlusion glaucoma model across four disease timepoints (baseline, early 2–4 wk, mid 6–8 wk, late 12+ wk), linking RNA velocity-based fate trajectories to spatial disease spread and druggable therapeutic targets.

```
Raw Slide-tag FASTQ
    ↓  STARsolo — spliced/unspliced/Velocyto
    ↓  Seurat v5 — QC, clustering, label transfer (MRCA reference)
    ↓  scVI — batch correction + 7-test validation suite
    ↓  scVelo — dynamical RNA velocity
    ↓  CellRank 2 — VelocityKernel + SpatialKernel (GPCCA)
    ↓
    ├── Spatial statistics — Moran's I + LISA → primary endpoint
    └── Lineage drivers → pathway enrichment → therapeutic targets
```

**Primary endpoint:** Moran's I of RGC apoptotic fate probability, mid vs baseline, 100 µm bandwidth, Mann-Whitney U (one-sided), BH-FDR, α = 0.05.

---

## Repository Structure

```
Spatial-transcriptomic/
├── glaucoma-spatial-pipeline/        # Pipeline scripts
│   ├── config/
│   │   └── preregistration.py        # SHA256-locked analysis contract
│   ├── scripts/
│   │   ├── 01_starsolo_alignment.sh  # STARsolo + spatial barcode registration
│   │   ├── 02_seurat_processing.R    # Seurat v5 QC, clustering, label transfer
│   │   ├── 03_scvi_integration.py    # scVI batch correction + 7-test validation
│   │   ├── 04_scvelo_velocity.py     # RNA velocity (dynamical mode)
│   │   ├── 05_cellrank_trajectories.py  # Fate probabilities (GPCCA)
│   │   ├── 06_spatial_statistics.py  # Moran's I, LISA, hotspot tracking
│   │   └── 07_lineage_drivers_and_targets.py  # Drivers, pathways, targets
│   ├── validation/
│   │   └── validate_pipeline.py      # 11-checkpoint validation orchestrator
│   └── nextflow/
│       └── main.nf                   # Nextflow DSL2 orchestration
├── glaucoma_quarto_docs/             # Documentation site source
│   ├── _quarto.yml
│   ├── index.qmd
│   ├── methods.qmd
│   ├── validation.qmd
│   ├── targets.qmd
│   └── assets/
├── HD/                               # Visium HD analysis
└── sequence-based/                   # Sequence-based spatial analysis
```

---

## Experimental Design

| Parameter | Value |
|---|---|
| **Model** | C57BL/6J microbead occlusion |
| **Platform** | Slide-tag (true single-nucleus + RNA velocity) |
| **Timepoints** | Baseline / Early (2–4 wk) / Mid (6–8 wk) / Late (12+ wk) |
| **Replicates** | 5 animals × 4 timepoints × 2 eyes = 20 runs |
| **Tissues** | Retina cross-sections + Optic Nerve Head |
| **Control** | Contralateral eye (within-animal) |

---

## Key Features

**Preregistration integrity** — All analysis parameters are SHA256-hashed before data collection. The pipeline refuses to run if the hash is violated.

**7-test scVI validation suite** — Posterior predictive check, batch silhouette, reconstruction accuracy, marker preservation, velocity compatibility, seed stability, and biological validation must all pass before proceeding.

**Combined spatial + velocity kernel** — CellRank 2 uses `w × VelocityKernel + (1−w) × SpatialKernel` with kernel weights locked before fate extraction via grid search.

**3-tier therapeutic target scoring** — Composite score weighting tractability + genetic evidence (GWAS/Mendelian) + driver characteristics, with OpenTargets API integration.

---

## Installation

```bash
# Python environment
conda create -n glaucoma_st python=3.10 -y
conda activate glaucoma_st
pip install scvi-tools==1.2.0 scvelo==0.3.2 cellrank==2.5.0 \
    scanpy anndata libpysal esda scipy statsmodels gseapy

# R packages
Rscript -e 'install.packages(c("Seurat", "SeuratDisk", "optparse", "jsonlite"))'

# Nextflow (for full pipeline orchestration)
curl -s https://get.nextflow.io | bash
```

---

## Quick Start (Test Run)

```bash
# Generate synthetic test data (no real data required)
python glaucoma-spatial-pipeline/scripts/make_test_data.py

# Run scVI integration + 7-test validation
python glaucoma-spatial-pipeline/scripts/03_scvi_integration.py \
    --h5ad_dir test_data \
    --sample_csv test_data/sample_sheet.csv \
    --outdir test_out/scvi_out \
    --skip_t6_stability

# Run full validation report
python glaucoma-spatial-pipeline/validation/validate_pipeline.py \
    --results_dir test_out \
    --outdir test_out/validation_report
```

---

## Reference Datasets

| Dataset | Access | Role |
|---|---|---|
| MRCA — Li et al. 2024 | [GSE243413](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243413) | Primary label transfer reference |
| Tran/Sanes et al. 2019 | [GSE133382](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133382) | RGC subtype atlas |
| Benhar et al. 2023 | [GSE199317](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199317) | Non-neuronal temporal atlas |
| Keuthan et al. 2023 | [GSE241782](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241782) | Microbead bulk RNA-seq benchmark |

---

## Citation

If you use this pipeline, please cite the key underlying tools:

- **Seurat v5** — Hao et al. *Nature Biotechnology* 2024
- **scVI** — Lopez et al. *Nature Methods* 2018
- **scVelo** — Bergen et al. *Nature Biotechnology* 2020
- **CellRank 2** — Weiler et al. *Nature Methods* 2024
- **STARsolo** — Kaminow et al. *bioRxiv* 2021

---

## License

MIT License — see [LICENSE](LICENSE) for details.
> © 2026 tjmb03. This project is provided for educational and methodological
demonstration purposes.
