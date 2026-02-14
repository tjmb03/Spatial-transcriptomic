# Spatial Transcriptomics Analysis (R / Seurat)

Reproducible spatial transcriptomics workflows built in R using Seurat.

This repository demonstrates structured pipelines for spatial data preprocessing, clustering, and spatially informed domain detection across high-resolution and sequencing-based platforms.

---

## Contents

### 📂 HD/
High-definition spatial workflow examples.

### 📂 sequence-based/
Sequencing-based spatial analysis pipeline.

### 📄 Seurat Spatial Workflow.R  
Core Seurat spatial pipeline:
- QC and filtering  
- Normalization  
- Dimensionality reduction  
- Clustering  
- Spatial visualization  

### 📄 BANKSY.R  
Spatially informed clustering using BANKSY:
- Neighborhood modeling  
- Spatial smoothing  
- Domain refinement  

---

## Capabilities Demonstrated

- Seurat spatial object construction  
- QC and normalization of spatial data  
- PCA / UMAP embedding  
- Tissue overlay visualization  
- Spatially aware clustering  

---

## Requirements

- R ≥ 4.2  
- Seurat  
- ggplot2  
- dplyr  

Install core packages:

```r
install.packages(c("Seurat", "ggplot2", "dplyr"))

