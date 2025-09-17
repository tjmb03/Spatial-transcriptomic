# Seurat Spatial Transcriptomics Workflow in R

# 1. Install and load Seurat
install.packages("Seurat")  # or use Satija Lab repo for latest version
library(Seurat)
library(ggplot2)

# 2. Load spatial data (10x Genomics Visium)
spatial_data <- Load10X_Spatial(data.dir = "path/to/your/data")

# 3. Quality control and filtering
spatial_data[["percent.mt"]] <- PercentageFeatureSet(spatial_data, pattern = "^MT-")
spatial_data <- subset(spatial_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 4. Normalize data using SCTransform
spatial_data <- SCTransform(spatial_data, assay = "Spatial", verbose = FALSE)

# 5. Dimensionality reduction
spatial_data <- RunPCA(spatial_data, assay = "SCT", verbose = FALSE)
spatial_data <- FindNeighbors(spatial_data, dims = 1:30)
spatial_data <- FindClusters(spatial_data, resolution = 0.5)
spatial_data <- RunUMAP(spatial_data, dims = 1:30)

# 6. Visualize spatial gene expression
SpatialDimPlot(spatial_data, label = TRUE)
SpatialFeaturePlot(spatial_data, features = "GeneName")

# 7. Identify spatially variable features
spatial_data <- FindSpatiallyVariableFeatures(spatial_data, assay = "SCT", selection.method = "markvariogram")
VariableFeatures(spatial_data)

# 8. Downstream analysis
# - Compare clusters to tissue anatomy
# - Perform differential expression between clusters
# - Integrate with single-cell RNA-seq data if available
