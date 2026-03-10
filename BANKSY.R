

# first install remotes if not installed
install.packages("remotes")
library(remotes)
# then install SeuratWrappers
remotes::install_github("satijalab/seurat-wrappers")

if (!requireNamespace("Banksy", quietly = TRUE)) {
  remotes::install_github("prabhakarlab/Banksy@devel")
}


library(SeuratWrappers)
library(Banksy)


# Make sure you have variable features first
Spatial.008um <- FindVariableFeatures(Spatial.008um, assay = "Spatial", nfeatures = 2000)

# Then run Banksy
object <- RunBanksy(
  object = Spatial.008um,
  lambda = 0.8,
  verbose = TRUE,
  assay = "Spatial",
  slot = "counts",   # <- use slot here (wrapper is old)
  k_geom = 20
)


# 1. Get variable features
var.features <- VariableFeatures(object, assay = "Spatial.008um")

# 2. Make sure they exist in the counts slot
counts.features <- rownames(GetAssayData(object, assay = "Spatial.008um", slot = "counts"))
common.features <- intersect(var.features, counts.features)

# 3. Take the top 10,000 (or all if fewer)
selected.features <- head(common.features, 10000)

set.seed(123)
cells.subset <- sample(Cells(object), 20000)  # e.g. 20k cells
Spatial.sub  <- subset(object, cells = cells.subset)


# 4. Run Banksy
object <- RunBanksy(
  object   = Spatial.sub,
  lambda   = 0.8,
  verbose  = TRUE,
  assay    = "Spatial.008um",
  slot    = "counts",
  k_geom   = 30
)

DefaultAssay(object) <- "BANKSY"
object <- FindVariableFeatures(object, assay = "BANKSY", nfeatures = 2000)

object <- RunPCA(
  object,
  assay = "BANKSY",
  features = VariableFeatures(object, assay = "BANKSY"),
  reduction.name = "banksy.pca"
)
# Use PCA from your Banksy assay
object <- RunPCA(object, assay = "BANKSY", reduction.name = "banksy.pca")

# Build neighbors from Banksy PCA
object <- FindNeighbors(object, reduction = "banksy.pca", dims = 1:30)

# Cluster
object <- FindClusters(object, resolution = 0.5)

# Save clusters into metadata with a clearer name
object$banksy_cluster <- object$seurat_clusters


# Example: cluster the Banksy graph
object <- FindNeighbors(object, graph.name = "banksy")  # graph from RunBanksy
object <- FindClusters(object, resolution = 0.5, graph.name = "banksy")

# This will create a column in metadata
head(object$seurat_clusters)  # default column

cortex <- NormalizeData(cortex, assay = "sketch", verbose = TRUE)
