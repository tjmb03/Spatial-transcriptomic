library(Seurat)

# Assume you already loaded each binning separately
Spatial.008um <- Load10X_Spatial(data.dir = "~/Documents/Seurat/spatial-mousebrain/square_008um", filename = "filtered_feature_bc_matrix.h5")
Spatial.016um <- Load10X_Spatial(data.dir = "~/Documents/Seurat/spatial-mousebrain/square_016um", filename = "filtered_feature_bc_matrix.h5")

# Merge the two Seurat objects
object <- merge(
  x = Spatial.008um,
  y = Spatial.016um,
  add.cell.ids = c("Spatial.008um", "Spatial.016um"),
  project = "SpatialMerged"
)

# Check assays
Assays(object)
DefaultAssay(object) <- "Spatial.008um"


#If you want both binnings as separate assays
# -----------------------------
# 1) Specify file paths
# -----------------------------
h5_8um  <- "~/Documents/Seurat/spatial-mousebrain/square_008um/filtered_feature_bc_matrix.h5"
h5_16um <- "~/Documents/Seurat/spatial-mousebrain/square_016um/filtered_feature_bc_matrix.h5"

# -----------------------------
# 2) Load each H5 as a Seurat object
# -----------------------------
obj16 <- CreateSeuratObject(counts = Read10X_h5(h5_16um),
                            assay = "Spatial.016um",
                            project = "VisiumHD")

obj8  <- CreateSeuratObject(counts = Read10X_h5(h5_8um),
                            assay = "Spatial.008um",
                            project = "VisiumHD")

# -----------------------------
# 3) Merge the objects
# -----------------------------
# Add a prefix to each barcode to keep them unique
obj <- merge(obj16, y = obj8, add.cell.ids = c("16um", "08um"))

# Check the assays
Assays(obj)
# Should show: "Spatial.016um" and "Spatial.008um"


library(Seurat)
library(SeuratData)  # optional, for example datasets
library(hdf5r)

# -----------------------------
# 1) Specify file paths
# -----------------------------
# 16 µm
h5_16um <- "~/Documents/Seurat/spatial-mousebrain/square_016um/filtered_feature_bc_matrix.h5"
img_16um_dir <- "~/Documents/Seurat/spatial-mousebrain/square_016um/spatial"  # folder with tissue_lowres_image.png and tissue_positions_list.csv

# 8 µm
h5_8um <- "~/Documents/Seurat/spatial-mousebrain/square_008um/filtered_feature_bc_matrix.h5"
img_8um_dir <- "~/Documents/Seurat/spatial-mousebrain/square_008um/spatial"

# -----------------------------
# 2) Load counts
# -----------------------------
counts16 <- Read10X_h5(h5_16um)
counts8  <- Read10X_h5(h5_8um)

# -----------------------------
# 3) Create Seurat objects
# -----------------------------
obj16 <- CreateSeuratObject(counts = counts16, assay = "Spatial.016um", project = "VisiumHD")
obj8  <- CreateSeuratObject(counts = counts8, assay = "Spatial.008um", project = "VisiumHD")

# -----------------------------
# 4) Attach spatial images
# -----------------------------
# Helper to read and attach a Visium image
attach_image <- function(obj, image_dir, assay_name){
  # Read the tissue image
  img <- Read10X_Image(image_dir)
  # Name the image slot same as assay
  img_name <- paste0(assay_name, "_image")
  obj[[img_name]] <- img
  # Set the default assay for this image
  DefaultAssay(obj[[img_name]]) <- assay_name
  return(obj)
}

obj16 <- attach_image(obj16, img_16um_dir, "Spatial.016um")
obj8  <- attach_image(obj8, img_8um_dir, "Spatial.008um")

# -----------------------------
# 5) Merge Seurat objects
# -----------------------------
obj <- merge(obj16, y = obj8, add.cell.ids = c("16um", "08um"))

# Check assays and images
Assays(obj)
Images(obj)

# -----------------------------
# 6) Preprocess each assay separately
# -----------------------------
DefaultAssay(obj) <- "Spatial.008um"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)

DefaultAssay(obj) <- "Spatial.016um"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)

# -----------------------------
# 7) Run PCA/UMAP per assay
# -----------------------------
DefaultAssay(obj) <- "Spatial.008um"
obj <- RunPCA(obj, features = VariableFeatures(obj), reduction.name = "pca.008um")
obj <- RunUMAP(obj, reduction = "pca.008um", dims = 1:30, reduction.name = "umap.008um")

DefaultAssay(obj) <- "Spatial.016um"
obj <- RunPCA(obj, features = VariableFeatures(obj), reduction.name = "pca.016um")
obj <- RunUMAP(obj, reduction = "pca.016um", dims = 1:30, reduction.name = "umap.016um")

# -----------------------------
# 8) Quick spatial plot example
# -----------------------------
DefaultAssay(obj) <- "Spatial.008um"
SpatialFeaturePlot(obj, features = c("Mbp", "Calm2"), images = "Spatial.008um_image")

DefaultAssay(obj) <- "Spatial.016um"
SpatialFeaturePlot(obj, features = c("GeneA", "GeneB"), images = "Spatial.016um_image")








#Tried but didn't work cause [[<- cannot add new cells that don’t exist in the original Seurat object.
# Rename assays first
Spatial.008um[["Spatial.008um"]] <- Spatial.008um@assays$Spatial

Assays(Spatial.008um)

# Setting default assay
DefaultAssay(Spatial.008um) <- "Spatial.008um"

# Delete one assay, e.g. "Spatial"
Spatial.008um[["Spatial"]] <- NULL

# First make sure barcodes are identical
all.equal(colnames(Spatial.008um), colnames(Spatial.016um))

# Combine assays into one Seurat object
Spatial.008um[["Spatial.016um"]] <- Spatial.016um@assays$Spatial

# Check assays
Assays(Spatial.008um)
