# First make sure you have devtools or remotes
install.packages("remotes")

# Then install SeuratData from GitHub
remotes::install_github("satijalab/seurat-data")

# Load the package
library(SeuratData)

# Install the brain dataset
InstallData("stxBrain")

# Load the dataset
data("stxBrain")

# Load the expression data
expr.url <- 'http://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Brain_Sagittal_Anterior/V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5'
curl::curl_download(url = expr.url, destfile = basename(path = expr.url))
expr.data <- Seurat::Read10X_h5(filename = basename(path = expr.url))
anterior1 <- Seurat::CreateSeuratObject(counts = expr.data, project = 'anterior1', assay = 'Spatial')
anterior1$slice <- 1
anterior1$region <- 'anterior'
# Load the image data
img.url <- 'http://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Brain_Sagittal_Anterior/V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz'
curl::curl_download(url = img.url, destfile = basename(path = img.url))
untar(tarfile = basename(path = img.url))
img <- Seurat::Read10X_Image(image.dir = 'spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = anterior1)]
anterior1[['anterior1']] <- img
# Clean up downloaded files
unlink(x = c(basename(path = c(expr.url, img.url)), 'spatial'), recursive = TRUE)
}

# Create directories if missing
dir.create("~/checks/output/images", recursive = TRUE, showWarnings = FALSE)

brain <- Load10X_Spatial(
  data.dir = "~/checks/spatial_vignette/",
  filename = "V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5"
)
colnames(brain@meta.data)

# Extract pixel coordinates
coords <- brain@images[["slice1"]]@boundaries[["centroids"]]@coords

# Build logical mask (invert = TRUE means we negate the condition)
cells_to_keep <- rownames(coords[
  !(coords[, "y"] > 400 | coords[, "x"] < 150) &   # first condition
    !(coords[, "y"] > 275 & coords[, "x"] > 370) &   # second condition
    !(coords[, "y"] > 250 & coords[, "x"] > 440)     # third condition
  , ])


# Subset the Seurat object
cortex <- subset(cortex, cells = cells_to_keep)
