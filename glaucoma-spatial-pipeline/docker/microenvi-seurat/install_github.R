
# Only install what conda-forge cannot provide
cat("Installing SeuratDisk from GitHub...\n")
remotes::install_github("mojaveazure/seurat-disk", dependencies = TRUE)
if (!requireNamespace("SeuratDisk", quietly = TRUE)) stop("FAILED: SeuratDisk")
cat("SeuratDisk OK\n")
cat("All GitHub packages installed.\n")
