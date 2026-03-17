#!/usr/bin/env Rscript
# install_r_packages.R
# Run AFTER conda env create -f glaucoma-r.yml
# conda activate glaucoma-r
# Rscript envs/install_r_packages.R

message("Installing R packages not available on conda-forge...")
message("This will take 10-20 minutes on first run.\n")

# ── CRAN ─────────────────────────────────────────────────────────────────────
cran_pkgs <- c("Seurat", "jsonlite", "optparse", "patchwork",
               "ggplot2", "dplyr", "Matrix", "scales")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    message("Already installed: ", pkg)
  }
}

# ── Bioconductor ──────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(version = "3.18", ask = FALSE, update = FALSE)

bioc_pkgs <- c(
  "spdep",            # Moran's I, LISA
  "clusterProfiler",  # pathway enrichment
  "org.Mm.eg.db",     # mouse gene annotation
  "ReactomePA",       # Reactome pathways
  "DESeq2",           # compositional analysis
  "limma",            # DE analysis
  "BiocParallel"      # parallel backend
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing (Bioconductor): ", pkg)
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    message("Already installed: ", pkg)
  }
}

# ── SeuratDisk (GitHub only — not on CRAN or conda-forge) ────────────────────
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  message("Installing SeuratDisk from GitHub...")
  remotes::install_github("mojaveazure/seurat-disk", quiet = TRUE)
} else {
  message("Already installed: SeuratDisk")
}

# ── Verification ──────────────────────────────────────────────────────────────
message("\nVerifying all packages load correctly...")

required <- c("Seurat", "SeuratDisk", "spdep", "clusterProfiler",
              "org.Mm.eg.db", "DESeq2", "optparse", "jsonlite")

failed <- c()
for (pkg in required) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message("  OK: ", pkg, " (", as.character(packageVersion(pkg)), ")")
  } else {
    message("  FAIL: ", pkg)
    failed <- c(failed, pkg)
  }
}

if (length(failed) > 0) {
  stop("The following packages failed to install: ",
       paste(failed, collapse = ", "))
} else {
  message("\nAll R packages installed successfully.")
}
