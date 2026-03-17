
# Force C++17 for packages like fgsea that require it
cat("Setting up Makevars for C++17...\n")
makevar_dir <- path.expand("~/.R")
dir.create(makevar_dir, showWarnings = FALSE)
writeLines(
    c("CXX17 = g++ -std=c++17",
      "CXX17FLAGS = -O2 -fPIC",
      "CXXFLAGS = -O2 -fPIC -std=c++17"),
    file.path(makevar_dir, "Makevars")
)
cat("Makevars written\n")


cran <- "https://cloud.r-project.org"

# ── Step 1: base dependencies first ───────────────────────────────────────
cat("Installing base deps...\n")
for (p in c("ggplot2", "scales", "Rcpp", "remotes",
            "optparse", "jsonlite", "RColorBrewer", "viridis")) {
    install.packages(p, repos = cran, dependencies = TRUE)
    if (!requireNamespace(p, quietly = TRUE))
        stop("FAILED: ", p)
    cat("OK:", p, "\n")
}

# ── Step 2: pinned packages for R 4.3.x ───────────────────────────────────
cat("Installing Matrix (pinned)...\n")
install.packages(
    "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz",
    repos = NULL, type = "source"
)
if (!requireNamespace("Matrix", quietly = TRUE)) stop("FAILED: Matrix")
cat("Matrix:", as.character(packageVersion("Matrix")), "\n")

cat("Installing ggrepel (pinned)...\n")
install.packages(
    "https://cran.r-project.org/src/contrib/Archive/ggrepel/ggrepel_0.9.5.tar.gz",
    repos = NULL, type = "source"
)
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("FAILED: ggrepel")
cat("ggrepel:", as.character(packageVersion("ggrepel")), "\n")

# ── Step 3: remaining CRAN deps ────────────────────────────────────────────
for (p in c("patchwork","dplyr","tibble","stringr","hdf5r",
            "future","progressr","cowplot","uwot","igraph",
            "leiden","irlba","fitdistrplus","reticulate",
            "spatstat.explore","spatstat.geom")) {
    cat("Installing:", p, "\n")
    install.packages(p, repos = cran, dependencies = TRUE)
    if (!requireNamespace(p, quietly = TRUE))
        cat("WARN:", p, "\n")
    else
        cat("OK:", p, "\n")
}

# ── Step 4: Seurat ─────────────────────────────────────────────────────────
cat("Installing Seurat...\n")
install.packages("Seurat", repos = cran, dependencies = TRUE)
if (!requireNamespace("Seurat", quietly = TRUE)) stop("FAILED: Seurat")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")

# ── Step 5: SeuratDisk ─────────────────────────────────────────────────────
cat("Installing SeuratDisk...\n")
remotes::install_github("mojaveazure/seurat-disk", dependencies = TRUE)
if (!requireNamespace("SeuratDisk", quietly = TRUE)) stop("FAILED: SeuratDisk")
cat("SeuratDisk OK\n")

# Bioconductor installed in separate Dockerfile layers