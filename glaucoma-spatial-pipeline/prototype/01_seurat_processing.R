#!/usr/bin/env Rscript
# =============================================================================
# 01_seurat_processing.R
# Seurat v5 QC, clustering, and h5ad export for the R→Python handoff prototype.
# Uses Tran 2019 ONC data (GSE133382) as a public stand-in for real GRAPE data.
#
# What this validates:
#   1. Seurat v5 loads 10x-format matrices correctly
#   2. QC filtering matches expected cell counts
#   3. SeuratDisk::SaveH5Seurat() produces valid HDF5
#   4. SeuratDisk::Convert() → .h5ad preserves all layers
#   5. obs metadata (timepoint, cell_type) survives conversion
#
# Usage:
#   conda activate glaucoma-r
#   Rscript prototype/01_seurat_processing.R --datadir prototype/data/tran2019
# =============================================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratDisk)
    library(optparse)
    library(jsonlite)
    library(dplyr)
    library(ggplot2)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
option_list <- list(
    make_option("--datadir", type="character",
                default="prototype/data/tran2019",
                help="Directory with Tran 2019 data files"),
    make_option("--outdir", type="character",
                default="prototype/output/seurat",
                help="Output directory"),
    make_option("--min_cells", type="integer", default=3,
                help="Minimum cells per gene"),
    make_option("--min_features", type="integer", default=200,
                help="Minimum features per cell"),
    make_option("--max_pct_mito", type="double", default=20.0,
                help="Maximum % mitochondrial reads"),
    make_option("--n_hvg", type="integer", default=2000,
                help="Number of highly variable genes")
)

opts <- parse_args(OptionParser(option_list=option_list))
dir.create(opts$outdir, recursive=TRUE, showWarnings=FALSE)

cat("\n════════════════════════════════════════\n")
cat(" Step 4 — R↔Python Handoff Prototype\n")
cat(" Seurat v5 Processing (Tran 2019 ONC)\n")
cat("════════════════════════════════════════\n\n")

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Loading Tran 2019 data...\n")

# The GEO supplementary files are typically per-sample count matrices.
# Look for either 10x MTX format or dense matrix files.
data_files <- list.files(opts$datadir,
                          pattern="\\.(mtx|mtx\\.gz|csv\\.gz|txt\\.gz)$",
                          full.names=TRUE, recursive=TRUE)

if (length(data_files) == 0) {
    # Try loading directly from a single matrix file
    matrix_files <- list.files(opts$datadir, pattern="matrix", full.names=TRUE)
    cat("No MTX files found. Available files:\n")
    cat(paste(list.files(opts$datadir), collapse="\n"), "\n")
    cat("\nAttempting to load from available files...\n")
}

# ── Flexible loading: handle multiple file formats ────────────────────────────
load_tran_data <- function(datadir) {

    # Format 1: 10x MTX directories
    mtx_dirs <- list.dirs(datadir, recursive=TRUE)
    mtx_dirs <- mtx_dirs[sapply(mtx_dirs, function(d) {
        any(grepl("matrix.mtx", list.files(d)))
    })]

    if (length(mtx_dirs) > 0) {
        cat(sprintf("Found %d 10x MTX directories\n", length(mtx_dirs)))
        sample_list <- lapply(mtx_dirs, function(d) {
            sample_name <- basename(d)
            cat(sprintf("  Loading: %s\n", sample_name))
            counts <- Read10X(d)
            CreateSeuratObject(counts, project=sample_name,
                               min.cells=3, min.features=200)
        })
        names(sample_list) <- basename(mtx_dirs)
        if (length(sample_list) > 1) {
            return(merge(sample_list[[1]],
                         y=sample_list[-1],
                         add.cell.ids=names(sample_list)))
        }
        return(sample_list[[1]])
    }

    # Format 2: Single dense matrix (rows=genes, cols=cells)
    csv_files <- list.files(datadir, pattern="\\.csv(\\.gz)?$",
                             full.names=TRUE)
    if (length(csv_files) > 0) {
        cat(sprintf("Loading dense matrix: %s\n", basename(csv_files[1])))
        mat <- read.csv(csv_files[1], row.names=1, check.names=FALSE)
        return(CreateSeuratObject(as.matrix(mat), min.cells=3,
                                   min.features=200))
    }

    # Format 3: RDS file (pre-processed Seurat object)
    rds_files <- list.files(datadir, pattern="\\.rds$",
                             full.names=TRUE, ignore.case=TRUE)
    if (length(rds_files) > 0) {
        cat(sprintf("Loading RDS: %s\n", basename(rds_files[1])))
        obj <- readRDS(rds_files[1])
        if (!inherits(obj, "Seurat")) {
            obj <- CreateSeuratObject(obj, min.cells=3, min.features=200)
        }
        return(obj)
    }

    stop("No supported data format found in ", datadir,
         "\nSupported: 10x MTX directories, CSV matrix, RDS Seurat object")
}

seurat_obj <- load_tran_data(opts$datadir)
cat(sprintf("\nLoaded: %d cells × %d genes\n",
            ncol(seurat_obj), nrow(seurat_obj)))

# ── Add metadata ──────────────────────────────────────────────────────────────
cat("\nAdding metadata...\n")

# Extract cell names AFTER Seurat object creation
# (R may mangle special characters in column names during CSV loading)
cell_names <- colnames(seurat_obj)
cat(sprintf("  Example cell name (post-Seurat): %s\n", cell_names[1]))

# Try to parse timepoint from cell name
timepoint_map <- c(
    "0h" = "baseline", "0H" = "baseline",
    "12h" = "early", "12H" = "early",
    "1d" = "early", "1D" = "early",
    "3d" = "mid", "3D" = "mid",
    "7d" = "late", "7D" = "late",
    "14d" = "late", "14D" = "late",
    "ctrl" = "baseline", "Ctrl" = "baseline", "control" = "baseline"
)

# Parse sample/timepoint from cell barcodes
seurat_obj$sample_id <- sapply(strsplit(cell_names, "_"), `[`, 1)
seurat_obj$timepoint_raw <- seurat_obj$sample_id

# Map to canonical timepoints
seurat_obj$timepoint <- timepoint_map[seurat_obj$timepoint_raw]
seurat_obj$timepoint[is.na(seurat_obj$timepoint)] <- "unknown"

cat(sprintf("  Timepoints: %s\n",
            paste(table(seurat_obj$timepoint), collapse=", ")))

# ── QC metrics ────────────────────────────────────────────────────────────────
cat("\nComputing QC metrics...\n")

seurat_obj[["pct_mito"]] <- PercentageFeatureSet(
    seurat_obj, pattern="^mt-|^MT-"
)
seurat_obj[["log_total_counts"]] <- log1p(seurat_obj$nCount_RNA)

# QC summary before filtering
qc_before <- data.frame(
    metric = c("n_cells", "median_genes", "median_counts", "median_pct_mito"),
    value  = c(
        ncol(seurat_obj),
        median(seurat_obj$nFeature_RNA),
        median(seurat_obj$nCount_RNA),
        median(seurat_obj$pct_mito)
    )
)
cat("  Before filtering:\n")
print(qc_before)

# ── QC filtering ──────────────────────────────────────────────────────────────
cat("\nApplying QC filters...\n")
cat(sprintf("  min_features: %d\n", opts$min_features))
cat(sprintf("  max_pct_mito: %.1f%%\n", opts$max_pct_mito))

seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA >= opts$min_features &
             pct_mito < opts$max_pct_mito
)

cat(sprintf("  After filtering: %d cells\n", ncol(seurat_obj)))

if (ncol(seurat_obj) < 100) {
    stop("Too few cells after QC filtering (", ncol(seurat_obj), "). ",
         "Check QC thresholds and input data format.")
}

# ── Save raw counts layer BEFORE normalization ────────────────────────────────
# Critical: scVI requires raw counts in a separate layer.
# This is the most common handoff failure — forgetting to save raw counts.
cat("\nSaving raw counts layer (required for scVI)...\n")
seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay5")
seurat_obj <- JoinLayers(seurat_obj)

# Store raw counts in 'counts' layer explicitly
counts_matrix <- GetAssayData(seurat_obj, assay="RNA", layer="counts")
cat(sprintf("  Raw counts matrix: %d genes × %d cells\n",
            nrow(counts_matrix), ncol(counts_matrix)))
cat(sprintf("  Count range: [%.0f, %.0f]\n",
            min(counts_matrix), max(counts_matrix)))

# Verify these are truly raw counts (integers, not normalized)
if (max(counts_matrix) > 0 && all(counts_matrix == floor(counts_matrix))) {
    cat("  ✓ Confirmed: counts are integers (raw, not normalized)\n")
} else {
    warning("Counts appear to be non-integer — scVI requires raw integer counts")
}

# ── Normalization and HVG selection ──────────────────────────────────────────
cat("\nNormalizing and selecting HVGs...\n")
seurat_obj <- NormalizeData(seurat_obj, verbose=FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures=opts$n_hvg,
                                    verbose=FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose=FALSE)

cat(sprintf("  HVGs selected: %d\n", length(VariableFeatures(seurat_obj))))

# ── Dimensionality reduction ──────────────────────────────────────────────────
cat("\nRunning PCA and UMAP...\n")
seurat_obj <- RunPCA(seurat_obj, npcs=30, verbose=FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims=1:20, verbose=FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims=1:20, verbose=FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution=0.5, verbose=FALSE)

cat(sprintf("  Clusters found: %d\n", length(unique(seurat_obj$seurat_clusters))))

# ── Export QC figures ──────────────────────────────────────────────────────────
cat("\nSaving QC figures...\n")

p1 <- VlnPlot(seurat_obj,
              features=c("nFeature_RNA", "nCount_RNA", "pct_mito"),
              ncol=3, pt.size=0) &
      theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(file.path(opts$outdir, "qc_violin.png"),
       p1, width=12, height=5, dpi=150)

p2 <- DimPlot(seurat_obj, reduction="umap",
              group.by=c("seurat_clusters", "timepoint"),
              ncol=2)
ggsave(file.path(opts$outdir, "umap_clusters.png"),
       p2, width=12, height=5, dpi=150)

# ── Save h5seurat ─────────────────────────────────────────────────────────────
cat("\nSaving h5seurat file...\n")
h5seurat_path <- file.path(opts$outdir, "tran2019.h5seurat")

SaveH5Seurat(seurat_obj, filename=h5seurat_path, overwrite=TRUE)
cat(sprintf("  Saved: %s (%.1f MB)\n",
            h5seurat_path,
            file.size(h5seurat_path) / 1e6))

# ── Convert to h5ad ───────────────────────────────────────────────────────────
cat("\nConverting to h5ad (Python-readable)...\n")
h5ad_path <- file.path(opts$outdir, "tran2019.h5ad")

Convert(h5seurat_path, dest="h5ad", overwrite=TRUE)
cat(sprintf("  Saved: %s (%.1f MB)\n",
            h5ad_path,
            file.size(h5ad_path) / 1e6))

# ── Write handoff manifest ────────────────────────────────────────────────────
# This JSON tells the Python side exactly what to expect in the h5ad file.
manifest <- list(
    created       = Sys.time(),
    seurat_version= as.character(packageVersion("Seurat")),
    h5ad_path     = normalizePath(h5ad_path),
    h5seurat_path = normalizePath(h5seurat_path),
    n_cells       = ncol(seurat_obj),
    n_genes       = nrow(seurat_obj),
    n_hvg         = length(VariableFeatures(seurat_obj)),
    n_clusters    = length(unique(seurat_obj$seurat_clusters)),
    layers_expected = list("counts", "data", "scale.data"),
    obs_cols_expected = c("sample_id", "timepoint", "seurat_clusters",
                          "nFeature_RNA", "nCount_RNA", "pct_mito",
                          "log_total_counts"),
    qc_thresholds = list(
        min_features = opts$min_features,
        max_pct_mito = opts$max_pct_mito
    ),
    count_range = list(
        min = min(counts_matrix),
        max = max(counts_matrix)
    ),
    handoff_checksum = digest::digest(h5ad_path, file=TRUE, algo="md5")
)

manifest_path <- file.path(opts$outdir, "handoff_manifest.json")
writeLines(toJSON(manifest, pretty=TRUE, auto_unbox=TRUE), manifest_path)
cat(sprintf("\nHandoff manifest: %s\n", manifest_path))

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n════════════════════════════════════════\n")
cat(" R Processing Complete\n")
cat("════════════════════════════════════════\n")
cat(sprintf("  Cells:    %d\n", ncol(seurat_obj)))
cat(sprintf("  Genes:    %d\n", nrow(seurat_obj)))
cat(sprintf("  HVGs:     %d\n", length(VariableFeatures(seurat_obj))))
cat(sprintf("  Clusters: %d\n", length(unique(seurat_obj$seurat_clusters))))
cat(sprintf("\n  h5ad → Python:  %s\n", h5ad_path))
cat(sprintf("  Manifest:       %s\n", manifest_path))
cat("\n  Next step:\n")
cat("    conda activate glaucoma-scvi\n")
cat(sprintf("    python3 prototype/02_python_validation.py \\\n"))
cat(sprintf("        --h5ad %s \\\n", h5ad_path))
cat(sprintf("        --manifest %s\n", manifest_path))
cat("════════════════════════════════════════\n\n")
