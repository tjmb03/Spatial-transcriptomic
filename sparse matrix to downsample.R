library(Matrix)
library(Seurat)

# Convert to sparse first
counts <- as(counts, "dgCMatrix")

# Downsample cells to, e.g., 20k
set.seed(123)
cells_subset <- sample(colnames(counts), 20000)
counts_subset <- counts[, cells_subset]
ncol(counts_subset)
# Now create Seurat object
temp_obj <- CreateSeuratObject(counts = counts_subset)


# Create a temporary Seurat object to compute variable genes
temp_obj <- CreateSeuratObject(counts = counts)

# Normalize & find variable features
temp_obj <- NormalizeData(temp_obj)
temp_obj <- FindVariableFeatures(temp_obj, selection.method = "vst", nfeatures = 5000)

# Keep only the top 5k variable genes
features <- VariableFeatures(temp_obj)
counts_subset <- counts[features, ]
ncol(counts_subset)
#Downsample cells (optional, recommended for huge datasets)
set.seed(123)
n_cells <- 20000  # adjust based on available RAM
cells_subset <- sample(colnames(counts_subset), n_cells)


counts_subset <- counts_subset[, cells_subset]


# Assuming you have cluster labels and nUMI
# cluster <- your_cluster_vector[cells_subset]
# nUMI    <- colSums(counts_subset)

# counts_subset: sparse matrix of genes × cells
# cell_types: vector of cell type labels for the same cells
# nUMI: optional vector of UMI counts per cell
# Suppose you have cluster/annotation info somewhere
# It must have the same order as the columns of counts_subset

head(cortex@meta.data)
colnames(cortex@meta.data)
# Example: use clustering info
cortex$cell_type <- Idents(cortex)  # assigns current identities to 'cell_type'

# Example: if you have a Seurat object with cell type info
cell_types <- cortex$cell_type  # or whatever metadata column

# Subset to match your counts_subset cells
cell_types_subset <- cells_subset
cell_types_subset <- factor(cell_types_subset)

# Assign cell barcodes as names
names(cell_types_subset) <- colnames(counts_subset)

head(cell_types_subset)


reference <- Reference(
  counts     = counts_subset,
  cell_types = cell_types_subset,   # must match counts_subset columns
  nUMI       = nUMI[cells_subset], # optional
  n_max_cells = 20000               # downsample per type if needed
)

#Remove types with <25 cells:
# Keep only cell types with >=25 cells
valid_types <- names(which(table(reference@cell_types) >= 25))
ncol(counts_subset)
# Subset the reference
keep_cells <- names(reference@cell_types)[reference@cell_types %in% valid_types]
counts_filtered <- reference@counts[, keep_cells]
cell_types_filtered <- reference@cell_types[keep_cells]
nUMI_filtered <- reference@nUMI[keep_cells]

# Rebuild reference
reference <- Reference(
  counts     = counts_filtered,
  cell_types = cell_types_filtered,
  nUMI       = nUMI_filtered,
  n_max_cells = 10000
)

#1. Use broader cell type labels
cluster <- as.factor(ref$class_label)
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

#2. Merge rare clusters into “Other”
tbl <- table(cluster)
rare_types <- names(tbl[tbl < 25])
cluster[cluster %in% rare_types] <- "Other"
cluster <- droplevels(cluster)

#3. Increase n_max_cells
reference <- Reference(
  counts     = counts_subset,
  cell_types = cluster[cells_subset],
  nUMI       = nUMI[cells_subset],
  n_max_cells = 20000   # instead of 10000
)
valid_types <- names(which(table(reference@cell_types) >= 25))
colnames(ref@meta.data)

#1. Restrict genes before calling create.RCTD

#Since you can’t pass gene_list, you need to subset your objects directly:
gene_list <- rownames(reference@counts)[1:5000]

# Subset reference
reference_small <- reference
reference_small@counts <- reference@counts[gene_list, ]

# Subset query
query_small <- query
query_small@counts <- query@counts[gene_list, ]

RCTD <- create.RCTD(
  spatialRNA = query_small,
  reference  = reference_small,
  max_cores  = 28
)

cortex <- NormalizeData(cortex, assay = "Spatial.008um")
cortex <- ScaleData(cortex, assay = "Spatial.008um")

SpatialDimPlot(
  cortex,
  group.by = "first_type",
  cols = c("Glutamatergic" = "#FFFF00",
           "GABAergic"     = "#00BFFF",
           "Non-Neuronal"  = "#7FFF00",
           "Unknown"       = "grey80")
)


#Optional: process in batches if still too large
# Split your 200k cells into 10k-cell chunks, create reference for each, then merge
chunks <- split(colnames(counts_subset), ceiling(seq_along(colnames(counts_subset))/10000))
reference_list <- lapply(chunks, function(cells) {
  Reference(counts = counts_subset[, cells],
            cluster = cluster[cells],
            nUMI = nUMI[cells])
})
# Merge references if your package supports it

