# Load libraries
library(hdf5r)
library(Seurat)

# -----------------------------
# Step 1: Combine the two H5 files
# -----------------------------
h5_8um <- "path/to/Spatial.008um.h5"
h5_16um <- "path/to/Spatial.016um.h5"
combined_h5 <- "path/to/Spatial_combined.h5"

# Open the original H5 files
h5_8 <- H5File$new(h5_8um, mode = "r")
h5_16 <- H5File$new(h5_16um, mode = "r")

# Create a new H5 file
h5_combined_file <- H5File$new(combined_h5, mode = "w")

# Copy everything except binned matrices from 8um
for (name in h5_8$ls()$name) {
  if (name != "binned_matrices") {
    h5_combined_file$copy(h5_8[[name]], name)
  }
}

# Create group for binned matrices
grp <- h5_combined_file$create_group("binned_matrices")

# Copy 8um and 16um matrices
grp$copy(h5_8[["binned_matrices"]][["8um"]], "8um")
grp$copy(h5_16[["binned_matrices"]][["16um"]], "16um")

# Close files
h5_8$close_all()
h5_16$close_all()
h5_combined_file$close_all()

cat("Combined H5 file created at:", combined_h5, "\n")

# -----------------------------
# Step 2: Load combined file in Seurat
# -----------------------------
object <- Load10X_Spatial(
  data.dir = dirname(combined_h5),   # directory containing the combined h5
  filename = basename(combined_h5),  # combined h5 file
  bin.size = c(8, 16)
)

# -----------------------------
# Step 3: Check assays
# -----------------------------
Assays(object)
