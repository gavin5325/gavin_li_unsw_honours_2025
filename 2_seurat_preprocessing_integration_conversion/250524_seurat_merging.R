library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(Matrix)
source("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Scripts/250520_seurat_preprocessing.R")

data_directory <- "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_subsetting/"
seurat_objects_individual <- c("KO52_D20", "KO_D60", "KO_D100_B4", "KO52_D100_B8",
                               "KO58_D100_B8", "WT_D20", "WT_D60", "WT_D100_B4", 
                               "WT_D100_B8")

for (object in seurat_objects_individual) {
  full_path <- paste0(data_directory, object, "_cells_filtered.h5Seurat")
  assign(object, get(load(full_path)))
}

seurat_obj_list <- list(WT_D20, KO52_D20, WT_D60, KO_D60, WT_D100_B4, KO_D100_B4, WT_D100_B8, KO52_D100_B8, KO58_D100_B8)
org_merge <- merge(seurat_obj_list[[1]], 
                   y = c(seurat_obj_list[[2]], seurat_obj_list[[3]], seurat_obj_list[[4]], seurat_obj_list[[5]],
                         seurat_obj_list[[6]], seurat_obj_list[[7]], seurat_obj_list[[8]], seurat_obj_list[[9]]),
                   add.cell.ids = c("D20.org.WT", "D20.org.KO52", "D60.org.WT", 
                                    "D60.org.KO58", "D100.org.WT", "D100.org.KO58",
                                    "D100.B8.org.WT", "D100.B8.org.KO52", "D100.B8.org.KO58"),
                   project = 'org.Merge')

org_merge$Sample_ID <- as.character(gsub("_.*$","",rownames(org_merge[[]])))

save(org_merge, file="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_merging/org.Merge_Preintegration.rda")

# Observation
load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_merging/org.Merge_Preintegration.rda")
org_merge
# An object of class Seurat 
# 242700 features across 94790 samples within 4 assays 
# Active assay: RNA (60675 features, 0 variable features)
# 2 layers present: counts, data
# 3 other assays present: spliced, unspliced, ambiguous


########## I need to extract spliced/unspliced count matrices for both WT and KO samples.
## Not written as a function as I only need the matrices and won't go back on it later

WT <- list(WT_D20, WT_D60, WT_D100_B4, WT_D100_B8)
KO <- list(KO52_D20, KO_D60, KO_D100_B4, KO52_D100_B8, KO58_D100_B8)

for (i in seq_along(WT)) {
  WT[[i]] <- WT[[i]][rownames(julis_seurat), ]
}

for (i in seq_along(KO)) {
  KO[[i]] <- KO[[i]][rownames(julis_seurat), ]
}



spliced_matrices <- lapply(seq_along(KO), function(i) {
  # Change the following accordingly
  # obj <- WT[[i]]
  obj <- KO[[i]]
  
  # Change the following accordingly
  # mat <- GetAssayData(obj, assay = "spliced", slot = "counts")
  mat <- GetAssayData(obj, assay = "unspliced", slot = "counts")
  
  # Make cell names unique (optional but recommended)
  colnames(mat) <- paste0("Sample", i, "_", colnames(mat))
  
  return(mat)
})

# Union of all genes (rows)
all_genes <- unique(unlist(lapply(spliced_matrices, rownames)))

# Function to align matrix to union gene set
align_genes <- function(mat, gene_set) {
  missing_genes <- setdiff(gene_set, rownames(mat))
  if (length(missing_genes) > 0) {
    # Add missing rows as 0
    zero_mat <- Matrix(0, nrow = length(missing_genes), ncol = ncol(mat), sparse = TRUE)
    rownames(zero_mat) <- missing_genes
    colnames(zero_mat) <- colnames(mat)
    mat <- rbind(mat, zero_mat)
  }
  # Reorder rows to match the unified gene set
  mat <- mat[gene_set, , drop = FALSE]
  return(mat)
}

aligned_matrices <- lapply(spliced_matrices, align_genes, gene_set = all_genes)
combined_spliced <- do.call(cbind, aligned_matrices)
dim(combined_spliced)

# Comment out which ever ones you don't need
# writeMM(combined_spliced, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/WT_spliced_Juli.mtx")
# writeMM(combined_spliced, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/WT_unspliced_Juli.mtx")
# writeMM(combined_spliced, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/KO_spliced_Juli.mtx")
writeMM(combined_spliced, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/KO_unspliced_Juli.mtx")

write.csv(rownames(combined_spliced),
          file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/KO_genes_Juli.csv",
          row.names = FALSE, quote = FALSE)

# Save cell names (columns)
write.csv(colnames(combined_spliced),
          file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/KO_cells_Juli.csv",
          row.names = FALSE, quote = FALSE)

# 39208 WT cells; 55582 KO cells

# Just testing
matrix_WT_spliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/WT_spliced.mtx")
matrix_WT_unspliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/WT_unspliced.mtx")
matrix_KO_spliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/KO_spliced.mtx")
matrix_KO_unspliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/KO_unspliced.mtx")

WT_genes <- read.csv("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/WT_genes.csv", header = FALSE, stringsAsFactors = FALSE)[,1]
WT_cells <- read.csv("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/WT_cells.csv", header = FALSE, stringsAsFactors = FALSE)[,1]
KO_genes <- read.csv("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/KO_genes.csv", header = FALSE, stringsAsFactors = FALSE)[,1]
KO_cells <- read.csv("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/KO_cells.csv", header = FALSE, stringsAsFactors = FALSE)[,1]

WT_genes <- WT_genes[-1]
WT_cells <- WT_cells[-1]
KO_genes <- KO_genes[-1]
KO_cells <- KO_cells[-1]

rownames(matrix_WT_spliced) <- WT_genes
colnames(matrix_WT_spliced) <- WT_cells
rownames(matrix_WT_unspliced) <- WT_genes
colnames(matrix_WT_unspliced) <- WT_cells

rownames(matrix_KO_spliced) <- KO_genes
colnames(matrix_KO_spliced) <- KO_cells
rownames(matrix_KO_unspliced) <- KO_genes
colnames(matrix_KO_unspliced) <- KO_cells

##########################################
##########################################
##########################################

# ----- STEP 1: Build spliced and unspliced matrices for KO and WT -----

get_aligned_matrices <- function(samples, assay_type = "spliced") {
  raw_matrices <- lapply(seq_along(samples), function(i) {
    obj <- samples[[i]]
    mat <- GetAssayData(obj, assay = assay_type, slot = "counts")
    colnames(mat) <- paste0("Sample", i, "_", colnames(mat))  # Unique colnames
    return(mat)
  })
  
  # Union of all genes across samples
  all_genes <- unique(unlist(lapply(raw_matrices, rownames)))
  
  align_genes <- function(mat, gene_set) {
    missing_genes <- setdiff(gene_set, rownames(mat))
    if (length(missing_genes) > 0) {
      zero_mat <- Matrix(0, nrow = length(missing_genes), ncol = ncol(mat), sparse = TRUE)
      rownames(zero_mat) <- missing_genes
      colnames(zero_mat) <- colnames(mat)
      mat <- rbind(mat, zero_mat)
    }
    mat <- mat[gene_set, , drop = FALSE]
    return(mat)
  }
  
  aligned <- lapply(raw_matrices, align_genes, gene_set = all_genes)
  combined <- do.call(cbind, aligned)
  return(combined)
}

# ----- STEP 2: Get combined spliced and unspliced matrices -----

combined_spliced <- cbind(
  get_aligned_matrices(WT, assay_type = "spliced"),
  get_aligned_matrices(KO, assay_type = "spliced")
)

combined_unspliced <- cbind(
  get_aligned_matrices(WT, assay_type = "unspliced"),
  get_aligned_matrices(KO, assay_type = "unspliced")
)

# Check dimensions
dim(combined_spliced)
dim(combined_unspliced)

write.csv(rownames(combined_spliced),
          file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/genes.csv",
          row.names = FALSE, quote = FALSE)

# Save cell names (columns)
write.csv(colnames(combined_spliced),
          file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/cells.csv",
          row.names = FALSE, quote = FALSE)

# ----- STEP 3: Write to Matrix Market format (.mtx) -----

writeMM(combined_spliced, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/spliced.mtx")
writeMM(combined_unspliced, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/unspliced.mtx")
