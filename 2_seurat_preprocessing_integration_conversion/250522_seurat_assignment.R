library(Seurat)
library(SeuratDisk)
source("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Scripts/250520_seurat_preprocessing.R")

data_directory <- "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/"

# This part of the script performs the pre-processing functions in 
# /srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Scripts/250520_seurat_preprocessing.R
seurat_objects_individual <- c("KO52_D20.rda", "KO_D60.rda", "KO_D100_B4.rda", "KO52_D100_B8.rda",
                               "KO58_D100_B8.rda", "WT_D20.rda", "WT_D60.rda", "WT_D100_B4.rda", 
                               "WT_D100_B8.rda")


for (object in seurat_objects_individual) {
  combined_path <- paste0(data_directory, object)
  Pre_Processing(combined_path)
}

results_directory <- "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250520_seurat_preprocessing/"

# Load the pre-processed seurat objects
for (object in seurat_objects_individual) {
  object_without_ext <- sub("\\.rda$", "", object)
  combined_path <- paste0(results_directory, object_without_ext, '_preprocessed', '.h5Seurat')
  seurat_obj <- LoadH5Seurat(combined_path)
  assign(object_without_ext, seurat_obj)
}

# Inspect each seurat object
WT_D20
# An object of class Seurat 
# 182025 features across 50791 samples within 3 assays 
# Active assay: spliced (60675 features, 0 variable features)
# 1 layer present: counts
# 2 other assays present: unspliced, ambiguous
WT_D60
WT_D100_B4
WT_D100_B8
KO52_D20
KO_D60
KO_D100_B4
KO52_D100_B8
KO58_D100_B8

