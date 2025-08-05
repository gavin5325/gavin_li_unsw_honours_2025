library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)


add_tags_to_cell_id <- function(seurat_obj_name) {
  if (seurat_obj_name == 'KO52_D20') {
    return(c("D20.org.KO52", "D20KO52", "KO_52_D20_B7"))
  } else if (seurat_obj_name == 'KO_D60') {
    return(c("D60.org.KO58", "D60KO58", "sample-S2_KO_Zmynd8"))
  } else if (seurat_obj_name == 'KO_D100_B4') {
    return(c("D100.org.KO58", "D100KO58", "sample-ZMYND8_KO_D102"))
  } else if (seurat_obj_name == 'KO52_D100_B8') {
    return(c("D100.B8.org.KO52", "D100B8KO52", "KO52_D100_B8"))
  } else if (seurat_obj_name == 'KO58_D100_B8') {
    return(c("D100.B8.org.KO58", "D100B8KO58", "KO58_D100_B8"))
  } else if (seurat_obj_name == 'WT_D20') {
    return(c("D20.org.WT", "D20WT", "WT_D20_B7"))
  } else if (seurat_obj_name == 'WT_D60') {
    return(c("D60.org.WT", "D60WT", "sample-S1_WT"))
  } else if (seurat_obj_name == 'WT_D100_B4') {
    return(c("D100.org.WT", "D100WT", "sample-ZMYND8_WT_D102"))
  } else {
    return(c("D100.B8.org.WT", "D100B8WT", "WT_D100_B8"))
  }
}

filter_MT_content <- function(seurat_obj, seurat_obj_merged) {
  test_seurat <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/KO52_D20.rda"))
  julis_mt_list <- grep("^MT-", rownames(seurat_obj_merged), value = TRUE)
  my_mt_list <- grep("^MT-", rownames(seurat_obj), value = TRUE) 
  mt_diff_list <- setdiff(my_mt_list, julis_mt_list)
  seurat_obj_mt_subsetted <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), mt_diff_list))
  
  spliced_counts <- GetAssayData(seurat_obj, assay = "spliced", slot = "counts")
  unspliced_counts <- GetAssayData(seurat_obj, assay = "unspliced", slot = "counts")
  ambiguous_counts <- GetAssayData(seurat_obj, assay = "ambiguous", slot = "counts")
  total_counts <- spliced_counts + unspliced_counts + ambiguous_counts
  # seurat_obj[["RNA"]] <- CreateAssayObject(counts = total_counts)
  seurat_obj[["RNA"]] <- CreateAssayObject(counts = spliced_counts)
  DefaultAssay(seurat_obj) <- "RNA"
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Pre-process the seurat object
  seurat_obj_plus_min_pc <- Standard_Workflow(seurat_obj)
  seurat_obj <- seurat_obj_plus_min_pc$obj
  
  a <- VlnPlot(seurat_obj, features = c("percent.mt"))
  print(a)
  
  return(seurat_obj)
}


keep_julis_cells <- function(seurat_obj_path) {
  
  # seurat_obj_path <- "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/KO52_D20.rda"
  # seurat_obj_merged_path <- load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/Juli_seurat/org.Merge_Preintegration.clustered.Integrated.rda")
  # seurat_obj_merged <- get(seurat_obj_merged_path)
  
  seurat_obj_static <- load(seurat_obj_path)
  seurat_obj <- get(seurat_obj_static)

  
  spliced_counts <- GetAssayData(seurat_obj[["spliced"]], layer = "counts")
  spliced_assay <- CreateAssayObject(counts = spliced_counts)
  seurat_obj[["spliced"]] <- spliced_assay
  seurat_obj[["RNA"]] <- seurat_obj[["spliced"]]
  
  DefaultAssay(seurat_obj) <- "spliced"
  seurat_obj_name <- seurat_obj_static[1]
  
  tags <- add_tags_to_cell_id(seurat_obj_name)
  
  seurat_tag_prefix <- tags[1]
  seurat_tag_suffix <- paste0("-", tags[2])
  seurat_tag_recorded_name <- tags[3]
  
  seurat_obj_name_colon <- paste0(seurat_tag_recorded_name, ":")
  seurat_obj <- RenameCells(
    seurat_obj,
    new.names = gsub("x", seurat_tag_suffix,
                     gsub(seurat_obj_name_colon, "", colnames(seurat_obj)))
  )
  
  cell_ids <- colnames(seurat_obj)
  prefixed_ids <- paste0(seurat_tag_prefix, "_", cell_ids)
  merged_ids <- colnames(seurat_obj_merged)
  cells_to_keep <- cell_ids[prefixed_ids %in% merged_ids]
  seurat_obj_filtered <- subset(seurat_obj, cells = cells_to_keep)
  seurat_obj_filtered <- filter_MT_content(seurat_obj_filtered, seurat_obj_merged)
  seurat_obj_file_path <- paste0("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_subsetting/", seurat_obj_name, "_cells_filtered.h5Seurat")
  save(seurat_obj_filtered, file = seurat_obj_file_path)
  print(seurat_obj_filtered)
  # return(seurat_obj_filtered)
}

source("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Scripts/250520_seurat_preprocessing.R")

seurat_obj_merged_path <- load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/Juli_seurat/org.Merge_Preintegration.clustered.Integrated.rda")
seurat_obj_merged <- get(seurat_obj_merged_path)

data_directory <- "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/"
seurat_objects_individual <- c("KO52_D20.rda", "KO_D60.rda", "KO_D100_B4.rda", "KO52_D100_B8.rda",
                               "KO58_D100_B8.rda", "WT_D20.rda", "WT_D60.rda", "WT_D100_B4.rda", 
                               "WT_D100_B8.rda")

for (object in seurat_objects_individual) {
  set.seed(521)
  seurat_obj_path <- paste0(data_directory, object)
  keep_julis_cells(seurat_obj_path)
  # Only discrepancy is KO52_D20:(Information about the seurat obj is below)
  # An object of class Seurat 
  # 182025 features across 6247 samples within 3 assays 
  # Active assay: spliced (60675 features, 0 variable features)
  # 1 layer present: counts
  # 2 other assays present: unspliced, ambiguous
}
