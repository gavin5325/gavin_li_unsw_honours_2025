library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(parallel)
library(SeuratDisk)
library(SeuratWrappers)
library(hdf5r)
library(sceasy)
library(reticulate)
use_virtualenv("~/scvelo_env", required = TRUE)
py_config()
loompy <- reticulate::import('loompy')
scv <- import("scvelo")
sc <- import("scanpy")

### This function performs standard pre-processing steps on Seurat objects.
Pre_Processing_Filter <- function(seurat_obj) {

  # Generate nFeature_RNA meta data column
  spliced_mat <- GetAssayData(seurat_obj, assay = "spliced", slot = "counts")
  unspliced_mat <- GetAssayData(seurat_obj, assay = "unspliced", slot = "counts")
  ambiguous_mat <- GetAssayData(seurat_obj, assay = "ambiguous", slot = "counts")
  combined_mat <- (spliced_mat > 0) | (unspliced_mat > 0) | (ambiguous_mat > 0)
  nFeature_RNA <- Matrix::colSums(combined_mat)
  seurat_obj$nFeature_RNA <- nFeature_RNA
  
  # Generate nCount_RNA meta data column
  spliced_counts <- GetAssayData(seurat_obj, assay = "spliced", slot = "counts")
  unspliced_counts <- GetAssayData(seurat_obj, assay = "unspliced", slot = "counts")
  ambiguous_counts <- GetAssayData(seurat_obj, assay = "ambiguous", slot = "counts")
  nCount_RNA <- Matrix::colSums(spliced_counts) +
    Matrix::colSums(unspliced_counts) +
    Matrix::colSums(ambiguous_counts)
  seurat_obj$nCount_RNA <- nCount_RNA
  
  # Define the % of mitochondrial, ribosomal small and ribosomal large genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rpl"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPL")
  seurat_obj[["percent.rps"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS")
  
  P1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rps", "percent.rpl"), ncol=5)
  
  # Get rid of cells with abnormally low/high RNA count (99.8 percentile)
  sum(seurat_obj$nCount_RNA)/length(seurat_obj$nCount_RNA) # 4679.9
  nCount_RNA_99.8pct <- quantile(seurat_obj$nCount_RNA, 0.998)
  
  # Subset seurat objects to filter out cells
  seurat_obj <- subset(seurat_obj, subset = nCount_RNA > 1000 & nCount_RNA < nCount_RNA_99.8pct & nFeature_RNA < 25000 & percent.mt <10 & percent.rps < 20 & percent.rpl < 20)
  return(seurat_obj)
}

# This function uses DoubletFinder to identify and remove doublets from the Seurat object
Remove_Doublet <- function(seurat_obj, num_pcs) {
  
  # Perform sweep on 10000 cells
  pn_pk_sweep <- paramSweep(seurat_obj, PCs = 1:num_pcs, sct = FALSE, num.cores = 1)
  pn_pk_sweep_stats <- summarizeSweep(pn_pk_sweep, GT = FALSE)
  
  # Maximise BCmvn coefficient
  bcmvn <- find.pK(pn_pk_sweep_stats)
  bcmvn_max <- bcmvn[which.max(bcmvn$BCmetric), ]
  optimal_pK <- bcmvn_max$pK
  optimal_pK <- as.numeric(levels(optimal_pK))[optimal_pK]
  
  # Find homotypic doublets
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic_prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.076*nrow(seurat_obj@meta.data))
  nExp_poi_adj <- round(nExp_poi*(1-homotypic_prop))
  
  # DoubletFinder
  seurat_obj <- doubletFinder(seurat_obj, 
                              PCs = 1:num_pcs, 
                              pN = 0.25, 
                              pK = optimal_pK, 
                              nExp = nExp_poi_adj, 
                              reuse.pANN = NULL, 
                              sct = FALSE)

  # Remove doublets
  df_class_col <- grep("^DF.classifications_0.25_", colnames(seurat_obj@meta.data), value = TRUE)
  DimPlot(seurat_obj, reduction = 'umap', group.by = df_class_col)
  singlet_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[df_class_col]] == "Singlet"]
  seurat_obj_singlets <- subset(seurat_obj, cells = singlet_cells)
  
  
  return(seurat_obj_singlets)
  
}

Standard_Workflow <- function(seurat_obj) {

  # Standard Seurat workflow filtering
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  
  # Identify the number of PCs to keep using standard deviation from PCA
  st_dev <- seurat_obj[["pca"]]@stdev
  st_dev_sum <- sum(seurat_obj[["pca"]]@stdev)
  st_dev_percent <- (st_dev / st_dev_sum) * 100
  st_dev_cum_sum <- cumsum(st_dev_percent)
  component_candidate1 <- which(st_dev_cum_sum > 90 & st_dev_percent < 5)[1]
  component_candidate2 <- sort(which((st_dev_percent[1:length(st_dev_percent) - 1] - st_dev_percent[2:length(st_dev_percent)]) > 0.1), decreasing = T)[1] + 1
  final_pc <- min(component_candidate1, component_candidate2)
  
  # Filtering continued using number of PCs to keep
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:final_pc, reduction = "pca", reduction.name = "umap.unintegrated")
  seurat_obj <- FindNeighbors(object = seurat_obj, dims = 1:final_pc, reduction = "pca")
  seurat_obj <- FindClusters(object = seurat_obj, resolution = 0.6, cluster.name = "unintegrated_clusters")
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  return(list(obj = seurat_obj, min_pcs = final_pc))
}

Pre_Processing <- function(seurat_obj_path) {

  # Set seed for reproducibility
  set.seed(521)
  
  # Load the seurat object within the path
  seurat_obj_static <- load(seurat_obj_path)
  seurat_obj_name <- seurat_obj_static[1]
  seurat_obj <- get(seurat_obj_static)
  
  # seurat_obj[["RNA"]] <- seurat_obj[["spliced"]]
  # assay_v3_spliced <- CreateAssayObject(counts = seurat_obj[["spliced"]]$counts)
  # seurat_obj[["spliced"]] <- assay_v3_spliced
  # assay_v3 <- CreateAssayObject(counts = seurat_obj[["RNA"]]$counts)
  # seurat_obj[["RNA3"]] <- assay_v3
  # DefaultAssay(seurat_obj) <- "RNA3"
  # seurat_obj[["RNA"]] <- NULL
  # seurat_obj <- RenameAssays(seurat_obj, RNA3 = 'RNA')
  # DefaultAssay(seurat_obj) <- "spliced"
  # 
  # seurat_obj <- Pre_Processing_Filter(seurat_obj)
  # 
  # # Create a RNA layer necessary for DoubletFinder
  # spliced_counts <- GetAssayData(seurat_obj, assay = "spliced", slot = "counts")
  # seurat_obj[["RNA"]] <- CreateAssayObject(counts = spliced_counts)
  
  spliced_counts <- GetAssayData(seurat_obj[["spliced"]], layer = "counts")
  spliced_assay <- CreateAssayObject(counts = spliced_counts)
  seurat_obj[["spliced"]] <- spliced_assay
  seurat_obj[["RNA"]] <- seurat_obj[["spliced"]]
  
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  seurat_obj <- Pre_Processing_Filter(seurat_obj)
  
  # Carry out standard workflow  
  seurat_obj_plus_min_pc <- Standard_Workflow(seurat_obj)
  seurat_obj <- seurat_obj_plus_min_pc$obj
  final_pc <- seurat_obj_plus_min_pc$min_pcs
  
  # Doublet removed
  seurat_obj <- Remove_Doublet(seurat_obj, final_pc)
  
  # Perform the standard workflow again for the singlet dataset
  seurat_obj_plus_min_pc <- Standard_Workflow(seurat_obj)
  seurat_obj <- seurat_obj_plus_min_pc$obj
  final_pc <- seurat_obj_plus_min_pc$min_pcs
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Below is necessary for conversion to h5ad
  spliced_counts <- GetAssayData(object = seurat_obj, assay = "spliced", slot = "counts")
  unspliced_counts <- GetAssayData(object = seurat_obj, assay = "unspliced", slot = "counts")
  ambiguous_counts <- GetAssayData(object = seurat_obj, assay = "ambiguous", slot = "counts")
  save_path <- paste0("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250520_seurat_preprocessing/",
                      seurat_obj_name, "_preprocessed", ".h5Seurat")
  SaveH5Seurat(seurat_obj, filename = save_path, overwrite = TRUE)
  
  # save_path_h5ad <- paste0("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250520_seurat_preprocessing/",
  #                          seurat_obj_name, ".h5ad")
  # sceasy::convertFormat(seurat_obj, from="seurat", to="anndata", outFile=save_path_h5ad)
  # 
  # # Load adata
  # adata <- sc$read(save_path_h5ad)
  # adata$layers['spliced'] <- t(spliced_counts)
  # adata$layers['unspliced'] <- t(unspliced_counts)
  # adata$layers['ambiguous'] <- t(ambiguous_counts)
  
  # return(adata)
}
