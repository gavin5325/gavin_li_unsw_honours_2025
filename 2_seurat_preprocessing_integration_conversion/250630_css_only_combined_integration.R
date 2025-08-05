library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(reticulate)
library(Matrix)
library(simspec)
use_virtualenv("~/scvelo_env", required = TRUE)
py_config()
scv <- import("scvelo")
sc <- import("scanpy")
ad <- import("anndata")
np <- import("numpy")
pd <- import("pandas")
scipy_sparse <- import("scipy.sparse")
source("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Scripts/250520_seurat_preprocessing.R")

org_merge_WT <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/Seurat.objects.PostIntegration_CSS_only_combined_WT.h5Seurat")

output <- Standard_Workflow(org_merge_WT)
org_merge_WT <- output$obj
final_pc <- output$min_pcs

org_merge_WT <- cluster_sim_spectrum(object = org_merge_WT, label_tag = "Sample_ID", dims_use = 1:final_pc)
org_merge_WT <- RunUMAP(org_merge_WT, reduction = "css", dims = 1:final_pc, reduction.name = "umap.css.sigPCs")
org_merge_WT <- FindNeighbors(org_merge_WT, reduction = "css", dims = 1:final_pc)
org_merge_WT <- FindClusters(org_merge_WT, resolution = 0.6, cluster.name = "css_cluster")
org_merge_WT <- RunTSNE(org_merge_WT, reduction = "css", dims = 1:final_pc, reduction.name = "tsne.css.sigPCs", check_duplicates = FALSE)

SaveH5Seurat(org_merge_WT, filename = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined_WT.h5Seurat", overwrite = TRUE)
Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined_WT.h5Seurat", dest = "h5ad", overwrite = TRUE)

############## KO

org_merge_KO <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/Seurat.objects.PostIntegration_CSS_only_combined_KO.h5Seurat")

output <- Standard_Workflow(org_merge_KO)
org_merge_KO <- output$obj
final_pc <- output$min_pcs

org_merge_KO <- cluster_sim_spectrum(object = org_merge_KO, label_tag = "Sample_ID", dims_use = 1:final_pc)
org_merge_KO <- RunUMAP(org_merge_KO, reduction = "css", dims = 1:final_pc, reduction.name = "umap.css.sigPCs")
org_merge_KO <- FindNeighbors(org_merge_KO, reduction = "css", dims = 1:final_pc)
org_merge_KO <- FindClusters(org_merge_KO, resolution = 0.6, cluster.name = "css_cluster")
org_merge_KO <- RunTSNE(org_merge_KO, reduction = "css", dims = 1:final_pc, reduction.name = "tsne.css.sigPCs", check_duplicates = FALSE)

SaveH5Seurat(org_merge_KO, filename = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined_KO.h5Seurat", overwrite = TRUE)
Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined_KO.h5Seurat", dest = "h5ad", overwrite = TRUE)


############ PLOTTING #########
palette <- c(
  "#E41A1C",  # red
  "#377EB8",  # blue
  "#4DAF4A",  # green
  "#00008B",  # purple
  "#FF7F00",  # orange
  "#A65628",  # brown
  "#FFD92F",
  "#999999",  # grey
  "#66C2A5",  # teal
  "#984EA3",
  "#E78AC3",  # rose
  "#A6D854",  # lime green
  "#B3B3B3",  # light grey
  "#8DD3C7",  # aqua
  "#FB8072",  # coral
  "#80B1D3"   # sky blue
)
cell_types <- sort(unique(org_merge_WT$HNOCA_annot_level_1_pruned))
names(palette) <- cell_types


pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/WT_cells_development.pdf", width = 30, height = 7)
DimPlot(org_merge_WT, reduction = "umap.css.sigPCs", group.by = "HNOCA_annot_level_1_pruned", split.by = "Sample_ID", cols = palette[cell_types])
dev.off()
