# Set Python path early (optional if already in ~/.Renviron)
Sys.setenv(RETICULATE_PYTHON = "~/scvelo_env/bin/python")

# Load reticulate FIRST
library(reticulate)

# Attach Python environment before anything else
use_virtualenv("~/scvelo_env", required = TRUE)
py_config()  # optional, just for verification

# Then import Python packages
scv <- import("scvelo")
sc <- import("scanpy")
ad <- import("anndata")
np <- import("numpy")
pd <- import("pandas")
scipy_sparse <- import("scipy.sparse")
plt <- import("matplotlib.pyplot")

# Now safely load R packages that may depend on Python
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(simspec)
library(Matrix)

output_dir <- "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/"
orgs_list <- list(org_merge_KO_changed)
names(orgs_list) <- c("org_merge_KO_changed")
org_merge_WT <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/Seurat.objects.PostIntegration_CSS_only_combined_WT.h5Seurat")
org_merge_KO <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/Seurat.objects.PostIntegration_CSS_only_combined_KO.h5Seurat")


ldat_day60_KO <- ReadVelocity(file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/loom_files/D60_KO_hg19.loom")
rownames(ldat_day60_KO$spliced) <- make.unique(rownames(ldat_day60_KO$spliced))
rownames(ldat_day60_KO$unspliced) <- make.unique(rownames(ldat_day60_KO$unspliced))
rownames(ldat_day60_KO$ambiguous) <- make.unique(rownames(ldat_day60_KO$ambiguous))
D60_KO <- as.Seurat(x = ldat_day60_KO)

spliced_counts <- GetAssayData(D60_KO[["spliced"]], layer = "counts")
spliced_assay <- CreateAssayObject(counts = spliced_counts)
D60_KO[["spliced"]] <- spliced_assay
D60_KO[["RNA"]] <- D60_KO[["spliced"]]
DefaultAssay(D60_KO) <- "RNA"

old_names <- colnames(D60_KO)
new_names <- sub(
  pattern = "sample-S2_KO_Zmynd8:([ACGT]+)x$",
  replacement = "D60.org.KO58_\\1-D60KO58",
  x = old_names
)
D60_KO <- RenameCells(D60_KO, new.names = new_names)

shared_cells <- intersect(colnames(D60_KO), colnames(org_merge_KO))
D60_KO_subset <- subset(D60_KO, cells = shared_cells)

shared_genes <- intersect(rownames(D60_KO_subset[["RNA"]]), rownames(org_merge_KO[["RNA"]]))
D60_KO_subset <- subset(D60_KO_subset, features = shared_genes)
save(D60_KO_subset, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/KO_D60_hg19.rda") 


org_merge_KO_changed <- subset(org_merge_KO, subset = Day != 60)
org_merge_KO_changed <- merge(org_merge_KO_changed, y = c(D60_KO_subset))
org_merge_KO_changed$Day[is.na(org_merge_KO_changed$Day)] <- 60
org_merge_KO_changed$Sample_ID[is.na(org_merge_KO_changed$Sample_ID)] <- "D60.org.KO.hg19"


getMatrix <- function(seurat_obj, name) {
  assays<- c("spliced", "unspliced", "ambiguous")
  for (assay in assays) {
    mat <- GetAssayData(seurat_obj, assay = assay, slot = "counts")
    Matrix::writeMM(mat, file = paste0(output_dir, name, "_", assay, ".mtx"))
    writeLines(rownames(mat), paste0(output_dir, name, "_", assay, "_genes.txt"))
    writeLines(colnames(mat), paste0(output_dir, name, "_", assay, "_cells.txt"))
  }
}

Map(getMatrix, orgs_list, names(orgs_list))

# SaveH5Seurat(D60_KO_subset, filename="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/D60_KO.h5Seurat", overwrite = TRUE)

SaveH5Seurat(org_merge_KO_changed, filename="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_KO_changed.h5Seurat", overwrite = TRUE)
# Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_KO_changed.h5Seurat", dest="h5ad", overwrite = TRUE)

#######
adata <- sc$read("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed.h5ad")

spliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed_spliced.mtx")
unspliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed_unspliced.mtx")
ambiguous <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed_ambiguous.mtx")

genes <- readLines("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed_spliced_genes.txt")
cells <- readLines("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed_spliced_cells.txt")

rownames(spliced) <- genes
colnames(spliced) <- cells
rownames(unspliced) <- genes
colnames(unspliced) <- cells
rownames(ambiguous) <- genes
colnames(ambiguous) <- cells

obs_r <- py_to_r(adata$obs)
obs_copy <- r_to_py(obs_r) 

adata_full <- ad$AnnData(
  X=t(as.matrix(spliced)),    
  obs=obs_copy,
  var=data.frame(gene_names = rownames(spliced))
)

adata_full$uns = reticulate::r_to_py(adata$uns)$copy()
adata_full$obsm = reticulate::r_to_py(adata$obsm)$copy()
adata_full$obsp = reticulate::r_to_py(adata$obsp)$copy()

spliced_csr <- as(t(spliced), "dgCMatrix")
adata_full$layers["spliced"] <- scipy_sparse$csr_matrix(spliced_csr)
unspliced_csr <- as(t(unspliced), "dgCMatrix")
adata_full$layers["unspliced"] <- scipy_sparse$csr_matrix(unspliced_csr)
ambiguous_csr <- as(t(ambiguous), "dgCMatrix")
adata_full$layers["ambiguous"] <- scipy_sparse$csr_matrix(ambiguous_csr)
adata_full$X <- scipy_sparse$csr_matrix(spliced_csr)
adata_full$write("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed.h5ad")

##### Normalise

source("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Scripts/250520_seurat_preprocessing.R")

lseurat_integration <- function(seurat_obj, path) {
  
  preprocessing_outputs <- Standard_Workflow(seurat_obj)
  seurat_obj <- preprocessing_outputs$obj
  final_pc <- preprocessing_outputs$min_pcs
  seurat_obj@misc$final_pc <- final_pc
  print(final_pc)
  
  seurat_obj <- cluster_sim_spectrum(seurat_obj, label_tag = "Sample_ID", dims_use = 1:final_pc)

  seurat_obj <- RunUMAP(seurat_obj, reduction = "css", dims = 1:final_pc, reduction.name = "umap.css.sigPCs")
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "css", dims = 1:final_pc)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.6, cluster.name = "css_cluster")

  save(seurat_obj, file = path)
  return(seurat_obj)
}

org_merge_WT <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_WT_changed.h5Seurat")
org_merge_KO <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_KO_changed.h5Seurat")
org_merge <- merge(org_merge_WT, y = c(org_merge_KO))
org_merge@meta.data$Genotype <- ifelse(grepl("WT", org_merge@meta.data$Sample_ID), "WT", "KO")
org_merge@meta.data$Day <- sub("D(\\d+).*", "\\1", org_merge@meta.data$Sample_ID)
JW_s_annot <- read.csv("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250618_further_CSS_integration/HNOCA_annot_level1n2.csv")
annot<-JW_s_annot$cell_ID
GL_s_barcode <- rownames(org_merge@meta.data)
JW_s_barcode <- annot
table(is.na(match(GL_s_barcode,JW_s_barcode)))
annot_merge<-JW_s_annot[match(GL_s_barcode,JW_s_barcode),]
dim(annot_merge)
org_merge$HNOCA_annot_level_1_pruned<- annot_merge$HNOCA_annot_level_1_pruned
org_merge$HNOCA_annot_level_2_pruned <- annot_merge$HNOCA_annot_level_2_pruned
SaveH5Seurat(org_merge, filename="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed.h5Seurat", overwrite=TRUE)

org_merge_WT <- subset(org_merge, subset = Genotype == "WT")
org_merge_WT <- lseurat_integration(org_merge_WT, "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_WT_changed.h5Seurat")
SaveH5Seurat(org_merge_WT, filename="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_WT_changed.h5Seurat", overwrite=TRUE)
Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_WT_changed.h5Seurat", dest="h5ad", overwrite = TRUE)

org_merge_KO <- subset(org_merge, subset = Genotype == "KO")
org_merge_KO <- lseurat_integration(org_merge_KO, "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_KO_changed.h5Seurat")
SaveH5Seurat(org_merge_KO, filename="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_KO_changed.h5Seurat", overwrite=TRUE)
Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_KO_changed.h5Seurat", dest="h5ad", overwrite = TRUE)

#######
org_merge_WT <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_WT_changed.h5Seurat")
org_merge_KO <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_KO_changed.h5Seurat")


r_palette <- c(
  "Dorsal Telencephalic IP" = "#E41A1C",
  "Dorsal Telencephalic NPC" = "#377EB8",
  "MC" = "#4DAF4A",
  "NA" = "#984EA3",
  "NC Derivatives" = "#FF7F00",
  "Neuroepithelium" = "#A65628",
  "Non-telencephalic NPC" = "#F781BF",
  "Non-telencephalic Neuron" = "#999999",
  "PSC" = "#66C2A5",
  "Ventral Telencephalic NPC" = "#FFD92F",
  "Ventral Telencephalic Neuron" = "#7F00FF",
  "Astrocyte" = "#4B0082",
  "CP" = "#056517",
  "Dorsal Telencephalic Neuron" = "#8DD3C7",
  "Glioblast" = "#FFA33F",
  "OPC" = "#80B1D3"
)

pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/UMAPs.pdf", width = 15, height = 10)
DimPlot(org_merge_WT, reduction = "umap.css.sigPCs", group.by = "HNOCA_annot_level_2_pruned", cols = r_palette)
DimPlot(org_merge_WT, reduction = "umap.css.sigPCs", group.by = "Day")
# DimPlot(org_merge_WT, reduction = "umap.css.sigPCs", group.by = "HNOCA_annot_level_2_pruned", split.by = "Day", cols = r_palette)
DimPlot(org_merge_KO, reduction = "umap.css.sigPCs", group.by = "HNOCA_annot_level_2_pruned", cols = r_palette)
DimPlot(org_merge_KO, reduction = "umap.css.sigPCs", group.by = "Day")
# DimPlot(org_merge_KO, reduction = "umap.css.sigPCs", group.by = "HNOCA_annot_level_2_pruned", split.by = "Day", cols = r_palette)
dev.off()

pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/UMAP_composition_day.pdf", width = 35, height = 10)
DimPlot(org_merge_WT, reduction = "umap.css.sigPCs", group.by = "HNOCA_annot_level_2_pruned", split.by = "Day", cols = r_palette)
DimPlot(org_merge_KO, reduction = "umap.css.sigPCs", group.by = "HNOCA_annot_level_2_pruned", split.by = "Day", cols = r_palette)
dev.off()

#### MERGE WT AND KO TOGETHER
org_merge_changed <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed.h5Seurat")
org_merge_changed <- lseurat_integration(org_merge_changed, "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed.h5Seurat")
orgs_list <- list(org_merge_changed)
names(orgs_list) <- c("org_merge_changed")
Map(getMatrix, orgs_list, names(orgs_list))
SaveH5Seurat(org_merge_changed, filename="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed.h5Seurat", overwrite=TRUE)
Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed.h5Seurat", dest="h5ad", overwrite = TRUE)



