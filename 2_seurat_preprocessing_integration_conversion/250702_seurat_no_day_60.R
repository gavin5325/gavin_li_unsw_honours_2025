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

source("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Scripts/250520_seurat_preprocessing.R")
org_merge <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_merging/org.Merge_Preintegration.rda"))

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

org_merge_day20 <- subset(org_merge, subset = Day == "20")
org_merge_no_day60 <- subset(org_merge, subset = Day != "60")
org_merge_no_day60_WT <- subset(org_merge_no_day60, subset = Genotype == "WT")

save(org_merge, file="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_merging/org.Merge_Preintegration.rda")
save(org_merge_day20, file="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_day20.rda")
save(org_merge_no_day60, file="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60.rda")
save(org_merge_no_day60_WT, file="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT.rda")

lseurat_integration <- function(seurat_obj, path) {
  
  preprocessing_outputs <- Standard_Workflow(seurat_obj)
  seurat_obj <- preprocessing_outputs$obj
  final_pc <- preprocessing_outputs$min_pcs
  
  seurat_obj <- cluster_sim_spectrum(seurat_obj, label_tag = "Sample_ID", dims_use = 1:final_pc)
  
  seurat_obj <- RunUMAP(seurat_obj, reduction = "css", dims = 1:final_pc, reduction.name = "umap.css.sigPCs")
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "css", dims = 1:final_pc)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.6, cluster.name = "css_cluster")
  
  save(seurat_obj, file = path)
}

org_merge_day20 <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_day20.rda"))
org_merge_no_day60 <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60.rda"))
org_merge_no_day60_WT <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT.rda"))

seurat_objects <- list(
  org_merge_day20 = org_merge_day20,
  org_merge_no_day60 = org_merge_no_day60,
  org_merge_no_day60_WT = org_merge_no_day60_WT
)

seurat_paths <- list(
  org_merge_day20 = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_day20.rda",
  org_merge_no_day60 = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60.rda",
  org_merge_no_day60_WT = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT.rda"
)

mapply(lseurat_integration, seurat_objects, seurat_paths, SIMPLIFY = FALSE)

orgs_list <- list(org_merge_day20, org_merge_no_day60, org_merge_no_day60_WT)
names(orgs_list) <- c("org_merge_day20", "org_merge_no_day60", "org_merge_no_day60_WT")
output_dir <- "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/"
group_bys <- c("Genotype", "Sample_ID")

for (name in names(orgs_list)) {
  org <- orgs_list[[name]]

  
  pdf(paste0(output_dir, name, ".pdf"), width = 10, height = 7)
  for (group_by in group_bys) {
    p <- DimPlot(org, reduction = "umap.css.sigPCs", group.by = group_by)
    print(p)
  }
  print(DimPlot(org, reduction = "umap.css.sigPCs", group.by = "HNOCA_annot_level_1_pruned", split.by = "Sample_ID"))
  print(DimPlot(org, reduction = "umap.css.sigPCs", group.by = "HNOCA_annot_level_1_pruned", split.by = "Day"))
  
  dev.off()
}

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

SaveH5Seurat(org_merge_day20, filename = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_day20.h5Seurat", overwrite = TRUE)
SaveH5Seurat(org_merge_no_day60, filename = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60.h5Seurat", overwrite = TRUE)
SaveH5Seurat(org_merge_no_day60_WT, filename = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT.h5Seurat", overwrite = TRUE)

Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_day20.h5Seurat", dest = "h5ad", overwrite = TRUE)
Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60.h5Seurat", dest = "h5ad", overwrite = TRUE)
Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT.h5Seurat", dest = "h5ad", overwrite = TRUE)

###########

org_merge_no_day60_WT <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT.h5Seurat")
adata <- sc$read("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT.h5ad")

spliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT_spliced.mtx")
unspliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT_unspliced.mtx")
ambiguous <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT_ambiguous.mtx")

genes <- readLines("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT_spliced_genes.txt")
cells <- readLines("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT_spliced_cells.txt")

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
adata_full$write("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_WT.h5ad")

#############
org_merge_day20 <- sc$read("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_day20.h5ad")
scv$pp$filter_genes_dispersion(org_merge_day20)
sc$pp$neighbors(org_merge_day20, n_neighbors=as.integer(20), use_rep='X_pca')
scv$pp$moments(org_merge_day20, n_pcs = 20, n_neighbors = 20)
scv$tl$recover_dynamics(org_merge_day20)
org_merge_day20$write("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_day20.h5ad")
scv$tl$velocity(org_merge_day20, mode='dynamical')
scv$tl$velocity_graph(org_merge_day20)

palette <- c(
  "#E41A1C",  # red
  "#377EB8",  # blue
  "#4DAF4A",  # green
  "#984EA3",  # purple
  "#FF7F00",  # orange
  "#A65628",  # brown
  "#F781BF",  # pink
  "#999999",  # grey
  "#66C2A5",  # teal
  "#FFD92F",  # yellow
  "#7F00FF",  # rose
  "#A6D854",  # lime green
  "#B3B3B3",  # light grey
  "#8DD3C7",  # aqua
  "#FB8072",  # coral
  "#80B1D3"   # sky blue
)

group_bys <- c("Genotype", "Sample_ID", "HNOCA_annot_level_1_pruned", "HNOCA_annot_level_2_pruned")

for (group_by in group_bys) {
  plt$figure(figsize = c(15, 7))  # Set figure size in inches
  
  scv$pl$velocity_embedding_grid(
    org_merge_day20,
    basis = 'X_umap.css.sigPCs',
    color = group_by,
    legend_loc = 'lower left',
    size = 50,
    palette = palette,
    alpha = 1.0,
    scale = 0.3,
    legend_fontsize=3,
    show = FALSE  # Don't show in notebook/console
  )
  
  plt$savefig(
    paste0("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/velocity_stream_org_merge_day20_", group_by, ".pdf"),
    dpi = 300 
  )
  
  plt$close()
}

###### MAKE THE PLOTS ########
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
pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250702_seurat_no_day_60/org_merge_no_day60_UMAP_HNOCA_annot_level_2_pruned.pdf", width = 10, height = 7)
p <- DimPlot(org_merge_no_day60, reduction = "umap.css.sigPCs", group.by = "HNOCA_annot_level_2_pruned", cols = r_palette)
print(p)
dev.off()
