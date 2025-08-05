library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(simspec)
library(harmony)

source("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Scripts/250520_seurat_preprocessing.R")
org_merge <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_merging/org.Merge_Preintegration.rda"))

#### TRY SUBSETTING
julis_seurat <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/julis_seurat_integrated.rda"))

diff_cells <- setdiff(colnames(julis_seurat), colnames(org_merge))
julis_seurat <- subset(julis_seurat, cells = setdiff(colnames(julis_seurat), diff_cells))
org_merge <- org_merge[rownames(julis_seurat), colnames(julis_seurat)]

common_genes <- intersect(rownames(julis_seurat), rownames(org_merge))
julis_seurat <- subset(julis_seurat, features = common_genes)

julis_seurat[["spliced"]] <- org_merge[["spliced"]]
julis_seurat[["unspliced"]] <- org_merge[["unspliced"]]
julis_seurat[["ambiguous"]] <- org_merge[["ambiguous"]]

org_merge <- julis_seurat
save(org_merge, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_CSS_only_combined.rda")


# julis_genes <- rownames(julis_seurat)
# common_genes <- intersect(julis_genes, rownames(org_merge))
# org_merge_subset <- org_merge[common_genes, ]
# org_merge <- org_merge_subset

# spliced_counts <- GetAssayData(org_merge[["spliced"]], slot = "counts")
# org_merge[["RNA"]] <- CreateAssayObject(counts = spliced_counts)
# DefaultAssay(org_merge) <- "RNA"

total_counts <- GetAssayData(org_merge[["spliced"]], slot = "counts") +
  GetAssayData(org_merge[["unspliced"]], slot = "counts") +
  GetAssayData(org_merge[["ambiguous"]], slot = "counts")

total_assay <- CreateAssayObject(counts = total_counts)
org_merge[["RNA"]] <- total_assay
DefaultAssay(org_merge) <- "RNA"

cell_ids <- colnames(org_merge)
batches <- sub(".*-", "", cell_ids)
org_merge$Batch <- batches

# days <- unique(org_merge$Day)
# hvg_list <- list()
# 
# for (day_i in days) {
#   org_merge <- subset(org_merge, subset = Day == day_i)
#   org_merge <- NormalizeData(org_merge)
#   org_merge <- FindVariableFeatures(org_merge, selection.method = "vst", nfeatures = 2000)
#   hvg_list[[day_i]] <- VariableFeatures(org_merge)
# }
# 
# # Take union of all HVGs across days
# union_hvgs <- unique(unlist(hvg_list))
# length(union_hvgs)
# 
# # Now set HVGs for the whole object to this union
# VariableFeatures(org_merge) <- union_hvgs

org_merge <- NormalizeData(org_merge)
org_merge <- FindVariableFeatures(org_merge)
org_merge <- ScaleData(org_merge)
org_merge <- RunPCA(org_merge)

# Identify the number of PCs to keep using standard deviation from PCA
st_dev <- org_merge[["pca"]]@stdev
st_dev_sum <- sum(org_merge[["pca"]]@stdev)
st_dev_percent <- (st_dev / st_dev_sum) * 100
st_dev_cum_sum <- cumsum(st_dev_percent)
component_candidate1 <- which(st_dev_cum_sum > 90 & st_dev_percent < 5)[1]
component_candidate2 <- sort(which((st_dev_percent[1:length(st_dev_percent) - 1] - st_dev_percent[2:length(st_dev_percent)]) > 0.1), decreasing = T)[1] + 1
final_pc <- min(component_candidate1, component_candidate2)

#####
#####

# final_pc <- 42
org_merge <- FindNeighbors(object = org_merge, dims = 1:final_pc, reduction = "pca")
org_merge <- FindClusters(object = org_merge, resolution = 0.6, cluster.name = "unintegrated_clusters")
# org_merge <- FindClusters(org_merge, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0))
org_merge <- RunUMAP(org_merge, dims = 1:final_pc, reduction = "pca", reduction.name = "umap.unintegrated")

save(org_merge, file="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/org.Merge_Preintegration.clustered.rda")
# save(org_merge, file="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/org.Merge_Preintegration.clustered_julis_genes.rda")

DimPlot(org_merge, reduction = "umap.unintegrated", group.by = "Sample_ID") +ggtitle("resolution = 0.2")
DimPlot(org_merge, reduction = "umap.unintegrated", group.by = "unintegrated_clusters")

pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Plot_1c.pdf", width = 30, height = 4)
DimPlot(org_merge, reduction = "umap.unintegrated", split.by = "Sample_ID")
dev.off()

org_merge <- cluster_sim_spectrum(org_merge, label_tag = "Sample_ID", dims_use = 1:final_pc)
# org_merge <- cluster_sim_spectrum(object = org_merge, label_tag = "Sample_ID", dims_use = 1:final_pc)

# https://github.com/quadbio/simspec/issues/10
css_mat <- org_merge[["css"]]@cell.embeddings
num_na <- sum(is.na(css_mat))
print(num_na)
css_mat[is.na(css_mat)] <- 0
org_merge[["css"]]@cell.embeddings <- css_mat

org_merge <- FindNeighbors(org_merge, reduction = "css", dims = 1:final_pc)
org_merge <- FindClusters(org_merge, resolution = 0.6, cluster.name = "clusters_CSS_SigPs_0.6")
org_merge <- RunUMAP(org_merge, dims = 1:final_pc, reduction.name = "umap.css.sigPCs", reduction = "css")
save(org_merge, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_CSS_only.rda")


DimPlot(org_merge, group.by = "Sample_ID")
# org_merge <- RunUMAP(org_merge, reduction = "css", dims = 1:ncol(Embeddings(org_merge, "css")))
# org_merge <- FindNeighbors(org_merge, reduction = "css", dims = 1:ncol(Embeddings(org_merge, "css")))
# org_merge <- FindClusters(org_merge, resolution = 0.2)
# UMAPPlot(org_merge, group.by = "Sample_ID") + UMAPPlot(org_merge)

# Harmony
# org_merge <- RunHarmony(org_merge, group.by.vars = "Sample_ID", theta=2)
org_merge <- RunHarmony(org_merge, group.by.vars = "Sample_ID", theta=4)
dim(org_merge@reductions[["harmony"]]@cell.embeddings) # 94790 50
sum(org_merge@reductions[["harmony"]]@cell.embeddings) # 56471.63
min(org_merge@reductions[["harmony"]]@cell.embeddings) # -44.23942
max(org_merge@reductions[["harmony"]]@cell.embeddings) # 29.53911

Harmony_Seurat_pcs<-org_merge@reductions[["harmony"]]@cell.embeddings
write.csv(Harmony_Seurat_pcs, file="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Harmony_Seurat_pcs.csv", row.names = T, col.names = T)
dim(org_merge@reductions[["harmony"]]@feature.loadings) #[1] 2000   50 (2000 variable genes, 50 harmony dimensions)
set.seed(521)

harmony_embed <- harmony::HarmonyMatrix(
  data_mat = org_merge@reductions[["pca"]]@cell.embeddings,  
  meta_data = org_merge@meta.data,
  vars_use = 'Sample_ID',
  do_pca = FALSE, 
  theta = 4, 
)

identical(rownames(harmony_embed),Cells(org_merge))

Harmony_Original_pcs <- harmony_embed
dim(harmony_embed); sum(harmony_embed); min(harmony_embed); max(harmony_embed)
write.csv(Harmony_Original_pcs, file="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Harmony_Original_pcs.csv", row.names = T, col.names = T)
# write.csv(Harmony_Original_pcs, file="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Harmony_Original_pcs_julis_genes.csv", row.names = T, col.names = T)
org_merge@reductions[["harmony"]]@cell.embeddings<-harmony_embed

L5L8<-function(obj, reduction, reduction.name, resolution, cluster.name, dims) {
  obj <- FindNeighbors(obj, reduction = reduction, dims = dims) 
  # L6
  obj <- FindClusters(obj, resolution = resolution, cluster.name = cluster.name) # cluster.name is where you get to define the name of the cluster
  # L7
  obj <- RunUMAP(obj, reduction = reduction, dims = dims, reduction.name = reduction.name)
  return(obj)
}

Resolution.options <- c(0.1, 0.2, 0.5)
names(Resolution.options) <- c("Point1","Point2","Point5")

for(i in c(1:length(Resolution.options)))  {
  resolution<-Resolution.options[i]
  res.sufix<-names(Resolution.options)[i]
  
  print(paste0("Trying resolution ", resolution))
  
  org_merge<-L5L8(org_merge, reduction = "css", reduction.name = "umap.css.AllEmbed", resolution = resolution, cluster.name= paste0("clusters_CSS_allEmbeddingDims_",res.sufix), dims = 1:ncol(Embeddings(org_merge, "css")))
  org_merge<-L5L8(org_merge, reduction = "css", reduction.name = "umap.css.sigPCs", resolution = resolution, cluster.name= paste0("clusters_CSS_SigPs_",res.sufix), dims = 1:final_pc)
  # org_merge<-L5L8(org_merge, reduction = "harmony",reduction.name = "umap.harmony", resolution = resolution, cluster.name=paste0("clusters_Harmony_",res.sufix), dims = 1:final_pc)
}

org_merge[["spliced"]] <- as(org_merge[["spliced"]], Class="Assay")

org_merge@meta.data$Genotype <- ifelse(grepl("WT", org_merge@meta.data$Sample_ID), "WT", "KO")
org_merge@meta.data$Day <- sub("D(\\d+).*", "\\1", org_merge@meta.data$Sample_ID)
# Add Juli's cell annotation
JW_s_annot <- read.csv("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250618_further_CSS_integration/HNOCA_annot_level1n2.csv")
annot<-JW_s_annot$cell_ID
GL_s_barcode <- rownames(org_merge@meta.data)
JW_s_barcode <- annot
table(is.na(match(GL_s_barcode,JW_s_barcode)))
annot_merge<-JW_s_annot[match(GL_s_barcode,JW_s_barcode),]
dim(annot_merge)
org_merge$HNOCA_annot_level_1_pruned<- annot_merge$HNOCA_annot_level_1_pruned
org_merge$HNOCA_annot_level_2_pruned <- annot_merge$HNOCA_annot_level_2_pruned
save(org_merge, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_CSS_only_Julis_genes.rda")

save(org_merge, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_CSS_only.rda")
# save(org_merge, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_harmony5_julis_genes.rda")
# save(org_merge_subset_WT, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_harmony5_julis_genes_WT.rda")
# save(org_merge_subset_KO, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_harmony5_julis_genes_KO.rda")
# 

org_merge <- subset(org_merge, features = rownames(juli_WT))
