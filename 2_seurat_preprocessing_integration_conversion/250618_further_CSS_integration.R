library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(Matrix)
library(patchwork)
library(simspec)

org_merge <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_CSS_only_combined.rda"))

# Create keys in metadata
org_merge@meta.data$Genotype <- ifelse(grepl("WT", org_merge@meta.data$Sample_ID), "WT", "KO")
org_merge@meta.data$Day <- sub("D(\\d+).*", "\\1", org_merge@meta.data$Sample_ID)
# save(org_merge, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_CSS_only.rda")
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
# save(org_merge, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_CSS_only.rda")


output_dir <- "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250618_further_CSS_integration/"
group_bys <- c("Genotype", "Day", "css_cluster", "HNOCA_annot_level_1_pruned")

org1 <- cluster_sim_spectrum(object = org_merge, label_tag = "Sample_ID")
org2 <-cluster_sim_spectrum(object = org_merge, label_tag = "Sample_ID", cluster_resolution = 0.8)
org3 <-cluster_sim_spectrum(object = org_merge, label_tag = "Sample_ID", cluster_resolution = 1)
org4 <-cluster_sim_spectrum(object = org_merge, label_tag = "Sample_ID", train_on = "pseudo")
org5 <-cluster_sim_spectrum(object = org_merge, label_tag = "Sample_ID", num_pcs_use = 19)
org6 <-cluster_sim_spectrum(object = org_merge, label_tag = "Sample_ID", num_pcs_use = 30)
org7 <-cluster_sim_spectrum(object = org_merge, label_tag = "Sample_ID", redo_pca = TRUE)
org8<-cluster_sim_spectrum(object = org_merge, label_tag = "Sample_ID",
                           cluster_resolution = 0.4,
                           corr_method = "pearson", 
                           spectrum_type = "corr_kernel")

orgs <- list(org1, org6)
names(orgs) <- c("org1", "org6")
for (i in seq_along(orgs)) {
  css_mat <- orgs[[i]][["css"]]@cell.embeddings
  num_na <- sum(is.na(css_mat))
  print(num_na)
  css_mat[is.na(css_mat)] <- 0
  orgs[[i]][["css"]]@cell.embeddings <- css_mat
  print(sum(is.na(orgs[[i]][["css"]]@cell.embeddings)))
}
# save(orgs, file = paste0(output_dir, 'org.rda'))

# a <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250618_further_CSS_integration/org.rda"))

org_merge <- RunUMAP(org_merge, reduction = "css", dims = 1:12, reduction.name = "umap.css.sigPCs")
org_merge <- FindNeighbors(org_merge, reduction = "css", dims = 1:12)
org_merge <- FindClusters(org_merge, resolution = 0.6, cluster.name = "clusters_CSS_SigPs_Point6")

css_mat <- org1[["css"]]@cell.embeddings
sum(is.na(css_mat)) 
css_mat[is.na(css_mat)] <- 0
org1[["css"]]@cell.embeddings <- css_mat
org1 <- RunUMAP(org1, reduction = "css", dims = 1:20, reduction.name = "umap.css.parameterised")
org1 <- FindNeighbors(org1, reduction = "css", dims = 1:20)
org1 <- FindClusters(org1, resolution = 0.6, cluster.name = "css_cluster")
org1 <- RunTSNE(org1, reduction = "css", dims = 1:20, reduction.name = "tsne.css.parameterised", check_duplicates = FALSE)

org2 <- RunUMAP(org2, reduction = "css", dims = 1:20, reduction.name = "umap.css.parameterised")
org2 <- FindNeighbors(org2, reduction = "css", dims = 1:20)
org2 <- FindClusters(org2, resolution = 0.8, cluster.name = "css_cluster")

org3 <- RunUMAP(org3, reduction = "css", dims = 1:20, reduction.name = "umap.css.parameterised")
org3 <- FindNeighbors(org3, reduction = "css", dims = 1:20)
org3 <- FindClusters(org3, resolution = 1, cluster.name = "css_cluster")

org4 <- RunUMAP(org4, reduction = "css", dims = 1:20, reduction.name = "umap.css.parameterised")
org4 <- FindNeighbors(org4, reduction = "css", dims = 1:20)
org4 <- FindClusters(org4, resolution = 0.6, cluster.name = "css_cluster")

org5 <- RunUMAP(org5, reduction = "css", dims = 1:19, reduction.name = "umap.css.parameterised")
org5 <- FindNeighbors(org5, reduction = "css", dims = 1:19)
org5 <- FindClusters(org5, resolution = 0.6, cluster.name = "css_cluster")

css_mat <- org6[["css"]]@cell.embeddings
sum(is.na(css_mat))  
css_mat[is.na(css_mat)] <- 0
org6[["css"]]@cell.embeddings <- css_mat
org6 <- RunUMAP(org6, reduction = "css", dims = 1:30, reduction.name = "umap.css.parameterised")
org6 <- FindNeighbors(org6, reduction = "css", dims = 1:30)
org6 <- FindClusters(org6, resolution = 0.6, cluster.name = "css_cluster")
org6 <- RunTSNE(org6, reduction = "css", dims = 1:30, reduction.name = "tsne.css.parameterised", check_duplicates = FALSE)


org7 <- RunUMAP(org7, reduction = "css", dims = 1:20, reduction.name = "umap.css.parameterised")
org7 <- FindNeighbors(org7, reduction = "css", dims = 1:20)
org7 <- FindClusters(org7, resolution = 0.6, cluster.name = "css_cluster")

org8 <- RunUMAP(org8, reduction = "css", dims = 1:20, reduction.name = "umap.css.parameterised")
org8 <- FindNeighbors(org8, reduction = "css", dims = 1:20)
org8 <- FindClusters(org8, resolution = 0.4, cluster.name = "css_cluster")

for (name in names(orgs)) {

  pdf(paste0(output_dir, '_tsne_', name, ".pdf"), width = 10, height = 7)
  for (group_by in group_bys) {
    p <- DimPlot(orgs[[name]], reduction = "umap.css.parameterised", group.by = group_by)
    print(p)
    t <- DimPlot(orgs[[name]], reduction = "tsne.css.parameterised", group.by = group_by)
    print(t)
  }
  dev.off()
}


########################### 
org1_WT <- subset(org1, subset = Genotype == "WT")
org1_KO <- subset(org1,  subset = Genotype == "KO")
org6_WT <- subset(org6,  subset = Genotype == "WT")
org6_KO <- subset(org6,  subset = Genotype == "KO")

orgs_list <- list(org1_WT, org1_KO, org6_WT, org6_KO)
names(orgs_list) <- c("org1_WT", "org1_KO", "org6_WT", "org6_KO")

for (name in names(orgs_list)) {

  org_merge <- orgs_list[[name]]
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

  org_merge <- FindNeighbors(object = org_merge, dims = 1:final_pc, reduction = "pca")
  org_merge <- FindClusters(object = org_merge, resolution = 0.6, cluster.name = "unintegrated_clusters")
  org_merge <- RunUMAP(org_merge, dims = 1:final_pc, reduction = "pca", reduction.name = "umap.unintegrated")

  if (name == "org6_WT" | name == "org6_KO") {
    org_merge <-cluster_sim_spectrum(object = org_merge, label_tag = "Sample_ID", num_pcs_use = 30)
  } else {
    org_merge <-cluster_sim_spectrum(object = org_merge, label_tag = "Sample_ID")
  }
  
  css_mat <- org_merge[["css"]]@cell.embeddings
  sum(is.na(css_mat)) 
  css_mat[is.na(css_mat)] <- 0
  org_merge[["css"]]@cell.embeddings <- css_mat

  if (name == "org6_WT" | name == "org6_KO") {
    org_merge <- RunUMAP(org_merge, reduction = "css", dims = 1:30, reduction.name = "umap.css.parameterised")
    org_merge <- FindNeighbors(org_merge, reduction = "css", dims = 1:30)
    org_merge <- FindClusters(org_merge, resolution = 0.6, cluster.name = "css_cluster")
    org_merge <- RunTSNE(org_merge, reduction = "css", dims = 1:30, reduction.name = "tsne.css.parameterised", check_duplicates = FALSE)
  } else {
    org_merge <- RunUMAP(org_merge, reduction = "css", dims = 1:20, reduction.name = "umap.css.parameterised")
    org_merge <- FindNeighbors(org_merge, reduction = "css", dims = 1:20)
    org_merge <- FindClusters(org_merge, resolution = 0.6, cluster.name = "css_cluster")
    org_merge <- RunTSNE(org_merge, reduction = "css", dims = 1:20, reduction.name = "tsne.css.parameterised", check_duplicates = FALSE)
  }

  pdf(paste0(output_dir, name, ".pdf"), width = 10, height = 7)
  for (group_by in group_bys) {
    p <- DimPlot(org_merge, reduction = "umap.css.parameterised", group.by = group_by)
    print(p)
    t <- DimPlot(org_merge, reduction = "tsne.css.parameterised", group.by = group_by)
    print(t)
  }
  dev.off()
  
}
