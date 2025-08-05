library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(reticulate)
library(Matrix)
use_virtualenv("~/scvelo_env", required = TRUE)
py_config()
scv <- import("scvelo")
sc <- import("scanpy")
ad <- import("anndata")
np <- import("numpy")
pd <- import("pandas")
scipy_sparse <- import("scipy.sparse")


org_merge <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_CSS_only_combined.rda"))

# org_merge <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_CSS_only_Julis_genes.rda"))
# org_merge <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration_CSS_only.rda"))
# org_merge <- get(load("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250524_seurat_integration/Seurat.objects.PostIntegration.rda"))
# org_merge_WT <- subset(org_merge, grepl("WT", Sample_ID))
# org_merge <- org_merge_WT

DefaultAssay(org_merge) <- "RNA"



# org_merge_clean <- org_merge
# org_merge_clean[["spliced"]] <- NULL
# org_merge_clean[["unspliced"]] <- NULL
# org_merge_clean[["ambiguous"]] <- NULL
# SaveH5Seurat(org_merge_clean, filename = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/org_merge_RNAonly_harmony5.h5Seurat")
# Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/org_merge_RNAonly_harmony5.h5Seurat", dest = "h5ad")

# SaveH5Seurat(org_merge_clean, filename = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/org_merge_RNAonly.h5Seurat")
# Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/org_merge_RNAonly.h5Seurat", dest = "h5ad")
# org_merge_h5ad <- sc$read("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/org_merge_RNAonly_harmony5.h5ad")

SaveH5Seurat(org_merge, filename = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined.h5Seurat", overwrite = TRUE)
# Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/Seurat.objects.PostIntegration_CSS_only_combined.h5Seurat", dest = "h5ad", overwrite = TRUE)
Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined.h5Seurat", dest = "h5ad", overwrite = TRUE)


org_merge_WT <- subset(org_merge, subset = Genotype == "WT")
SaveH5Seurat(org_merge_WT, filename = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined_WT.h5Seurat", overwrite = TRUE)
# Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/Seurat.objects.PostIntegration_CSS_only_combined_WT.h5Seurat", dest = "h5ad", overwrite = TRUE)
Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined_WT.h5Seurat", dest = "h5ad", overwrite = TRUE)

org_merge_KO <- subset(org_merge, subset = Genotype == "KO")
SaveH5Seurat(org_merge_KO, filename = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined_KO.h5Seurat", overwrite = TRUE)
# Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250528_seurat_WTKO_conversion/Seurat.objects.PostIntegration_CSS_only_combined_KO.h5Seurat", dest = "h5ad", overwrite = TRUE)
Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined_KO.h5Seurat", dest = "h5ad", overwrite = TRUE)



###############
org_merge <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined.h5Seurat")
adata <- sc$read("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined.h5ad")
adata_WT <- sc$read("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined_WT.h5ad")

WT_spliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/WT_spliced_Juli.mtx")
WT_unspliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/WT_unspliced_Juli.mtx")
WT_genes <- read.csv("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/WT_genes_Juli.csv", header = FALSE, stringsAsFactors = FALSE)[,1]
WT_cells <- read.csv("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250522_seurat_merging/WT_cells_Juli.csv", header = FALSE, stringsAsFactors = FALSE)[,1]

WT_genes <- WT_genes[-1]
WT_cells <- WT_cells[-1]

rownames(WT_spliced) <- WT_genes
colnames(WT_spliced) <- WT_cells
rownames(WT_unspliced) <- WT_genes
colnames(WT_unspliced) <- WT_cells

obs_r <- py_to_r(adata_WT$obs)
obs_copy <- r_to_py(obs_r) 

adata_WT_full <- ad$AnnData(
  X=t(as.matrix(WT_spliced)),    
  obs=obs_copy,
  var=data.frame(gene_names = rownames(WT_spliced))
)

adata_WT_full$uns = reticulate::r_to_py(adata_WT$uns)$copy()
adata_WT_full$obsm = reticulate::r_to_py(adata_WT$obsm)$copy()
# adata_full$obsp = adata$obsp$copy()
adata_WT_full$obsp = reticulate::r_to_py(adata_WT$obsp)$copy()

spliced_csr <- as(t(WT_spliced), "dgCMatrix")
adata_WT_full$layers["spliced"] <- scipy_sparse$csr_matrix(spliced_csr)
unspliced_csr <- as(t(WT_unspliced), "dgCMatrix")
adata_WT_full$layers["unspliced"] <- scipy_sparse$csr_matrix(unspliced_csr)
adata_WT_full$X <- scipy_sparse$csr_matrix(spliced_csr)

emb <- adata_WT_full$obsm['X_umap.css.sigPCs']
clusters <- adata_WT_full$obs$HNOCA_annot_level_1_pruned
clusters <- factor(clusters)
rownames(emb) <- names(clusters) <- adata_WT_full$obs_names$values
col <- c(
  "#E6194B",  # red
  "#3CB44B",  # green
  "#0082C8",  # blue
  "#F58231",  # orange
  "#911EB4",  # purple
  "#46F0F0",  # cyan
  "#F032E6",  # pink/magenta
  "#D2F53C",  # lime
  "#FABEBE",  # light pink
  "#008080",  # teal
  "#E6BEFF",  # lavender
  "#AA6E28"   # brown
)
cell.cols <- col[clusters]
names(cell.cols) <- names(clusters)

plot(emb, col=cell.cols, pch=20, xlab='UMAP X', ylab='UMAP Y')
legend("bottomright", legend=levels(clusters), col=col, pch=16, cex=0.4)


# scv$pp$filter_and_normalize(adata_WT_full, n_top_genes=2000L)
# # Numbers were from the upstream org_merge seurat object
# sc$tl$pca(adata_WT_full, n_comps = 12L)
# sc$pp$neighbors(adata_WT_full, n_pcs = 12L, n_neighbors = 20L)
# scv$pp$moments(adata_WT_full, n_pcs = 12L, n_neighbors = 20L)
# scv$tl$recover_dynamics(adata_WT_full)

scv$pp$filter_genes_dispersion(adata_WT_full, flavor="seurat")
adata_WT_full$obsm["X_pca"] = adata_WT_full$obsm["X_css"]
scv$pp$neighbors(adata_WT_full, use_rep="X_pca")
scv$pp$moments(adata_WT_full, n_pcs = 12L, n_neighbors = 20L)
scv$tl$recover_dynamics(adata_WT_full)

scv$tl$velocity(adata_WT_full, mode='dynamical')
scv$tl$velocity_graph(adata_WT_full)
pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250618_further_CSS_integration/combined_WT_dynamical_dispersion_Juli.pdf", width = 11, height = 7)
scv$pl$velocity_embedding_grid(
  adata_WT_full,
  basis = 'X_umap.css.sigPCs',
  color = 'HNOCA_annot_level_1_pruned',
  legend_fontsize = 7,
  legend_loc = 'lower left',
  size = 50,
  alpha = 1.0,
  scale = 0.6,
  title = 'WT velocity (Sample ID)'
)

scv$pl$velocity_embedding_grid(
  adata_WT_full,
  basis = 'X_umap.css.sigPCs',
  color = 'Day',
  legend_fontsize = 7,
  legend_loc = 'lower left',
  size = 50,
  alpha = 1.0,
  scale = 0.6,
  title = 'WT velocity (Harmony clusters)'
)
dev.off()

a <- adata_WT_full$layers["velocity"]
adata_WT_full$layers["velocity_dynamical"] <- a

b <- adata_WT_full$uns["velocity_graph"]
adata_WT_full$uns["velocity_graph_dynamical"] <- b

c <- adata_WT_full$uns["velocity_params"]
adata_WT_full$uns["velocity_params_dynamical"] <- c

scv$tl$velocity(adata_WT_full, mode="stochastic")
scv$tl$velocity_graph(adata_WT_full)

pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250618_further_CSS_integration/combined_WT_stochastic_dispersion_Juli.pdf", width = 11, height = 7)
scv$pl$velocity_embedding_grid(
  adata_WT_full,
  basis = 'X_umap.css.sigPCs',
  color = 'HNOCA_annot_level_1_pruned',
  legend_fontsize = 7,
  legend_loc = 'lower left',
  size = 50,
  alpha = 1.0,
  scale = 0.6,
  title = 'WT velocity (Sample ID)'
)

scv$pl$velocity_embedding_grid(
  adata_WT_full,
  basis = 'X_umap.css.sigPCs',
  color = 'Day',
  legend_fontsize = 7,
  legend_loc = 'lower left',
  size = 50,
  alpha = 1.0,
  scale = 0.6,
  title = 'WT velocity (Harmony clusters)'
)
dev.off()

d <- adata_WT_full$layers["velocity"]
adata_WT_full$layers["velocity_stochastic"] <- d

e <- adata_WT_full$uns["velocity_graph"]
adata_WT_full$uns["velocity_graph_stochastic"] <- e

f <- adata_WT_full$uns["velocity_params"]
adata_WT_full$uns["velocity_params_stochastic"] <- f

scv$tl$latent_time(adata_WT_full)
pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250618_further_CSS_integration/latent_time_Juli.pdf", width = 10, height = 10)
scv$pl$scatter(adata_WT_full, basis="X_umap.css.sigPCs", color='latent_time', color_map='gnuplot', size=30)
dev.off()

pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250618_further_CSS_integration/WT_root_end_cels_umap_refined.pdf", width = 20, height = 7)
scv$pl$velocity_embedding_grid(adata_WT_full,
                                 basis = 'X_umap.css.sigPCs',
                                 color = c("root_cells", "end_points"),
                                 legend_fontsize = 7,
                                 legend_loc = "lower left",
                                 size = 50,
                                 alpha = 1.0,
                                 arrow_size = 1.5,
                                 density = 1.5,
                                 title = c('root', 'end'))
dev.off()

scv$tl$rank_velocity_genes(adata_WT_full, groupby='Sample_ID', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250618_further_CSS_integration/MIR4688.pdf", width = 20, height = 7)
scv$pl$scatter(adata_WT_full, basis="X_umap.css.sigPCs", color=c("latent_time", "MIR4688"))
dev.off()

pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250618_further_CSS_integration/MIR4688.pdf", width = 40, height = 10)
scv$pl$velocity(adata_WT_full, var_names=c('9922'), color='HNOCA_annot_level_1_pruned', basis = "X_umap.css.sigPCs", size = 10, legend_loc = 'upper left')
dev.off()

topgenes <- adata_WT_full$var["fit_likelihood"]
topgenes_vals <- topgenes[,1]
names(topgenes_vals) <- rownames(topgenes)
topgenes_vals <- sort(topgenes_vals, decreasing=TRUE)
head(topgenes_vals)

pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250604_WT_retry/phase.pdf", width = 30, height = 7)
scv$pl$scatter(adata_WT_full, basis=names(topgenes_vals)[1:5], ncols=5, frameon=FALSE, color="Sample_ID")
dev.off()

scv$tl$velocity_confidence(adata_WT_full)
keys = c('velocity_length', 'velocity_confidence')
pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250618_further_CSS_integration/phase_refined.pdf", width = 20, height = 7)
scv$pl$scatter(adata_WT_full, basis="X_umap.css.sigPCs", color=keys, cmap='coolwarm', perc=c(5, 95))
dev.off()

adata1 <- adata_WT_full
# adata1$raw <- NULL
adata1$write("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250630_css_only_combined_integration/PostIntegration_CSS_only_combined_WT.h5ad")


################
# adata_WT <- sc$read("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250604_WT_retry/combined_WT_modified.h5ad")
