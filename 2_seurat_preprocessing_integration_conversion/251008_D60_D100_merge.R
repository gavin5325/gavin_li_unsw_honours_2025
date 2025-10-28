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

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(simspec)
library(harmony)
library(Matrix)

source("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Scripts/250520_seurat_preprocessing.R")

output_dir <- "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/251008_D60_D100_merge/"

org_merge <- LoadH5Seurat("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/org_merge_changed.h5Seurat")
org_merge_D60 <- subset(org_merge, subset = Day == 60)
org_merge_D100 <- subset(org_merge, subset = Day == 100)
org_merge_D20 <- subset(org_merge, subset = Day == 20)

org_merge_by_D60 <- subset(org_merge, subset = Day %in% c(20, 60))

org_merge_by_D60_WT <- subset(org_merge_by_D60, subset = Genotype == "WT")
org_merge_by_D60_KO <- subset(org_merge_by_D60, subset = Genotype == "KO")

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

org_merge_by_D60 <- lseurat_integration(org_merge_by_D60, "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/251008_D60_D100_merge/org_merge_by_D60.h5Seurat")
orgs_list <- list(org_merge_by_D60)
names(orgs_list) <- c("org_merge_by_D60")

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

SaveH5Seurat(org_merge_by_D60, filename="/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/251008_D60_D100_merge/org_merge_by_D60.h5Seurat", overwrite = TRUE)
Convert("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/251008_D60_D100_merge/org_merge_by_D60.h5Seurat", dest="h5ad", overwrite = TRUE)


adata <- sc$read("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/251008_D60_D100_merge/org_merge_by_D60.h5ad")

spliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/251008_D60_D100_merge/org_merge_by_D60_spliced.mtx")
unspliced <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/251008_D60_D100_merge/org_merge_by_D60_unspliced.mtx")
ambiguous <- readMM("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/251008_D60_D100_merge/org_merge_by_D60_ambiguous.mtx")

genes <- readLines("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/251008_D60_D100_merge/org_merge_by_D60_spliced_genes.txt")
cells <- readLines("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/251008_D60_D100_merge/org_merge_by_D60_spliced_cells.txt")

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
adata_full$write("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/251008_D60_D100_merge/org_merge_by_D60.h5ad")

#######

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# ---------------------------
# 1️⃣ Genes of interest
# ---------------------------
# genes_of_interest <- c("PAX6", "NEUROD1", "NEUROD6")
genes_of_interest <- c("BLC11B", "MAP2", "ENO2")

# ---------------------------
# 2️⃣ Extract counts with genotype and day
# ---------------------------
extract_counts <- function(assay_name, genes) {
  counts_matrix <- GetAssayData(org_merge, assay = assay_name, slot = "counts")
  df <- as.data.frame(t(as.matrix(counts_matrix[genes, ])))
  df$Cell <- rownames(df)
  
  # Include only Day 20, 60, 100
  metadata <- org_merge@meta.data %>%
    rownames_to_column(var = "Cell") %>%
    filter(Day %in% c(20, 60, 100)) %>%
    select(Cell, Day, Genotype)
  
  df <- df %>%
    left_join(metadata, by = "Cell") %>%
    pivot_longer(cols = all_of(genes),
                 names_to = "Gene",
                 values_to = "Expression") %>%
    mutate(Type = assay_name)  # spliced/unspliced
  
  return(df)
}

# Extract spliced and unspliced counts
df_spliced <- extract_counts("spliced", genes_of_interest)
df_unspliced <- extract_counts("unspliced", genes_of_interest)

# Combine both assays
df_all <- bind_rows(df_spliced, df_unspliced)

# ---------------------------
# 3️⃣ Apply log2 transformation with offset
# ---------------------------
df_all <- df_all %>%
  mutate(LogExpression = log2(Expression + 0.5))

# ---------------------------
# 4️⃣ Summarize mean log expression per Day per Gene per Type per Genotype
# ---------------------------
df_summary <- df_all %>%
  group_by(Gene, Day, Type, Genotype) %>%
  summarize(MeanExpression = mean(LogExpression), .groups = "drop") %>%
  mutate(LineGroup = paste(Genotype, Type, sep = "-"))

# Ensure Day is a factor in correct order
df_summary$Day <- factor(df_summary$Day, levels = c(20, 60, 100))

# ---------------------------
# 5️⃣ Plot (log2(exp + 0.5)) expression per gene
# ---------------------------
pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/neuronal_misc_log2.pdf",
    width = 8, height = 5)

for (gene in genes_of_interest) {
  
  df_gene <- df_summary %>% filter(Gene == gene)
  
  p <- ggplot(df_gene, aes(x = Day, y = MeanExpression, group = LineGroup, color = LineGroup)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(title = paste0(gene, " (log2 expression + 0.5)"),
         x = "Day",
         y = "Mean log2(Expression + 0.5)",
         color = "Genotype-RNA Type") +
    scale_color_manual(values = c("WT-spliced" = "forestgreen",
                                  "WT-unspliced" = "green",
                                  "KO-spliced" = "red",
                                  "KO-unspliced" = "orange"))
  
  print(p)  # Each print creates a new page in the PDF
}

dev.off()

#####
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# ---------------------------
# Genes of interest
# ---------------------------
genes_of_interest <- c("PAX6", "NEUROD1", "NEUROD6")
# genes_of_interest <- c("DLX1", "DLX2", "DLX6")

# ---------------------------
# Helper: Extract counts from an assay
# ---------------------------
extract_counts <- function(seurat_obj, assay_name, genes) {
  counts <- GetAssayData(seurat_obj, assay = assay_name, slot = "counts")
  genes_present <- intersect(genes, rownames(counts))
  if (length(genes_present) == 0) {
    stop(paste("None of the genes found in assay:", assay_name))
  }
  
  # Transpose counts to cells x genes
  df <- as.data.frame(t(as.matrix(counts[genes_present, , drop = FALSE])))
  df$cell <- rownames(df)
  
  # Extract metadata
  meta_df <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, Day, Genotype)
  
  # Join counts with metadata
  df <- left_join(df, meta_df, by = "cell")
  df$Assay <- assay_name
  return(df)
}

# ---------------------------
# Extract and combine assays
# ---------------------------
df_spliced <- extract_counts(org_merge, "spliced", genes_of_interest)
df_unspliced <- extract_counts(org_merge, "unspliced", genes_of_interest)
combined_df <- bind_rows(df_spliced, df_unspliced)

# Convert to long format
combined_long <- combined_df %>%
  pivot_longer(all_of(genes_of_interest),
               names_to = "Gene",
               values_to = "Expression") %>%
  # Log2 transform and factor ordering
  mutate(
    LogExpression = log2(Expression + 0.5),
    Day = factor(Day, levels = c(20, 60, 100)),
    Genotype = factor(Genotype, levels = c("WT", "KO"))  # WT first
  )

# ---------------------------
# Plot violins with clear outlines
# ---------------------------
pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/dorsal_violin_plots_log2_outlines.pdf",
    width = 10, height = 6)

for (g in genes_of_interest) {
  
  df_gene <- combined_long %>% filter(Gene == g)
  
  p <- ggplot(df_gene, aes(x = Day, y = LogExpression, fill = Genotype)) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      alpha = 0.7,
      width = 0.8,
      color = "black",  # outline
      size = 0.3        # line thickness
    ) +
    facet_wrap(~ Assay, ncol = 2) +
    scale_fill_manual(values = c("WT" = "#6aa84f", "KO" = "#f1c232")) +
    labs(title = paste0(g, " (log2(Expression + 0.5))"),
         x = "Day",
         y = "log2(Expression + 0.5)",
         fill = "Genotype") +
    theme_bw(base_size = 13) +
    theme(
      panel.background = element_blank(),   # removes the grey background
      panel.grid = element_blank(),         # removes grid lines
      panel.border = element_blank(),       # removes the box around the plot
      strip.background = element_rect(fill = "#f0f0f0"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )
  
  print(p)
}

dev.off()


####
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# ---------------------------
# Genes of interest
# ---------------------------
genes_of_interest <- c("BCL11B", "ENO2", "MAP2", "HES1")
# genes_of_interest <- c("DLX1", "DLX2", "DLX6")

# ---------------------------
# Helper: Extract counts from the counts layer
# ---------------------------
extract_counts <- function(seurat_obj, genes) {
  counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")  # counts layer
  genes_present <- intersect(genes, rownames(counts))
  if (length(genes_present) == 0) {
    stop(paste("None of the genes found in counts layer"))
  }
  
  # Transpose counts to cells x genes
  df <- as.data.frame(t(as.matrix(counts[genes_present, , drop = FALSE])))
  df$cell <- rownames(df)
  
  # Extract metadata
  meta_df <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, Day, Genotype)
  
  # Join counts with metadata
  df <- left_join(df, meta_df, by = "cell")
  return(df)
}

# ---------------------------
# Extract counts
# ---------------------------
combined_df <- extract_counts(org_merge, genes_of_interest)

# Convert to long format
combined_long <- combined_df %>%
  pivot_longer(all_of(genes_of_interest),
               names_to = "Gene",
               values_to = "Expression") %>%
  # Log2 transform and factor ordering
  mutate(
    LogExpression = log2(Expression + 0.5),
    Day = as.numeric(Day),  # convert to numeric for filtering
    Genotype = factor(Genotype, levels = c("WT", "KO"))
  )

# ---------------------------
# Create cumulative day grouping
# ---------------------------
combined_long <- combined_long %>%
  mutate(
    DayGroup = case_when(
      Day == 20 ~ 20,
      Day == 60 ~ 60,
      Day == 100 ~ 100
    )
  )

# Define a helper function for cumulative filtering
get_cumulative_data <- function(df, day_cutoff) {
  df %>% filter(Day <= day_cutoff)
}

# ---------------------------
# Plot violins with cumulative data
# ---------------------------
pdf("/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Results/250726_D60_hg19/MISC_cumulative.pdf",
    width = 8, height = 6)

for (g in genes_of_interest) {
  plots_list <- list()
  
  for (day_cutoff in c(20, 60, 100)) {
    df_gene <- combined_long %>%
      filter(Gene == g) %>%
      get_cumulative_data(day_cutoff) %>%
      mutate(DayGroup = factor(day_cutoff, levels = c(20, 60, 100)))
    
    plots_list[[as.character(day_cutoff)]] <- df_gene
  }
  
  df_gene_cum <- bind_rows(plots_list)
  
  p <- ggplot(df_gene_cum, aes(x = DayGroup, y = LogExpression, fill = Genotype)) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      alpha = 0.7,
      width = 0.8,
      color = "black",
      size = 0.3
    ) +
    scale_fill_manual(values = c("WT" = "#6aa84f", "KO" = "#f1c232")) +
    labs(
      title = paste0(g, " (log2(Expression + 0.5))"),
      x = "Day (cumulative)",
      y = "log2(Expression + 0.5)",
      fill = "Genotype"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      strip.background = element_rect(fill = "#f0f0f0"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )
  
  print(p)
}

dev.off()
