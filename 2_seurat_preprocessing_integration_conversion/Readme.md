# About this directory

These scripts contain the steps I took to pre-process, integrate and convert Seurat objects containing the loom count matrices. Analysis were performed on UNSW’s HPC cluster Katana (DOI: 10.26190/669X-A286). The chronological order of all scripts in this directory are listed below.


`250511_seurat_generation.R`: Generates seurat object from sample-specific looom count matrices, aligned to hg38.

`250520_seurat_preprocessing.R`: A series of functions that performs pre-processing on seurat objects created from previous script.

`250522_seurat_assignment.R`: Assigns each seurat object sample and genotype specific names.

`250524_seurat_merging.R`: Merge each integrated seurat object (regardless of timepoint and genotype) into one merged object.

`250524_seurat_subsetting.R`: Subset the merged Seurat object to just include cells and genes in Juli's gene-by-cell seurat object.

`250524_seurat_integration.R`: Integrate the subsetted, merged Seurat object. Both CSS and Harmony methods used. Also attached Juli's HNOCA cell annotations to the Seurat object.

`250528_seurat_WTKO_conversion.R`: Subset the Seurat object into WT/KO, convert it into annData (h5ad) objects and perform initial round of RNA velocity using scVelo on R, via the reticulate package.

`250618_further_CSS_integration.R`: Applies a variety of CSS parameters to integrate the Seurat object. CSS has been chosen as the to-go integration method.

`250630_css_only_combined_integration.R`: All WT and KO samples integrated, however as separate Seurat objects.

`250702_seurat_no_day_60.R`: Subsets the genotype combined Seurat object to include D20, no D60 (WT). Performs pre-processing and integration again.  

`250726_D60_hg19.R`: Creates a WT/KO seurat object for D60 samples aligned to hg19. Repleaced D60 cells in the genotype-combined Seurat object with these cells. Pre-processed and integrated the resultant Seurat object and converted it to annData. **Edit**: As of 25/08/13, D60 samples have now been aligned to hg38, just like D20 and D100. For convinience, the name is kept as such.
