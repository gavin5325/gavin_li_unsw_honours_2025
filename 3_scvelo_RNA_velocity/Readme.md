# About this directory

This directory contains scripts for performing RNA velocity analysis on the genotype (ZMYND8 WT/KO) specific annData objects generated from the previous step. Analysis was performed on UNSW’s HPC cluster Katana (DOI: 10.26190/669X-A286).

`250621_RNA_velocity.ipynb`: Initial RNA velocity with all the subsetted objects from 250702_seurat_no_day_60.R. 

`250702_seurat_no_day_60.ipynb`: Performs RNA velocity on WT, KO and combined objects once D60 has been aligned to hg19. 

`250813_differential_kinetics.ipynb`: Performs differential kinetics testing on day 100 samples to determine whether or not splicing and degradation parameters are different for different celltypes for a particular gene. Subsequently, velocity transition graph is updated. Performs a bunch of miscellaneous operations as well. Including latent time plots, as well as Wilcoxon rank-sum tests.

`251008_D60_D100.ipynb`: Performs RNA velocity on the cumulative day 60 and 100 annData objects.
