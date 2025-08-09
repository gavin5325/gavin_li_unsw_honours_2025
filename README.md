# Gavin Li Honours 2025 project
Project title: Uncovering developmental trajectories in brain organoids using single-cell RNA sequencing

![Project flowchart](https://raw.githubusercontent.com/gavin5325/gavin_li_unsw_honours_2025/main/finalflowchart_README.png)

As per the project flowchart, each folder in this repository represents a specific stage of the project, starting from scRNA-seq reads obtained from the Cell Ranger pipeline.

### 1_loom_file_generation
Assigns the aforementioned scRNA-seq reads into spliced, unspliced and ambiguous count matrices. Occurred on Voineagu lab’s RNA2 server via X2Go client.

### 2_seurat_preprocessing_integration_conversion
Pre-processes sample-wide Seurat objects and removes batch effects due to sequencing depths sequences. Occurred on RStudio with R-based packages. Also converts these Seurat objects into annData objects for downstream RNA velocity analysis.

### 3_scvelo_RNA_velocity
RNA velocity analysis of annData object containing genotype-specific (ZMYND8 WT/KO) cells. Occurred on JupyterLab.
