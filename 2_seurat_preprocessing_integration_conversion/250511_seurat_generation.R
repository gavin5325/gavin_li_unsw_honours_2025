# THIS FILE WILL BE USED FOR GENERATING SEURAT OBJECTS FROM ALL SAMPLES

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(hdf5r)


# KO
ldat_day20 <- ReadVelocity(file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/loom_files/KO_52_D20_B7.loom")
rownames(ldat_day20$spliced) <- make.unique(rownames(ldat_day20$spliced))
rownames(ldat_day20$unspliced) <- make.unique(rownames(ldat_day20$unspliced))
rownames(ldat_day20$ambiguous) <- make.unique(rownames(ldat_day20$ambiguous))
KO52_D20 <- as.Seurat(x = ldat_day20)
save(KO52_D20, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/KO52_D20.rda") 

ldat_day60 <- ReadVelocity(file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/loom_files/sample-S2_KO_Zmynd8.loom")
rownames(ldat_day60$spliced) <- make.unique(rownames(ldat_day60$spliced))
rownames(ldat_day60$unspliced) <- make.unique(rownames(ldat_day60$unspliced))
rownames(ldat_day60$ambiguous) <- make.unique(rownames(ldat_day60$ambiguous))
KO_D60 <- as.Seurat(x = ldat_day60)
save(KO_D60, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/KO_D60.rda") 


ldat_day100_b4 <- ReadVelocity(file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/loom_files/sample-ZMYND8_KO_D102.loom")
rownames(ldat_day100_b4$spliced) <- make.unique(rownames(ldat_day100_b4$spliced))
rownames(ldat_day100_b4$unspliced) <- make.unique(rownames(ldat_day100_b4$unspliced))
rownames(ldat_day100_b4$ambiguous) <- make.unique(rownames(ldat_day100_b4$ambiguous))
KO_D100_B4 <- as.Seurat(x = ldat_day100_b4)
save(KO_D100_B4, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/KO_D100_B4.rda") 

ldat_day100_b8_KO52 <- ReadVelocity(file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/loom_files/KO52_D100_B8.loom")
rownames(ldat_day100_b8_KO52$spliced) <- make.unique(rownames(ldat_day100_b8_KO52$spliced))
rownames(ldat_day100_b8_KO52$unspliced) <- make.unique(rownames(ldat_day100_b8_KO52$unspliced))
rownames(ldat_day100_b8_KO52$ambiguous) <- make.unique(rownames(ldat_day100_b8_KO52$ambiguous))
KO52_D100_B8 <- as.Seurat(x = ldat_day100_b8_KO52)
save(KO52_D100_B8, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/KO52_D100_B8.rda") 

ldat_day100_b8_KO58 <- ReadVelocity(file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/loom_files/KO58_D100_B8.loom")
rownames(ldat_day100_b8_KO58$spliced) <- make.unique(rownames(ldat_day100_b8_KO58$spliced))
rownames(ldat_day100_b8_KO58$unspliced) <- make.unique(rownames(ldat_day100_b8_KO58$unspliced))
rownames(ldat_day100_b8_KO58$ambiguous) <- make.unique(rownames(ldat_day100_b8_KO58$ambiguous))
KO58_D100_B8 <- as.Seurat(x = ldat_day100_b8_KO58)
save(KO58_D100_B8, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/KO58_D100_B8.rda")

# Day 20
ldat_day20 <- ReadVelocity(file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/loom_files/WT_D20_B7.loom")
rownames(ldat_day20$spliced) <- make.unique(rownames(ldat_day20$spliced))
rownames(ldat_day20$unspliced) <- make.unique(rownames(ldat_day20$unspliced))
rownames(ldat_day20$ambiguous) <- make.unique(rownames(ldat_day20$ambiguous))
WT_D20 <- as.Seurat(x = ldat_day20)
save(WT_D20, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/WT_D20.rda") 

# Day 60
ldat_day60 <- ReadVelocity(file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/loom_files/sample-S1_WT.loom")
rownames(ldat_day60$spliced) <- make.unique(rownames(ldat_day60$spliced))
rownames(ldat_day60$unspliced) <- make.unique(rownames(ldat_day60$unspliced))
rownames(ldat_day60$ambiguous) <- make.unique(rownames(ldat_day60$ambiguous))
WT_D60 <- as.Seurat(x = ldat_day60)
save(WT_D60, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/WT_D60.rda") 

# Day 100, batch 4
ldat_day100_b4 <- ReadVelocity(file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/loom_files/sample-ZMYND8_WT_D102.loom")
rownames(ldat_day100_b4$spliced) <- make.unique(rownames(ldat_day100_b4$spliced))
rownames(ldat_day100_b4$unspliced) <- make.unique(rownames(ldat_day100_b4$unspliced))
rownames(ldat_day100_b4$ambiguous) <- make.unique(rownames(ldat_day100_b4$ambiguous))
WT_D100_B4 <- as.Seurat(x = ldat_day100_b4)
save(WT_D100_B4, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/WT_D100_B4.rda") 

# Day 100, batch 8
ldat_day100_b8 <- ReadVelocity(file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/loom_files/WT_D100_B8.loom")
rownames(ldat_day100_b8$spliced) <- make.unique(rownames(ldat_day100_b8$spliced))
rownames(ldat_day100_b8$unspliced) <- make.unique(rownames(ldat_day100_b8$unspliced))
rownames(ldat_day100_b8$ambiguous) <- make.unique(rownames(ldat_day100_b8$ambiguous))
WT_D100_B8 <- as.Seurat(x = ldat_day100_b8)
save(WT_D100_B8, file = "/srv/scratch/voineagu/PROJECTS/GavinLi/Combined/Data/WT_D100_B8.rda") 
