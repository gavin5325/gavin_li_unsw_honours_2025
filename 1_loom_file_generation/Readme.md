# About this directory

This directory contains scripts that I used to align scRNA-seq reads from the Cell Ranger (v7.1.0) pipeline to generate spliced, unspliced and ambiguous count matrices for D20, D60 and D100 ZMYND8 WT/KO brain organoid samples. D20 and D100 were aligned to hg38 whilst D60 was aligned to hg19.

`250424_Velocyto_organoidB2478.sh` performs the first round of read alignment using velocyto to generate spliced, unspliced and ambiguous loom count matrices across D20, 60 and 100 WT and KO samples using hg38.

`250725_D60_hg19.sh` is an additional round of read alignment just for aligning D60 WT and KO samples to hg19, which is the version used in the upstream CellRanger pipeline.
## Note:
These results overwrite the D60 results from the previous script.

All loom count matrices were transferred from RNA2 server to Katana using `rsync -avh`.
