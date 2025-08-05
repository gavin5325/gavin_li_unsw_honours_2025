#!/bin/bash

#PBS -l select=1:ncpus=16:mem=200gb
#PBS -l walltime=12:00:00
#PBS -o job_output.log
#PBS -e job_error.log

# source ~/.bashrc
# source /home/rna2/PROGRAMS/miniconda3/etc/profile.d/conda.sh
# conda activate velocyto

# Download and install dependencies for velocyto if not done
dependencies=("numpy" "scipy" "cython" "numba" "matplotlib" "scikit-learn" "h5py" "click")

for dependency in "${dependencies[@]}"; do

    if ! pip show "$dependency" &>/dev/null 2>&1; then
        pip install "$dependency"
    else
        echo "$dependency already installed"
    fi 
done

# Install velocyto if not already installed
if ! pip show velocyto &>/dev/null 2>&1; then
    pip install velocyto
else
    echo "velocyto already installed"
fi

# Load Samtools (I dont think its used?)
# module load /usr/bin/samtools

# Global variable for reference genome
REFERENCE_GENOME_PATH="/mnt/Data0/REFERENCE/HUMAN/GRCh38_hg38/GRCh38/GTF/Homo_sapiens.GRCh38.84.gtf"
EXPECTED_VELOCYTO_OUTPUT="/mnt/Scratch/PROJECTS/ZMYND8/RNA_velocity/Gavin_Li/RESULTS/250424_Velocyto_organoidB2478_retry"
# MASK_FILE_PATH="/mnt/Scratch/PROJECTS/ZMYND8/RNA_velocity/Gavin_Li/SCRIPTS/hg38_rmsk.gtf"

cellranger_to_loom_1() {
    sample_path=$1
    reference_genome_path=$2

    # velocyto run10x -m "$MASK_FILE_PATH" "$sample_path" "$reference_genome_path"
    velocyto run10x "$sample_path" "$reference_genome_path"

    default_velocyto_output_dir="${sample_path}/velocyto"

    if [ -d "$default_velocyto_output_dir" ]; then
        mv "$default_velocyto_output_dir"/*.loom "$EXPECTED_VELOCYTO_OUTPUT/"
        ### Remove the default_velocyto_output_dir variable
    else
        echo "Warning: velocyto directory for ${sample_path} not found, skipping file move."
    fi
}

samples=(
    # Day 20, Batch 7
    "/mnt/Data1/PROJECTS/ZMYND8/JW/Katana230202/Zmynd8_Pub_D20_231215/LAUNCH_240110/results/cellranger/count/WT_D20_B7/"
    "/mnt/Data1/PROJECTS/ZMYND8/JW/Katana230202/Zmynd8_Pub_D20_231215/LAUNCH_240110/results/cellranger/count/KO_52_D20_B7/"

    #Day 60, Batch 2
    "/mnt/Data1/PROJECTS/ZMYND8/JW/Katana230202/Zmynd8_Pub_D60/DATA/WAN11267/scrnaseq/cellranger/count/sample-S1_WT"
    "/mnt/Data1/PROJECTS/ZMYND8/JW/Katana230202/Zmynd8_Pub_D60/DATA/WAN11267/scrnaseq/cellranger/count/sample-S2_KO_Zmynd8"

    #Day 100, Batch 4
    "/mnt/Data1/PROJECTS/ZMYND8/JW/Katana230202/Zmynd8_Pub_D100/DATA/GOK12012/launch/results/cellranger/count/sample-ZMYND8_WT_D102/"
    "/mnt/Data1/PROJECTS/ZMYND8/JW/Katana230202/Zmynd8_Pub_D100/DATA/GOK12012/launch/results/cellranger/count/sample-ZMYND8_KO_D102/"

    #Day 100, Batch 8
    "/mnt/Data1/PROJECTS/ZMYND8/JW/Katana230202/Zmynd8_Pub_D100_Batch8/LAUNCH_240820/results/cellranger/count/WT_D100_B8/"
    "/mnt/Data1/PROJECTS/ZMYND8/JW/Katana230202/Zmynd8_Pub_D100_Batch8/LAUNCH_240820/results/cellranger/count/KO52_D100_B8/"
    "/mnt/Data1/PROJECTS/ZMYND8/JW/Katana230202/Zmynd8_Pub_D100_Batch8/LAUNCH_240820/results/cellranger/count/KO58_D100_B8/"
)

for sample in "${samples[@]}"; do
    cellranger_to_loom_1 "$sample" "$REFERENCE_GENOME_PATH"
done
