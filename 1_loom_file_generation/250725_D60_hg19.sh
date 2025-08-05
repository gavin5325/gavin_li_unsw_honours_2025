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
REFERENCE_GENOME_PATH="/mnt/Data0/REFERENCE/HUMAN/GRCh37_hg19/GRCh37/Annotation/Genes/genes.gtf"
EXPECTED_VELOCYTO_OUTPUT="/mnt/Scratch/PROJECTS/ZMYND8/RNA_velocity/Gavin_Li/RESULTS/250725_Velocyto_organoidB2478_retry"
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

    #Day 60, Batch 2
    "/mnt/Data1/PROJECTS/ZMYND8/JW/Katana230202/Zmynd8_Pub_D60/DATA/WAN11267/scrnaseq/cellranger/count/sample-S1_WT"
    "/mnt/Data1/PROJECTS/ZMYND8/JW/Katana230202/Zmynd8_Pub_D60/DATA/WAN11267/scrnaseq/cellranger/count/sample-S2_KO_Zmynd8"
)

for sample in "${samples[@]}"; do
    cellranger_to_loom_1 "$sample" "$REFERENCE_GENOME_PATH"
done
