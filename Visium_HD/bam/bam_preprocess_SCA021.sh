#!/bin/bash
#SBATCH --job-name=SCA021_bam
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00

source ~/.bashrc
conda activate bamtools

BAM_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA021/SpaceRanger/outs/possorted_genome_bam.bam"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA021/harpy/barcodes_in_tissue.txt"
OUTPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA021/bam_output/possorted_genome_bam_preprocessed.csv"

python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/2024-07-UTBenchmark-SpC/Visium_HD/bam/bam_preprocess_VisiumHD_3prime.py \
  --bam_path $BAM_PATH \
  --tissue_barcodes_path $TISSUE_BARCODES_PATH \
  --output_path $OUTPUT_PATH
