#!/bin/bash
#SBATCH --job-name=SCA024_bam
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=400G
#SBATCH --time=30:00:00

source ~/.bashrc
conda activate bamtools

INPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA024/bam_output/possorted_genome_bam_preprocessed.csv"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA024/harpy/barcodes_in_tissue.txt"
OUTPUT_DIR="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA024/bam_output/output_metrics" # directory will be created
TISSUE_AREA=13.03
STEP_SIZE=10000000 # put to 0 for no sampling

python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/2024-07-UTBenchmark-SpC/Visium_HD/bam/calculate_metrics_VisiumHD_3prime.py \
  --input_path $INPUT_PATH \
  --output_dir $OUTPUT_DIR \
  --tissue_barcodes_path $TISSUE_BARCODES_PATH \
  --tissue_area $TISSUE_AREA \
  --step_size $STEP_SIZE
