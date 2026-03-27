#!/bin/bash
#SBATCH --job-name=SCA015_bam
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100G
#SBATCH --time=03:00:00
#SBATCH --output=SCA015_metrics_%j.out
#SBATCH --error=SCA015_metrics_%j.err



INPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA015/bam_output/possorted_genome_bam_preprocessed.csv"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA015/harpy/barcodes_in_tissue.txt"
OUTPUT_DIR="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA015/bam_output/output_metrics" # directory will be created
TISSUE_AREA=34.46
STEP_SIZE=10000000 # put to 0 for no sampling

python bam_calculate_metrics.py \
  --input_path $INPUT_PATH \
  --output_dir $OUTPUT_DIR \
  --tissue_barcodes_path $TISSUE_BARCODES_PATH \
  --tissue_area $TISSUE_AREA \
  --step_size $STEP_SIZE
