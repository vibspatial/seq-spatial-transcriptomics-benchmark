#!/bin/bash
#SBATCH --job-name=spc022_metrics
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256G
#SBATCH --time=36:00:00

source ~/.bashrc
conda activate bamtools

INPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC022/bam_output/preprocessed_reads.csv"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC022/harpy/barcodes_in_tissue.txt"
OUTPUT_DIR="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC022/bam_output/output_metrics"
TISSUE_AREA=47.07
STEP_SIZE=10000000

python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/2024-07-UTBenchmark-SpC/Nova-ST/bam/novast_calculate_metrics.py \
  --input_path $INPUT_PATH \
  --output_dir $OUTPUT_DIR \
  --tissue_barcodes_path $TISSUE_BARCODES_PATH \
  --tissue_area $TISSUE_AREA \
  --step_size $STEP_SIZE