#!/bin/bash
#SBATCH --job-name=SC13_preprocess
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=400G
#SBATCH --time=12:00:00

source ~/.bashrc
conda activate bamtools

PARQUET_DIR="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC13/curioseeker/work/42/04ecac911243b6a098838885b06b1d/SC13-r1-r2-merged/"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC13/harpy/barcodes_in_tissue.txt"
OUTPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC13/bam_output/preprocessed_reads.csv"

python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/Curio_seeker/bam/curio_preprocess.py \
  --parquet_dir "$PARQUET_DIR" \
  --tissue_barcodes_path "$TISSUE_BARCODES_PATH" \
  --output_path "$OUTPUT_PATH" \