#!/bin/bash
#SBATCH --job-name=SC11_bam
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256G
#SBATCH --time=12:00:00

ml SAMtools
source ~/.bashrc
conda activate bamtools

BAM_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/bam_output/SC11_tagged.bam"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/harpy/barcodes_in_tissue.txt"
OUTPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/bam_output"
TISSUE_AREA_MM2=8.57
STEP_SIZE=10_000_000

python3 /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/Curio_seeker/bam/bam_metrics_curio.py \
  --bam_path "$BAM_PATH" \
  --tissue_barcodes_path "$TISSUE_BARCODES_PATH" \
  --output_path "$OUTPUT_PATH" \
  --tissue_area_mm2 "$TISSUE_AREA_MM2" \
  --step_size "$STEP_SIZE"