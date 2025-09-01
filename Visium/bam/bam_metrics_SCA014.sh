#!/bin/bash
#SBATCH --job-name=SCA014_bam
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256G
#SBATCH --time=24:00:00

BAM_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium/SCA014/COUNTS/SCA014/SCA014_reseq/outs/possorted_genome_bam.bam"
ALL_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/harpy/all_barcodes.txt"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/harpy/barcodes_in_tissue.txt"
PROBE_SET_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/raw/Visium/SCA014/COUNTS/SCA014/SCA014_reseq/outs/probe_set.csv"
OUTPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/bam_output"
TISSUE_AREA_MM2=28.03
STEP_SIZE=10_000_000

python3 /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/Visium/bam/bam_metrics_visium.py \
  --bam_path "$BAM_PATH" \
  --all_barcodes_path "$ALL_BARCODES_PATH" \
  --tissue_barcodes_path "$TISSUE_BARCODES_PATH" \
  --probe_set_path "$PROBE_SET_PATH" \
  --output_path "$OUTPUT_PATH" \
  --tissue_area_mm2 "$TISSUE_AREA_MM2" \
  --step_size "$STEP_SIZE"