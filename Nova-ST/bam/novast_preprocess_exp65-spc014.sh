#!/bin/bash
#SBATCH --job-name=spc014_preprocess
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=350G
#SBATCH --time=12:00:00

ml SAMtools
source ~/.bashrc
conda activate bamtools

BAM="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp65-sampleSPC014/NovaScope/novascope_output/align/NovaSeq6000_20241107/03D1/Exp65-sampleSPC014/bam/sttoolsAligned.sortedByCoord.out.bam"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp65-sampleSPC014/harpy/barcodes_in_tissue.txt"
OUTPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp65-sampleSPC014/bam_output/preprocessed_reads.csv"
FASTQ_READ_COUNT=

python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/2024-07-UTBenchmark-SpC/Nova-ST/bam/novast_preprocess.py \
  --bam_path "$BAM" \
  --tissue_barcodes_path "$TISSUE_BARCODES_PATH" \
  --output_path "$OUTPUT_PATH" \
  --fastq_read_count "$FASTQ_READ_COUNT"