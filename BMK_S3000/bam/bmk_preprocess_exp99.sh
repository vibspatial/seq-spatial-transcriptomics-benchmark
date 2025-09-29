#!/bin/bash
#SBATCH --job-name=exp99_preprocess
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=350G
#SBATCH --time=12:00:00

source ~/.bashrc
conda activate bamtools

UMI_GENE="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/BSTMAtrix/exp99_28Feb2025/umi_gene.tsv.gz"
FEATURES="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/BSTMAtrix/exp99_28Feb2025/features.tsv"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/harpy/barcodes_in_tissue.txt"
OUTPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/bam_output/preprocessed_reads.csv"
FASTQ_READ_COUNT=522855800

python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/BMK_S1000/bam/bmk_preprocess.py \
  --umi_gene_path "$UMI_GENE" \
  --features_path "$FEATURES" \
  --tissue_barcodes_path "$TISSUE_BARCODES_PATH" \
  --output_path "$OUTPUT_PATH" \
  --fastq_read_count "$FASTQ_READ_COUNT"