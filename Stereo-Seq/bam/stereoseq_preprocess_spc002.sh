#!/bin/bash
#SBATCH --job-name=spc002_preprocess
#SBATCH --partition=gp_64C_128T_512GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=350G
#SBATCH --time=12:00:00

source ~/.bashrc
conda activate bamtools

ID=spc002
SN=

RAW_BARCODE_GENE_EXP="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/SAW/outs/feature_expression/${SN}_raw_barcode_gene_exp.txt"
GEM="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/SAW/outs/feature_expression/${SN}.raw.gem"
TISSUE_BARCODES_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/harpy/barcodes_in_tissue.txt"
OUTPUT_PATH="/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/${ID}/bam_output/preprocessed_reads.csv"
FASTQ_READ_COUNT=

python /data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/scripts/2024-07-UTBenchmark-SpC/Stereo-Seq/bam/stereoseq_preprocess.py \
  --raw_barcode_gene_exp_path "$RAW_BARCODE_GENE_EXP" \
  --gem_path "$GEM" \
  --tissue_barcodes_path "$TISSUE_BARCODES_PATH" \
  --output_path "$OUTPUT_PATH" \
  --fastq_read_count "$FASTQ_READ_COUNT"